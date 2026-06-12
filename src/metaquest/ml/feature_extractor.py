#!/usr/bin/env python3
"""
MetaQuest Pathogenic Functional Fingerprint (PFF) v4.0

A novel multi-scale feature architecture for metagenomic pathogen protein
classification. Combines:

1. Positional Decomposition — N-terminal (secretion), Core (catalytic),
   C-terminal (sorting) zone-specific features
2. Quantitative PWM Scoring — position weight matrices for virulence
   mechanism active sites (not binary regex)
3. Host-Mimicry Detection — eukaryotic-like motifs bacterial effectors
   use to hijack host signaling
4. Evolutionary Signal — compositional heterogeneity and GC-deviation
   as horizontal gene transfer proxy
5. Fourier-Based Repeat Analysis — spectral features capturing
   beta-solenoids, RTX repeats, TM topology without alignment

Total: ~95 features. No external tool dependencies at inference.

Reference: Novel approach — no existing tool combines positional decomposition
with quantitative PWM scoring and host-mimicry detection in a single feature
vector for metagenomic pathogen classification.
"""

import re
import logging
import math
from collections import Counter
from typing import Dict, List

import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from .pwm_data import (
    AA_ORDER, AA_INDEX, BG_FREQ,
    T3SS_NTERM_PROFILE, SEC_SIGNAL_PROFILE, TAT_MOTIF_PROFILE,
    METALLOPROTEASE_CONTEXT, BLACTAMASE_SERINE, BLACTAMASE_METALLO,
    ELM_PATTERNS, DISORDER_WEIGHTS, KD_HYDROPHOBICITY,
    EISENBERG_HYDROPHOBICITY, AGGREGATION_PROPENSITY,
    AA_GC3_CONTENT, EXPECTED_BACTERIAL_GC3,
    HELIX_PROPENSITY, SHEET_PROPENSITY, TURN_PROPENSITY,
)

logger = logging.getLogger(__name__)

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, **kwargs):
        return iterable


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def _seq_to_indices(sequence: str) -> np.ndarray:
    """Convert sequence to integer index array (20 AAs + unknown=20)."""
    return np.array([AA_INDEX.get(aa, 20) for aa in sequence])


def _disorder_score(sequence: str) -> float:
    """IUPred-style disorder propensity for a sequence stretch."""
    if not sequence:
        return 0.0
    return np.mean([DISORDER_WEIGHTS.get(aa, 0.0) for aa in sequence])


def _hydrophobic_moment(sequence: str, window: int = 11, angle: float = 100.0) -> float:
    """
    Eisenberg hydrophobic moment — measures amphipathicity.
    High values indicate amphipathic helix (pore-forming toxins, membrane-active).
    """
    if len(sequence) < window:
        return 0.0

    angle_rad = math.radians(angle)
    max_moment = 0.0

    for i in range(len(sequence) - window + 1):
        sin_sum = 0.0
        cos_sum = 0.0
        for j in range(window):
            h = EISENBERG_HYDROPHOBICITY.get(sequence[i + j], 0.0)
            theta = j * angle_rad
            sin_sum += h * math.sin(theta)
            cos_sum += h * math.cos(theta)
        moment = math.sqrt(sin_sum ** 2 + cos_sum ** 2) / window
        if moment > max_moment:
            max_moment = moment

    return max_moment


def _score_pwm(sequence: str, pwm: np.ndarray) -> float:
    """
    Score a sequence against a PWM. Returns best log-odds score
    over all positions in the sequence.
    """
    pwm_len = pwm.shape[0]
    if len(sequence) < pwm_len:
        return 0.0

    best_score = -999.0
    for i in range(len(sequence) - pwm_len + 1):
        score = 0.0
        for j in range(pwm_len):
            aa_idx = AA_INDEX.get(sequence[i + j], -1)
            if aa_idx >= 0 and aa_idx < 20:
                score += pwm[j, aa_idx]
        if score > best_score:
            best_score = score

    return max(0.0, best_score)


def _charge_density(sequence: str) -> tuple:
    """Return (positive, negative, net) charge density."""
    if not sequence:
        return 0.0, 0.0, 0.0
    pos = sum(1 for aa in sequence if aa in "RKH") / len(sequence)
    neg = sum(1 for aa in sequence if aa in "DE") / len(sequence)
    return pos, neg, pos - neg


def _polar_nonpolar_ratio(sequence: str) -> float:
    """Ratio of polar to nonpolar residues."""
    if not sequence:
        return 0.0
    polar = sum(1 for aa in sequence if aa in "STNQCYHKRDE")
    nonpolar = sum(1 for aa in sequence if aa in "GALMFWVIP")
    return polar / max(1, nonpolar)


def _sequence_entropy(sequence: str) -> float:
    """Shannon entropy of amino acid composition (0-4.32 bits)."""
    if not sequence:
        return 0.0
    counts = Counter(sequence)
    total = len(sequence)
    entropy = 0.0
    for count in counts.values():
        p = count / total
        if p > 0:
            entropy -= p * math.log2(p)
    return entropy


# ============================================================================
# A. POSITIONAL DECOMPOSITION (30 features)
# ============================================================================

def extract_positional_features(sequence: str) -> Dict[str, float]:
    """
    Tripartite positional decomposition: N-terminal (secretion zone),
    Core (catalytic zone), C-terminal (sorting/anchoring zone).

    Novel: No other pathogen classifier decomposes features by protein region.
    Captures the structured positional bias of virulence factors:
    - Effectors: disordered N-termini for secretion injection
    - Toxins: catalytic cores with specific active sites
    - Surface proteins: C-terminal anchoring signals
    """
    features = {}
    n = len(sequence)

    if n < 60:
        # Protein too short for meaningful decomposition
        for prefix in ("nterm_", "core_", "cterm_"):
            features[f"{prefix}disorder"] = 0.0
            features[f"{prefix}charge_pos"] = 0.0
            features[f"{prefix}charge_neg"] = 0.0
            features[f"{prefix}charge_net"] = 0.0
            features[f"{prefix}amphipathic"] = 0.0
            features[f"{prefix}polar_ratio"] = 0.0
            features[f"{prefix}entropy"] = 0.0
            features[f"{prefix}enrichment1"] = 0.0
            features[f"{prefix}enrichment2"] = 0.0
            features[f"{prefix}enrichment3"] = 0.0
        return features

    # Define zones
    nterm = sequence[:50]
    cterm = sequence[-50:]
    core = sequence[50:-50] if n > 100 else sequence[25:-25]

    for zone, prefix in [(nterm, "nterm_"), (core, "core_"), (cterm, "cterm_")]:
        if not zone:
            zone = sequence  # fallback

        # Common features per zone
        features[f"{prefix}disorder"] = _disorder_score(zone)
        pos, neg, net = _charge_density(zone)
        features[f"{prefix}charge_pos"] = pos
        features[f"{prefix}charge_neg"] = neg
        features[f"{prefix}charge_net"] = net
        features[f"{prefix}amphipathic"] = _hydrophobic_moment(zone)
        features[f"{prefix}polar_ratio"] = _polar_nonpolar_ratio(zone)
        features[f"{prefix}entropy"] = _sequence_entropy(zone)

    # Zone-specific enrichments
    # N-terminal: Ser+Thr (T3SS), Arg diplets (Tat), Met-start fraction
    st_frac = sum(1 for aa in nterm if aa in "ST") / len(nterm)
    rr_count = nterm.count("RR") + nterm.count("KR") + nterm.count("RK")
    pro_frac = nterm.count("P") / len(nterm)
    features["nterm_enrichment1"] = st_frac  # Ser/Thr enrichment (T3SS signal)
    features["nterm_enrichment2"] = rr_count / len(nterm) * 100  # Twin-Arg (Tat)
    features["nterm_enrichment3"] = pro_frac  # Pro enrichment (disorder)

    # Core: Cys pairs (disulfide), metal-coordinating (H/D/E clusters), catalytic
    if core:
        cys_pairs = sum(1 for i in range(len(core) - 1)
                        if core[i] == "C" and core[i+1:i+6].count("C") > 0) / max(1, len(core))
        metal_coord = sum(1 for aa in core if aa in "HDE") / len(core)
        # Catalytic potential: consecutive His-Xn-Asp/Glu patterns
        cat_pattern = len(re.findall(r"H.{1,8}[DE]", core)) / max(1, len(core)) * 100
    else:
        cys_pairs = metal_coord = cat_pattern = 0.0
    features["core_enrichment1"] = cys_pairs  # Disulfide potential
    features["core_enrichment2"] = metal_coord  # Metal coordination
    features["core_enrichment3"] = cat_pattern  # Catalytic site density

    # C-terminal: LPXTG-like, hydrophobic stretch (GPI), charged tail
    lpxtg_score = 1.0 if re.search(r"LP.T[GAS]", cterm) else 0.0
    # Hydrophobic stretch in last 20 aa (GPI-anchor signal)
    last20 = sequence[-20:] if n >= 20 else sequence
    hydro_stretch = max(
        sum(1 for aa in last20[i:i+12] if aa in "AILMFWV")
        for i in range(max(1, len(last20) - 12))
    ) / 12.0 if len(last20) >= 12 else 0.0
    # Charged tail (common in secreted toxins)
    charged_tail = sum(1 for aa in cterm[-15:] if aa in "RKDE") / 15.0 if n >= 15 else 0.0
    features["cterm_enrichment1"] = lpxtg_score  # Sortase substrate
    features["cterm_enrichment2"] = hydro_stretch  # GPI/membrane anchor
    features["cterm_enrichment3"] = charged_tail  # Charged tail

    return features


# ============================================================================
# B. HOST-MIMICRY & EVOLUTIONARY FEATURES (15 features)
# ============================================================================

def extract_mimicry_features(sequence: str) -> Dict[str, float]:
    """
    Detect eukaryotic-like motifs (ELMs) that bacterial effectors use to
    hijack host cell signaling pathways.

    Novel: No pathogen classifier explicitly scores host-mimicry motifs.
    This directly captures the biological mechanism of effector evolution.
    """
    features = {}
    n = len(sequence)
    if n < 30:
        return {
            "mimicry_elm_density": 0.0, "mimicry_lrr_periodicity": 0.0,
            "mimicry_coiled_coil": 0.0, "mimicry_prion_like": 0.0,
            "mimicry_acidic_domain": 0.0, "mimicry_gc_deviation": 0.0,
            "mimicry_compositional_jsd": 0.0, "mimicry_rare_dipeptide": 0.0,
            "mimicry_nls_signal": 0.0, "mimicry_sumo_motif": 0.0,
            "mimicry_death_domain": 0.0, "mimicry_14_3_3": 0.0,
            "mimicry_cyclin_binding": 0.0, "mimicry_ankyrin_repeat": 0.0,
            "mimicry_elm_total": 0.0,
        }

    # ELM density — count eukaryotic-like motifs per 100 residues
    elm_counts = {}
    total_elm = 0
    for elm_name, pattern in ELM_PATTERNS.items():
        try:
            count = len(re.findall(pattern, sequence))
        except re.error:
            count = 0
        elm_counts[elm_name] = count
        total_elm += count

    features["mimicry_elm_density"] = total_elm / n * 100
    features["mimicry_elm_total"] = float(total_elm)

    # Individual important ELMs
    features["mimicry_nls_signal"] = elm_counts.get("elm_nls", 0) / n * 100
    features["mimicry_sumo_motif"] = elm_counts.get("elm_sim", 0) / n * 100
    features["mimicry_death_domain"] = elm_counts.get("elm_death", 0) / n * 100
    features["mimicry_14_3_3"] = elm_counts.get("elm_1433", 0) / n * 100
    features["mimicry_cyclin_binding"] = elm_counts.get("elm_cyclin", 0) / n * 100
    features["mimicry_ankyrin_repeat"] = elm_counts.get("elm_ankyrin", 0) / n * 100

    # LRR periodicity via FFT on hydrophobicity signal
    hydro_signal = np.array([KD_HYDROPHOBICITY.get(aa, 0.0) for aa in sequence])
    if n > 50:
        fft_vals = np.abs(np.fft.rfft(hydro_signal - hydro_signal.mean()))
        freqs = np.fft.rfftfreq(n)
        # LRR period ~24-28 residues → frequency 0.036-0.042
        lrr_band = np.where((freqs > 0.033) & (freqs < 0.045))[0]
        features["mimicry_lrr_periodicity"] = float(np.max(fft_vals[lrr_band])) / n if len(lrr_band) > 0 else 0.0
    else:
        features["mimicry_lrr_periodicity"] = 0.0

    # Coiled-coil heptad score (period ~3.5 residues in hydrophobic moment)
    if n > 28:
        # Use hydrophobic moment with 100° angle (alpha helix)
        features["mimicry_coiled_coil"] = _hydrophobic_moment(sequence, window=28, angle=100.0)
    else:
        features["mimicry_coiled_coil"] = 0.0

    # Low-complexity prion-like region (Q/N-rich windows)
    window = 40
    max_qn = 0.0
    for i in range(max(1, n - window)):
        qn = sum(1 for aa in sequence[i:i+window] if aa in "QN") / window
        if qn > max_qn:
            max_qn = qn
    features["mimicry_prion_like"] = max_qn

    # Acidic activation domain mimic (D/E-rich stretches)
    max_acidic = 0.0
    for i in range(max(1, n - 20)):
        acidic = sum(1 for aa in sequence[i:i+20] if aa in "DE") / 20
        if acidic > max_acidic:
            max_acidic = acidic
    features["mimicry_acidic_domain"] = max_acidic

    # GC-content deviation proxy (HGT indicator)
    # Back-translate to estimate GC3 from amino acid composition
    estimated_gc3 = np.mean([AA_GC3_CONTENT.get(aa, 0.5) for aa in sequence])
    features["mimicry_gc_deviation"] = abs(estimated_gc3 - EXPECTED_BACTERIAL_GC3)

    # Compositional heterogeneity (JSD between halves — multi-domain indicator)
    mid = n // 2
    half1_counts = Counter(sequence[:mid])
    half2_counts = Counter(sequence[mid:])
    p = np.array([half1_counts.get(aa, 0) for aa in AA_ORDER], dtype=float)
    q = np.array([half2_counts.get(aa, 0) for aa in AA_ORDER], dtype=float)
    p = p / max(1, p.sum())
    q = q / max(1, q.sum())
    m = (p + q) / 2
    # JSD = (KLD(p||m) + KLD(q||m)) / 2
    jsd = 0.0
    for i in range(20):
        if p[i] > 0 and m[i] > 0:
            jsd += p[i] * math.log2(p[i] / m[i])
        if q[i] > 0 and m[i] > 0:
            jsd += q[i] * math.log2(q[i] / m[i])
    features["mimicry_compositional_jsd"] = jsd / 2

    # Rare dipeptide enrichment
    # Dipeptides enriched in effectors vs housekeeping (pre-computed from literature)
    effector_dipeptides = {"SS", "NS", "SN", "NN", "TS", "ST", "PS", "SP", "QS", "SQ"}
    dp_count = sum(1 for i in range(n - 1) if sequence[i:i+2] in effector_dipeptides)
    features["mimicry_rare_dipeptide"] = dp_count / max(1, n - 1) * 100

    return features


# ============================================================================
# C. VIRULENCE MECHANISM PWM SCORES (25 features)
# ============================================================================

def extract_virulence_pwm_features(sequence: str) -> Dict[str, float]:
    """
    Quantitative PWM scoring against virulence mechanism active sites.

    Novel: Instead of binary regex (0/1), we use log-odds scoring against
    position weight matrices. A near-miss active site still contributes
    signal, allowing the model to detect degenerate/divergent virulence domains.
    """
    features = {}
    n = len(sequence)

    if n < 30:
        keys = [
            "pwm_t3ss_nterm", "pwm_t4ss_cterm", "pwm_t6ss_spike",
            "pwm_sec_signal", "pwm_tat_signal",
            "pwm_adprt", "pwm_pore_forming", "pwm_metalloprotease",
            "pwm_cys_protease", "pwm_phospholipase", "pwm_ab_toxin",
            "pwm_superantigen", "pwm_rtx_repeat", "pwm_dnase", "pwm_glycosyltransferase",
            "pwm_blactamase_ser", "pwm_blactamase_met", "pwm_aminoglyc",
            "pwm_efflux", "pwm_ribosomal_prot",
            "pwm_integrin_rgd", "pwm_fibronectin", "pwm_collagen_bind",
            "pwm_invasin", "pwm_biofilm_amyloid",
        ]
        return {k: 0.0 for k in keys}

    # --- Secretion system compatibility (5) ---

    # T3SS: N-terminal 25aa profile score
    nterm25 = sequence[1:26] if n > 26 else sequence[:25]  # Skip initiator Met
    features["pwm_t3ss_nterm"] = _score_pwm(nterm25, T3SS_NTERM_PROFILE[:len(nterm25)])

    # T4SS: C-terminal positively charged pattern
    cterm30 = sequence[-30:]
    t4ss_score = sum(1 for aa in cterm30[-15:] if aa in "RK") / 15.0
    # Bonus for R-X(5-9)-R-X-R type pattern
    if re.search(r"R.{5,9}R.R", cterm30):
        t4ss_score += 0.5
    features["pwm_t4ss_cterm"] = t4ss_score

    # T6SS: VgrG/PAAR spike-like composition (enriched in V,G,R at N-term)
    nterm50 = sequence[:50]
    vgrg_aa = sum(1 for aa in nterm50 if aa in "VGR") / len(nterm50)
    features["pwm_t6ss_spike"] = vgrg_aa * 3  # Scale factor

    # Sec signal peptide: von Heijne matrix score in first 30 aa
    if n > 30:
        # Find best cleavage site position (15-30)
        best_sec = 0.0
        for pos in range(15, min(35, n - 16)):
            window = sequence[pos - 13: pos + 3]
            if len(window) == 16:
                score = _score_pwm(window, SEC_SIGNAL_PROFILE)
                if score > best_sec:
                    best_sec = score
        features["pwm_sec_signal"] = best_sec
    else:
        features["pwm_sec_signal"] = 0.0

    # Tat twin-arginine: scan first 50 aa for RR motif with proper context
    features["pwm_tat_signal"] = _score_pwm(sequence[:50], TAT_MOTIF_PROFILE)

    # --- Toxin mechanism scores (10) ---

    # ADP-ribosyltransferase: catalytic Glu in hydrophobic context
    # Pattern: [LIVMF]-X-[LIVMF]-X-E-X-E (ExoA-like)
    adprt_score = 0.0
    for m in re.finditer(r"[LIVMF].[LIVMF].E.E", sequence):
        adprt_score += 1.0
    features["pwm_adprt"] = min(adprt_score, 3.0)

    # Pore-forming: max amphipathic helix potential (18-residue window, 100° angle)
    features["pwm_pore_forming"] = _hydrophobic_moment(sequence, window=18, angle=100.0)

    # Metalloprotease: HEXXH + downstream Glu (PWM-scored context)
    mp_score = 0.0
    for m in re.finditer(r"HE..H", sequence):
        pos = m.start()
        # Check for downstream Glu (20-40 aa later)
        downstream = sequence[pos + 5: pos + 45]
        if "E" in downstream:
            # Score the upstream context
            if pos >= 5:
                upstream = sequence[pos - 5: pos]
                context_score = sum(
                    METALLOPROTEASE_CONTEXT["upstream_preferences"][j, AA_INDEX.get(upstream[j], 20)]
                    for j in range(min(5, len(upstream)))
                    if AA_INDEX.get(upstream[j], 20) < 20
                )
                mp_score = max(mp_score, 1.0 + context_score)
            else:
                mp_score = max(mp_score, 1.0)
    features["pwm_metalloprotease"] = mp_score

    # Cysteine protease: C-H-D/N catalytic triad (C...H with 50-200 spacing, then D/N)
    cp_score = 0.0
    for i, aa in enumerate(sequence):
        if aa == "C":
            for j in range(i + 50, min(i + 200, n)):
                if sequence[j] == "H":
                    for k in range(j + 5, min(j + 50, n)):
                        if sequence[k] in "DN":
                            cp_score = 1.0
                            break
                    if cp_score > 0:
                        break
        if cp_score > 0:
            break
    features["pwm_cys_protease"] = cp_score

    # Phospholipase: GDSL/SGNH hydrolase motif
    pla_score = 0.0
    if re.search(r"G[DE]S[LM]", sequence):
        pla_score += 1.0
    if re.search(r"[SG]G[NH].{20,80}[DE]", sequence):
        pla_score += 0.5
    features["pwm_phospholipase"] = pla_score

    # AB toxin: beta-trefoil lectin fold (QXW repeats)
    ab_score = len(re.findall(r"Q.W", sequence)) / max(1, n) * 100
    features["pwm_ab_toxin"] = ab_score

    # Superantigen: MHC-II groove composition (enriched N,Y,Q in exposed loops)
    # Approximation via loop-forming + aromatic residue density
    super_aa = sum(1 for aa in sequence if aa in "NYQFW") / n
    features["pwm_superantigen"] = super_aa * 5 if super_aa > 0.15 else 0.0

    # RTX calcium-binding: nonapeptide repeats via spectral analysis
    if n > 50:
        hydro = np.array([KD_HYDROPHOBICITY.get(aa, 0.0) for aa in sequence])
        fft_vals = np.abs(np.fft.rfft(hydro - hydro.mean()))
        freqs = np.fft.rfftfreq(n)
        # RTX nonapeptide: period ~9 residues → freq 0.10-0.12
        rtx_band = np.where((freqs > 0.09) & (freqs < 0.13))[0]
        features["pwm_rtx_repeat"] = float(np.max(fft_vals[rtx_band])) / n if len(rtx_band) > 0 else 0.0
    else:
        features["pwm_rtx_repeat"] = 0.0

    # DNase/RNase: conserved His-Asp pair in catalytic context
    dnase_score = len(re.findall(r"H.{5,15}D.{5,15}H", sequence)) / max(1, n) * 100
    features["pwm_dnase"] = min(dnase_score, 2.0)

    # Glycosyltransferase: DXD motif in proper context (flanked by hydrophobic)
    gt_score = 0.0
    for m in re.finditer(r"D.D", sequence):
        pos = m.start()
        # Check flanking (should be in structured context)
        if pos > 3 and pos < n - 3:
            flank = sequence[pos-3:pos] + sequence[pos+3:pos+6]
            hydro_frac = sum(1 for aa in flank if aa in "LIVMFW") / len(flank)
            if hydro_frac > 0.4:
                gt_score += 1.0
    features["pwm_glycosyltransferase"] = min(gt_score, 3.0)

    # --- AMR mechanism scores (5) ---

    # Serine beta-lactamase: SDN + KTG with proper spacing
    sbl_score = 0.0
    sdn_positions = [m.start() for m in re.finditer(r"SDN", sequence)]
    ktg_positions = [m.start() for m in re.finditer(r"K[TS]G", sequence)]
    for sdn_pos in sdn_positions:
        for ktg_pos in ktg_positions:
            spacing = ktg_pos - sdn_pos
            if BLACTAMASE_SERINE["spacing"][0] <= spacing <= BLACTAMASE_SERINE["spacing"][1]:
                sbl_score = 2.0
                break
    features["pwm_blactamase_ser"] = sbl_score

    # Metallo beta-lactamase: HXHXD pattern
    features["pwm_blactamase_met"] = _score_pwm(sequence, BLACTAMASE_METALLO["hxhxd_profile"])

    # Aminoglycoside modifier: GNAT fold composition (enriched in G,N,A with specific topology)
    gnat_aa = sum(1 for aa in sequence if aa in "GNA") / n
    has_gxxg = 1.0 if re.search(r"G..G", sequence) else 0.0
    features["pwm_aminoglyc"] = (gnat_aa * 3 + has_gxxg) if gnat_aa > 0.12 else 0.0

    # Efflux pump: 12+ TM helix topology (many hydrophobic stretches)
    hydrophobic = set("AILMFWV")
    tm_count = 0
    in_tm = False
    hydro_run = 0
    for aa in sequence:
        if aa in hydrophobic:
            hydro_run += 1
            if hydro_run >= 18 and not in_tm:
                tm_count += 1
                in_tm = True
        else:
            if hydro_run < 10:
                in_tm = False
            hydro_run = 0
    features["pwm_efflux"] = min(tm_count / 12.0, 1.5) if tm_count >= 8 else 0.0

    # Ribosomal protection: GTPase domain (DXXG + NKXD motif pair)
    rp_score = 0.0
    if re.search(r"D..G", sequence) and re.search(r"NK.D", sequence):
        rp_score = 1.0
    features["pwm_ribosomal_prot"] = rp_score

    # --- Adhesion/invasion scores (5) ---

    # Integrin-binding: RGD in exposed (disordered) context
    rgd_score = 0.0
    for m in re.finditer(r"RGD", sequence):
        pos = m.start()
        # Check if in disordered region (flanking disorder)
        flank = sequence[max(0, pos-5): pos+8]
        if _disorder_score(flank) > 0.1:
            rgd_score += 1.0
    features["pwm_integrin_rgd"] = min(rgd_score, 2.0)

    # Fibronectin-binding: tandem repeats (period 35-45)
    if n > 100:
        fft_vals = np.abs(np.fft.rfft(hydro_signal := np.array(
            [KD_HYDROPHOBICITY.get(aa, 0.0) for aa in sequence]) - np.mean(
            [KD_HYDROPHOBICITY.get(aa, 0.0) for aa in sequence])))
        freqs = np.fft.rfftfreq(n)
        fn_band = np.where((freqs > 0.022) & (freqs < 0.030))[0]  # period 33-45
        features["pwm_fibronectin"] = float(np.max(fft_vals[fn_band])) / n if len(fn_band) > 0 else 0.0
    else:
        features["pwm_fibronectin"] = 0.0

    # Collagen-binding: GXX repeats (Cna-type)
    gxx_count = len(re.findall(r"G..", sequence))
    # Normalize — collagen-binding proteins have G at every 3rd position
    gxx_expected = n / 3 * BG_FREQ[AA_INDEX["G"]]
    features["pwm_collagen_bind"] = max(0, (gxx_count - gxx_expected) / max(1, n) * 10)

    # Invasin: beta-barrel indicator (alternating polar/hydrophobic pattern)
    alt_score = 0.0
    for i in range(n - 10):
        window = sequence[i:i+10]
        alternating = sum(1 for j in range(9) if
                         (window[j] in "AILMFWV") != (window[j+1] in "AILMFWV"))
        if alternating >= 7:
            alt_score += 1
    features["pwm_invasin"] = alt_score / max(1, n - 10) * 10

    # Biofilm/amyloid: aggregation propensity (Waltz-like)
    if n > 6:
        agg_scores = []
        for i in range(n - 6):
            window = sequence[i:i+6]
            score = sum(AGGREGATION_PROPENSITY.get(aa, 0.0) for aa in window) / 6
            agg_scores.append(score)
        features["pwm_biofilm_amyloid"] = max(agg_scores) if agg_scores else 0.0
    else:
        features["pwm_biofilm_amyloid"] = 0.0

    return features


# ============================================================================
# D. SEQUENCE INTELLIGENCE FEATURES (10 features)
# ============================================================================

def extract_sequence_intelligence(sequence: str) -> Dict[str, float]:
    """
    Spectral and statistical features capturing higher-order sequence patterns.

    Novel: FFT-based repeat detection and disorder transition counting for
    pathogen classification.
    """
    features = {}
    n = len(sequence)

    if n < 30:
        return {
            "intel_fft_peak1": 0.0, "intel_fft_peak2": 0.0, "intel_fft_peak3": 0.0,
            "intel_charged_cluster": 0.0, "intel_amphipathic_max": 0.0,
            "intel_disorder_transitions": 0.0, "intel_surface_access": 0.0,
            "intel_repeat_density": 0.0, "intel_hydro_asymmetry": 0.0,
            "intel_complexity_gradient": 0.0,
        }

    # Top 3 FFT peaks in hydrophobicity signal
    hydro = np.array([KD_HYDROPHOBICITY.get(aa, 0.0) for aa in sequence])
    fft_vals = np.abs(np.fft.rfft(hydro - hydro.mean()))
    # Skip DC component and very low frequencies
    fft_vals[:3] = 0
    freqs = np.fft.rfftfreq(n)

    # Get top 3 peaks (by magnitude, normalized by length)
    top_indices = np.argsort(fft_vals)[-3:][::-1]
    for i, idx in enumerate(top_indices):
        features[f"intel_fft_peak{i+1}"] = float(fft_vals[idx]) / n

    # Charged cluster score: longest same-charge run / length
    max_run = 0
    current_run = 0
    current_type = None
    for aa in sequence:
        if aa in "RKH":
            if current_type == "pos":
                current_run += 1
            else:
                current_run = 1
                current_type = "pos"
        elif aa in "DE":
            if current_type == "neg":
                current_run += 1
            else:
                current_run = 1
                current_type = "neg"
        else:
            max_run = max(max_run, current_run)
            current_run = 0
            current_type = None
    max_run = max(max_run, current_run)
    features["intel_charged_cluster"] = max_run / n

    # Max amphipathic helix potential (18-residue window)
    features["intel_amphipathic_max"] = _hydrophobic_moment(sequence, window=18, angle=100.0)

    # Intrinsic disorder transitions (order↔disorder boundaries)
    # Sliding window disorder score, count zero-crossings
    window = 15
    transitions = 0
    prev_disordered = None
    for i in range(0, n - window, 5):  # Step by 5 for speed
        d_score = _disorder_score(sequence[i:i+window])
        is_disordered = d_score > 0.1
        if prev_disordered is not None and is_disordered != prev_disordered:
            transitions += 1
        prev_disordered = is_disordered
    features["intel_disorder_transitions"] = transitions / max(1, n / 50)

    # Surface accessibility proxy (fraction of polar/charged residues in
    # turn-prone contexts)
    turn_prone = sum(1 for i in range(n - 3) if
                     TURN_PROPENSITY.get(sequence[i], 1.0) > 1.2 and
                     sequence[i] in "STNQDEKR")
    features["intel_surface_access"] = turn_prone / max(1, n)

    # Repeat density (exact 3-8mer repeats)
    repeat_score = 0
    for k in range(3, 9):
        kmers = Counter(sequence[i:i+k] for i in range(n - k + 1))
        repeat_score += sum(c - 1 for c in kmers.values() if c > 1)
    features["intel_repeat_density"] = repeat_score / max(1, n)

    # Hydrophobicity asymmetry (difference between N-term and C-term halves)
    mid = n // 2
    h1 = np.mean([KD_HYDROPHOBICITY.get(aa, 0.0) for aa in sequence[:mid]])
    h2 = np.mean([KD_HYDROPHOBICITY.get(aa, 0.0) for aa in sequence[mid:]])
    features["intel_hydro_asymmetry"] = abs(h1 - h2)

    # Complexity gradient (is one end more complex than the other?)
    e1 = _sequence_entropy(sequence[:mid])
    e2 = _sequence_entropy(sequence[mid:])
    features["intel_complexity_gradient"] = abs(e1 - e2)

    return features


# ============================================================================
# E. PHYSICOCHEMICAL FEATURES (10 features) — kept from v3
# ============================================================================

def extract_physicochemical(sequence: str) -> Dict[str, float]:
    """Category E: Global physicochemical properties (10 features)."""
    defaults = {
        "physchem_mw": 0.0, "physchem_pi": 0.0, "physchem_gravy": 0.0,
        "physchem_instability": 0.0, "physchem_aromaticity": 0.0,
        "physchem_pos_charge": 0.0, "physchem_neg_charge": 0.0,
        "physchem_net_charge": 0.0, "physchem_hydrophob_mean": 0.0,
        "physchem_hydrophob_std": 0.0,
    }

    clean = "".join(c for c in sequence if c in "ACDEFGHIKLMNPQRSTVWY")
    if len(clean) < 10:
        return defaults

    try:
        analysis = ProteinAnalysis(clean)
        aa_pct = analysis.get_amino_acids_percent()

        hydrophob = [KD_HYDROPHOBICITY.get(aa, 0.0) for aa in clean]
        pos = aa_pct.get("R", 0) + aa_pct.get("K", 0) + aa_pct.get("H", 0)
        neg = aa_pct.get("D", 0) + aa_pct.get("E", 0)

        return {
            "physchem_mw": analysis.molecular_weight() / 1000,
            "physchem_pi": analysis.isoelectric_point(),
            "physchem_gravy": analysis.gravy(),
            "physchem_instability": analysis.instability_index(),
            "physchem_aromaticity": analysis.aromaticity(),
            "physchem_pos_charge": pos,
            "physchem_neg_charge": neg,
            "physchem_net_charge": pos - neg,
            "physchem_hydrophob_mean": float(np.mean(hydrophob)),
            "physchem_hydrophob_std": float(np.std(hydrophob)),
        }
    except Exception:
        return defaults


# ============================================================================
# F. QUALITY & NEGATIVE CONTROL (5 features) — kept from v3
# ============================================================================

def extract_quality_features(sequence: str) -> Dict[str, float]:
    """Category F: Quality metrics and housekeeping signatures."""
    n = len(sequence)

    # Sequence entropy (complexity)
    entropy = _sequence_entropy(sequence)

    # Completeness indicator (starts with M, reasonable length)
    completeness = 1.0 if sequence.startswith("M") and n > 100 else 0.5

    # Length bin (log scale)
    length_bin = math.log2(max(50, n)) - math.log2(50)

    # Housekeeping protein signatures (negative control)
    housekeeping_score = 0.0

    # Ribosomal protein: small, basic, Lys/Arg rich
    if n < 300 and sum(1 for aa in sequence if aa in "KR") / n > 0.15:
        housekeeping_score += 0.3

    # DNA/RNA polymerase: large, contains DXXLYP-like, highly conserved
    if n > 500 and re.search(r"D..D", sequence) and re.search(r"Y.D.F", sequence):
        housekeeping_score += 0.3

    # ATP-binding (P-loop): GX{4}GK[TS]
    if re.search(r"G.{4}GK[TS]", sequence):
        housekeeping_score += 0.2

    # NAD-binding (Rossmann fold): GXGXXG
    if re.search(r"G.G..G", sequence):
        housekeeping_score += 0.2

    return {
        "quality_entropy": entropy,
        "quality_completeness": completeness,
        "quality_length_bin": length_bin,
        "quality_housekeeping": housekeeping_score,
        "quality_ribosomal": 1.0 if housekeeping_score >= 0.5 else 0.0,
    }


# ============================================================================
# MAIN EXTRACTOR CLASS
# ============================================================================

class ProteinFeatureExtractor:
    """
    MetaQuest Pathogenic Functional Fingerprint (PFF) v4.0

    Extracts ~95 biologically-informed features for pathogen classification.
    No external tool dependencies at inference time.
    """

    def __init__(self):
        self.feature_categories = [
            ("positional", extract_positional_features),
            ("mimicry", extract_mimicry_features),
            ("virulence_pwm", extract_virulence_pwm_features),
            ("intelligence", extract_sequence_intelligence),
            ("physicochemical", extract_physicochemical),
            ("quality", extract_quality_features),
        ]

    def extract_single(self, sequence: str) -> Dict[str, float]:
        """Extract all features for a single protein sequence."""
        clean = "".join(c for c in sequence.upper() if c in "ACDEFGHIKLMNPQRSTVWY")
        if len(clean) < 20:
            # Return zeros for all features
            features = {}
            for _, func in self.feature_categories:
                features.update(func("A" * 20))  # Get key names
            return {k: 0.0 for k in features}

        features = {}
        for _, func in self.feature_categories:
            features.update(func(clean))

        return features

    def extract_batch(self, sequences, desc: str = "Extracting features", show_progress: bool = True):
        """
        Extract features for a batch of sequences.

        Args:
            sequences: List of strings OR list of (accession, sequence) tuples.

        Returns:
            pandas DataFrame with features (and 'accession' column if tuples provided).
        """
        import pandas as pd

        results = []
        iterator = tqdm(sequences, desc=desc) if show_progress else sequences

        for item in iterator:
            if isinstance(item, tuple):
                acc, seq = item
                feats = self.extract_single(seq)
                feats["accession"] = acc
            else:
                feats = self.extract_single(item)
            results.append(feats)

        return pd.DataFrame(results)

    def feature_names(self) -> List[str]:
        """Return ordered list of all feature names."""
        # Extract from a dummy sequence to get all keys
        dummy = "M" + "A" * 200 + "K" * 50
        features = self.extract_single(dummy)
        return sorted(features.keys())


# Backward compatibility alias
MetaQuestProteinFeatureExtractor = ProteinFeatureExtractor
