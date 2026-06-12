#!/usr/bin/env python3
"""
MetaQuest SHAP Explainability Module

Generates per-prediction explanations showing WHY a protein
was classified as pathogenic or non-pathogenic.

Outputs:
- Text table of top contributing features
- HTML waterfall chart for dashboards
- Structured dict for JSON reports
"""

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)

# Human-readable feature descriptions
FEATURE_DESCRIPTIONS = {
    # --- Positional Decomposition ---
    "nterm_disorder": "N-terminal disorder (secretion signal)",
    "nterm_charge_pos": "N-terminal positive charge",
    "nterm_charge_neg": "N-terminal negative charge",
    "nterm_charge_net": "N-terminal net charge",
    "nterm_amphipathic": "N-terminal amphipathicity",
    "nterm_polar_ratio": "N-terminal polar/nonpolar ratio",
    "nterm_entropy": "N-terminal sequence complexity",
    "nterm_enrichment1": "N-terminal Ser/Thr (T3SS signal)",
    "nterm_enrichment2": "N-terminal twin-Arg (Tat signal)",
    "nterm_enrichment3": "N-terminal Pro (disorder signal)",
    "core_disorder": "Catalytic core disorder",
    "core_charge_pos": "Core positive charge",
    "core_charge_neg": "Core negative charge",
    "core_charge_net": "Core net charge",
    "core_amphipathic": "Core amphipathicity",
    "core_polar_ratio": "Core polar/nonpolar ratio",
    "core_entropy": "Core sequence complexity",
    "core_enrichment1": "Disulfide bond potential",
    "core_enrichment2": "Metal coordination residues",
    "core_enrichment3": "Catalytic site density",
    "cterm_disorder": "C-terminal disorder",
    "cterm_charge_pos": "C-terminal positive charge",
    "cterm_charge_neg": "C-terminal negative charge",
    "cterm_charge_net": "C-terminal net charge",
    "cterm_amphipathic": "C-terminal amphipathicity",
    "cterm_polar_ratio": "C-terminal polar/nonpolar ratio",
    "cterm_entropy": "C-terminal complexity",
    "cterm_enrichment1": "LPXTG sortase signal",
    "cterm_enrichment2": "GPI/membrane anchor",
    "cterm_enrichment3": "Charged tail (secreted toxin)",
    # --- Host-Mimicry & Evolutionary ---
    "mimicry_elm_density": "Eukaryotic-like motif density",
    "mimicry_elm_total": "Total host-mimicry motifs",
    "mimicry_nls_signal": "Nuclear localization signal",
    "mimicry_sumo_motif": "SUMO interaction motif",
    "mimicry_death_domain": "Death domain interaction",
    "mimicry_14_3_3": "14-3-3 binding motif",
    "mimicry_cyclin_binding": "Cyclin-binding (cell cycle hijack)",
    "mimicry_ankyrin_repeat": "Ankyrin repeat (host mimicry)",
    "mimicry_lrr_periodicity": "Leucine-rich repeat periodicity",
    "mimicry_coiled_coil": "Coiled-coil heptad score",
    "mimicry_prion_like": "Prion-like low-complexity region",
    "mimicry_acidic_domain": "Acidic activation domain mimic",
    "mimicry_gc_deviation": "GC-content deviation (HGT signal)",
    "mimicry_compositional_jsd": "Compositional heterogeneity (multi-domain)",
    "mimicry_rare_dipeptide": "Effector-enriched dipeptides",
    # --- Virulence PWM Scores ---
    "pwm_t3ss_nterm": "Type III secretion signal (PWM)",
    "pwm_t4ss_cterm": "Type IV secretion signal (PWM)",
    "pwm_t6ss_spike": "Type VI secretion spike",
    "pwm_sec_signal": "Sec signal peptide (von Heijne)",
    "pwm_tat_signal": "Tat twin-arginine signal (PWM)",
    "pwm_adprt": "ADP-ribosyltransferase",
    "pwm_pore_forming": "Pore-forming toxin (amphipathic)",
    "pwm_metalloprotease": "Metalloprotease (HEXXH+context)",
    "pwm_cys_protease": "Cysteine protease catalytic triad",
    "pwm_phospholipase": "Phospholipase (GDSL hydrolase)",
    "pwm_ab_toxin": "AB toxin lectin domain",
    "pwm_superantigen": "Superantigen MHC-binding",
    "pwm_rtx_repeat": "RTX calcium-binding repeats",
    "pwm_dnase": "DNase/RNase active site",
    "pwm_glycosyltransferase": "Glycosyltransferase (DXD)",
    "pwm_blactamase_ser": "Serine beta-lactamase (SDN+KTG)",
    "pwm_blactamase_met": "Metallo-beta-lactamase (HXHXD)",
    "pwm_aminoglyc": "Aminoglycoside modifier (GNAT)",
    "pwm_efflux": "Efflux pump (multi-TM topology)",
    "pwm_ribosomal_prot": "Ribosomal protection (GTPase)",
    "pwm_integrin_rgd": "Integrin-binding RGD (exposed)",
    "pwm_fibronectin": "Fibronectin-binding repeats",
    "pwm_collagen_bind": "Collagen-binding domain",
    "pwm_invasin": "Invasin beta-barrel pattern",
    "pwm_biofilm_amyloid": "Biofilm amyloid propensity",
    # --- Sequence Intelligence ---
    "intel_fft_peak1": "Spectral repeat (dominant period)",
    "intel_fft_peak2": "Spectral repeat (2nd period)",
    "intel_fft_peak3": "Spectral repeat (3rd period)",
    "intel_charged_cluster": "Charged residue clustering",
    "intel_amphipathic_max": "Max amphipathic helix potential",
    "intel_disorder_transitions": "Disorder-order transitions",
    "intel_surface_access": "Surface accessibility proxy",
    "intel_repeat_density": "Internal repeat density",
    "intel_hydro_asymmetry": "Hydrophobicity asymmetry (N/C)",
    "intel_complexity_gradient": "Complexity gradient (N/C)",
    # --- Physicochemical ---
    "physchem_mw": "Molecular weight",
    "physchem_pi": "Isoelectric point",
    "physchem_gravy": "Hydrophobicity (GRAVY)",
    "physchem_instability": "Instability index",
    "physchem_aromaticity": "Aromaticity",
    "physchem_pos_charge": "Positive charge density",
    "physchem_neg_charge": "Negative charge density",
    "physchem_net_charge": "Net charge",
    "physchem_hydrophob_mean": "Mean hydrophobicity",
    "physchem_hydrophob_std": "Hydrophobicity variation",
    # --- Quality & Negative Control ---
    "quality_entropy": "Sequence complexity",
    "quality_completeness": "Sequence completeness",
    "quality_length_bin": "Protein length",
    "quality_housekeeping": "Housekeeping protein score",
    "quality_ribosomal": "Ribosomal protein signature",
}


def explain_prediction(
    shap_values: np.ndarray,
    feature_names: List[str],
    feature_values: np.ndarray,
    top_n: int = 5,
) -> Dict:
    """
    Generate a human-readable explanation for a single prediction.

    Args:
        shap_values: SHAP values for this prediction (1D array).
        feature_names: List of feature names.
        feature_values: Actual feature values for this sample.
        top_n: Number of top features to include.

    Returns:
        Dict with explanation data.
    """
    # Pair features with their SHAP values
    contributions = list(zip(feature_names, shap_values, feature_values))

    # Sort by absolute SHAP value
    contributions.sort(key=lambda x: abs(x[1]), reverse=True)

    # Top positive contributors (pushing toward pathogenic)
    positive = [
        {
            "feature": name,
            "description": FEATURE_DESCRIPTIONS.get(name, name),
            "shap_value": float(sv),
            "feature_value": float(fv),
            "direction": "pathogenic",
        }
        for name, sv, fv in contributions if sv > 0
    ][:top_n]

    # Top negative contributors (pushing toward non-pathogenic)
    negative = [
        {
            "feature": name,
            "description": FEATURE_DESCRIPTIONS.get(name, name),
            "shap_value": float(sv),
            "feature_value": float(fv),
            "direction": "non-pathogenic",
        }
        for name, sv, fv in contributions if sv < 0
    ][:top_n]

    # Generate summary text
    summary_parts = []
    if positive:
        top_pos = positive[0]
        summary_parts.append(
            f"Strongest pathogenic signal: {top_pos['description']} "
            f"(value={top_pos['feature_value']:.3f})"
        )
    if negative:
        top_neg = negative[0]
        summary_parts.append(
            f"Strongest safe signal: {top_neg['description']} "
            f"(value={top_neg['feature_value']:.3f})"
        )

    return {
        "top_pathogenic_features": positive,
        "top_safe_features": negative,
        "summary": "; ".join(summary_parts),
        "total_shap_positive": float(sum(sv for _, sv, _ in contributions if sv > 0)),
        "total_shap_negative": float(sum(sv for _, sv, _ in contributions if sv < 0)),
    }


def render_shap_text(explanation: Dict, protein_id: str = "") -> str:
    """
    Render SHAP explanation as a text table for CLI/txt reports.

    Returns:
        Formatted text string.
    """
    lines = []
    header = f"SHAP Explanation: {protein_id}" if protein_id else "SHAP Explanation"
    lines.append(f"  {header}")
    lines.append(f"  {'─' * 60}")

    # Pathogenic drivers
    if explanation["top_pathogenic_features"]:
        lines.append("  Pathogenic signals:")
        for feat in explanation["top_pathogenic_features"]:
            bar = "+" * min(int(abs(feat["shap_value"]) * 20), 20)
            lines.append(
                f"    {bar:20s} {feat['description'][:35]:35s} "
                f"({feat['feature_value']:.3f})"
            )

    # Safe signals
    if explanation["top_safe_features"]:
        lines.append("  Non-pathogenic signals:")
        for feat in explanation["top_safe_features"]:
            bar = "-" * min(int(abs(feat["shap_value"]) * 20), 20)
            lines.append(
                f"    {bar:20s} {feat['description'][:35]:35s} "
                f"({feat['feature_value']:.3f})"
            )

    lines.append(f"  {'─' * 60}")
    return "\n".join(lines)


def render_shap_html(
    explanation: Dict, protein_id: str = "", confidence: float = 0.0
) -> str:
    """
    Render SHAP explanation as inline HTML for the dashboard.

    Returns:
        HTML string with a waterfall-style visualization.
    """
    all_features = (
        explanation["top_pathogenic_features"]
        + explanation["top_safe_features"]
    )
    # Sort by absolute SHAP for display
    all_features.sort(key=lambda x: abs(x["shap_value"]), reverse=True)
    all_features = all_features[:8]  # Limit to 8

    # Find max for scaling bars
    max_shap = max((abs(f["shap_value"]) for f in all_features), default=1.0)

    rows_html = ""
    for feat in all_features:
        pct = abs(feat["shap_value"]) / max_shap * 100
        color = "#e74c3c" if feat["shap_value"] > 0 else "#27ae60"
        sign = "+" if feat["shap_value"] > 0 else ""
        label = feat["description"][:40]

        rows_html += f"""
        <div style="display:flex;align-items:center;margin:2px 0;font-size:12px;">
          <span style="width:200px;text-overflow:ellipsis;overflow:hidden;white-space:nowrap;">{label}</span>
          <div style="flex:1;height:16px;background:#f0f0f0;margin:0 8px;position:relative;">
            <div style="width:{pct:.0f}%;height:100%;background:{color};"></div>
          </div>
          <span style="width:60px;text-align:right;color:{color};">{sign}{feat['shap_value']:.3f}</span>
        </div>"""

    conf_color = "#e74c3c" if confidence > 0.7 else "#f39c12" if confidence > 0.5 else "#27ae60"

    html = f"""
    <div style="border:1px solid #ddd;border-radius:4px;padding:12px;margin:8px 0;font-family:monospace;">
      <div style="font-weight:bold;margin-bottom:8px;">
        {protein_id} — Confidence: <span style="color:{conf_color};">{confidence:.1%}</span>
      </div>
      {rows_html}
      <div style="font-size:11px;color:#666;margin-top:8px;">
        Red = pathogenic signal, Green = safe signal
      </div>
    </div>"""

    return html
