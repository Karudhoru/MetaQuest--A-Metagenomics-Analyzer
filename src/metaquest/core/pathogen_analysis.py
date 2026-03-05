"""
Pathogen Analysis Module - Cleaned Version
==========================================
Fragment-aware pathogen detection with minimal output
"""

import subprocess
import math
import json
import pandas as pd
import re
from pathlib import Path
from typing import List, Dict, Tuple, Optional
from collections import defaultdict
from Bio import SeqIO

from ..config import PATHOGEN_DB_CUSTOM, CRITICAL_MOTIFS
from ..io.data_loaders import load_protein_sequences_streaming
from ..io.text_parsers import extract_organism_name
from ..io.output_formatter import get_formatter


def build_abundance_lookup(taxonomy_df: pd.DataFrame) -> Dict[str, float]:
    """
    Build fast O(1) lookup dictionary for organism abundance.
    
    Returns:
        Dict mapping genus/species names to abundance values
    """
    if taxonomy_df is None or taxonomy_df.empty:
        return {}
    
    abundance_lookup = {}
    
    for _, row in taxonomy_df.iterrows():
        species = row['name'].lower()
        abundance = row.get('fraction_total_reads', 0.0)
        
        # Add full species name
        abundance_lookup[species] = abundance
        
        # Also add genus name (for partial matches)
        genus = species.split()[0] if ' ' in species else species
        # Keep max abundance if genus appears multiple times
        abundance_lookup[genus] = max(abundance_lookup.get(genus, 0), abundance)
    
    return abundance_lookup

def get_organism_abundance_fast(organism: str, abundance_lookup: Dict[str, float]) -> float:
    """Fast O(1) organism abundance lookup."""
    if not abundance_lookup or organism == "Unknown":
        return 0.0
    
    organism_lower = organism.lower()
    
    # Try exact match first
    if organism_lower in abundance_lookup:
        return abundance_lookup[organism_lower]
    
    # Try genus match
    genus = organism_lower.split()[0] if ' ' in organism_lower else organism_lower
    if genus in abundance_lookup:
        return abundance_lookup[genus]
    
    return 0.00

# ============================================================================
# CONFIDENCE SCORING
# ============================================================================

def calculate_confidence_score(
    hit: Dict,
    abundance_lookup: Dict[str, float],
    contig_coverage: float,
    motifs_found: List[str],
    prokka_annotation: str
) -> Tuple[float, str, List[str]]:
    """
    Multi-factor confidence scoring for fragment hits.
    
    This function is intentionally quiet - just calculates the score.
    """
    evidence = []
    score = 0.0
    
    # FACTOR 1: Fragment Completeness (0-35 points) - CONTINUOUS SCORING
    query_len = hit['qlen']
    reference_len = hit.get('reference_length', hit['slen'])
    completeness = query_len / reference_len

    # Continuous sigmoid scoring (no sharp cliffs)
    # Formula: score = max_points * (1 / (1 + e^(-k*(completeness - 0.5))))
    # This gives smooth scaling from 0-35 points
    k = 8  # Steepness parameter (higher = steeper curve)
    completeness_score = 35 * (1 / (1 + math.exp(-k * (completeness - 0.5))))

    # Adjust for absolute length (short genes need less coverage)
    # FACTOR 7: Signal Peptide Check (SOFT PENALTY)
    # Only penalize if BOTH: very short AND high hydrophobic content
    if query_len < 60 and is_likely_signal_peptide(hit['sequence']):
        # Check if this is a known short gene (toxins, small AMR genes)
        gene_name = hit.get('sseqid', '').lower()
        known_short_genes = ['mec', 'van', 'mcr', 'ctx', 'oxa', 'kpc', 'ndm']
        
        if any(short_gene in gene_name for short_gene in known_short_genes):
            # Known short gene - no penalty
            evidence.append("Short gene (known)")
        else:
            # Likely just signal peptide - soft penalty
            score -= 10  # Reduced from 20
            evidence.append("⚠️ Possible signal peptide fragment")
    elif query_len < 300:  # Medium gene
        length_bonus = 1.0
        evidence.append(f"Medium gene fragment ({query_len}aa, {completeness*100:.0f}% complete)")
    else:  # Long gene (structural proteins)
        # Long genes need higher completeness to be confident
        length_bonus = 0.8 if completeness < 0.7 else 1.0
        evidence.append(f"Large gene fragment ({query_len}aa, {completeness*100:.0f}% complete)")

    score += completeness_score * length_bonus

    # Add descriptive interpretation
    if completeness >= 0.9:
        evidence.append("Near-complete gene")
    elif completeness >= 0.7:
        evidence.append("Substantial coverage")
    elif completeness >= 0.5:
        evidence.append("Moderate coverage")
    else:
        evidence.append("Partial fragment")
    
    # FACTOR 2: Sequence Identity with Coverage Context (0-25 points)
    identity = hit['pident']
    query_coverage = hit.get('query_coverage', (hit['qend'] - hit['qstart'] + 1) / hit['qlen'])

    # CRITICAL: Low identity OR low coverage = not reliable
    if identity < 85 or query_coverage < 0.5:
        # Below 85% identity or <50% coverage = minimal confidence
        identity_score = 0
        evidence.append(f"Low confidence: {identity:.1f}% identity, {query_coverage*100:.0f}% coverage")
    elif identity < 90:
        # 85-90% identity = moderate confidence (5-15 points)
        identity_score = 5 + (identity - 85) * 2  # Linear from 5 to 15
        evidence.append(f"Moderate identity ({identity:.1f}%)")
    elif identity < 95:
        # 90-95% identity = high confidence (15-20 points)
        identity_score = 15 + (identity - 90)  # Linear from 15 to 20
        evidence.append(f"High identity ({identity:.1f}%)")
    else:
        # 95%+ identity = very high confidence (20-25 points)
        identity_score = 20 + min((identity - 95), 5)  # Up to 25 at 100%
        evidence.append(f"Very high identity ({identity:.1f}%)")

    # Penalize short alignments even with high identity
    if hit['length'] < 50:
        identity_score *= 0.5
        evidence.append("⚠️ Short alignment (<50aa)")

    score += identity_score
    
    # FACTOR 3: Organism Abundance (0-20 points) - OPTIMIZED
    organism = extract_organism_name(hit.get('stitle', '')) or "Unknown"
    abundance = get_organism_abundance_fast(organism, abundance_lookup)

    # Continuous scoring based on log10 of abundance
    # Log scale because abundances span many orders of magnitude
    if abundance > 0:
        # Log scale: 0.1% = -3, 1% = -2, 10% = -1, 100% = 0
        log_abundance = math.log10(abundance * 100)  # Convert to percentage first
        
        # Map to 0-20 scale
        # -4 (0.01%) → 0 points
        # -3 (0.1%) → 5 points
        # -2 (1%) → 10 points
        # -1 (10%) → 15 points
        # 0 (100%) → 20 points (max)
        
        abundance_score = max(0, min(20, 5 * (log_abundance + 4)))
        
        if abundance >= 0.05:
            evidence.append(f"Organism abundant ({abundance*100:.1f}%)")
        elif abundance >= 0.01:
            evidence.append(f"Organism present ({abundance*100:.2f}%)")
        elif abundance >= 0.001:
            evidence.append(f"Organism detected ({abundance*100:.3f}%)")
        else:
            evidence.append(f"Organism rare ({abundance*100:.4f}%)")
    else:
        abundance_score = 0
        evidence.append("Organism not detected in taxonomy")

    score += abundance_score
    
    # FACTOR 4: Contig Quality (0-10 points)
    if contig_coverage >= 20:
        score += 10
        evidence.append(f"High coverage contig ({contig_coverage:.0f}x)")
    elif contig_coverage >= 10:
        score += 7
        evidence.append(f"Good coverage contig ({contig_coverage:.0f}x)")
    elif contig_coverage >= 5:
        score += 4
    else:
        score += 1
        evidence.append(f"Low coverage contig ({contig_coverage:.1f}x)")
    
    # FACTOR 5: Critical Motif Detection (0-15 points)
    if motifs_found:
        score += len(motifs_found) * 5
        evidence.append(f"Critical motifs: {', '.join(motifs_found)}")
    
    # FACTOR 6: Prokka Cross-Validation
    if prokka_annotation and prokka_annotation != "hypothetical protein":
        gene_name = hit.get('sseqid', '').split('_')[0].lower()
        
        if gene_name in prokka_annotation.lower():
            score += 10
            evidence.append("Prokka confirms gene identity")
        elif "hypothetical" not in prokka_annotation.lower():
            score -= 15
            evidence.append(f"⚠️ Prokka disagrees: {prokka_annotation[:50]}")
    
    # FACTOR 7: Signal Peptide Check (PENALTY)
    if query_len < 80 and is_likely_signal_peptide(hit['sequence']):
        score -= 20
        evidence.append("⚠️ Likely signal peptide only")
    
    # Determine Confidence Tier
    score = max(0, min(score, 100))
    
    if score >= 70:
        tier = "HIGH"
    elif score >= 45:
        tier = "MEDIUM"
    else:
        tier = "LOW"
    
    return score, tier, evidence


def is_likely_signal_peptide(sequence: str) -> bool:
    """Detect if sequence is likely just a signal peptide"""
    if len(sequence) > 100:
        return False
    
    hydrophobic_count = sum(sequence.count(aa) for aa in 'LIVAFM')
    hydrophobic_ratio = hydrophobic_count / len(sequence)
    has_cleavage_site = bool(re.search(r'[AV][ASTG][AS]', sequence[:30]))
    
    return (hydrophobic_ratio > 0.4 and len(sequence) < 50) or \
           (hydrophobic_ratio > 0.35 and has_cleavage_site and len(sequence) < 80)


def check_critical_motifs(sequence: str, gene_name: str) -> Tuple[float, List[str]]:
    """Fast motif-based validation for critical genes"""
    base_gene = re.split(r'[_\-]', gene_name)[0].upper()
    
    motif_data = None
    for key in CRITICAL_MOTIFS:
        if key.upper() in base_gene or base_gene in key.upper():
            motif_data = CRITICAL_MOTIFS[key]
            break
    
    if not motif_data:
        return 0.0, []
    
    found_motifs = []
    for motif in motif_data['motifs']:
        regex_pattern = motif.replace('X', '[A-Z]')
        if re.search(regex_pattern, sequence):
            found_motifs.append(motif)
    
    if len(found_motifs) >= motif_data['min_motifs']:
        motif_score = 1.0
    elif found_motifs:
        motif_score = 0.5
    else:
        motif_score = 0.0
    
    return motif_score, found_motifs


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def get_contig_coverage(contig_name: str, contigs_file: Path) -> float:
    """Extract coverage from contig header"""
    try:
        with open(contigs_file) as f:
            for line in f:
                if line.startswith('>') and contig_name in line:
                    match = re.search(r'cov_([\d.]+)', line)
                    if match:
                        return float(match.group(1))
    except Exception:
        pass
    return 0.0


def get_prokka_annotation(gene_id: str, prokka_files: List[Path]) -> str:
    """Get Prokka annotation for a gene"""
    for prokka_file in prokka_files:
        if prokka_file.suffix == '.tsv':
            try:
                df = pd.read_csv(prokka_file, sep='\t')
                match = df[df['locus_tag'] == gene_id]
                if not match.empty:
                    return match.iloc[0].get('product', 'hypothetical protein')
            except Exception:
                pass
    return "hypothetical protein"

def export_legacy_tsv(results_dict: Dict, output_file: Path) -> Path:
    """
    Export new JSON format to legacy TSV format for backward compatibility.
    
    This ensures the risk scoring module can still read the results.
    """
    # Combine HIGH and MEDIUM confidence hits
    all_hits = results_dict['high_confidence_hits'] + results_dict['medium_confidence_hits']
    
    if not all_hits:
        # Create empty file
        pd.DataFrame(columns=[
            'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
            'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
            'qlen', 'slen', 'stitle'
        ]).to_csv(output_file, sep='\t', index=False, header=False)
        return output_file
    
    # Convert to flat dataframe format
    rows = []
    for hit in all_hits:
        rows.append({
            'qseqid': hit['qseqid'],
            'sseqid': hit['sseqid'],
            'pident': hit['pident'],
            'length': hit['length'],
            'mismatch': 0,
            'gapopen': 0,
            'qstart': hit['qstart'],
            'qend': hit['qend'],
            'sstart': hit['sstart'],
            'send': hit['send'],
            'evalue': hit['evalue'],
            'bitscore': hit['bitscore'],
            'qlen': hit['qlen'],
            'slen': hit['slen'],
            'stitle': hit['stitle']
        })
    
    df = pd.DataFrame(rows)
    df.to_csv(output_file, sep='\t', index=False, header=False)
    
    return output_file


# ============================================================================
# MAIN PATHOGEN SCANNING FUNCTION
# ============================================================================

def run_pathogen_scan_v2(
    protein_file: Path,
    output_dir: Path,
    contigs_file: Path,
    prokka_dir: Path,
    bracken_results: Optional[Path] = None,
    min_fragment_length: int = 80,
    min_identity: float = 80.0,
    min_query_coverage: float = 0.7,
    min_subject_coverage: float = 0.7
) -> Path:
    """
    Enhanced pathogen scan with fragment-aware validation.
    
    Output is controlled by caller - this function only logs at DEBUG level.
    """
    formatter = get_formatter()
    
    formatter.debug("Loading taxonomy and protein data...")
    
    # Load taxonomy data
    taxonomy_df = None
    abundance_lookup = {}
    if bracken_results and bracken_results.exists():
        try:
            taxonomy_df = pd.read_csv(bracken_results, sep='\t')
            abundance_lookup = build_abundance_lookup(taxonomy_df)
            formatter.debug(f"Loaded {len(taxonomy_df)} species from Bracken")
        except Exception as e:
            formatter.debug(f"Could not load Bracken: {e}")
    
    # Get Prokka annotation files
    prokka_files = list(prokka_dir.glob("*.tsv")) + list(prokka_dir.glob("*.gff"))
    
    # Run DIAMOND
    selected_db = PATHOGEN_DB_CUSTOM
    blast_out = output_dir / "pathogen_results_raw.txt"
    
    if not selected_db or not selected_db.exists():
        formatter.error(f"Pathogen database not found: {selected_db}")
        return None
    
    db_base = selected_db.with_suffix('') if selected_db.suffix == '.dmnd' else selected_db
    
    outfmt_cols = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle"
    
    cmd = [
        "diamond", "blastp",
        "-d", str(db_base),
        "-q", str(protein_file),
        "-o", str(blast_out),
        "--outfmt", "6", *outfmt_cols.split(),
        "--top", "1",
        "--evalue", "1e-5",
        "--id", str(min_identity),
        "--threads", "4",
        "--more-sensitive"
    ]
    
    formatter.debug(f"DIAMOND command: {' '.join(cmd)}")
    
    # Run DIAMOND - suppress output
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        formatter.error(f"DIAMOND failed: {result.stderr}")
        return None
    
    if not blast_out.exists() or blast_out.stat().st_size == 0:
        formatter.debug("No pathogen hits found")
        
        # Create empty results
        empty_results = {
            'summary': {
                'total_hits': 0,
                'high_confidence': 0,
                'medium_confidence': 0,
                'low_confidence': 0,
                'filters_applied': {}
            },
            'high_confidence_hits': [],
            'medium_confidence_hits': [],
            'low_confidence_hits': []
        }
        
        filtered_file = output_dir / "pathogen_detections_validated.json"
        with open(filtered_file, 'w') as f:
            json.dump(empty_results, f, indent=2)
        
        legacy_file = output_dir / "pathogen_results.tsv"
        export_legacy_tsv(empty_results, legacy_file)
        
        return legacy_file
    
    # STEP 1: Collect needed protein IDs from BLAST hits (memory efficient)
    formatter.debug("Collecting protein IDs from pathogen hits...")
    needed_ids = set()
    
    with open(blast_out) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 1:
                needed_ids.add(parts[0])  # query ID
    
    # STEP 2: Load ONLY the sequences we need (99% memory reduction)
    formatter.debug(f"Loading {len(needed_ids)} protein sequences (streaming mode)...")
    protein_sequences = load_protein_sequences_streaming(protein_file, needed_ids)
    formatter.debug(f"Loaded {len(protein_sequences)} sequences")
    
    # STEP 3: Parse and filter results
    formatter.debug("Processing pathogen hits with confidence scoring...")
    
    hits_high_conf = []
    hits_medium_conf = []
    hits_low_conf = []
    
    with open(blast_out) as f:
        for line_num, line in enumerate(f, 1):
            if not line.strip():
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 15:
                continue
            
            try:
                hit = {
                    'qseqid': parts[0],
                    'sseqid': parts[1],
                    'pident': float(parts[2]),
                    'length': int(parts[3]),
                    'qstart': int(parts[6]),
                    'qend': int(parts[7]),
                    'sstart': int(parts[8]),
                    'send': int(parts[9]),
                    'evalue': float(parts[10]),
                    'bitscore': float(parts[11]),
                    'qlen': int(parts[12]),
                    'slen': int(parts[13]),
                    'stitle': parts[14] if len(parts) > 14 else ""
                }
                
                # Calculate coverages
                hit['query_coverage'] = (hit['qend'] - hit['qstart'] + 1) / hit['qlen']
                hit['subject_coverage'] = abs(hit['send'] - hit['sstart'] + 1) / hit['slen']
                
                # Filter
                if hit['qlen'] < min_fragment_length:
                    continue
                if hit['query_coverage'] < min_query_coverage:
                    continue
                if hit['subject_coverage'] < min_subject_coverage:
                    continue
                
                # Get additional context
                hit['sequence'] = protein_sequences.get(hit['qseqid'], '')
                hit['contig_id'] = hit['qseqid'].rsplit('_', 1)[0]
                hit['contig_coverage'] = get_contig_coverage(hit['contig_id'], contigs_file)
                hit['prokka_annotation'] = get_prokka_annotation(hit['qseqid'], prokka_files)
                
                # Check for critical motifs
                motif_score, motifs_found = check_critical_motifs(hit['sequence'], hit['sseqid'])
                hit['motifs_found'] = motifs_found
                hit['has_critical_motif'] = motif_score > 0.5
                
                # Calculate confidence score
                confidence_score, confidence_tier, evidence = calculate_confidence_score(
                    hit, abundance_lookup, hit['contig_coverage'], motifs_found, hit['prokka_annotation']
                )
                
                hit['confidence_score'] = confidence_score
                hit['confidence_tier'] = confidence_tier
                hit['evidence'] = evidence
                
                # Categorize by confidence
                if confidence_tier == "HIGH":
                    hits_high_conf.append(hit)
                elif confidence_tier == "MEDIUM":
                    hits_medium_conf.append(hit)
                else:
                    hits_low_conf.append(hit)
                
            except (ValueError, IndexError) as e:
                formatter.debug(f"Parse error line {line_num}: {e}")
                continue
    
    # Summary
    total_hits = len(hits_high_conf) + len(hits_medium_conf) + len(hits_low_conf)
    
    formatter.debug(
        f"Pathogen scan complete: {total_hits} hits "
        f"(HIGH: {len(hits_high_conf)}, MEDIUM: {len(hits_medium_conf)}, LOW: {len(hits_low_conf)})"
    )
    
    # Save filtered results (JSON)
    filtered_file = output_dir / "pathogen_detections_validated.json"
    
    results = {
        'summary': {
            'total_hits': total_hits,
            'high_confidence': len(hits_high_conf),
            'medium_confidence': len(hits_medium_conf),
            'low_confidence': len(hits_low_conf),
            'filters_applied': {
                'min_fragment_length': min_fragment_length,
                'min_identity': min_identity,
                'min_query_coverage': min_query_coverage,
                'min_subject_coverage': min_subject_coverage
            }
        },
        'high_confidence_hits': hits_high_conf,
        'medium_confidence_hits': hits_medium_conf,
        'low_confidence_hits': hits_low_conf
    }
    
    with open(filtered_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    # Export legacy TSV format
    legacy_file = output_dir / "pathogen_results.tsv"
    export_legacy_tsv(results, legacy_file)
    
    formatter.debug(f"Results saved: {legacy_file.name}")
    
    return legacy_file


# ============================================================================
# BACKWARD COMPATIBILITY WRAPPER
# ============================================================================

def run_pathogen_scan(protein_file, output_dir, bracken_results=None, taxonomy_results=None):
    """Backward compatibility wrapper for existing code"""
    # Find required files
    contigs_file = output_dir.parent / "megahit_assembly" / "contigs.fasta"
    prokka_dir = output_dir.parent / "prokka_annotation"
    
    if not contigs_file.exists():
        contigs_file = output_dir / "contigs.fasta"
    
    if not prokka_dir.exists():
        prokka_dir = output_dir
    
    return run_pathogen_scan_v2(
        protein_file=Path(protein_file),
        output_dir=Path(output_dir),
        contigs_file=contigs_file,
        prokka_dir=prokka_dir,
        bracken_results=Path(bracken_results) if bracken_results else None
    )