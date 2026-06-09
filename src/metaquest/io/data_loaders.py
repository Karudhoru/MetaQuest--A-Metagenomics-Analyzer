"""
MetaQuest Data Loading Utilities
=================================
Centralized file loading with validation and consistent column naming.
"""

import pandas as pd
import json
from pathlib import Path
from Bio import SeqIO
from typing import Dict, List, Optional, Literal
from ..io.output_formatter import get_formatter

def load_bracken_report(bracken_file: Path) -> pd.DataFrame:
    """
    Load Bracken species abundance report.
    
    Args:
        bracken_file: Path to bracken_report.tsv
        
    Returns:
        DataFrame with standardized columns
        
    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If required columns missing
        
    Example:
        >>> df = load_bracken_report(Path("bracken_report.tsv"))
        >>> print(df.columns)
        ['name', 'taxonomy_id', 'taxonomy_lvl', 'kraken_assigned_reads', 
         'added_reads', 'new_est_reads', 'fraction_total_reads']
    """
    if not bracken_file.exists():
        raise FileNotFoundError(f"Bracken report not found: {bracken_file}")
    
    df = pd.read_csv(bracken_file, sep='\t')
    
    # Validate required columns
    required = ['name', 'taxonomy_id', 'new_est_reads', 'fraction_total_reads']
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"Bracken file missing required columns: {missing}")
    
    return df


def load_annotation_file(file_path: Path) -> pd.DataFrame:
    """Load annotation file safely."""
    if not file_path.exists():
        return pd.DataFrame()
        
    try:
        df = pd.read_csv(file_path, sep='\t', header=None)
            
        standard_cols = ['query_id', 'subject_id', 'identity', 'length', 'mismatches',
                        'gap_opens', 'q_start', 'q_end', 's_start', 's_end', 'evalue',
                        'bit_score', 'description']
            
        if len(df.columns) > len(standard_cols):
            extra_count = len(df.columns) - len(standard_cols)
            ann_cols = standard_cols + [f'extra_{i}' for i in range(extra_count)]
        else:
            ann_cols = standard_cols[:len(df.columns)]
        
        df.columns = ann_cols
            
        return df
    except Exception as e:
        get_formatter().debug(f"Error loading annotation file: {e}")
        return pd.DataFrame()
    
def load_ml_predictions(ml_file: Path) -> List[Dict]:
    """Load ML predictions from JSON."""
    if not ml_file.exists():
        return []
        
    try:
        with open(ml_file, 'r') as f:
            data = json.load(f)
            return data.get('predictions', [])
    except Exception as e:
        get_formatter().debug(f"Could not load ML predictions: {e}")
        return []


def load_pathogen_hits(pathogen_file: Path) -> pd.DataFrame:
    """
    Load pathogen database hits (supports both JSON and TSV formats).
    
    Args:
        pathogen_file: Path to pathogen detection results
                      (pathogen_results.tsv or pathogen_detections_validated.json)
        
    Returns:
        DataFrame with pathogen hits (empty if no hits)
        
    Example:
        >>> df = load_pathogen_hits(Path("pathogen_results.tsv"))
        >>> print(len(df))  # Number of pathogen hits
    """
    if not pathogen_file.exists():
        return pd.DataFrame()
    
    try:
        # Try JSON format first (new format from pathogen_analysis.py)
        if pathogen_file.suffix == '.json':
            with open(pathogen_file) as f:
                json_data = json.load(f)
            
            # Combine HIGH and MEDIUM confidence hits
            all_hits = (json_data.get('high_confidence_hits', []) + 
                       json_data.get('medium_confidence_hits', []))
            
            if all_hits:
                return pd.DataFrame(all_hits)
            else:
                return pd.DataFrame()
        
        # Try TSV format (legacy format)
        else:
            df = pd.read_csv(pathogen_file, sep='\t', header=None)
            
            # Assign standard column names
            tsv_cols = [
                'query_id', 'subject_id', 'identity', 'length', 'mismatch', 
                'gapopen', 'q_start', 'q_end', 's_start', 's_end', 'evalue', 
                'bitscore', 'qlen', 'slen', 'description'
            ]
            
            if len(df.columns) == 15:
                df.columns = tsv_cols
            elif len(df.columns) < 15:
                df.columns = tsv_cols[:len(df.columns)]
            else:
                extra = len(df.columns) - 15
                df.columns = tsv_cols + [f'extra_{i}' for i in range(extra)]
            
            return df
            
    except Exception as e:
        get_formatter().debug(f"Could not load pathogen hits from {pathogen_file}: {e}")
        return pd.DataFrame()


def load_prokka_stats(sample_txt: Path) -> Dict[str, int]:
    """
    Parse Prokka sample.txt statistics file.
    
    Args:
        sample_txt: Path to prokka_annotation/sample.txt
        
    Returns:
        Dict with keys: contigs, bases, CDS, rRNA, tRNA, etc.
        
    Example:
        >>> stats = load_prokka_stats(Path("prokka_annotation/sample.txt"))
        >>> print(stats['CDS'])  # Number of coding sequences
    """
    stats = {}
    
    if not sample_txt.exists():
        return stats
    
    try:
        with open(sample_txt, 'r') as f:
            for line in f:
                if ':' in line:
                    key, value = line.strip().split(':', 1)
                    key = key.strip()
                    value = value.strip()
                    
                    # Try to convert to int
                    try:
                        stats[key] = int(value)
                    except ValueError:
                        stats[key] = value
    except Exception as e:
        get_formatter().debug(f"Could not parse Prokka stats: {e}")
    
    return stats

def load_protein_sequences_streaming(faa_file: Path, needed_ids: set) -> Dict[str, str]:
    """
    Load only needed protein sequences (99% memory reduction).
    
    FIXED: Properly handles duplicate IDs in FASTA files by tracking
    which IDs have been found rather than counting sequences.
    
    Args:
        faa_file: Path to protein FASTA file (.faa)
        needed_ids: Set of protein IDs to load
        
    Returns:
        Dict mapping protein IDs to sequences (only requested IDs)
        
    Example:
        >>> needed = {'gene_001', 'gene_042', 'gene_123'}
        >>> seqs = load_protein_sequences_streaming(Path("sample.faa"), needed)
        >>> len(seqs)  # Will be ≤3 (some IDs may not be in file)
        3
    """
    sequences = {}
    needed_ids_remaining = needed_ids.copy()
    
    for record in SeqIO.parse(faa_file, "fasta"):
        if record.id in needed_ids_remaining:
            sequences[record.id] = str(record.seq)
            needed_ids_remaining.remove(record.id)
            

            if not needed_ids_remaining:
                break
    

    if needed_ids_remaining:
        from ..io.output_formatter import get_formatter
        fmt = get_formatter()
        fmt.debug(f"Warning: {len(needed_ids_remaining)} protein IDs not found in {faa_file.name}")
        fmt.debug(f"Missing IDs (first 5): {list(needed_ids_remaining)[:5]}")
    
    return sequences
