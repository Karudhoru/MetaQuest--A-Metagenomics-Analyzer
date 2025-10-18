import subprocess
import os
from pathlib import Path
import json
from Bio import SeqIO
from typing import List, Optional
import pandas as pd
import re
import numpy as np
from collections import Counter
from ..config import *
import importlib
from .output_formatter import get_formatter

formatter = get_formatter()

def _check_command(name, version_cmd="--version"):
    """Helper to check for a single command-line tool."""
    try:
        # Some tools like seqtk don't have a version flag and exit with an error
        if name == 'seqtk':
            subprocess.run([name], capture_output=True, text=True, check=False)
            return True
        
        cmd = name.split() + [version_cmd] if version_cmd else name.split()
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False

def run_system_check(formatter):
    """
    Runs a comprehensive check of all dependencies using the formatter for output.
    """
    errors = []
    warnings = []
    
    formatter.info("Performing system-wide checks...")

    # --- 1. Command-Line Tools ---
    formatter.substep("Checking command-line tools...")
    tools = {
        'kraken2': '--version', 'bracken': '--version', 'diamond': 'version',
        'prokka': '--version', 'seqkit': 'version', 'ktImportText': None,
        'seqtk': None, 'spades.py': '--version', 'reformat.sh': None
    }
    for tool, cmd in tools.items():
        if not _check_command(tool, cmd):
            errors.append(f"Tool not found: '{tool}'. Please install it, e.g., via 'conda install -c bioconda {tool}'.")

    # --- 2. Python Packages ---
    formatter.substep("Checking Python packages...")
    packages = {
        'Bio': 'biopython', 'pandas': 'pandas', 'numpy': 'numpy', 'plotly': 'plotly',
        'sklearn': 'scikit-learn', 'joblib': 'joblib', 'requests': 'requests',
        'tqdm': 'tqdm', 'scipy': 'scipy', 'matplotlib': 'matplotlib', 'xgboost': 'xgboost',
        'lightgbm': 'lightgbm', 'catboost': 'catboost'
    }
    for pkg, install_name in packages.items():
        try:
            importlib.import_module(pkg)
        except ImportError:
            errors.append(f"Python package not found: '{install_name}'. Please install it, e.g., 'pip install {install_name}'.")

    # --- 3. Machine Learning Artifacts ---
    formatter.substep("Checking ML model artifacts...")
    if not MODEL_ARTIFACTS_DIR.exists():
        errors.append(f"Model artifacts directory not found: {MODEL_ARTIFACTS_DIR}")
    else:
        expected_models = ['best_model.pkl', 'scaler.pkl', 'feature_selector.pkl', 'selected_feature_names.json']
        for model_file in expected_models:
            if not (MODEL_ARTIFACTS_DIR / model_file).exists():
                warnings.append(f"ML artifact missing: '{model_file}'. ML predictions may fail.")

    # --- 4. Databases ---
    formatter.substep("Checking database files...")
    required_dbs = {
        "Kraken2 DB": KRAKEN_DB / "hash.k2d", "SwissProt DB": SWISSPROT_DB,
        "Pathogen DB": PATHOGEN_DB_CUSTOM
    }
    for name, path in required_dbs.items():
        if not path.exists():
            errors.append(f"Required database not found: {name} (expected at {path}).")

    # --- Final Report ---
    if warnings:
        for warning in warnings:
            formatter.warning(warning)
    
    if errors:
        formatter.error(
            "System check failed due to missing critical dependencies.",
            solutions=errors + ["Please resolve the issues above before running an analysis."]
        )
        raise SystemExit()

def get_tool_version(tool_name: str) -> str:
    """
    Dynamically gets the version of a command-line tool.
    
    Args:
        tool_name: The name of the tool (e.g., 'prokka', 'diamond').
        
    Returns:
        The version string or 'N/A' if not found.
    """
    # Dictionary mapping tool names to their version commands
    version_commands = {
        'prokka': ['prokka', '--version'],
        'kraken2': ['kraken2', '--version'],
        'diamond': ['diamond', '--version'],
        'spades': ['spades.py', '--version'],
        'megahit': ['megahit', '--version']
    }
    
    cmd = version_commands.get(tool_name.lower())
    if not cmd:
        return "N/A"
        
    try:
        # Many tools print version info to stderr, so we capture both
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=10)
        output = result.stdout + result.stderr
        
        # Use regex to find version numbers like X.Y.Z or X.Y
        match = re.search(r'(\d+\.\d+(\.\d+)?)', output)
        if match:
            return match.group(1)
        return "Unknown"
    except (FileNotFoundError, subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return "Not Found"

def assemble_reads_to_fasta(reads: List[str], output_dir: Path, formatter, threads: int = 8) -> Path:
    """
    Performs metagenomic assembly on FASTQ reads using SPAdes and the output formatter.
    """
    spades_outdir = output_dir / "spades_assembly"
    final_contigs = spades_outdir / "contigs.fasta"
    
    cmd = ["spades.py", "--meta", "--only-assembler"]
    if len(reads) == 1:
        cmd.extend(["-s", reads[0]])
    elif len(reads) == 2:
        cmd.extend(["-1", reads[0], "-2", reads[1]])
    
    cmd.extend(["-o", str(spades_outdir), "-t", str(threads)])
    
    returncode, _, stderr = formatter.run_subprocess(cmd, "Running SPAdes assembler", show_command=True)

    if returncode != 0 or not final_contigs.exists():
        formatter.error(
            "SPAdes assembly failed.",
            solutions=[f"Check the log file for details: {spades_outdir}/spades.log"]
        )
        raise RuntimeError(f"SPAdes failed. Stderr: {stderr}")

    return final_contigs

def split_interleaved(interleaved_fastq: str, output_dir: Path, formatter) -> list:
    """
    Splits an interleaved FASTQ into two files using reformat.sh and the output formatter.
    """
    r1 = output_dir / "split_R1.fastq.gz"
    r2 = output_dir / "split_R2.fastq.gz"

    cmd = [
        "reformat.sh", f"in={interleaved_fastq}",
        f"out1={r1}", f"out2={r2}", "overwrite=true",
    ]
    
    returncode, _, stderr = formatter.run_subprocess(cmd, "Splitting interleaved file", show_command=True)

    if returncode != 0:
        formatter.error(
            "reformat.sh failed to split the file.",
            solutions=["Ensure BBTools is installed and the input file is a valid FASTQ."]
        )
        raise RuntimeError(f"reformat.sh failed. Stderr: {stderr}")
    
    return [str(r1), str(r2)]

def parse_prokka_gff(gff_file):
    """Parse Prokka GFF file to count features"""
    feature_counts = Counter()
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                feature_type = parts[2]
                feature_counts[feature_type] += 1
    
    return feature_counts

def has_taxonomy_info(pathogen_db_path):
    """Check if Diamond database has taxonomy information"""
    if pathogen_db_path is None:
        return False
    
    try:
        db_path = Path(pathogen_db_path)
        if not db_path.exists():
            print(f"Warning: Pathogen database not found at {db_path}")
            return False
        
        # Check if database was built with taxonomy
        # Try a test command to see if taxonomy fields work
        test_cmd = f"diamond help | grep -q taxon"
        result = subprocess.run(test_cmd, shell=True, capture_output=True)
        
        # For now, assume no taxonomy (safer default)
        return False
    
    except Exception as e:
        print(f"Error checking taxonomy info: {e}")
        return False

def should_use_ml_prediction(prokka_dir, formatter):
    """Determine if ML prediction is appropriate, logging with the formatter."""
    try:
        protein_files = list(Path(prokka_dir).glob("*.faa"))
        if not protein_files: return False
        
        sequences = list(SeqIO.parse(protein_files[0], "fasta"))
        if not sequences: return False
        
        # **FIX: Check if ANY proteins meet threshold, not average**
        lengths = [len(seq.seq) for seq in sequences]
        min_length_threshold = 200
        
        # Count proteins >= threshold
        eligible_proteins = [l for l in lengths if l >= min_length_threshold]
        avg_length = np.mean(lengths)
        
        formatter.info("Protein length analysis for ML suitability:")
        formatter.result({
            "Total proteins": len(sequences),
            "Average protein length": f"{avg_length:.1f} aa",
            "Proteins ≥ 200 aa": len(eligible_proteins),
            "ML Threshold": f"≥ {min_length_threshold} aa",
            "Eligible for ML": f"{len(eligible_proteins)/len(sequences)*100:.1f}%"
        }, indent=2)
        
        # **NEW LOGIC: Run ML if at least 10% of proteins are long enough**
        min_eligible_fraction = 0.10  # 10% threshold
        
        if len(eligible_proteins) / len(sequences) >= min_eligible_fraction:
            formatter.info(f"✓ Sufficient long proteins for ML prediction", indent=2)
            return True
        else:
            formatter.info(f"✗ Too few long proteins ({len(eligible_proteins)}) for reliable ML prediction", indent=2)
            return False
            
    except Exception as e:
        formatter.warning(f"Could not analyze protein lengths for ML: {e}")
        return False

def extract_pathogens_from_bracken(bracken_file, formatter):
    """Extract pathogenic organisms from Bracken output, logging with the formatter."""
    # (Function logic remains the same, but now uses formatter for errors)
    pathogenic_keywords = [
        'salmonella', 'escherichia coli', 'staphylococcus aureus', 'clostridium tetani',
        'klebsiella pneumoniae', 'yersinia', 'streptococcus', 'pseudomonas aeruginosa',
        'acinetobacter baumannii', 'vibrio', 'brucella', 'mycobacterium',
        'bacillus anthracis', 'listeria', 'clostridium difficile', 'shigella'
    ]
    pathogens = []
    try:
        df = pd.read_csv(bracken_file, sep='\t')
        for _, row in df.iterrows():
            if any(keyword in row['name'].lower() for keyword in pathogenic_keywords):
                pathogens.append({
                    'organism': row['name'], 'abundance': row['fraction_total_reads'],
                    'reads': int(row['new_est_reads']), 'taxonomy_id': row.get('taxonomy_id', None)
                })
    except Exception as e:
        formatter.warning(f"Error extracting pathogens from Bracken: {e}")
    return pathogens

def extract_pathogens_from_blast_taxonomy(blast_results, formatter):
    """
    Extract pathogenic organisms from BLAST taxonomy results using the formatter.
    """
    pathogenic_keywords = [
        'salmonella', 'escherichia coli', 'staphylococcus aureus', 'clostridium tetani',
        'klebsiella pneumoniae', 'yersinia enterocolitica', 'yersinia pestis',
        'streptococcus', 'enterococcus faecalis', 'pseudomonas aeruginosa',
        'acinetobacter baumannii', 'vibrio', 'brucella', 'mycobacterium tuberculosis',
        'bacillus anthracis', 'listeria monocytogenes', 'clostridium difficile',
        'francisella tularensis', 'burkholderia', 'rickettsia', 'coxiella',
        'campylobacter', 'helicobacter pylori', 'neisseria gonorrhoeae',
        'shigella', 'enterobacter', 'serratia', 'proteus', 'providencia', 'morganella'
    ]
    
    pathogens = []
    
    try:
        if not blast_results:
            return pathogens
        
        for result in blast_results:
            # Improvement: handle cases where 'organism' might be missing in a result
            organism = result.get('organism', '')
            if not organism:
                continue

            # Check if organism matches any pathogenic keyword
            if any(keyword in organism.lower() for keyword in pathogenic_keywords):
                pathogens.append({
                    'organism': result.get('organism', 'Unknown'),
                    'identity': result.get('identity', 0),
                    'evalue': result.get('evalue', 1.0),
                    'query_id': result.get('query_id', ''),
                    'accession': result.get('accession', '')
                })
    
    except Exception as e:
        # Fix: Use formatter instead of print()
        formatter.warning(f"Error extracting pathogens from BLAST results: {e}")
    
    return pathogens

# Add this to your utils.py or wherever parse_diamond_progress is defined

import re
from typing import Optional

def parse_diamond_progress(line: str) -> Optional[int]:
    """
    Parse DIAMOND output to extract current progress.
    
    DIAMOND outputs progress in several formats:
    - "Processed 12000 queries"
    - "12000 queries aligned"
    - Progress percentage indicators
    - "Queries: 12000/50000"
    
    Returns:
        The absolute number of queries processed, or None if no match found
    """
    # Pattern 1: "Processed X queries" or "X queries processed"
    match = re.search(r'(?:Processed|processed)\s+(\d+)\s+(?:queries|sequences)', line)
    if match:
        return int(match.group(1))
    
    # Pattern 2: "X queries aligned"
    match = re.search(r'(\d+)\s+queries\s+aligned', line)
    if match:
        return int(match.group(1))
    
    # Pattern 3: "Queries: X/Y" format
    match = re.search(r'Queries:\s*(\d+)/\d+', line)
    if match:
        return int(match.group(1))
    
    # Pattern 4: Progress indicator with numbers
    match = re.search(r'(\d+)\s*/\s*\d+\s+queries', line)
    if match:
        return int(match.group(1))
    
    # Pattern 5: Simple number followed by query/sequence indicators
    match = re.search(r'^[>\s]*(\d+)\s+(?:query|queries|sequence|sequences)', line, re.IGNORECASE)
    if match:
        return int(match.group(1))
    
    return None


def parse_diamond_progress_verbose(line: str) -> Optional[int]:
    """
    More verbose version for debugging - prints what it's trying to match.
    Use this temporarily to see what DIAMOND is actually outputting.
    """
    patterns = [
        (r'(?:Processed|processed)\s+(\d+)\s+(?:queries|sequences)', "Processed X queries"),
        (r'(\d+)\s+queries\s+aligned', "X queries aligned"),
        (r'Queries:\s*(\d+)/\d+', "Queries: X/Y"),
        (r'(\d+)\s*/\s*\d+\s+queries', "X / Y queries"),
        (r'^[>\s]*(\d+)\s+(?:query|queries|sequence|sequences)', "Line starts with X query/queries"),
    ]
    
    for pattern, description in patterns:
        match = re.search(pattern, line, re.IGNORECASE)
        if match:
            value = int(match.group(1))
            # Uncomment for debugging:
            print(f"  [DEBUG] Matched '{description}': {value} from line: {line.strip()}")
            return value
    
    return None