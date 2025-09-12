import subprocess
import os
from pathlib import Path
import json
from Bio import SeqIO
import pandas as pd
from collections import Counter
from ..config import *
import importlib

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

def run_system_check():
    """
    Runs a comprehensive check of all dependencies: tools, packages, models, and databases.
    Provides a single, clean report of the system's status.
    """
    errors = []
    warnings = []
    
    print("Performing system-wide checks...")

    # --- 1. Command-Line Tools ---
    print("  -> Checking command-line tools...")
    tools = {
        'kraken2': '--version',
        'bracken': '--version',
        'diamond': 'version',
        'prokka': '--version',
        'seqkit': 'version',
        'ktImportText': None, # ktImportText just needs to be runnable
        'seqtk': None
    }
    for tool, cmd in tools.items():
        if not _check_command(tool, cmd):
            errors.append(f"Tool not found: '{tool}'. Please install it, e.g., via 'conda install -c bioconda {tool}'.")

    # --- 2. Python Packages ---
    print("  -> Checking Python packages...")
    packages = {
        'Bio': 'biopython',
        'pandas': 'pandas',
        'numpy': 'numpy',
        'plotly': 'plotly',
        'sklearn': 'scikit-learn',
        'joblib': 'joblib',
        'requests': 'requests'
    }
    for pkg, install_name in packages.items():
        try:
            importlib.import_module(pkg)
        except ImportError:
            errors.append(f"Python package not found: '{install_name}'. Please install it, e.g., 'pip install {install_name}'.")

    # --- 3. Machine Learning Artifacts ---
    print("  -> Checking ML model artifacts...")
    if not MODEL_ARTIFACTS_DIR.exists():
        errors.append(f"Model artifacts directory not found: {MODEL_ARTIFACTS_DIR}")
    else:
        expected_models = [
            'best_model.pkl', 'scaler.pkl', 'feature_selector.pkl',
            'all_feature_names.pkl', 'feature_names.pkl'
        ]
        for model_file in expected_models:
            if not (MODEL_ARTIFACTS_DIR / model_file).exists():
                warnings.append(f"ML artifact missing: '{model_file}'. ML predictions may fail.")

    # --- 4. Databases ---
    print("  -> Checking database files...")
    required_dbs = {
        "Kraken2 DB (hash)": KRAKEN_DB / "hash.k2d",
        "SwissProt DB": SWISSPROT_DB,
        "Taxonomy Nodes": TAXDUMP_NODES,
        "Taxonomy Names": TAXDUMP_NAMES,
    }
    optional_dbs = {
        "Pathogen Screening DB": PATHOGEN_DB,
        "CARD Protein DB": CARD_PROTEIN_DB,
        "VFDB": VFDB_DB,
    }
    for name, path in required_dbs.items():
        if not path.exists():
            errors.append(f"Required database not found: {name} (expected at {path}).")
    for name, path in optional_dbs.items():
        if not path.exists():
            warnings.append(f"Optional database not found: {name} (expected at {path}). Some features may be disabled.")

    # --- Final Report ---
    print("\n" + "="*50)
    print("          SYSTEM CHECK COMPLETE")
    print("="*50)

    if not errors and not warnings:
        print("\n✅  Success! All dependencies and databases are correctly configured.")
        return True

    if warnings:
        print("\n⚠️  Warnings Found:")
        for warning in warnings:
            print(f"  - {warning}")
    
    if errors:
        print("\n❌ CRITICAL ERRORS FOUND:")
        for error in errors:
            print(f"  - {error}")
        print("\nPlease resolve the critical errors above before running an analysis.")
        raise SystemExit("System check failed due to missing critical dependencies.")
    
    print("\nSystem check passed with warnings. Some functionality may be limited.")
    return True


def convert_fastq_to_fasta(fastq_path, output_dir):
    """Convert FASTQ to FASTA for downstream analysis, handling duplicate IDs"""
    fasta_path = output_dir/"converted.fasta"
    
    # Method 1: Use seqtk with deduplication by adding sequence counter
    temp_fasta = output_dir/"temp_converted.fasta"
    # fastq_path may be str or list of two strings
    if isinstance(fastq_path, (list,tuple)):
        # merge both reads
        fastq_str = " ".join(fastq_path)
    else:
        fastq_str = fastq_path
    cmd = f"seqtk seq -a {fastq_str} > {temp_fasta}"
    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    
    # Remove duplicates and add unique identifiers
    seen_ids = set()
    counter = 0
    
    with open(temp_fasta, 'r') as infile, open(fasta_path, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                original_id = line.strip()[1:]  # Remove '>'
                base_id = original_id.split()[0]  # Take first part before any spaces
                
                # Create unique ID
                if base_id in seen_ids:
                    unique_id = f"{base_id}_{counter}"
                    counter += 1
                else:
                    unique_id = base_id
                    seen_ids.add(base_id)
                
                outfile.write(f">{unique_id}\n")
            else:
                outfile.write(line)
    
    # Clean up temp file
    temp_fasta.unlink()
    
    print(f"✓ Converted FASTQ to FASTA with unique IDs: {fasta_path}")
    return fasta_path


def split_interleaved(interleaved_fastq: str, output_dir: Path) -> list:
    """
    Splits a single interleaved FASTQ into two gzipped files (R1, R2)
    using the reformat.sh utility from the BBTools suite.
    """
    
    r1 = output_dir / "split_R1.fastq.gz"
    r2 = output_dir / "split_R2.fastq.gz"

    cmd = [
        "reformat.sh",
        f"in={interleaved_fastq}",
        f"out1={r1}",
        f"out2={r2}"
    ]
    
    try:
        print(f"Running: {' '.join(cmd)}")
        # Using a list of args is safer than shell=True
        subprocess.run(cmd, check=True, capture_output=True)
        return [str(r1), str(r2)]
    except FileNotFoundError:
        print("❌ ERROR: 'reformat.sh' not found. Please install BBTools via 'conda install -c bioconda bbmap'.")
        raise
    except subprocess.CalledProcessError as e:
        print(f"❌ reformat.sh failed. Please ensure your input file is a valid FASTQ file.")
        print(f"Error details: {e.stderr.decode()}")
        raise

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
