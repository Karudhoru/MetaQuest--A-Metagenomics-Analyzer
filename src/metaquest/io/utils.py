import subprocess
import os
from pathlib import Path
import json
from Bio import SeqIO
from typing import List
import pandas as pd
import re
import numpy as np
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
        'seqtk': None,
        'spades.py': '--version',
        'reformat.sh': None  # from BBTools
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
        "Annotated DB": SWISSPROT_DB,
        "Pathogen Screening DB": PATHOGEN_DB_CUSTOM,

    }

    for name, path in required_dbs.items():
        if not path.exists():
            errors.append(f"Required database not found: {name} (expected at {path}).")

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

def assemble_reads_to_fasta(reads: List[str], output_dir: Path, threads: int = 8) -> Path:
    """
    Performs metagenomic assembly on FASTQ reads using SPAdes.
    Handles both single-end and paired-end reads.
    """
    print("--- Assembling reads into contigs with SPAdes ---")
    spades_outdir = output_dir / "spades_assembly"
    final_contigs = spades_outdir / "contigs.fasta"  # SPAdes uses this output name
    
    # Build the SPAdes command
    cmd = ["spades.py", "--meta", "--only-assembler"] # Use '--meta' for metagenomic assembly
    if len(reads) == 1:
        cmd.extend(["-s", reads[0]])
    elif len(reads) == 2:
        cmd.extend(["-1", reads[0], "-2", reads[1]])
    else:
        raise ValueError(f"Invalid number of read files for SPAdes: {len(reads)}. Expected 1 or 2.")

    cmd.extend(["-o", str(spades_outdir), "-t", str(threads)])
    
    try:
        print(f"Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        
        if final_contigs.exists():
            print(f"✅ Assembly complete. Contigs saved to: {final_contigs}")
            return final_contigs
        else:
            raise FileNotFoundError("SPAdes finished, but the final contig file was not found.")
            
    except FileNotFoundError:
        print("❌ ERROR: 'spades.py' command not found. Please ensure it is installed.")
        raise
    except subprocess.CalledProcessError as e:
        print(f"❌ SPAdes assembly failed. See error details below:")
        # SPAdes often writes detailed errors to a log file
        print(f"Check the log file for more details: {spades_outdir}/spades.log")
        print(e.stderr)
        raise

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
        f"out2={r2}",
        "overwrite=true",
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

def should_use_ml_prediction(prokka_dir):
    """Determine if ML prediction is appropriate based on sequence lengths"""
    try:
        prokka_path = Path(prokka_dir)
        protein_files = list(prokka_path.glob("*.faa"))
        
        if not protein_files:
            return False
        
        from Bio import SeqIO
        sequences = list(SeqIO.parse(protein_files[0], "fasta"))
        
        if not sequences:
            return False
        
        lengths = [len(seq.seq) for seq in sequences]
        avg_length = np.mean(lengths)
        min_length_threshold = 200  
        
        print(f"🔍 Protein length analysis:")
        print(f"   • Average protein length: {avg_length:.1f} amino acids")
        print(f"   • Total proteins: {len(sequences)}")
        
        if avg_length >= min_length_threshold:
            print(f"   ✅ Suitable for ML prediction (avg ≥ {min_length_threshold} aa)")
            return True
        else:
            print(f"   ⚠️ Too short for ML prediction (avg < {min_length_threshold} aa)")
            print(f"   📏 ML model trained on 500-1000 aa proteins")
            return False
            
    except Exception as e:
        print(f"⚠️ Could not analyze protein lengths: {e}")
        return False


def extract_pathogens_from_bracken(bracken_file):
    """
    Extract pathogenic organisms from Bracken output file.
    
    Args:
        bracken_file: Path to Bracken output file
        
    Returns:
        List of dictionaries containing pathogen information
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
        df = pd.read_csv(bracken_file, sep='\t')
        
        for _, row in df.iterrows():
            organism_name = row['name'].lower()
            
            # Check if organism matches any pathogenic keyword
            if any(keyword in organism_name for keyword in pathogenic_keywords):
                pathogens.append({
                    'organism': row['name'],
                    'abundance': row['fraction_total_reads'],
                    'reads': int(row['new_est_reads']),
                    'taxonomy_id': row.get('taxonomy_id', None)
                })
    
    except Exception as e:
        print(f"⚠️ Error extracting pathogens from Bracken: {e}")
    
    return pathogens


def extract_pathogens_from_blast_taxonomy(blast_results):
    """
    Extract pathogenic organisms from BLAST taxonomy results.
    
    Args:
        blast_results: List of BLAST result dictionaries or JSON data
        
    Returns:
        List of dictionaries containing pathogen information
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
            organism = result.get('organism', '').lower()
            
            # Check if organism matches any pathogenic keyword
            if any(keyword in organism for keyword in pathogenic_keywords):
                pathogens.append({
                    'organism': result.get('organism', 'Unknown'),
                    'identity': result.get('identity', 0),
                    'evalue': result.get('evalue', 1.0),
                    'query_id': result.get('query_id', ''),
                    'accession': result.get('accession', '')
                })
    
    except Exception as e:
        print(f"⚠️ Error extracting pathogens from BLAST: {e}")
    
    return pathogens
