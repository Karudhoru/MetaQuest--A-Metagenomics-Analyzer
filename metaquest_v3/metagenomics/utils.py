import subprocess
import os
from pathlib import Path
import json
from Bio import SeqIO
import pandas as pd
from collections import Counter
from .config import *
import requests
import time
from datetime import datetime, timedelta

def check_blast_dependencies():
    """Check if required Python packages for BLAST API are available"""
    required_packages = ['Bio', 'requests']
    missing = []
    
    for package in required_packages:
        try:
            __import__(package)
            print(f"✓ {package} found")
        except ImportError:
            missing.append(package)
            print(f"✗ {package} not found")
    
    if missing:
        print(f"\nMissing Python packages: {', '.join(missing)}")
        print("Install with:")
        if 'Bio' in missing:
            print("conda install -c bioconda biopython")
        if 'requests' in missing:
            print("pip install requests")
        return False
    return True

def check_ncbi_api_status():
    """Check if NCBI BLAST API is accessible"""
    try:
        print("Checking NCBI BLAST API status...")
        response = requests.get('https://blast.ncbi.nlm.nih.gov/Blast.cgi', timeout=10)
        if response.status_code == 200:
            print("✓ NCBI BLAST API is accessible")
            return True
        else:
            print(f"⚠️ NCBI BLAST API returned status code: {response.status_code}")
            return False
    except requests.RequestException as e:
        print(f"✗ Cannot access NCBI BLAST API: {e}")
        return False

def cleanup_blast_cache(cache_dir, max_age_days=30):
    """Clean up old cache entries"""
    if not cache_dir.exists():
        return
    
    cache_file = cache_dir / "blast_cache.json"
    if not cache_file.exists():
        return
    
    try:
        with open(cache_file, 'r') as f:
            cache = json.load(f)
        
        current_time = time.time()
        cutoff_time = current_time - (max_age_days * 24 * 60 * 60)
        
        # Remove old entries
        cleaned_cache = {}
        removed_count = 0
        
        for key, value in cache.items():
            if isinstance(value, dict) and 'timestamp' in value:
                if value['timestamp'] > cutoff_time:
                    cleaned_cache[key] = value
                else:
                    removed_count += 1
            else:
                # Keep entries without timestamp (backwards compatibility)
                cleaned_cache[key] = value
        
        if removed_count > 0:
            with open(cache_file, 'w') as f:
                json.dump(cleaned_cache, f, indent=2)
            print(f"✓ Cleaned {removed_count} old cache entries")
        
    except Exception as e:
        print(f"Warning: Could not clean cache: {e}")

def estimate_blast_time(num_sequences, avg_length=1000):
    """Estimate how long BLAST analysis will take"""
    # Rough estimates based on API rate limits and processing time
    seconds_per_sequence = 2  # Conservative estimate including rate limiting
    total_seconds = num_sequences * seconds_per_sequence
    
    if total_seconds < 60:
        return f"~{total_seconds} seconds"
    elif total_seconds < 3600:
        return f"~{total_seconds/60:.1f} minutes"
    else:
        return f"~{total_seconds/3600:.1f} hours"

def validate_fasta_for_blast(fasta_path, min_length=50, max_sequences=1000):
    """Validate FASTA file for BLAST analysis"""
    if not fasta_path.exists():
        return False, "FASTA file does not exist"
    
    sequences = list(SeqIO.parse(fasta_path, "fasta"))
    
    if not sequences:
        return False, "No sequences found in FASTA file"
    
    valid_sequences = [seq for seq in sequences if len(seq.seq) >= min_length]
    
    if not valid_sequences:
        return False, f"No sequences longer than {min_length} bp found"
    
    if len(sequences) > max_sequences:
        return True, f"Warning: {len(sequences)} sequences found, will process first {max_sequences}"
    
    return True, f"Ready to process {len(valid_sequences)} sequences"

def check_dependencies():
    """Verify required tools are installed"""
    required = {
        'kraken2': 'kraken2 --version',
        'bracken': 'bracken --version', 
        'diamond': 'diamond version',
        'prokka': 'prokka --version',
        'taxonkit': 'taxonkit version',
        'ktImportText': 'ktImportText --version',
        'seqtk': 'seqtk'
    }
    
    missing = []
    for cmd, ver_cmd in required.items():
        try:
            if cmd == 'seqtk':
                result = subprocess.run([cmd], capture_output=True, text=True)
                # seqtk returns error code when run without args, but that's expected
            else:
                subprocess.run(ver_cmd.split(), check=True, 
                             stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            print(f"✓ {cmd} found")
        except (subprocess.CalledProcessError, FileNotFoundError):
            missing.append(cmd)
            print(f"✗ {cmd} not found")
    
    # Check BLAST API dependencies
    blast_ready = check_blast_dependencies()
    api_accessible = check_ncbi_api_status()
    
    if missing:
        print(f"\nMissing tools: {', '.join(missing)}")
        print("Install with conda:")
        print(f"conda install -c bioconda {' '.join(missing)}")
        print("\nOr with mamba (faster):")
        print(f"mamba install -c bioconda {' '.join(missing)}")
        raise SystemExit("Please install missing dependencies")
    
    if not blast_ready:
        print("\nBLAST API dependencies missing - FASTA analysis will be limited")
    
    if not api_accessible:
        print("\nWarning: NCBI API not accessible - FASTA taxonomic analysis may fail")
    
    return True

def check_database_status():
    """Check the status of all required databases"""
    print("=== Database Status Check ===")
    
    # Required files
    required_files = {
        'Kraken DB': KRAKEN_DB / "hash.k2d",
        'SwissProt DB': SWISSPROT_DB,
        'Pathogen FASTA': CAT_FASTA_SOURCE,
        'Current Pathogen DB': CAT_DB,
        'Taxonomy Nodes': TAXDUMP_NODES,
        'Taxonomy Names': TAXDUMP_NAMES,
    }
    
    # Optional files
    optional_files = {
        'Fixed Pathogen DB': PATHOGEN_DB_WITH_TAX,
        'Protein Accession2TaxID': PROT_ACCESSION2TAXID,
        'SeqID2TaxID Map': SEQID2TAXID_MAP,
        'Pathogenic TaxIDs': PATHOGEN_TAXIDS,
        'CARD DB': CARD_PROTEIN_DB,
        'VFDB DB': VFDB_DB,
    }
    
    missing_required = []
    missing_optional = []
    
    for name, path in required_files.items():
        if path.exists():
            print(f"✓ {name}: {path}")
        else:
            print(f"✗ {name}: {path} (MISSING - REQUIRED)")
            missing_required.append(name)
    
    for name, path in optional_files.items():
        if path.exists():
            print(f"✓ {name}: {path}")
        else:
            print(f"- {name}: {path} (missing - optional)")
            missing_optional.append(name)
    
    if missing_required:
        print(f"\n⚠️  CRITICAL: {len(missing_required)} required files are missing!")
        print("The pipeline will fail without these files.")
    
    if 'Fixed Pathogen DB' in missing_optional:
        print(f"\n⚠️  IMPORTANT: Pathogen database lacks taxonomy information.")
        print("This will cause pathogen scanning to fail (as noted in the document).")
        print("Run the database rebuild command to fix this.")
    
    return len(missing_required) == 0

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
    Split a single interleaved FASTQ into two files (R1, R2) using seqkit.
    Requires seqkit in your PATH.
    """
    r1 = output_dir / "split_R1.fastq"
    r2 = output_dir / "split_R2.fastq"
    cmd = (
        f"seqkit split2 -1 {r1} -2 {r2} "
        f"--by-pair {interleaved_fastq}"
    )
    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)
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