import subprocess
import os
from pathlib import Path
import json
from Bio import SeqIO
import pandas as pd
from collections import Counter
from .config import *

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
    
    if missing:
        print(f"\nMissing tools: {', '.join(missing)}")
        print("Install with conda:")
        print(f"conda install -c bioconda {' '.join(missing)}")
        print("\nOr with mamba (faster):")
        print(f"mamba install -c bioconda {' '.join(missing)}")
        raise SystemExit("Please install missing dependencies")

def check_database_status():
    """Check the status of all required databases"""
    print("=== Database Status Check ===")
    
    # Required files
    required_files = {
        'Kraken DB': KRAKEN_DB / "hash.k2d",
        'SwissProt DB': SWISSPROT_DB,
        'Pathogen FASTA': PATHOGEN_FASTA_SOURCE,
        'Current Pathogen DB': PATHOGEN_DB,
        'Taxonomy Nodes': TAXDUMP_NODES,
        'Taxonomy Names': TAXDUMP_NAMES,
    }
    
    # Optional files
    optional_files = {
        'Fixed Pathogen DB': PATHOGEN_DB_WITH_TAX,
        'Protein Accession2TaxID': PROT_ACCESSION2TAXID,
        'SeqID2TaxID Map': SEQID2TAXID_MAP,
        'Pathogenic TaxIDs': PATHOGEN_TAXIDS,
        'CARD DB': CARD_DB,
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
        
def create_pathogenic_taxids():
    """Create or update pathogenic taxids file"""
    # Use the pathogen_db directory location
    pathogen_taxids_path = DB_DIR / "pathogen_db" / "pathogenic_taxids.txt"
    
    if not pathogen_taxids_path.exists():
        print(f"Creating pathogenic taxids file: {pathogen_taxids_path}")
        
        # Ensure directory exists
        pathogen_taxids_path.parent.mkdir(exist_ok=True)
        
        # Common pathogenic bacteria and viruses with more complete list
        pathogenic_data = [
            ("1280", "Staphylococcus aureus"),
            ("1392", "Bacillus anthracis"),
            ("620", "Shigella dysenteriae"),
            ("590", "Salmonella enterica"),
            ("287", "Pseudomonas aeruginosa"),
            ("1313", "Streptococcus pneumoniae"),
            ("1350", "Enterococcus faecalis"),
            ("573", "Klebsiella pneumoniae"),
            ("562", "Escherichia coli"),
            ("1763", "Mycobacterium tuberculosis"),
            ("11234", "Human immunodeficiency virus 1"),
            ("11103", "Hepatitis C virus"),
            ("10376", "Epstein-Barr virus"),
            ("10298", "Human herpesvirus 1"),
            ("1747", "Cutibacterium acnes"),
            ("1351", "Enterococcus faecium"),
            ("1423", "Bacillus subtilis"),
            ("1491", "Clostridium botulinum"),
            ("1496", "Clostridioides difficile"),
        ]
        
        with open(pathogen_taxids_path, 'w') as f:
            f.write("TaxID\tPathogen_Name\n")
            for taxid, name in pathogenic_data:
                f.write(f"{taxid}\t{name}\n")
        
        print(f"✓ Created pathogenic taxids file with {len(pathogenic_data)} entries")
    
    # Update global variable
    globals()['PATHOGEN_TAXIDS'] = pathogen_taxids_path

def convert_fastq_to_fasta(fastq_path, output_dir):
    """Convert FASTQ to FASTA for downstream analysis, handling duplicate IDs"""
    fasta_path = output_dir/"converted.fasta"
    
    # Method 1: Use seqtk with deduplication by adding sequence counter
    temp_fasta = output_dir/"temp_converted.fasta"
    cmd = f"seqtk seq -a {fastq_path} > {temp_fasta}"
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

def compute_sequence_statistics(fasta_path, output_dir):
    """Compute comprehensive statistics for FASTA sequences"""
    from Bio import SeqIO
    import json
    
    stats = {
        'total_sequences': 0,
        'total_length': 0,
        'lengths': [],
        'gc_contents': [],
        'n_contents': []
    }
    
    sequences = list(SeqIO.parse(fasta_path, "fasta"))
    
    for seq in sequences:
        seq_str = str(seq.seq).upper()
        length = len(seq_str)
        
        stats['total_sequences'] += 1
        stats['total_length'] += length
        stats['lengths'].append(length)
        
        # GC content
        gc_count = seq_str.count('G') + seq_str.count('C')
        gc_content = (gc_count / length * 100) if length > 0 else 0
        stats['gc_contents'].append(gc_content)
        
        # N content
        n_count = seq_str.count('N')
        n_content = (n_count / length * 100) if length > 0 else 0
        stats['n_contents'].append(n_content)
    
    # Calculate summary statistics
    stats['mean_length'] = stats['total_length'] / stats['total_sequences'] if stats['total_sequences'] > 0 else 0
    stats['mean_gc_content'] = sum(stats['gc_contents']) / len(stats['gc_contents']) if stats['gc_contents'] else 0
    stats['mean_n_content'] = sum(stats['n_contents']) / len(stats['n_contents']) if stats['n_contents'] else 0
    
    # N50 calculation
    sorted_lengths = sorted(stats['lengths'], reverse=True)
    total_len = sum(sorted_lengths)
    cumulative = 0
    stats['n50'] = 0
    
    for length in sorted_lengths:
        cumulative += length
        if cumulative >= total_len * 0.5:
            stats['n50'] = length
            break
    
    # Save statistics
    stats_file = output_dir / "sequence_statistics.json"
    with open(stats_file, 'w') as f:
        json.dump(stats, f, indent=2)
    
    print(f"✓ Computed statistics for {stats['total_sequences']} sequences")
    print(f"  - Total length: {stats['total_length']:,} bp")
    print(f"  - N50: {stats['n50']:,} bp")
    print(f"  - Mean GC content: {stats['mean_gc_content']:.1f}%")
    
    return stats

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