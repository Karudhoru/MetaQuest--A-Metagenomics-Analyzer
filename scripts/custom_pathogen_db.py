#!/usr/bin/env python3
"""
MetaQuest Custom Pathogen Database Builder
Downloads sequences for common pathogens and builds a DIAMOND database
"""

import os
import sys
import subprocess
import requests
import gzip
from pathlib import Path
from Bio import Entrez, SeqIO
import time
import json

# Set your email for NCBI (required)
Entrez.email = "your@gmail.com"  # CHANGE THIS TO YOUR EMAIL

# Common pathogenic organisms
PATHOGEN_LIST = {


    # High priority bacterial pathogens
    "Brucella melitensis": {"taxid": "29459", "priority": "high"},
    "Brucella abortus": {"taxid": "29450", "priority": "high"},
    "Brucella suis": {"taxid": "29452", "priority": "high"},
    "Salmonella enterica": {"taxid": "28901", "priority": "high"},
    "Escherichia coli": {"taxid": "511145", "priority": "high"},  # K-12 strain
    "Clostridium difficile": {"taxid": "272563", "priority": "high"},
    "Staphylococcus aureus": {"taxid": "93061", "priority": "high"},
    "Streptococcus pneumoniae": {"taxid": "170187", "priority": "high"},
    "Mycobacterium tuberculosis": {"taxid": "83332", "priority": "high"},
    "Pseudomonas aeruginosa": {"taxid": "287", "priority": "high"},
    "Vibrio cholerae": {"taxid": "243277", "priority": "high"},
    "Yersinia pestis": {"taxid": "214092", "priority": "high"},
    "Francisella tularensis": {"taxid": "177416", "priority": "high"},
    "Bacillus anthracis": {"taxid": "198094", "priority": "high"},
    "Listeria monocytogenes": {"taxid": "169963", "priority": "high"},
    
    # Medium priority pathogens
    "Klebsiella pneumoniae": {"taxid": "272620", "priority": "medium"},
    "Acinetobacter baumannii": {"taxid": "470", "priority": "medium"},
    "Enterococcus faecium": {"taxid": "333849", "priority": "medium"},
    "Campylobacter jejuni": {"taxid": "197", "priority": "medium"},
    "Helicobacter pylori": {"taxid": "85962", "priority": "medium"},
    "Shigella flexneri": {"taxid": "198214", "priority": "medium"},
    "Neisseria meningitidis": {"taxid": "122586", "priority": "medium"},
    "Haemophilus influenzae": {"taxid": "71421", "priority": "medium"},
    "Burkholderia pseudomallei": {"taxid": "320372", "priority": "medium"},
    "Legionella pneumophila": {"taxid": "272624", "priority": "medium"},
    "Streptococcus pyogenes": {"taxid": "160490", "priority": "medium"},
    "Chlamydia trachomatis": {"taxid": "272561", "priority": "medium"},
    "Rickettsia rickettsii": {"taxid": "272947", "priority": "medium"},
    "Coxiella burnetii": {"taxid": "227377", "priority": "medium"},
    
    # Additional important pathogens
    "Borrelia burgdorferi": {"taxid": "224326", "priority": "low"},
    "Treponema pallidum": {"taxid": "243276", "priority": "low"},
    "Mycoplasma pneumoniae": {"taxid": "272634", "priority": "low"},
    "Bartonella henselae": {"taxid": "38323", "priority": "low"},
}

def setup_directories():
    """Create necessary directories"""
    dirs = ['pathogen_sequences', 'pathogen_db_custom', 'logs']
    for d in dirs:
        Path(d).mkdir(exist_ok=True)
    return Path('pathogen_sequences'), Path('pathogen_db_custom'), Path('logs')

def check_dependencies():
    """Check if required tools are available"""
    required_tools = ['diamond']
    missing = []
    
    for tool in required_tools:
        try:
            subprocess.run([tool, '--version'], capture_output=True, check=True)
            print(f"✓ {tool} found")
        except (subprocess.CalledProcessError, FileNotFoundError):
            missing.append(tool)
    
    if missing:
        print(f"❌ Missing required tools: {', '.join(missing)}")
        print("Install with: conda install -c bioconda diamond")
        return False
    return True

def download_pathogen_proteome(organism, taxid, seq_dir, max_retries=3):
    """Download proteome for a specific pathogen"""
    print(f"Downloading {organism} (taxid: {taxid})...")
    
    output_file = seq_dir / f"{organism.replace(' ', '_')}.faa"
    
    # Skip if already exists
    if output_file.exists() and output_file.stat().st_size > 1000:
        print(f"  ✓ Already exists: {output_file}")
        return output_file
    
    for attempt in range(max_retries):
        try:
            # Search for protein sequences
            search_handle = Entrez.esearch(
                db="protein",
                term=f"txid{taxid}[Organism] AND refseq[filter]",
                retmax=2000,  # Limit to avoid huge downloads
                sort="relevance"
            )
            search_results = Entrez.read(search_handle)
            search_handle.close()
            
            if not search_results["IdList"]:
                print(f"  ⚠️ No sequences found for {organism}")
                return None
            
            print(f"  Found {len(search_results['IdList'])} sequences")
            
            # Download sequences in batches
            sequences = []
            batch_size = 200
            id_list = search_results["IdList"]
            
            for i in range(0, len(id_list), batch_size):
                batch_ids = id_list[i:i+batch_size]
                
                fetch_handle = Entrez.efetch(
                    db="protein",
                    id=batch_ids,
                    rettype="fasta",
                    retmode="text"
                )
                
                batch_sequences = list(SeqIO.parse(fetch_handle, "fasta"))
                sequences.extend(batch_sequences)
                fetch_handle.close()
                
                print(f"  Downloaded {len(sequences)} sequences...")
                time.sleep(0.5)  # Be nice to NCBI
            
            # Write sequences to file
            if sequences:
                with open(output_file, 'w') as f:
                    for seq in sequences:
                        # Add organism info to description
                        seq.description = f"{seq.description} [{organism}]"
                        SeqIO.write(seq, f, "fasta")
                
                print(f"  ✓ Saved {len(sequences)} sequences to {output_file}")
                return output_file
            else:
                print(f"  ⚠️ No sequences retrieved for {organism}")
                return None
                
        except Exception as e:
            print(f"  ❌ Attempt {attempt + 1} failed: {e}")
            if attempt < max_retries - 1:
                print(f"  Retrying in 10 seconds...")
                time.sleep(10)
            else:
                print(f"  ❌ Failed to download {organism} after {max_retries} attempts")
                return None

def download_vfdb():
    """Download VFDB (Virulence Factor Database)"""
    print("Downloading VFDB virulence factors...")
    
    vfdb_url = "http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz"
    vfdb_file = Path("pathogen_sequences/VFDB_proteins.faa")
    
    if vfdb_file.exists():
        print("  ✓ VFDB already exists")
        return vfdb_file
    
    try:
        response = requests.get(vfdb_url, stream=True)
        response.raise_for_status()
        
        with open("pathogen_sequences/VFDB_proteins.faa.gz", 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        
        # Decompress
        with gzip.open("pathogen_sequences/VFDB_proteins.faa.gz", 'rt') as f_in:
            with open(vfdb_file, 'w') as f_out:
                f_out.write(f_in.read())
        
        # Clean up
        os.remove("pathogen_sequences/VFDB_proteins.faa.gz")
        
        print(f"  ✓ Downloaded VFDB to {vfdb_file}")
        return vfdb_file
        
    except Exception as e:
        print(f"  ❌ Failed to download VFDB: {e}")
        return None

def combine_sequences(seq_dir, output_file):
    """Combine all downloaded sequences into one file"""
    print("Combining all pathogen sequences...")
    
    combined_count = 0
    with open(output_file, 'w') as outf:
        for fasta_file in seq_dir.glob("*.faa"):
            if fasta_file.name.startswith("combined"):
                continue
                
            print(f"  Adding {fasta_file.name}...")
            with open(fasta_file, 'r') as inf:
                sequences = list(SeqIO.parse(inf, "fasta"))
                for seq in sequences:
                    SeqIO.write(seq, outf, "fasta")
                    combined_count += 1
    
    print(f"✓ Combined {combined_count} sequences into {output_file}")
    return combined_count

def build_diamond_database(fasta_file, db_dir):
    """Build DIAMOND database from combined FASTA file"""
    print("Building DIAMOND database...")
    
    db_path = db_dir / "pathogen_markers"
    
    cmd = [
        "diamond", "makedb",
        "--in", str(fasta_file),
        "--db", str(db_path),
        "--threads", "4"
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print(f"✓ DIAMOND database created: {db_path}.dmnd")
        return f"{db_path}.dmnd"
        
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed to build DIAMOND database: {e}")
        print(f"STDOUT: {e.stdout}")
        print(f"STDERR: {e.stderr}")
        return None

def create_pathogen_list_files(db_dir):
    """Create pathogen taxonomy and name files"""
    
    # Create pathogen taxids file
    taxids_file = db_dir / "pathogen_taxids.txt"
    with open(taxids_file, 'w') as f:
        f.write("# Pathogenic organism taxonomy IDs\n")
        for organism, info in PATHOGEN_LIST.items():
            f.write(f"{info['taxid']}\t{organism}\t{info['priority']}\n")
    
    # Create pathogen names file
    names_file = db_dir / "pathogen_names.txt"
    with open(names_file, 'w') as f:
        f.write("# Pathogenic organism names\n")
        for organism, info in PATHOGEN_LIST.items():
            f.write(f"{organism}\t{info['priority']}\n")
    
    # Create JSON metadata
    metadata_file = db_dir / "pathogen_metadata.json"
    metadata = {
        "total_pathogens": len(PATHOGEN_LIST),
        "high_priority": len([p for p in PATHOGEN_LIST.values() if p["priority"] == "high"]),
        "medium_priority": len([p for p in PATHOGEN_LIST.values() if p["priority"] == "medium"]),
        "low_priority": len([p for p in PATHOGEN_LIST.values() if p["priority"] == "low"]),
        "pathogens": PATHOGEN_LIST
    }
    
    with open(metadata_file, 'w') as f:
        json.dump(metadata, f, indent=2)
    
    print(f"✓ Created pathogen reference files in {db_dir}")

def main():
    """Main function"""
    print("=== Pathogen Database Builder ===")
    print("This script will download sequences for common pathogens and build a DIAMOND database")
    print()
    
    # Check email
    if Entrez.email == "your.email@example.com":
        print("❌ Please set your email address in the script (line 13)")
        print("   This is required by NCBI for API access")
        return
    
    # Check dependencies
    if not check_dependencies():
        return
    
    # Setup directories
    seq_dir, db_dir, log_dir = setup_directories()
    
    # Download sequences for each pathogen
    downloaded_files = []
    
    # Prioritize high-priority pathogens
    high_priority = {k: v for k, v in PATHOGEN_LIST.items() if v["priority"] == "high"}
    medium_priority = {k: v for k, v in PATHOGEN_LIST.items() if v["priority"] == "medium"}
    low_priority = {k: v for k, v in PATHOGEN_LIST.items() if v["priority"] == "low"}
    
    print(f"\n=== Downloading High Priority Pathogens ({len(high_priority)}) ===")
    for organism, info in high_priority.items():
        result = download_pathogen_proteome(organism, info["taxid"], seq_dir)
        if result:
            downloaded_files.append(result)
    
    print(f"\n=== Downloading Medium Priority Pathogens ({len(medium_priority)}) ===")
    for organism, info in medium_priority.items():
        result = download_pathogen_proteome(organism, info["taxid"], seq_dir)
        if result:
            downloaded_files.append(result)
    
    # Ask user if they want low priority pathogens
    print(f"\n=== Low Priority Pathogens Available ({len(low_priority)}) ===")
    download_low = input("Download low priority pathogens? (y/n): ").lower().startswith('y')
    
    if download_low:
        for organism, info in low_priority.items():
            result = download_pathogen_proteome(organism, info["taxid"], seq_dir)
            if result:
                downloaded_files.append(result)
    
    # Download VFDB
    print("\n=== Downloading Virulence Factors ===")
    vfdb_file = download_vfdb()
    if vfdb_file:
        downloaded_files.append(vfdb_file)
    
    if not downloaded_files:
        print("❌ No sequences downloaded. Exiting.")
        return
    
    # Combine sequences
    print(f"\n=== Processing {len(downloaded_files)} sequence files ===")
    combined_file = seq_dir / "combined_pathogens.faa"
    total_sequences = combine_sequences(seq_dir, combined_file)
    
    # Build DIAMOND database
    print("\n=== Building DIAMOND Database ===")
    db_file = build_diamond_database(combined_file, db_dir)
    
    if db_file:
        # Create reference files
        create_pathogen_list_files(db_dir)
        
        print("\n=== SUCCESS ===")
        print(f"✓ Downloaded sequences for {len(downloaded_files)} pathogen sets")
        print(f"✓ Combined {total_sequences} protein sequences")
        print(f"✓ Built DIAMOND database: {db_file}")
        print(f"✓ Created reference files in: {db_dir}")
        print()
        print("Usage in your pipeline:")
        print(f"  Database path: {db_file}")
        print(f"  Reference files: {db_dir}/pathogen_*.txt")
        print()
        print("Update your config.py to use:")
        print(f"  PATHOGEN_DB = Path('{db_file}')")
    else:
        print("❌ Failed to build DIAMOND database")

if __name__ == "__main__":
    main()
