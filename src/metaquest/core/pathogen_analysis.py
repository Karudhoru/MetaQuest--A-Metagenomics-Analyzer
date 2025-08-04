import subprocess
import os
import json
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import re
from pathlib import Path
from collections import Counter
from ..config import *
from ..io.utils import check_dependencies

def run_virulence_factor_scan(fasta_path, output_dir):
    """Scan for virulence factors using VFDB"""
        
    if not VFDB_DB.exists():
        print("VFDB database not found, skipping virulence scan")
        return None
    
    vf_out = output_dir / "virulence_hits.txt"
    
    # Create temporary DIAMOND database
    temp_dmnd = output_dir / "vfdb_temp.dmnd"
    makedb_cmd = f"diamond makedb --in {VFDB_DB} --db {temp_dmnd.with_suffix('')}"
    
    try:
        subprocess.run(makedb_cmd, shell=True, check=True)
        
        cmd = f"diamond blastx -d {temp_dmnd.with_suffix('')} -q {fasta_path} -o {vf_out} " \
              "--outfmt 6 qseqid sseqid pident length evalue bitscore stitle " \
              "--top 5 --evalue 1e-5 --threads 4"
        
        subprocess.run(cmd, shell=True, check=True)
        print("✓ Virulence factor scan completed")
        
        # Clean up
        temp_dmnd.unlink()
        
        return vf_out
        
    except subprocess.CalledProcessError as e:
        print(f"Virulence scan failed: {e}")
        return None

# Enhanced wrapper functions to maintain compatibility
def run_antimicrobial_resistance_scan(fasta_path, output_dir):
    """Enhanced AMR scan with improved reporting"""
    
    if not CARD_PROTEIN_DB.exists():
        print("CARD database not found, skipping AMR scan")
        return None
    
    amr_out = output_dir / "amr_hits.txt"
    
    # Create temporary DIAMOND database if needed
    temp_dmnd = output_dir / "card_temp.dmnd"
    makedb_cmd = f"diamond makedb --in {CARD_PROTEIN_DB} --db {temp_dmnd.with_suffix('')}"
    
    try:
        subprocess.run(makedb_cmd, shell=True, check=True)
        
        # Run DIAMOND search
        cmd = f"diamond blastx -d {temp_dmnd.with_suffix('')} -q {fasta_path} -o {amr_out} " \
              "--outfmt 6 qseqid sseqid pident length evalue bitscore stitle " \
              "--top 5 --evalue 1e-5 --threads 4"
        
        subprocess.run(cmd, shell=True, check=True)
        print("✓ AMR scan completed")

        # Clean up
        for ext in ['.dmnd']:
            temp_file = temp_dmnd.with_suffix(ext)
            if temp_file.exists():
                temp_file.unlink()
        
        return amr_out
        
    except subprocess.CalledProcessError as e:
        print(f"AMR scan failed: {e}")
        return None

def run_pathogen_scan(fasta_path, output_dir, bracken_results=None, taxonomy_results=None):
    """Enhanced pathogen scan using both Bracken (FASTQ) and taxonomy (BLAST) data"""
    
    # First, check for pathogens in Bracken results (FASTQ taxonomic classification)
    bracken_pathogens = []
    if bracken_results and bracken_results.exists():
        try:
            import pandas as pd
            df_bracken = pd.read_csv(bracken_results, sep='\t')
            
            # Known pathogen taxonomy IDs and names
            pathogen_taxa = {
                '29459': 'Brucella melitensis',
                '28901': 'Salmonella enterica', 
                '562': 'Escherichia coli',
                '1280': 'Staphylococcus aureus',
                '1396': 'Bacillus cereus',
                '1428': 'Bacillus thuringiensis',
                '573': 'Klebsiella pneumoniae',
                '470': 'Acinetobacter baumannii',
                '1513': 'Clostridium tetani',
                '630': 'Yersinia enterocolitica',
                '1076': 'Rhodopseudomonas palustris',
                '1719': 'Corynebacterium pseudotuberculosis'
            }
            
            # Process Bracken results
            for _, row in df_bracken.iterrows():
                taxid = str(row['taxonomy_id'])
                name = row['name']
                abundance = row['fraction_total_reads']
                reads = row['new_est_reads']
                
                # Direct taxid match
                if taxid in pathogen_taxa:
                    bracken_pathogens.append({
                        'organism': pathogen_taxa[taxid],
                        'abundance': abundance,
                        'reads': reads,
                        'source': 'bracken_taxonomic',
                        'confidence': 'high',
                        'taxonomy_id': taxid
                    })
                # Name-based matching for common pathogens
                elif any(pathogen in name.lower() for pathogen in [
                    'brucella', 'salmonella', 'escherichia', 'staphylococcus aureus',
                    'bacillus cereus', 'klebsiella', 'acinetobacter', 'clostridium',
                    'yersinia', 'pseudomonas aeruginosa', 'mycobacterium tuberculosis',
                    'listeria', 'vibrio', 'campylobacter', 'shigella'
                ]):
                    bracken_pathogens.append({
                        'organism': name,
                        'abundance': abundance,
                        'reads': reads,
                        'source': 'bracken_taxonomic',
                        'confidence': 'medium',
                        'taxonomy_id': taxid
                    })
            
            print(f"Found {len(bracken_pathogens)} pathogens in Bracken classification")
            
        except Exception as e:
            print(f"Warning: Could not parse Bracken results: {e}")
    
    # Second, check for pathogens in taxonomy results (BLAST-based)
    taxonomy_pathogens = []
    if taxonomy_results and taxonomy_results.exists():
        try:
            # Read taxonomy results (assuming it's tab-separated like pathogen_results.txt)
            with open(taxonomy_results, 'r') as f:
                for line in f:
                    if line.strip() and not line.startswith('#'):
                        parts = line.strip().split('\t')
                        if len(parts) >= 7:  # Expected format from your pathogen_results.txt
                            query_id = parts[0]
                            subject_id = parts[1]
                            identity = float(parts[2])
                            length = int(parts[3])
                            evalue = float(parts[4])
                            bitscore = float(parts[5])
                            description = parts[6]
                            
                            # Extract organism from description
                            organism_matches = re.findall(r'\[([^\]]+)\]', description)
                            if organism_matches:
                                organism = organism_matches[-1]  # Take the last match (usually species)
                                
                                # Check if it's a known pathogen
                                org_lower = organism.lower()
                                if any(pathogen in org_lower for pathogen in [
                                    'coxiella', 'mycobacterium', 'shigella', 'brucella', 
                                    'salmonella', 'escherichia', 'staphylococcus', 'bacillus',
                                    'klebsiella', 'acinetobacter', 'clostridium', 'yersinia',
                                    'pseudomonas', 'listeria', 'vibrio', 'campylobacter'
                                ]):
                                    taxonomy_pathogens.append({
                                        'organism': organism,
                                        'identity': identity,
                                        'evalue': evalue,
                                        'bitscore': bitscore,
                                        'length': length,
                                        'query_id': query_id,
                                        'source': 'taxonomy_blast',
                                        'confidence': 'high' if identity > 80 else 'medium'
                                    })
            
            print(f"Found {len(taxonomy_pathogens)} pathogens in taxonomy BLAST results")
            
        except Exception as e:
            print(f"Warning: Could not parse taxonomy results: {e}")
    
    # Now run additional sequence-based pathogen database search if available
    selected_db = PATHOGEN_DB_CUSTOM if 'PATHOGEN_DB_CUSTOM' in globals() else None
    blast_out = output_dir / "pathogen_results.txt"
    sequence_pathogens = []
    
    if selected_db and selected_db.exists():
        print(f"Running additional pathogen screening against custom database: {selected_db}")
        
        # Remove .dmnd extension for diamond command
        db_base = selected_db.with_suffix('') if selected_db.suffix == '.dmnd' else selected_db
        
        cmd = f"diamond blastx -d {db_base} -q {fasta_path} -o {blast_out} " \
              "--outfmt 6 qseqid sseqid pident length evalue bitscore stitle " \
              "--top 5 --evalue 1e-3 --threads 4 --more-sensitive"
        
        try:
            print("Running DIAMOND search...")
            result = subprocess.run(cmd, shell=True, check=True, timeout=1800, 
                                  capture_output=True, text=True)
            
            if result.stderr:
                print(f"DIAMOND warnings: {result.stderr}")
            
            # Parse additional sequence-based results
            if blast_out.exists() and blast_out.stat().st_size > 0:
                try:
                    df_seq = pd.read_csv(blast_out, sep='\t', header=None,
                                       names=['qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore', 'stitle'])
                    
                    for _, row in df_seq.iterrows():
                        # Extract organism from title
                        organism_match = re.search(r'\[([^\]]+)\]', row['stitle'])
                        if organism_match:
                            organism = organism_match.group(1)
                            sequence_pathogens.append({
                                'organism': organism,
                                'identity': row['pident'],
                                'evalue': row['evalue'],
                                'bitscore': row['bitscore'],
                                'source': 'sequence_database',
                                'confidence': 'high' if row['pident'] > 80 else 'medium'
                            })
                except Exception as e:
                    print(f"Warning: Could not parse additional sequence results: {e}")
            
            print("✓ Additional pathogen screening completed")
            
        except subprocess.TimeoutExpired:
            print("Warning: Additional pathogen screening timed out after 30 minutes")
        except subprocess.CalledProcessError as e:
            print(f"Warning: Additional pathogen screening failed: {e}")
        except FileNotFoundError:
            print("Error: DIAMOND not found. Please install DIAMOND.")
    else:
        print("No additional custom pathogen database found")
    
    return blast_out