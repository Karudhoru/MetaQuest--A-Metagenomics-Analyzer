import subprocess
import os
import json
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from Bio import SeqIO
from typing import List, Optional
import re
from pathlib import Path
from collections import Counter
from ..config import *
from ..io.utils import parse_diamond_progress
from ..io.output_formatter import get_formatter


def run_pathogen_scan(protein_file, output_dir, bracken_results=None, taxonomy_results=None):
    """
    Enhanced pathogen scan using multiple data sources with professional output.
    
    *** REWRITTEN WITH TIGHTER RESTRICTIONS TO IMPROVE PRECISION ***
    
    Sources:
    1. Bracken results (FASTQ taxonomic classification)
    2. Taxonomy results (BLAST-based)
    3. Custom pathogen database (sequence-based with strict filters)
    """
    formatter = get_formatter()
    
    formatter.section_header("Pathogen Detection")
    
    # Initialize pathogen lists
    bracken_pathogens = []
    taxonomy_pathogens = []
    sequence_pathogens = []
    
    # ===== SOURCE 1: Bracken Results =====
    if bracken_results and bracken_results.exists():
        formatter.operation("Scanning Bracken taxonomic results", show_in_standard=True)
        try:
            df_bracken = pd.read_csv(bracken_results, sep='\t')
            
            pathogen_taxa = {
            # Existing entries
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
            '1719': 'Corynebacterium pseudotuberculosis',
            
            # *** CRITICAL FIX: Add MISSING pathogens ***
            '1496': 'Clostridioides difficile',  # THIS WAS THE KEY MISSING ONE!
            '1308': 'Streptococcus thermophilus',
            '1304': 'Streptococcus salivarius',
            '824': 'Campylobacter gracilis',
            '1301': 'Streptococcus',  # Genus level
            
            # Additional common pathogens for completeness
            '1313': 'Streptococcus pneumoniae',
            '1314': 'Streptococcus pyogenes',
            '1351': 'Enterococcus faecalis',
            '1352': 'Enterococcus faecium',
            '287': 'Pseudomonas aeruginosa',
            '1763': 'Mycobacterium',  # Genus
            '1773': 'Mycobacterium tuberculosis',
            '620': 'Shigella',  # Genus
            '666': 'Vibrio',  # Genus
            '197': 'Campylobacter',  # Genus
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
                # Name-based matching
                elif any(pathogen in name.lower() for pathogen in [
                    'brucella', 'salmonella', 'escherichia', 
                    'staphylococcus aureus', 'staphylococcus',
                    'bacillus cereus', 'klebsiella', 'acinetobacter', 
                    'clostridium', 'clostridioides',  # Both C. difficile names
                    'yersinia', 'pseudomonas aeruginosa', 
                    'mycobacterium tuberculosis', 'mycobacterium',
                    'listeria', 'vibrio', 'campylobacter', 'shigella',
                    # *** FIX: Make sure all Streptococcus variants are caught ***
                    'streptococcus thermophilus', 'streptococcus salivarius',
                    'streptococcus pneumoniae', 'streptococcus pyogenes',
                    'streptococcus',  # Catch-all for genus
                ]):
                    bracken_pathogens.append({
                        'organism': name,
                        'abundance': abundance,
                        'reads': reads,
                        'source': 'bracken_taxonomic',
                        'confidence': 'medium',
                        'taxonomy_id': taxid
                    })
            
            if bracken_pathogens:
                formatter.success(f"Found {len(bracken_pathogens)} pathogens in Bracken classification")
            else:
                formatter.info("No pathogens detected in Bracken results")
            
        except Exception as e:
            formatter.warning(f"Could not parse Bracken results: {e}")
    
    # ===== SOURCE 2: Taxonomy Results (BLAST) =====
    if taxonomy_results and taxonomy_results.exists():
        formatter.operation("Scanning BLAST taxonomy results", show_in_standard=False)
        try:
            with open(taxonomy_results, 'r') as f:
                for line in f:
                    if line.strip() and not line.startswith('#'):
                        parts = line.strip().split('\t')
                        if len(parts) >= 7:
                            query_id = parts[0]
                            subject_id = parts[1]
                            identity = float(parts[2])
                            length = int(parts[3])
                            evalue = float(parts[4])
                            bitscore = float(parts[5])
                            description = parts[6]
                            
                            # Extract organism
                            organism_matches = re.findall(r'\[([^\]]+)\]', description)
                            if organism_matches:
                                organism = organism_matches[-1]
                                
                                # Check if pathogen
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
            
            if taxonomy_pathogens:
                formatter.success(f"Found {len(taxonomy_pathogens)} pathogens in BLAST taxonomy results")
            else:
                formatter.info("No pathogens detected in BLAST results")
            
        except Exception as e:
            formatter.warning(f"Could not parse taxonomy results: {e}")
    
    # ===== SOURCE 3: Custom Pathogen Database =====
    
    # --- NEW: Define strict filter thresholds ---
    MIN_IDENTITY = 80.0  # Require 80% identity
    MIN_COVERAGE = 0.8   # Require 80% coverage (alignment_length / query_length)
    # ---
    
    selected_db = PATHOGEN_DB_CUSTOM
    blast_out = output_dir / "pathogen_results.txt"

    if selected_db and selected_db.exists():
        try:
            query_count = len(list(SeqIO.parse(protein_file, "fasta")))
            formatter.info(f"Scanning {query_count:,} sequences against pathogen database...")
            formatter.info(f"Applying strict filters: Identity >= {MIN_IDENTITY}%, Coverage >= {MIN_COVERAGE*100}%")
        except Exception:
            formatter.warning("Could not count input sequences for progress bar.")
            query_count = 1

        db_base = selected_db.with_suffix('') if selected_db.suffix == '.dmnd' else selected_db
        
        # --- MODIFIED: Add qlen to output format ---
        outfmt_cols = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle qlen"
        
        cmd = [
            "diamond", "blastp",
            "-d", str(db_base),
            "-q", str(protein_file),
            "-o", str(blast_out),
            "--outfmt", "6", *outfmt_cols.split(),
            "--top", "5",        # Keep top 5 to allow Python to do the fine-filtering
            "--evalue", "1e-3", # Lenient e-value, Python will filter
            "--threads", "4",
            "--more-sensitive",
            "--log"  # CRITICAL: Enable progress output
        ]

        return_code = formatter.run_subprocess_with_progress(
            cmd=cmd,
            operation_name="Scanning pathogen database",
            total=query_count,
            unit="sequences",
            parser_func=parse_diamond_progress
        )

        # Parse results after the run is complete
        if return_code == 0 and blast_out.exists() and blast_out.stat().st_size > 0:
            sequence_pathogens = []
            hits_by_gene = {}  # Key: query_id, Value: best hit info
            
            try:
                with open(blast_out, 'r') as f:
                    for line in f:
                        if line.strip() and not line.startswith('#'):
                            parts = line.strip().split('\t')
                            
                            # --- MODIFIED: Check for 14 columns ---
                            if len(parts) >= 14:
                                query_id = parts[0]
                                subject_id = parts[1]
                                identity = float(parts[2])
                                align_length = int(parts[3])
                                evalue = float(parts[10])
                                bitscore = float(parts[11])
                                description = parts[12]
                                query_len = int(parts[13]) # <-- NEW: Get query length
                                
                                # --- NEW: Calculate coverage ---
                                if query_len == 0:
                                    continue # Avoid division by zero
                                coverage = align_length / query_len
                                
                                # --- NEW: Apply strict filters ---
                                if identity < MIN_IDENTITY or coverage < MIN_COVERAGE:
                                    continue # Skip this low-quality hit
                                
                                # Extract organism from description [brackets]
                                organism_matches = re.findall(r'\[([^\]]+)\]', description)
                                organism = organism_matches[-1] if organism_matches else "Unknown"
                                
                                # Keep only BEST hit per gene (lowest e-value) that passes filters
                                if query_id not in hits_by_gene:
                                    hits_by_gene[query_id] = {
                                        'organism': organism,
                                        'identity': identity,
                                        'evalue': evalue,
                                        'bitscore': bitscore,
                                        'length': align_length,
                                        'coverage': coverage, # Store for reference
                                        'query_id': query_id,
                                        'subject_id': subject_id,
                                        'description': description,
                                        'source': 'pathogen_database',
                                        'confidence': 'high' # Passed strict filter
                                    }
                                else:
                                    # Keep hit with better e-value (lower is better)
                                    if evalue < hits_by_gene[query_id]['evalue']:
                                        hits_by_gene[query_id] = {
                                            'organism': organism,
                                            'identity': identity,
                                            'evalue': evalue,
                                            'bitscore': bitscore,
                                            'length': align_length,
                                            'coverage': coverage,
                                            'query_id': query_id,
                                            'subject_id': subject_id,
                                            'description': description,
                                            'source': 'pathogen_database',
                                            'confidence': 'high'
                                        }
                
                # Convert to list (now deduplicated - one hit per gene)
                sequence_pathogens = list(hits_by_gene.values())
                
                if sequence_pathogens:
                    formatter.success(
                        f"Found {len(sequence_pathogens)} unique genes with high-confidence pathogen hits "
                        f"(passed {MIN_IDENTITY}% ident, {MIN_COVERAGE*100}% cov filters)"
                    )
                else:
                    formatter.info("No high-confidence pathogen hits found in database search")
                    
            except Exception as e:
                formatter.warning(f"Could not parse pathogen scan results: {e}")
    
    # ===== Summary =====
    all_pathogens = bracken_pathogens + taxonomy_pathogens + sequence_pathogens
    total_pathogens = len(all_pathogens)
    results_file = output_dir / "pathogen_detections.json"

    formatter.summary("Pathogen Detection Summary", {
        'Bracken detections': len(bracken_pathogens),
        'BLAST detections': len(taxonomy_pathogens),
        'Database detections (strict)': len(sequence_pathogens),
        'Total detections': total_pathogens
    })

    with open(results_file, 'w') as f:
        json.dump({
            'summary': {
                'total_detections': total_pathogens,
                'bracken': len(bracken_pathogens),
                'blast': len(taxonomy_pathogens),
                'database_strict': len(sequence_pathogens)
            },
            'detections': all_pathogens,
            'bracken_pathogens': bracken_pathogens,
            'taxonomy_pathogens': taxonomy_pathogens,
            'sequence_pathogens_strict': sequence_pathogens
        }, f, indent=2)

    if total_pathogens > 0:
        formatter.info(f"Detailed results saved: {results_file.name}")
    else:
        formatter.info(f"No pathogens detected - empty results saved: {results_file.name}")

    # === INTEGRATION: Merge sequence-based detections with taxonomy ===
    if bracken_results and bracken_results.exists() and blast_out.exists():
        formatter.operation("Integrating pathogen hits with taxonomic classification")
        
        integrated_file = output_dir / "integrated_pathogen_taxonomy.tsv"
        integrated_df = integrate_pathogen_hits_with_taxonomy(
            bracken_results,
            blast_out,  # Pass the path to the raw DIAMOND output
            integrated_file
        )
        
        formatter.info(f"Integrated taxonomy saved: {integrated_file.name}")


    return blast_out

def integrate_pathogen_hits_with_taxonomy(
    bracken_file: Path,
    pathogen_hits_file: Path,
    output_file: Path
) -> pd.DataFrame:
    """
    Integrate sequence-based pathogen hits with Bracken taxonomic classification.
    
    *** MODIFIED to parse 14 columns from DIAMOND output ***
    
    Args:
        bracken_file: Path to bracken_report.tsv
        pathogen_hits_file: Path to pathogen_results.txt (raw DIAMOND output)
        output_file: Path to save integrated results
    
    Returns:
        DataFrame with integrated pathogen detections
    """
    formatter = get_formatter()
    
    # Load Bracken results
    df_bracken = pd.read_csv(bracken_file, sep='\t')
    
    # Parse pathogen hits to extract organisms
    pathogen_organisms = {}
    hit_counts = Counter()
    
    if pathogen_hits_file.exists():
        with open(pathogen_hits_file, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    parts = line.strip().split('\t')
                    
                    # --- MODIFIED: Check for 14 columns ---
                    if len(parts) >= 14: 
                        query_id = parts[0]
                        description = parts[12]
                        
                        # Extract organism name from [brackets]
                        organism_matches = re.findall(r'\[([^\]]+)\]', description)
                        if organism_matches:
                            organism = organism_matches[-1]
                            hit_counts[organism] += 1
                            
                            # Store best hit info (not strictly needed here, but good practice)
                            if organism not in pathogen_organisms:
                                pathogen_organisms[organism] = {
                                    'identity': float(parts[2]),
                                    'evalue': float(parts[10]),
                                    'hit_count': 1
                                }
    
    # Find organisms detected by pathogen DB but NOT in Bracken
    bracken_organisms = set(df_bracken['name'].values)
    new_pathogens = []
    
    for organism, info in pathogen_organisms.items():
        # Check if organism (or genus) already in Bracken
        already_present = False
        for bracken_org in bracken_organisms:
            if organism.lower() in bracken_org.lower() or bracken_org.lower() in organism.lower():
                already_present = True
                break
        
        # --- MODIFIED: Use hit_counts (which is not filtered) ---
        if not already_present and hit_counts[organism] >= 3:  # At least 3 hits
            # Add as low-abundance detection
            # Estimate abundance based on hit count (rough approximation)
            estimated_abundance = min(0.01, hit_counts[organism] / 1000)  # Max 1%
            
            new_row = {
                'name': organism,
                'taxonomy_id': '0',  # Unknown taxonomy ID
                'taxonomy_lvl': 'S',
                'kraken_assigned_reads': 0,
                'added_reads': hit_counts[organism],
                'new_est_reads': hit_counts[organism],
                'fraction_total_reads': estimated_abundance,
                'source': 'pathogen_database'
            }
            new_pathogens.append(new_row)
            formatter.info(
                f"Adding pathogen detected by sequence search: {organism} "
                f"({hit_counts[organism]} hits, ~{estimated_abundance*100:.2f}%)"
            )
    
    # Append new pathogens to Bracken results
    if new_pathogens:
        df_integrated = pd.concat([
            df_bracken,
            pd.DataFrame(new_pathogens)
        ], ignore_index=True)
        
        formatter.success(
            f"Integrated {len(new_pathogens)} additional pathogens from sequence search"
        )
    else:
        df_integrated = df_bracken
        formatter.info("No additional pathogens to integrate")
    
    # Save integrated results
    df_integrated.to_csv(output_file, sep='\t', index=False)
    
    return df_integrated