import subprocess
import os
import pandas as pd
from pathlib import Path
from collections import Counter
import time
import json
import hashlib
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
import requests
from ..config import *

def run_kraken(input_files, output_dir):
    """Run Kraken2 classification for FASTQ files"""
    report = output_dir/"kraken_report.txt"
    classified = output_dir/"kraken_classified.txt"
    # input_files might be a list of one (single-end) or two paths (paired-end)
    if isinstance(input_files, (list,tuple)) and len(input_files)==2:
        reads_flags = f"--paired {input_files[0]} {input_files[1]}"
    else:
        reads_flags = input_files[0] if isinstance(input_files,(list,tuple)) else input_files
    cmd = f"kraken2 --db {KRAKEN_DB} --threads 8 --report {report} --output {classified} {reads_flags}"

    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    return report

def run_bracken(report_path, output_dir, is_fasta=False):
    """Estimate abundances with Bracken, with FASTA mode handling"""
    bracken_out = output_dir / "bracken_report.tsv"
    
    # Skip Bracken for FASTA files as it's read-based, not contig-based
    if is_fasta:
        print("Skipping Bracken for FASTA input (using BLAST API for taxonomic analysis)")
        return report_path
    
    # Original Bracken code for FASTQ files
    cmd = f"bracken -d {KRAKEN_DB} -i {report_path} -o {bracken_out} -w {output_dir}/bracken_report.txt -r 150 -l S -t 10"
    print(f"Running: {cmd}")
    
    try:
        subprocess.run(cmd, shell=True, check=True)
        print("✓ Bracken abundance estimation completed")
        return bracken_out
    except subprocess.CalledProcessError as e:
        print(f"Warning: Bracken failed: {e}")
        print("Continuing with Kraken2 report...")
        return report_path

def get_sequence_cache_key(sequence):
    """Generate cache key for sequence"""
    return hashlib.md5(sequence.encode()).hexdigest()

def load_blast_cache(cache_dir):
    """Load existing BLAST cache"""
    cache_file = cache_dir / "blast_cache.json"
    if cache_file.exists():
        try:
            with open(cache_file, 'r') as f:
                return json.load(f)
        except Exception as e:
            print(f"Warning: Could not load cache: {e}")
    return {}

def save_blast_cache(cache_dir, cache_data):
    """Save BLAST cache to disk"""
    cache_file = cache_dir / "blast_cache.json"
    try:
        with open(cache_file, 'w') as f:
            json.dump(cache_data, f, indent=2)
    except Exception as e:
        print(f"Warning: Could not save cache: {e}")

def blast_sequence_online(sequence, sequence_id, database="nt", cache=None, cache_key=None):
    """
    BLAST a single sequence against NCBI database with caching and rate limiting
    """
    # Check cache first
    if cache and cache_key and cache_key in cache:
        print(f"  Using cached result for {sequence_id}")
        return cache[cache_key]
    
    # Rate limiting: NCBI recommends no more than 3 requests per second
    time.sleep(0.4)  # 400ms delay between requests
    
    try:
        print(f"  BLASTing {sequence_id} against {database}...")
        
        # Submit BLAST job
        result_handle = NCBIWWW.qblast(
            program="blastn" if database == "nt" else "blastx",
            database=database,
            sequence=sequence,
            hitlist_size=10,  # Top 10 hits
            expect=1e-5,
            word_size=28 if database == "nt" else 6
        )
        
        # Parse results
        blast_records = NCBIXML.parse(result_handle)
        record = next(blast_records)
        
        # Extract taxonomy information
        taxonomy_results = []
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect <= 1e-5:  # Only significant hits
                    # Extract taxonomy info from hit description
                    hit_info = {
                        'hit_id': alignment.hit_id,
                        'hit_def': alignment.hit_def,
                        'length': alignment.length,
                        'e_value': hsp.expect,
                        'bit_score': hsp.bits,
                        'identity': hsp.identities / hsp.align_length * 100,
                        'query_cover': (hsp.query_end - hsp.query_start + 1) / len(sequence) * 100
                    }
                    
                    # Try to extract organism name from description
                    organism = extract_organism_from_description(alignment.hit_def)
                    hit_info['organism'] = organism
                    
                    taxonomy_results.append(hit_info)
        
        result_data = {
            'query_id': sequence_id,
            'query_length': len(sequence),
            'hits': taxonomy_results,
            'timestamp': time.time()
        }
        
        # Cache the result
        if cache is not None and cache_key:
            cache[cache_key] = result_data
        
        return result_data
        
    except Exception as e:
        print(f"  BLAST failed for {sequence_id}: {e}")
        return {
            'query_id': sequence_id,
            'query_length': len(sequence),
            'hits': [],
            'error': str(e),
            'timestamp': time.time()
        }

def extract_organism_from_description(description):
    """Extract organism name from BLAST hit description"""
    # Common patterns in NCBI descriptions
    if '[' in description and ']' in description:
        # Extract text between last brackets (usually organism)
        organism = description.split('[')[-1].replace(']', '').strip()
        return organism
    
    # If no brackets, try to extract from common patterns
    desc_lower = description.lower()
    
    # Look for species-like patterns
    words = description.split()
    for i, word in enumerate(words):
        if i < len(words) - 1:
            # Look for genus species pattern (capitalized followed by lowercase)
            if (word[0].isupper() and word[1:].islower() and 
                words[i+1][0].islower() and words[i+1].isalpha()):
                return f"{word} {words[i+1]}"
    
    # Fallback: return first few words
    return ' '.join(description.split()[:3])

def run_fasta_blast_taxonomy(fasta_path, output_dir, database="nt", max_sequences=100):
    """
    Run BLAST taxonomy classification for FASTA files using NCBI API
    """
    print(f"Running BLAST taxonomic classification on {fasta_path}")
    print(f"Using database: {database}")
    
    # Setup cache
    cache_dir = output_dir / "blast_cache"
    cache_dir.mkdir(exist_ok=True)
    cache = load_blast_cache(cache_dir)
    
    # Parse FASTA file
    sequences = list(SeqIO.parse(fasta_path, "fasta"))
    total_sequences = len(sequences)
    
    print(f"Found {total_sequences} sequences")
    
    # Limit sequences to avoid overwhelming the API
    if total_sequences > max_sequences:
        print(f"Limiting to first {max_sequences} sequences to avoid API overload")
        sequences = sequences[:max_sequences]
    
    blast_results = []
    
    print("Starting BLAST analysis (this may take several minutes)...")
    
    for i, seq_record in enumerate(sequences):
        print(f"Processing sequence {i+1}/{len(sequences)}: {seq_record.id}")
        
        sequence_str = str(seq_record.seq)
        cache_key = get_sequence_cache_key(sequence_str)
        
        # Skip very short sequences
        if len(sequence_str) < 50:
            print(f"  Skipping {seq_record.id} (too short: {len(sequence_str)} bp)")
            continue
        
        # BLAST the sequence
        result = blast_sequence_online(
            sequence_str, 
            seq_record.id, 
            database=database,
            cache=cache,
            cache_key=cache_key
        )
        
        blast_results.append(result)
        
        # Save cache periodically
        if i % 10 == 0:
            save_blast_cache(cache_dir, cache)
    
    # Final cache save
    save_blast_cache(cache_dir, cache)
    
    # Save raw results
    results_file = output_dir / "blast_taxonomy_results.json"
    with open(results_file, 'w') as f:
        json.dump(blast_results, f, indent=2)
    
    return results_file
