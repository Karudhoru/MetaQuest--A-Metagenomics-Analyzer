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
from .config import *
from .visualization import create_krona_plot
from .utils import check_dependencies

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
    
    # Create summary report
    summary_report = create_blast_taxonomy_summary(blast_results, output_dir)
    
    print(f"✓ BLAST taxonomy analysis completed")
    print(f"  Results saved to: {results_file}")
    print(f"  Summary saved to: {summary_report}")
    
    return results_file

def create_blast_taxonomy_summary(blast_results, output_dir):
    """Create summary report from BLAST taxonomy results"""
    
    summary_file = output_dir / "blast_taxonomy_summary.txt"
    organism_counts = Counter()
    total_hits = 0
    sequences_with_hits = 0
    
    # Analyze results
    for result in blast_results:
        if 'error' in result:
            continue
            
        if result['hits']:
            sequences_with_hits += 1
            for hit in result['hits']:
                organism = hit.get('organism', 'Unknown')
                organism_counts[organism] += 1
                total_hits += 1
    
    # Write summary
    with open(summary_file, 'w') as f:
        f.write("BLAST TAXONOMIC CLASSIFICATION SUMMARY\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total sequences analyzed: {len(blast_results)}\n")
        f.write(f"Sequences with hits: {sequences_with_hits}\n")
        f.write(f"Total BLAST hits: {total_hits}\n")
        f.write(f"Unique organisms identified: {len(organism_counts)}\n\n")
        
        if organism_counts:
            f.write("TOP ORGANISMS IDENTIFIED:\n")
            f.write("-" * 30 + "\n")
            
            for organism, count in organism_counts.most_common(20):
                f.write(f"{organism}: {count} hits\n")
        else:
            f.write("No significant taxonomic matches found.\n")
    
    # Create Kraken-style report for compatibility
    kraken_style_report = output_dir / "blast_kraken_style_report.txt"
    create_kraken_style_report_from_blast(blast_results, kraken_style_report)
    
    return summary_file


def create_kraken_style_report_from_blast(blast_results, output_file, organism_counts=None):
    """
    Create a Kraken-style report from BLAST results using consistent organism counts
    """
    
    # If organism_counts not provided, calculate them (for backward compatibility)
    if organism_counts is None:
        organism_counts = Counter()
        organism_sequences = {}
        
        for result in blast_results:
            if 'error' in result or not result['hits']:
                continue
                
            for hit in result['hits']:
                organism = hit.get('organism', 'Unknown')
                organism_counts[organism] += 1
                
                if organism not in organism_sequences:
                    organism_sequences[organism] = set()
                organism_sequences[organism].add(result['query_id'])
    else:
        # Recalculate organism_sequences for consistency
        organism_sequences = {}
        for result in blast_results:
            if 'error' in result or not result['hits']:
                continue
                
            for hit in result['hits']:
                organism = hit.get('organism', 'Unknown')
                if organism not in organism_sequences:
                    organism_sequences[organism] = set()
                organism_sequences[organism].add(result['query_id'])
    
    total_sequences = len([r for r in blast_results if 'error' not in r])
    classified_sequences = len([r for r in blast_results if 'error' not in r and r['hits']])
    total_hits = sum(organism_counts.values())
    
    # Write Kraken-style report using the same organism counts as summary
    with open(output_file, 'w') as f:
        # Header line for unclassified
        unclassified = total_sequences - classified_sequences
        unclassified_pct = (unclassified / total_sequences) * 100 if total_sequences > 0 else 0
        f.write(f"{unclassified_pct:.2f}\t{unclassified}\t{unclassified}\tU\t0\tunclassified\n")
        
        # Organism lines sorted by hit count (same as summary)
        for organism, hit_count in organism_counts.most_common():
            sequences_with_hits = len(organism_sequences.get(organism, set()))
            
            # Calculate percentage based on hit count relative to total hits
            percentage = (hit_count / total_hits) * 100 if total_hits > 0 else 0
            
            # Use hit count as the main count
            f.write(f"{percentage:.2f}\t{hit_count}\t{sequences_with_hits}\tS\t0\t{organism}\n")

def create_blast_taxonomy_summary(blast_results, output_dir):
    """Create comprehensive summary report from BLAST taxonomy results"""
    
    summary_file = output_dir / "blast_taxonomy_summary.txt"
    organism_counts = Counter()  # This will count total hits per organism
    organism_sequences = {}
    total_hits = 0
    sequences_with_hits = 0
    
    # Analyze results - count ALL hits to match Kraken-style report
    for result in blast_results:
        if 'error' in result:
            continue
            
        if result['hits']:
            sequences_with_hits += 1
            
            for hit in result['hits']:
                organism = hit.get('organism', 'Unknown')
                organism_counts[organism] += 1  # Count every hit
                total_hits += 1
                
                # Track which sequences hit each organism
                if organism not in organism_sequences:
                    organism_sequences[organism] = set()
                organism_sequences[organism].add(result['query_id'])
    
    # Write comprehensive summary
    with open(summary_file, 'w') as f:
        f.write("BLAST TAXONOMIC CLASSIFICATION SUMMARY\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total sequences analyzed: {len(blast_results)}\n")
        f.write(f"Sequences with hits: {sequences_with_hits}\n")
        f.write(f"Total BLAST hits: {total_hits}\n")
        f.write(f"Unique organisms identified: {len(organism_counts)}\n\n")
        
        if organism_counts:
            f.write("TOP ORGANISMS BY TOTAL HITS:\n")
            f.write("-" * 40 + "\n")
            f.write(f"{'Organism':<30} {'Total Hits':<10} {'Sequences':<10} {'Avg Hits/Seq':<12}\n")
            f.write("-" * 40 + "\n")
            
            for organism, hit_count in organism_counts.most_common(20):
                seq_count = len(organism_sequences.get(organism, set()))
                avg_hits = hit_count / seq_count if seq_count > 0 else 0
                f.write(f"{organism[:29]:<30} {hit_count:<10} {seq_count:<10} {avg_hits:<12.1f}\n")
        else:
            f.write("No significant taxonomic matches found.\n")
        
        # Add quality metrics
        f.write(f"\nQUALITY METRICS:\n")
        f.write("-" * 20 + "\n")
        f.write(f"Classification rate: {sequences_with_hits/len(blast_results)*100:.1f}%\n")
        f.write(f"Average hits per classified sequence: {total_hits/sequences_with_hits:.1f}\n" if sequences_with_hits > 0 else "")
    
    # Create Kraken-style report for compatibility - pass the same organism_counts
    kraken_style_report = output_dir / "blast_report.txt"
    create_kraken_style_report_from_blast(blast_results, kraken_style_report, organism_counts)
    
    # Create hit comparison data for visualization - use same counts
    create_hit_comparison_data(organism_counts, organism_sequences, output_dir)
    
    # Create Krona plot
    create_krona_plot(output_dir)
    
    return summary_file

def create_hit_comparison_data(organism_counts, organism_sequences, output_dir):
    """Create data file for organism hit comparison visualization using consistent counts"""
    
    comparison_data = []
    
    # Use the same organism_counts ranking as summary and Kraken report
    for organism, hit_count in organism_counts.most_common(10):  # Top 10
        seq_count = len(organism_sequences.get(organism, set()))
        avg_hits = hit_count / seq_count if seq_count > 0 else 0
        
        comparison_data.append({
            'organism': organism,
            'total_hits': hit_count,
            'sequences_with_hits': seq_count,
            'avg_hits_per_sequence': avg_hits
        })
    
    # Save as JSON for easy reading by visualization
    with open(output_dir / "organism_comparison_data.json", 'w') as f:
        json.dump(comparison_data, f, indent=2)
    
    # Also save as CSV
    df = pd.DataFrame(comparison_data)
    df.to_csv(output_dir / "organism_comparison_data.csv", index=False)
    
    return comparison_data