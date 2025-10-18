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
from ..io.output_formatter import get_formatter


def run_kraken(input_files, output_dir):
    """Run Kraken2 classification for FASTQ files with professional output"""
    formatter = get_formatter()
    
    report = output_dir / "kraken_report.txt"
    classified = output_dir / "kraken_classified.txt"
    
    # Handle single-end vs paired-end
    if isinstance(input_files, (list, tuple)) and len(input_files) == 2:
        reads_flags = f"--paired {input_files[0]} {input_files[1]}"
        formatter.info("Running Kraken2 on paired-end reads")
    else:
        reads_flags = input_files[0] if isinstance(input_files, (list, tuple)) else input_files
        formatter.info("Running Kraken2 on single-end reads")
    
    cmd = [
        "kraken2",
        "--db", str(KRAKEN_DB),
        "--threads", "8",
        "--report", str(report),
        "--output", str(classified)
    ] + reads_flags.split()
    
    formatter.debug(f"Kraken2 database: {KRAKEN_DB}")
    
    # Run with smart output handling
    returncode, stdout, stderr = formatter.run_subprocess(
        cmd,
        operation_name="Kraken2 taxonomic classification",
        show_command=True
    )
    
    if returncode != 0:
        formatter.error(
            "Kraken2 classification failed",
            solutions=[
                "Check that Kraken2 database exists and is properly built",
                f"Verify database path: {KRAKEN_DB}",
                "Ensure sufficient memory (Kraken2 requires ~8GB+ RAM)",
                "Check input file format and integrity"
            ]
        )
        raise RuntimeError(f"Kraken2 failed with exit code {returncode}")
    
    # Parse and display key metrics
    if report.exists():
        try:
            # Count classified sequences
            with open(classified, 'r') as f:
                total_seqs = sum(1 for _ in f)
            
            formatter.success("Kraken2 classification complete")
            formatter.result({
                'Classified sequences': f"{total_seqs:,}",
                'Report file': str(report.name),
                'Output file': str(classified.name)
            })
        except Exception as e:
            formatter.debug(f"Could not parse metrics: {e}")
    
    return report


def run_bracken(report_path, output_dir, is_fasta=False):
    """Estimate abundances with Bracken with professional output"""
    formatter = get_formatter()
    
    bracken_out = output_dir / "bracken_report.tsv"
    
    # Skip Bracken for FASTA files
    if is_fasta:
        formatter.info("Skipping Bracken for FASTA input (read-based tool)")
        formatter.info("Using BLAST-based taxonomic analysis instead")
        return report_path
    
    # Wait for Kraken report file
    report_file = Path(report_path)
    formatter.operation("Waiting for Kraken2 report file", show_in_standard=True)
    
    for attempt in range(5):
        if report_file.exists() and report_file.stat().st_size > 0:
            formatter.debug(f"Found Kraken2 report: {report_file}")
            break
        time.sleep(1)
    else:
        formatter.error(
            f"Kraken2 report file not found: {report_file}",
            solutions=[
                "Ensure Kraken2 completed successfully",
                "Check output directory permissions",
                "Verify sufficient disk space"
            ]
        )
        raise FileNotFoundError(f"Kraken2 report missing: {report_file}")
    
    # Build Bracken command
    cmd = [
        "bracken",
        "-d", str(KRAKEN_DB),
        "-i", str(report_path),
        "-o", str(bracken_out),
        "-w", str(output_dir / "bracken_report.txt"),
        "-r", "150",
        "-l", "S",
        "-t", "10"
    ]
    
    formatter.debug(f"Bracken parameters: read_len=150, level=Species, threshold=10")
    
    # Run with smart output handling
    returncode, stdout, stderr = formatter.run_subprocess(
        cmd,
        operation_name="Bracken abundance estimation",
        show_command=True
    )
    
    if returncode != 0:
        formatter.warning("Bracken failed, continuing with Kraken2 report")
        formatter.debug(f"Bracken stderr: {stderr}")
        return report_path
    
    # Parse and display results
    if bracken_out.exists():
        try:
            df = pd.read_csv(bracken_out, sep='\t')
            species_count = len(df)
            total_reads = df['new_est_reads'].sum() if 'new_est_reads' in df.columns else 0
            
            formatter.success("Bracken abundance estimation complete")
            formatter.result({
                'Species detected': f"{species_count:,}",
                'Total reads classified': f"{int(total_reads):,}",
                'Output file': str(bracken_out.name)
            })
        except Exception as e:
            formatter.debug(f"Could not parse Bracken metrics: {e}")
    
    return bracken_out


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
            formatter = get_formatter()
            formatter.debug(f"Could not load cache: {e}")
    return {}


def save_blast_cache(cache_dir, cache_data):
    """Save BLAST cache to disk"""
    cache_file = cache_dir / "blast_cache.json"
    try:
        with open(cache_file, 'w') as f:
            json.dump(cache_data, f, indent=2)
    except Exception as e:
        formatter = get_formatter()
        formatter.debug(f"Could not save cache: {e}")


def blast_sequence_online(sequence, sequence_id, database="nt", cache=None, cache_key=None):
    """BLAST a single sequence against NCBI database with caching and rate limiting"""
    formatter = get_formatter()
    
    # Check cache first
    if cache and cache_key and cache_key in cache:
        formatter.debug(f"Using cached result for {sequence_id}")
        return cache[cache_key]
    
    # Rate limiting
    time.sleep(0.4)
    
    try:
        formatter.debug(f"BLASTing {sequence_id} against {database}")
        
        # Submit BLAST job
        result_handle = NCBIWWW.qblast(
            program="blastn" if database == "nt" else "blastx",
            database=database,
            sequence=sequence,
            hitlist_size=10,
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
                if hsp.expect <= 1e-5:
                    organism = extract_organism_from_description(alignment.hit_def)
                    
                    hit_info = {
                        'hit_id': alignment.hit_id,
                        'hit_def': alignment.hit_def,
                        'length': alignment.length,
                        'e_value': hsp.expect,
                        'bit_score': hsp.bits,
                        'identity': hsp.identities / hsp.align_length * 100,
                        'query_cover': (hsp.query_end - hsp.query_start + 1) / len(sequence) * 100,
                        'organism': organism
                    }
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
        formatter.warning(f"BLAST failed for {sequence_id}: {e}")
        return {
            'query_id': sequence_id,
            'query_length': len(sequence),
            'hits': [],
            'error': str(e),
            'timestamp': time.time()
        }


def extract_organism_from_description(description):
    """Extract organism name from BLAST hit description"""
    # Extract text between brackets
    if '[' in description and ']' in description:
        organism = description.split('[')[-1].replace(']', '').strip()
        return organism
    
    # Look for genus species pattern
    words = description.split()
    for i, word in enumerate(words):
        if i < len(words) - 1:
            if (word[0].isupper() and word[1:].islower() and 
                words[i+1][0].islower() and words[i+1].isalpha()):
                return f"{word} {words[i+1]}"
    
    # Fallback
    return ' '.join(description.split()[:3])


def run_fasta_blast_taxonomy(fasta_path, output_dir, database="nt", max_sequences=100):
    """Run BLAST taxonomy classification for FASTA files with professional output"""
    formatter = get_formatter()
    
    formatter.section_header("BLAST Taxonomic Classification")
    formatter.info(f"Database: {database}")
    formatter.info(f"Input: {fasta_path.name}")
    
    # Setup cache
    cache_dir = output_dir / "blast_cache"
    cache_dir.mkdir(exist_ok=True)
    cache = load_blast_cache(cache_dir)
    
    # Parse FASTA file
    sequences = list(SeqIO.parse(fasta_path, "fasta"))
    total_sequences = len(sequences)
    
    formatter.result({
        'Total sequences': f"{total_sequences:,}",
        'Cache directory': str(cache_dir.name)
    })
    
    # Limit sequences
    if total_sequences > max_sequences:
        formatter.warning(f"Limiting to first {max_sequences} sequences to avoid API overload")
        sequences = sequences[:max_sequences]
    
    blast_results = []
    
    formatter.operation("Starting BLAST analysis (this may take several minutes)")
    
    # Process sequences with progress
    for i, seq_record in enumerate(sequences):
        # Show progress
        if formatter.verbosity >= formatter.MINIMAL:
            formatter.progress_bar(i, len(sequences), prefix="BLAST Progress", suffix=f"{seq_record.id}")
        
        sequence_str = str(seq_record.seq)
        cache_key = get_sequence_cache_key(sequence_str)
        
        # Skip very short sequences
        if len(sequence_str) < 50:
            formatter.debug(f"Skipping {seq_record.id} (too short: {len(sequence_str)} bp)")
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
    
    # Final progress update
    if formatter.verbosity >= formatter.MINIMAL:
        formatter.progress_bar(len(sequences), len(sequences), prefix="BLAST Progress", suffix="Complete")
    
    # Final cache save
    save_blast_cache(cache_dir, cache)
    
    # Save results
    results_file = output_dir / "blast_taxonomy_results.json"
    with open(results_file, 'w') as f:
        json.dump(blast_results, f, indent=2)
    
    # Summary
    successful_blasts = len([r for r in blast_results if r.get('hits')])
    total_hits = sum(len(r.get('hits', [])) for r in blast_results)
    
    formatter.success("BLAST taxonomic classification complete")
    formatter.result({
        'Sequences analyzed': f"{len(blast_results):,}",
        'Successful BLASTs': f"{successful_blasts:,}",
        'Total hits found': f"{total_hits:,}",
        'Results file': str(results_file.name)
    })
    
    return results_file