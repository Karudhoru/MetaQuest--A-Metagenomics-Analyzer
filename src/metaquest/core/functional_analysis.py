# functional_analysis.py

import subprocess
import os
import pandas as pd
from pathlib import Path
from Bio import SeqIO
import threading
from typing import Optional
import time
import re
import json
from collections import Counter
import psutil
from ..io.utils import parse_diamond_progress
from ..config import SWISSPROT_DB
from ..io.output_formatter import get_formatter

def parse_prokka_sample_txt(sample_file: Path) -> dict:
    """
    Parse Prokka's sample.txt file for accurate gene counts.
    
    This reads the stats directly from Prokka rather than recounting,
    ensuring consistency across the pipeline.
    """
    stats = {
        'organism': '',
        'contigs': 0,
        'bases': 0,
        'CDS': 0,      # ← Uppercase for display
        'cds': 0,      # ← Lowercase for compatibility
        'gene': 0,
        'rRNA': 0,
        'tRNA': 0,
        'repeat_region': 0
    }
    
    if not sample_file.exists():
        return stats
    
    with open(sample_file) as f:
        for line in f:
            line = line.strip()
            if ':' in line:
                key, value = line.split(':', 1)
                key = key.strip()
                value = value.strip()
                
                # Case-insensitive matching
                key_lower = key.lower()
                
                if key_lower == 'organism':
                    stats['organism'] = value
                elif key_lower in ['cds', 'gene', 'rrna', 'trna', 'contigs', 'bases', 'repeat_region']:
                    try:
                        # *** FIX: Store with proper capitalization ***
                        if key_lower == 'cds':
                            stats['CDS'] = int(value)
                            stats['cds'] = int(value)
                        elif key_lower == 'rrna':
                            stats['rRNA'] = int(value)
                        elif key_lower == 'trna':
                            stats['tRNA'] = int(value)
                        elif key_lower == 'gene':
                            stats['gene'] = int(value)
                        else:
                            stats[key_lower] = int(value)
                    except ValueError:
                        pass
    
    return stats

def filter_contigs(input_fasta: Path, output_fasta: Path, min_length: int = 1000) -> dict:
    formatter = get_formatter()
    formatter.operation(f"Filtering contigs (minimum length: {min_length} bp)...")
    total_contigs = 0
    filtered_contigs = 0
    total_bases_before = 0
    total_bases_after = 0

    with open(output_fasta, 'w') as out_handle:
        for record in SeqIO.parse(input_fasta, 'fasta'):
            total_contigs += 1
            total_bases_before += len(record.seq)
            if len(record.seq) >= min_length:
                SeqIO.write(record, out_handle, 'fasta')
                filtered_contigs += 1
                total_bases_after += len(record.seq)

    stats = {
        'total_contigs': total_contigs,
        'filtered_contigs': filtered_contigs,
        'removed_contigs': total_contigs - filtered_contigs,
        'total_bases_before': total_bases_before,
        'total_bases_after': total_bases_after,
        'filter_percentage': (filtered_contigs / total_contigs * 100) if total_contigs > 0 else 0
    }
    formatter.result(stats)
    if filtered_contigs > 0:
        formatter.success(f"Filtering complete: {filtered_contigs:,} contigs retained, {stats['removed_contigs']:,} removed")
    elif stats['removed_contigs'] > 0:
        formatter.warning(f"All contigs were below {min_length} bp, using original file")
    else:
        formatter.info(f"No contigs needed filtering (all ≥{min_length} bp)")
    
    # Log filtering details for report consistency
    formatter.info(
        f"Original assembly: {total_contigs:,} contigs, {total_bases_before:,} bp"
    )
    if filtered_contigs < total_contigs:
        formatter.info(
            f"After filtering: {filtered_contigs:,} contigs, {total_bases_after:,} bp"
        )
    
    stats = {
        'total_contigs': total_contigs,
        'filtered_contigs': filtered_contigs,
        'removed_contigs': total_contigs - filtered_contigs,
        'total_bases_before': total_bases_before,
        'total_bases_after': total_bases_after,
        'filter_percentage': (filtered_contigs / total_contigs * 100) if total_contigs > 0 else 0
    }
    
    return stats

def kill_tbl2asn_monitor(max_wait_seconds: int = 30):
    formatter = get_formatter()
    start_time = time.time()
    while time.time() - start_time < max_wait_seconds:
        for proc in psutil.process_iter(['pid', 'name', 'cmdline', 'create_time']):
            try:
                proc_info = proc.info
                if proc_info['name'] and 'tbl2asn' in proc_info['name'].lower():
                    proc_age = time.time() - proc_info['create_time']
                    if proc_age > max_wait_seconds:
                        formatter.warning(f"Terminating tbl2asn (PID: {proc_info['pid']}) - exceeded {max_wait_seconds}s timeout")
                        proc.terminate()
                        try:
                            proc.wait(timeout=5)
                            formatter.success(f"tbl2asn terminated gracefully")
                        except psutil.TimeoutExpired:
                            formatter.warning("Force killing tbl2asn...")
                            proc.kill()
                            formatter.success("tbl2asn force killed")
                        return True
            except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
                pass
        time.sleep(2)
    return False

def run_prokka(
    fasta_path: Path,
    output_dir: Path,
    filter_contigs_flag: bool = True,
    min_contig_length: int = 500,
    kill_tbl2asn: bool = True,
    tbl2asn_timeout: int = 30
) -> tuple[Path, dict]:
    formatter = get_formatter()
    prokka_dir = output_dir / "prokka_annotation"
    working_fasta = fasta_path

    if filter_contigs_flag:
        filtered_fasta = output_dir / f"{fasta_path.stem}_filtered.fasta"
        filter_stats = filter_contigs(fasta_path, filtered_fasta, min_contig_length)
        sample_txt = output_dir / "sample.txt"
        if sample_txt.exists():
        # Update bases count in sample.txt
            with open(sample_txt, 'a') as f:
                f.write(f"filtered_bases: {filter_stats['total_bases_after']}\n")
                f.write(f"removed_contigs: {filter_stats['removed_contigs']}\n")
        
        # FIX: Store filter stats for later reference
        filter_log = output_dir / "contig_filtering.log"
        with open(filter_log, 'w') as f:
            f.write(f"Original contigs: {filter_stats['total_contigs']}\n")
            f.write(f"Filtered contigs: {filter_stats['filtered_contigs']}\n")
            f.write(f"Removed contigs: {filter_stats['removed_contigs']}\n")
            f.write(f"Original bases: {filter_stats['total_bases_before']}\n")
            f.write(f"Final bases: {filter_stats['total_bases_after']}\n")
        
        if filter_stats["filtered_contigs"] > 0:
            working_fasta = filtered_fasta
            formatter.info(f"Using filtered assembly: {filter_stats['filtered_contigs']} contigs")
        else:
            working_fasta = fasta_path
            formatter.warning("All contigs below threshold, using original file")
    else:
        working_fasta = fasta_path

    cmd = [
        "prokka",
        "--outdir", str(prokka_dir),
        "--prefix", "sample",
        "--cpus", "12",
        "--force",
        "--metagenome",
        "--centre", "X",
        "--compliant",
        "--fast", str(working_fasta)
    ]
    formatter.section_header("Gene Prediction & Annotation")
    formatter.operation("Running Prokka", show_in_standard=True)
    formatter.debug(f"Command: {' '.join(cmd)}")

    killer_thread = None
    if kill_tbl2asn:
        killer_thread = threading.Thread(target=kill_tbl2asn_monitor, args=(tbl2asn_timeout,), daemon=True)
        killer_thread.start()

    returncode, stdout, stderr = formatter.run_subprocess(cmd, "Prokka gene prediction", capture_output=True, show_command=True)

    essential_files = [
        prokka_dir / "sample.faa",
        prokka_dir / "sample.ffn",
        prokka_dir / "sample.gff"
    ]
    all_essential_exist = all(f.exists() for f in essential_files)
    stats = {
        'gene_count': 0,
        'protein_count': 0,
        'cds_count': 0,
        'success': False
    }
    
    if all_essential_exist:
        formatter.success("Prokka annotation complete (essential files generated)")
        
        # FIX: Read gene counts from Prokka's sample.txt instead of counting .faa
        sample_txt = prokka_dir /  "sample.txt"
        prokka_stats = parse_prokka_sample_txt(sample_txt)
        
        formatter.result({
            "Total genes": prokka_stats.get('gene', 0),
            "CDS": prokka_stats.get('CDS', 0),
            "rRNA": prokka_stats.get('rRNA', 0),
            "tRNA": prokka_stats.get('tRNA', 0)
        })
        
        stats.update({
            'gene_count': prokka_stats.get('gene', 0),
            'protein_count': prokka_stats.get('CDS', 0),
            'cds_count': prokka_stats.get('CDS', 0),
            'rrna_count': prokka_stats.get('rRNA', 0),
            'trna_count': prokka_stats.get('tRNA', 0),
            'success': True
        })

        sample_txt = prokka_dir / "sample.txt"
    filter_log = output_dir / "contig_filtering.log"
    main_sample_txt = output_dir / "sample.txt"
    
    # Merge Prokka stats with filter stats
    if sample_txt.exists():
        # Read Prokka's sample.txt (has gene counts)
        prokka_stats = {}
        with open(sample_txt, 'r') as f:
            for line in f:
                if ':' in line:
                    key, value = line.strip().split(':', 1)
                    prokka_stats[key.strip()] = value.strip()
        
        # Read filter stats if available
        filter_stats = {}
        if filter_log.exists():
            with open(filter_log, 'r') as f:
                for line in f:
                    if 'Final bases:' in line:
                        filter_stats['bases'] = int(line.split(':')[1].strip())
                    elif 'Filtered contigs:' in line:
                        filter_stats['contigs'] = int(line.split(':')[1].strip())
        
        # Write merged sample.txt to main output directory
        with open(main_sample_txt, 'w') as f:
            # Write Prokka stats first (they have gene counts)
            for key, value in prokka_stats.items():
                # Override bases/contigs with filtered values if available
                if key == 'bases' and 'bases' in filter_stats:
                    f.write(f"{key}: {filter_stats['bases']}\n")
                elif key == 'contigs' and 'contigs' in filter_stats:
                    f.write(f"{key}: {filter_stats['contigs']}\n")
                else:
                    f.write(f"{key}: {value}\n")
        
        formatter.info(f"Merged Prokka and filter stats into {main_sample_txt.name}")

    else:
        missing = [f.name for f in essential_files if not f.exists()]
        formatter.warning(f"Missing essential files after Prokka: {missing}")
        err_file = prokka_dir / "sample.err"
        if err_file.exists() and err_file.stat().st_size > 0:
            formatter.info(f"tbl2asn errors logged to: {err_file.name}")
    return prokka_dir

def run_functional_annotation(prokka_dir: Path, output_dir: Path, threads: int = 8) -> Optional[Path]:
    formatter = get_formatter()
    formatter.section_header("Functional Annotation")
    protein_fasta = prokka_dir / "sample.faa"
    diamond_output = output_dir / "functional_annotations.tsv"

    if not protein_fasta.exists() or protein_fasta.stat().st_size == 0:
        formatter.warning("Protein FASTA file is missing or empty. Skipping functional annotation.")
        return None
    if not SWISSPROT_DB.exists():
        formatter.error(f"SwissProt database not found at {SWISSPROT_DB}", 
                       solutions=["Ensure the database path in your config is correct."])
        raise FileNotFoundError(f"SwissProt DB not found: {SWISSPROT_DB}")

    protein_count = len(list(SeqIO.parse(protein_fasta, 'fasta')))
    formatter.info(f"Annotating {protein_count:,} proteins against SwissProt database...")

    outfmt_cols = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"
    cmd = [
        "diamond", "blastp",
        "--db", str(SWISSPROT_DB),
        "--query", str(protein_fasta),
        "--out", str(diamond_output),
        "--outfmt", "6", *outfmt_cols.split(),
        "--top", "1",
        "--evalue", "1e-5",
        "--threads", str(threads),
        "--sensitive",
        "--block-size", "4.0",
        "--index-chunks", "1",
        "--log"  # CRITICAL: This makes DIAMOND output progress to stderr/stdout
    ]

    return_code = formatter.run_subprocess_with_progress(
        cmd=cmd,
        operation_name="DIAMOND Annotation (SwissProt)",
        total=protein_count,
        unit="proteins",
        parser_func=parse_diamond_progress
    )

    if return_code == 0 and diamond_output.exists() and diamond_output.stat().st_size > 0:
        # Count UNIQUE annotated proteins, not total lines
        unique_proteins = set()
        with open(diamond_output) as f:
            for line in f:
                if line.strip():
                    query_id = line.split('\t')[0]
                    unique_proteins.add(query_id)
        
        annotation_count = len(unique_proteins)
        annotation_rate = (annotation_count / protein_count * 100) if protein_count > 0 else 0
        
        formatter.success("Functional annotation complete")
        formatter.result({
            "Annotated": f"{annotation_count:,} / {protein_count:,} proteins ({annotation_rate:.1f}%)"
        }, indent=2)
        return diamond_output
    else:
        formatter.warning("DIAMOND produced no results or failed.")
        return None
    
def run_pathway_analysis(annotation_file: Path, output_dir: Path):
    """
    Performs a lightweight metabolic pathway analysis by parsing
    EC numbers from an existing annotation file and saving the raw counts.
    """
    formatter = get_formatter()
    formatter.section_header("Metabolic Pathway Analysis")
    
    output_json = output_dir / "pathway_ec_counts.json"
    
    if not annotation_file or not annotation_file.exists() or annotation_file.stat().st_size == 0:
        formatter.warning("Annotation file missing or empty. Skipping pathway analysis.")
        return None

    # Regex to find EC numbers like [EC:1.1.1.1] or [EC:1.1.1.-]
    ec_regex = re.compile(r'\[EC:([\d\.-]+)\]')
    
    ec_counts = Counter()
    total_annotations_parsed = 0
    
    # 1. Read the annotation file and find all EC numbers
    try:
        with open(annotation_file, 'r') as f:
            for line in f:
                total_annotations_parsed += 1
                matches = ec_regex.findall(line)
                for ec in matches:
                    if ec != '-' and not ec.endswith('.-'): # Ignore incomplete EC numbers
                        ec_counts[ec] += 1
    except Exception as e:
        formatter.error(f"Failed to parse annotation file for EC numbers: {e}")
        return None
    
    if not ec_counts:
        formatter.info("No EC numbers found in annotations. Cannot generate pathway report.")
        return None

    formatter.success(f"Found {len(ec_counts)} unique EC numbers from {sum(ec_counts.values())} total hits.")
    
    # 2. Save the raw EC counts to a JSON file
    try:
        with open(output_json, 'w') as f:
            json.dump(ec_counts, f, indent=2)
        
        formatter.success(f"Metabolic pathway raw data saved to {output_json.name}")
        return output_json
        
    except Exception as e:
        formatter.error(f"Failed to save EC counts JSON: {e}")
        return None