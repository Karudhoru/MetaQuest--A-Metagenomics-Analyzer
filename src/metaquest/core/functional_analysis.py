"""
Functional Analysis Module - CLEANED v2.0
==========================================
Minimal console output - caller controls verbosity.
Progress parsing removed - using spinners instead.
"""

import subprocess
import os
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from typing import Optional
import time
import psutil

from ..config import SWISSPROT_DB
from ..io.output_formatter import get_formatter


def parse_prokka_sample_txt(sample_file: Path) -> dict:
    """
    Parse Prokka's sample.txt file for accurate gene counts.
    This function is intentionally quiet - just returns the data.
    """
    stats = {
        'organism': '',
        'contigs': 0,
        'bases': 0,
        'CDS': 0,
        'cds': 0,
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
                key_lower = key.lower()
                
                if key_lower == 'organism':
                    stats['organism'] = value
                elif key_lower in ['cds', 'gene', 'rrna', 'trna', 'contigs', 'bases', 'repeat_region']:
                    try:
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


def run_prokka(
    fasta_path: Path,
    output_dir: Path,
    kill_tbl2asn: bool = True,
    tbl2asn_timeout: int = 5
) -> Path:
    """
    Run Prokka with proper subprocess pipe handling (no deadlock).
    FIXED: Reads stdout/stderr in background thread to prevent pipe blocking.
    """
    import subprocess
    import threading
    import queue
    
    formatter = get_formatter()
    prokka_dir = output_dir / "prokka_annotation"
    
    cmd = [
        "prokka",
        "--outdir", str(prokka_dir),
        "--prefix", "sample",
        "--cpus", "12",
        "--force",
        "--metagenome",
        "--centre", "X",
        "--compliant",
        "--fast",
        str(fasta_path)
    ]
    
    formatter.debug(f"Prokka command: {' '.join(cmd)}")
    

    prokka_process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    
    formatter.debug(f"Prokka started (PID: {prokka_process.pid})")
    

    stdout_queue = queue.Queue()
    stderr_queue = queue.Queue()
    
    def read_pipe(pipe, q):
        """Read pipe in background thread."""
        try:
            for line in iter(pipe.readline, ''):
                if line:
                    q.put(line)
            pipe.close()
        except:
            pass
    
    stdout_thread = threading.Thread(target=read_pipe, args=(prokka_process.stdout, stdout_queue), daemon=True)
    stderr_thread = threading.Thread(target=read_pipe, args=(prokka_process.stderr, stderr_queue), daemon=True)
    
    stdout_thread.start()
    stderr_thread.start()
    

    tbl2asn_first_seen = None
    killed_count = 0
    
    if kill_tbl2asn:
        while prokka_process.poll() is None:  # Now safe - pipes are being drained!
            try:
                # Find tbl2asn processes
                for proc in psutil.process_iter(['pid', 'name', 'create_time']):
                    try:
                        proc_name = proc.info.get('name', '').lower()
                        
                        if 'tbl2asn' in proc_name:
                            if tbl2asn_first_seen is None:
                                tbl2asn_first_seen = time.time()
                                formatter.debug(f"Detected {proc.info['name']} (PID: {proc.info['pid']})")
                            
                            time_since_first_seen = time.time() - tbl2asn_first_seen
                            
                            if time_since_first_seen > tbl2asn_timeout:
                                formatter.warning(
                                    f"🔥 Killing {proc.info['name']} (PID: {proc.info['pid']}, "
                                    f"runtime: {time_since_first_seen:.1f}s)"
                                )
                                try:
                                    proc.kill()
                                    proc.wait(timeout=2)
                                    killed_count += 1
                                except:
                                    pass
                    
                    except (psutil.NoSuchProcess, psutil.AccessDenied, KeyError):
                        continue
            
            except Exception as e:
                formatter.debug(f"Monitor error: {e}")
            
            time.sleep(0.5)
    else:
        prokka_process.wait()
    

    stdout_thread.join(timeout=5)
    stderr_thread.join(timeout=5)
    

    stdout_lines = []
    stderr_lines = []
    
    while not stdout_queue.empty():
        stdout_lines.append(stdout_queue.get())
    
    while not stderr_queue.empty():
        stderr_lines.append(stderr_queue.get())
    
    stdout = ''.join(stdout_lines)
    stderr = ''.join(stderr_lines)
    returncode = prokka_process.returncode
    
    if killed_count > 0:
        formatter.warning(f"Killed {killed_count} tbl2asn process(es) during Prokka run")
    

    if kill_tbl2asn:
        try:
            subprocess.run(['pkill', '-9', 'tbl2asn'], stderr=subprocess.DEVNULL, timeout=2)
            subprocess.run(['pkill', '-9', 'real-tbl2asn'], stderr=subprocess.DEVNULL, timeout=2)
        except:
            pass
    
    # Check essential files
    essential_files = [
        prokka_dir / "sample.faa",
        prokka_dir / "sample.ffn",
        prokka_dir / "sample.gff"
    ]
    
    all_essential_exist = all(f.exists() for f in essential_files)
    
    if not all_essential_exist:
        missing = [f.name for f in essential_files if not f.exists()]
        formatter.warning(f"Missing Prokka files: {', '.join(missing)}")
    else:
        formatter.debug("All essential Prokka files generated")
    
    return prokka_dir

def run_functional_annotation(
    prokka_dir: Path,
    output_dir: Path,
    threads: int = 8
) -> Optional[Path]:
    """
    Run DIAMOND functional annotation against SwissProt.
    
    Returns path to annotation file or None if failed.
    """
    formatter = get_formatter()
    protein_fasta = prokka_dir / "sample.faa"
    diamond_output = output_dir / "functional_annotations.tsv"
    
    # Validation
    if not protein_fasta.exists() or protein_fasta.stat().st_size == 0:
        formatter.warning("Protein FASTA file is missing or empty")
        return None
    
    if not SWISSPROT_DB.exists():
        formatter.error(
            f"SwissProt database not found at {SWISSPROT_DB}",
            solutions=["Ensure the database path in your config is correct"]
        )
        raise FileNotFoundError(f"SwissProt DB not found: {SWISSPROT_DB}")
    
    # Count proteins
    protein_count = sum(1 for _ in SeqIO.parse(protein_fasta, 'fasta'))
    formatter.debug(f"Annotating {protein_count:,} proteins against SwissProt")
    
    # Build DIAMOND command
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
        "--log"
    ]
    
    # Run with spinner
    return_code, stdout, stderr = formatter.run_subprocess(
        cmd,
        operation_name="DIAMOND Annotation (SwissProt)",
        capture_output=True,
        show_command=False
    )
    
    # Validation
    if return_code != 0:
        formatter.warning(f"DIAMOND failed with exit code {return_code}")
        return None
    
    if not diamond_output.exists():
        formatter.warning("DIAMOND did not produce output file")
        return None
    
    if diamond_output.stat().st_size == 0:
        formatter.warning("DIAMOND produced empty output (no annotations found)")
        return None
    
    # Validate output content
    try:
        annotation_count = 0
        unique_proteins = set()
        
        with open(diamond_output) as f:
            for line in f:
                if line.strip():
                    annotation_count += 1
                    query_id = line.split('\t')[0]
                    unique_proteins.add(query_id)
        
        if annotation_count == 0:
            formatter.warning("DIAMOND output file is empty or invalid")
            return None
        
        annotation_pct = (len(unique_proteins) / protein_count * 100) if protein_count > 0 else 0
        formatter.debug(
            f"Annotation complete: {len(unique_proteins)}/{protein_count} proteins "
            f"({annotation_pct:.1f}%), {annotation_count} total hits"
        )
        
        return diamond_output
    
    except Exception as e:
        formatter.warning(f"Could not parse DIAMOND output: {e}")
        return None
