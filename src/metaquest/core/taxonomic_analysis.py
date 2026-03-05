"""
Taxonomic Analysis Module - Kraken2/Bracken Only v2.0.0
========================================================

Minimal console output - caller controls verbosity via formatter context managers.
BLAST functionality removed - FASTQ-only pipeline with Kraken2/Bracken.
"""

import subprocess
import os
import pandas as pd
from pathlib import Path
from collections import Counter
import time
import json
from ..config import KRAKEN_DB
from ..io.output_formatter import get_formatter


def run_kraken(input_files, output_dir):
    """
    Run Kraken2 classification for FASTQ files.
    
    Args:
        input_files: FASTQ file path(s) - single file or list of paired files
        output_dir: Output directory for Kraken2 results
        
    Returns:
        Path: Path to Kraken2 report file
        
    Output is controlled by the caller - this function only logs to debug.
    Wrap with spinner or suppressed_output() from calling code.
    """
    formatter = get_formatter()
    
    report = output_dir / "kraken_report.txt"
    classified = output_dir / "kraken_classified.txt"
    
    # Build base command
    cmd = [
        "kraken2",
        "--db", str(KRAKEN_DB),
        "--threads", "8",
        "--report", str(report),
        "--output", str(classified)
    ]
    
    # Handle single-end vs paired-end
    if isinstance(input_files, (list, tuple)) and len(input_files) == 2:
        cmd.extend(["--paired", str(input_files[0]), str(input_files[1])])
        formatter.debug(f"Kraken2: paired-end mode")
    else:
        input_file = input_files[0] if isinstance(input_files, (list, tuple)) else input_files
        cmd.append(str(input_file))
        formatter.debug(f"Kraken2: single-end mode")
    
    formatter.debug(f"Kraken2 database: {KRAKEN_DB}")
    formatter.debug(f"Kraken2 command: {' '.join(cmd)}")
    
    # Run subprocess - caller controls output via spinner/suppression
    returncode, stdout, stderr = formatter.run_subprocess(
        cmd,
        operation_name="Kraken2 classification",
        capture_output=True,
        show_command=False  # Already logged above
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
    
    # Verify report was created
    if not report.exists() or report.stat().st_size == 0:
        formatter.error(
            "Kraken2 report file was not created or is empty",
            solutions=[
                "Check if Kraken2 completed successfully",
                "Verify database integrity",
                "Check available disk space"
            ]
        )
        raise RuntimeError(f"Kraken2 report missing or empty: {report}")
    
    # Log metrics at debug level only
    try:
        with open(report, 'r') as f:
            lines = f.readlines()
            if lines:
                total_reads = sum(int(line.split('\t')[1]) for line in lines if line.strip())
                formatter.debug(f"Kraken2 classified {total_reads:,} sequences")
    except Exception as e:
        formatter.debug(f"Could not parse Kraken2 metrics: {e}")
    
    return report


def run_bracken(report_path, output_dir):
    """
    Estimate abundances with Bracken.
    
    Args:
        report_path: Path to Kraken2 report file
        output_dir: Output directory for Bracken results
        
    Returns:
        Path: Path to Bracken output file (or original Kraken2 report if Bracken fails)
        
    Output is controlled by the caller - this function only logs to debug.
    """
    formatter = get_formatter()
    
    bracken_out = output_dir / "bracken_report.tsv"
    
    # Wait for Kraken report file
    report_file = Path(report_path)
    formatter.debug(f"Waiting for Kraken2 report: {report_file}")
    
    for attempt in range(5):
        if report_file.exists() and report_file.stat().st_size > 0:
            formatter.debug(f"Found Kraken2 report")
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
        "-r", "150",  # Read length
        "-l", "S",    # Taxonomic level (Species)
        "-t", "10"    # Threshold
    ]
    
    formatter.debug(f"Bracken parameters: read_len=150, level=Species, threshold=10")
    
    # Run subprocess - caller controls output
    returncode, stdout, stderr = formatter.run_subprocess(
        cmd,
        operation_name="Bracken abundance estimation",
        capture_output=True,
        show_command=False
    )
    
    if returncode != 0:
        formatter.warning("Bracken failed, using Kraken2 report instead")
        formatter.debug(f"Bracken stderr: {stderr}")
        return report_path
    
    # Log metrics at debug level
    if bracken_out.exists():
        try:
            df = pd.read_csv(bracken_out, sep='\t')
            species_count = len(df)
            total_reads = df['new_est_reads'].sum() if 'new_est_reads' in df.columns else 0
            formatter.debug(f"Bracken found {species_count:,} species, {int(total_reads):,} reads")
        except Exception as e:
            formatter.debug(f"Could not parse Bracken metrics: {e}")
    
    return bracken_out
