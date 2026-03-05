# src/metaquest/core/assembly.py

"""
Metagenomic Assembly Module
Handles SPAdes and MEGAHIT assembly orchestration.
"""
from pathlib import Path
from typing import List
import subprocess
import os
from Bio import SeqIO
from ..io.output_formatter import get_formatter


def assemble_reads_to_fasta(reads: List[str], output_dir: Path, formatter, threads: int = 8) -> Path:
    """
    Performs metagenomic assembly using MEGAHIT.
    
    CORRECTED: Validates output file contents
    """
    megahit_outdir = output_dir / "megahit_assembly"
    final_contigs = megahit_outdir / "final.contigs.fa"
    
    cmd = ["megahit", "-t", str(threads)]
    if len(reads) == 1:
        cmd.extend(["-r", reads[0]])
    elif len(reads) == 2:
        cmd.extend(["-1", reads[0], "-2", reads[1]])
    
    cmd.extend(["-o", str(megahit_outdir), "--force"])
    
    returncode, _, stderr = formatter.run_subprocess(cmd, "Running MEGAHIT assembler", show_command=True)

    # CORRECTED: Comprehensive validation of assembly output
    if returncode != 0:
        formatter.error(
            "MEGAHIT assembly failed (non-zero exit code).",
            solutions=[f"Check log: {megahit_outdir}/log"]
        )
        raise RuntimeError(f"MEGAHIT failed with exit code {returncode}")
    
    if not final_contigs.exists():
        formatter.error("MEGAHIT did not produce final.contigs.fa")
        raise RuntimeError(f"Expected output missing: {final_contigs}")
    
    # NEW: Validate file has content
    if final_contigs.stat().st_size == 0:
        formatter.error("MEGAHIT produced empty contigs file")
        raise RuntimeError(f"Assembly file is empty: {final_contigs}")
    
    # NEW: Validate FASTA format and count contigs
    try:
        contig_count = sum(1 for _ in SeqIO.parse(final_contigs, "fasta"))
        if contig_count == 0:
            formatter.error("MEGAHIT produced no valid contigs")
            raise RuntimeError(f"No sequences in assembly: {final_contigs}")
        
        formatter.debug(f"Assembly validation: {contig_count} contigs produced")
    except Exception as e:
        formatter.error(f"Assembly file is corrupted or invalid FASTA")
        raise RuntimeError(f"Invalid assembly output: {e}")

    return final_contigs