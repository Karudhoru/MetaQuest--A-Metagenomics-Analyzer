import subprocess
import os
import pandas as pd
from pathlib import Path
from ..config import *


def run_prokka(fasta_path, output_dir):
    """Run Prokka for gene prediction and annotation"""
    prokka_dir = output_dir/"prokka_annotation"
    cmd = f"prokka --outdir {prokka_dir} --prefix sample --cpus 12 --force --metagenome --centre X --compliant --noanno --fast {fasta_path}"
    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    return prokka_dir

def run_functional_annotation(prokka_dir: Path, output_dir: Path, threads: int = 8) -> Path:
    """
    Annotates proteins against the custom Swiss-Prot+COG database using DIAMOND.
    """
    print("--- Running Functional Annotation (DIAMOND + Custom DB) ---")
    protein_fasta = prokka_dir / "sample.faa"
    diamond_output = output_dir / "functional_annotations.tsv"
    
    if not SWISSPROT_DB.exists():
        print(f"❌ ERROR: Custom database not found at {SWISSPROT_DB}. Please create it first.")
        raise FileNotFoundError(f"Custom database not found: {SWISSPROT_DB}")

    cmd = (f"diamond blastp --db {SWISSPROT_DB} --query {protein_fasta} --out {diamond_output} "
           f"--outfmt 6 qseqid stitle --top 1 --evalue 1e-5 --threads {threads} --sensitive")

    try:
        print(f"Running: {cmd}")
        subprocess.run(cmd, shell=True, check=True)
        print(f"✅ Custom functional annotation complete. Results saved to {diamond_output}")
        return diamond_output
    except Exception as e:
        print(f"❌ DIAMOND search against custom DB failed: {e}")
        raise

