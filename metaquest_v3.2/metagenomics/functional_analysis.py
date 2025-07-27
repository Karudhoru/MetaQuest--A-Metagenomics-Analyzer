import subprocess
import os
import pandas as pd
from pathlib import Path
from .config import *
from .utils import check_dependencies, parse_prokka_gff

def run_prokka(fasta_path, output_dir):
    """Run Prokka for gene prediction and annotation"""
    prokka_dir = output_dir/"prokka_annotation"
    cmd = f"prokka --outdir {prokka_dir} --prefix sample --cpus 12 --force {fasta_path}"
    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    return prokka_dir

def run_swissprot_annotation(prokka_dir, output_dir):
    """Annotate proteins against SwissProt database"""
    swissprot_out = output_dir/"swissprot_annotation.tsv"
    protein_file = prokka_dir/"sample.faa"
    
    if not protein_file.exists():
        print(f"Warning: Protein file {protein_file} not found")
        return None
    
    cmd = f"diamond blastp -d {SWISSPROT_DB} -q {protein_file} -o {swissprot_out} " \
          "--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle " \
          "--top 1 --evalue 1e-5 --threads 4"
    
    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    return swissprot_out