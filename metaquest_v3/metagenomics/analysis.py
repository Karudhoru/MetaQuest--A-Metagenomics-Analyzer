import subprocess
import os
import pandas as pd
from pathlib import Path
from Bio import SeqIO
import numpy as np
import json
from .config import *
from .utils import check_dependencies, convert_fastq_to_fasta, split_interleaved
from .taxonomic_analysis import run_kraken, run_bracken, run_fasta_blast_taxonomy, create_blast_taxonomy_summary
from .pathogen_analysis import *
from .functional_analysis import *
from .visualization import create_visualizations, create_functional_plots, create_pathogen_visualization

def run_analysis(input_file, file_type, output_dir):
    """Main analysis controller"""
    try:
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True)

        print(f"Analyzing {input_file} as {file_type}")
        print(f"Output directory: {output_dir}")
        
        if file_type == 'fastq':
            analyze_fastq(input_file, output_dir)
        else:
            analyze_fasta(input_file[0], output_dir)

        print(f"\n🎉 Analysis complete! Open {output_dir}/analysis_dashboard.html to explore results")
    except Exception as e:
        print(f"❌ Analysis failed: {str(e)}")
        raise

def analyze_fastq(reads, output_dir):
    """Process FASTQ files"""
    print("\n=== FASTQ Analysis Pipeline ===")
    
    try:
        # Handle interleaved files
        if len(reads) == 1 and reads[0].endswith(('.fq', '.fastq')):
            print("Splitting interleaved FASTQ...")
            reads = split_interleaved(reads[0], output_dir)

        # Taxonomic classification
        print("1. Running Kraken2 classification...")
        kraken_report = run_kraken(reads, output_dir)

        print("2. Running Bracken abundance estimation...")
        bracken_report = run_bracken(kraken_report, output_dir)
        
        print("3. Creating taxonomy visualizations...")
        create_visualizations(bracken_report, output_dir)
        
    except Exception as e:
        print(f"⚠️ Taxonomic analysis failed: {str(e)}")
    
    try:
        # Functional annotation
        print("4. Converting to FASTA...")
        fasta_path = convert_fastq_to_fasta(reads if len(reads) > 1 else reads[0], output_dir)
        
        print("5. Running gene prediction...")
        prokka_dir = run_prokka(fasta_path, output_dir)
        
        # Check if proteins were found before proceeding
        protein_files = list(Path(prokka_dir).glob("*.faa"))
        if not protein_files or all(os.path.getsize(pf) == 0 for pf in protein_files):
            print("⚠️ Warning: No protein sequences found. Skipping functional annotation.")
            return
        
        print("6. Running functional annotation...")
        swissprot_results = run_swissprot_annotation(prokka_dir, output_dir)
        
        print("7. Creating functional plots...")
        create_functional_plots(prokka_dir, swissprot_results, output_dir)
        
    except Exception as e:
        print(f"⚠️ Functional analysis failed: {str(e)}")

def analyze_fasta(fasta_path, output_dir):
    """Process FASTA files"""
    print("\n=== FASTA Analysis Pipeline ===")
    
    blast_results_data = None  # Store the actual parsed data
    
    try:
        # Taxonomic classification using BLAST for FASTA files
        print("1. Running BLAST taxonomic classification...")
        blast_results_file = run_fasta_blast_taxonomy(fasta_path, output_dir, database="nt", max_sequences=50)
        
        # Load the actual BLAST results data
        if blast_results_file and Path(blast_results_file).exists():
            try:
                with open(blast_results_file, 'r') as f:
                    blast_results_data = json.load(f)
                print(f"✓ Loaded {len(blast_results_data)} BLAST results")
            except Exception as e:
                print(f"⚠️ Could not load BLAST results: {e}")
        
        print("2. Creating taxonomy visualizations...")
        # Use the Kraken-style report created by BLAST analysis
        kraken_style_report = output_dir / "blast_kraken_style_report.txt"
        if kraken_style_report.exists():
            create_visualizations(kraken_style_report, output_dir)
        
    except Exception as e:
        print(f"⚠️ BLAST taxonomic analysis failed: {str(e)}")
        print("   This may be due to network issues or API rate limits.")
    
    try:
        # Pathogen screening
        print("3. Screening for pathogens...")
        blast_report = run_pathogen_scan(fasta_path, output_dir)
        
        print("4. Scanning for AMR genes...")
        amr_results = run_antimicrobial_resistance_scan(fasta_path, output_dir)
        
        print("5. Scanning for virulence factors...")
        vf_results = run_virulence_factor_scan(fasta_path, output_dir)
        
        # Generate reports
        print("6. Generating reports...")
        
        # Only create BLAST taxonomy summary if we have the actual data
        if blast_results_data:
            try:
                create_blast_taxonomy_summary(blast_results_data, output_dir)
                print("✓ BLAST taxonomy summary created")
            except Exception as e:
                print(f"⚠️ Could not create BLAST taxonomy summary: {e}")
        
        if amr_results:
            try:
                generate_amr_report(amr_results, output_dir)
                print("✓ AMR report generated")
            except Exception as e:
                print(f"⚠️ Could not generate AMR report: {e}")
        
        if vf_results:
            try:
                generate_vf_report(vf_results, output_dir)  
                print("✓ Virulence factor report generated")
            except Exception as e:
                print(f"⚠️ Could not generate VF report: {e}")
        
        print("7. Creating pathogen visualization...")
        if blast_report and Path(blast_report).exists():
            try:
                create_pathogen_visualization(blast_report, output_dir)
                print("✓ Pathogen visualization created")
            except Exception as e:
                print(f"⚠️ Could not create pathogen visualization: {e}")
            
    except Exception as e:
        print(f"⚠️ Pathogen analysis failed: {str(e)}")
    
    try:
        # Functional annotation
        print("8. Running gene prediction...")
        prokka_dir = run_prokka(fasta_path, output_dir)
        
        # Check if proteins were found before proceeding
        protein_files = list(Path(prokka_dir).glob("*.faa"))
        if not protein_files or all(os.path.getsize(pf) == 0 for pf in protein_files):
            print("⚠️ Warning: No protein sequences found. Skipping functional annotation.")
            print("   This is likely because the input sequence is too short (<300 bp) for gene prediction.")
            return
        
        print("9. Running functional annotation...")
        swissprot_results = run_swissprot_annotation(prokka_dir, output_dir)
        
        print("10. Creating functional plots...")
        create_functional_plots(prokka_dir, swissprot_results, output_dir)
        
    except Exception as e:
        print(f"⚠️ Functional analysis failed: {str(e)}")
        print("   This may be due to short sequences or annotation issues.")