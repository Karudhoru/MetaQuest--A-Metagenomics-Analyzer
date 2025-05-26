import subprocess
import os
import pandas as pd
from pathlib import Path
from collections import Counter
from .config import *
from .utils import check_dependencies

def run_kraken(input_file, output_dir):
    """Run Kraken2 classification"""
    report = output_dir/"kraken_report.txt"
    classified = output_dir/"kraken_classified.txt"
    cmd = f"kraken2 --db {KRAKEN_DB} --threads 8 --report {report} --output {classified} {input_file}"
    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    return report

def run_bracken(report_path, output_dir, is_fasta=False):
    """Estimate abundances with Bracken, with FASTA mode handling"""
    bracken_out = output_dir / "bracken_report.tsv"
    
    # Skip Bracken for FASTA files as it's read-based, not contig-based
    if is_fasta:
        print("Skipping Bracken for FASTA input (contig-based analysis)")
        print("Using Kraken2 report directly for taxonomic analysis")
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
    
def create_krona_plot(df, abundance_col, name_col, output_dir):
    """Create Krona hierarchical plot"""
    try:
        krona_input = output_dir/"krona_input.txt"
        with open(krona_input, 'w') as f:
            for _, row in df.iterrows():
                if row[abundance_col] > 0.001:
                    count = int(row[abundance_col] * 1000000)
                    f.write(f"{count}\t{row[name_col]}\n")
        
        krona_output = output_dir/"taxonomy_krona.html"
        cmd = f"ktImportText {krona_input} -o {krona_output}"
        subprocess.run(cmd, shell=True, check=True)
    except Exception as e:
        print(f"Krona plot error: {e}") 
        
def run_fasta_kraken(fasta_path, output_dir):
    """Run Kraken2 classification specifically for FASTA files"""
    report = output_dir/"fasta_kraken_report.txt"
    classified = output_dir/"fasta_kraken_classified.txt"
    cmd = f"kraken2 --db {KRAKEN_DB} --threads 8 --report {report} --output {classified} {fasta_path}"
    print(f"Running: {cmd}")
    try:
        subprocess.run(cmd, shell=True, check=True)
        return report
    except subprocess.CalledProcessError as e:
        print(f"Kraken2 classification failed: {e}")
        # Create empty report file
        with open(report, 'w') as f:
            f.write("0.00\t0\t0\tU\t0\tunclassified\n")
        return report
    
def generate_taxonomy_report(blast_results_path, output_dir):
    """Generate a summary report from pathogen blast results"""
    
    if not blast_results_path or not Path(blast_results_path).exists():
        print("No blast results to process for taxonomy report")
        return
    
    try:
        report_path = output_dir / "pathogen_summary.txt"
        pathogen_hits = {}
        total_hits = 0
        
        with open(blast_results_path, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) >= 7:
                    query_id = parts[0]
                    subject_id = parts[1]
                    identity = float(parts[2])
                    evalue = float(parts[4])
                    description = parts[6] if len(parts) > 6 else subject_id
                    
                    total_hits += 1
                    
                    # Extract organism name from description
                    organism = description.split('[')[-1].replace(']', '') if '[' in description else subject_id
                    
                    if organism not in pathogen_hits:
                        pathogen_hits[organism] = {
                            'count': 0,
                            'best_identity': 0,
                            'best_evalue': 1.0,
                            'queries': set()
                        }
                    
                    pathogen_hits[organism]['count'] += 1
                    pathogen_hits[organism]['queries'].add(query_id)
                    pathogen_hits[organism]['best_identity'] = max(pathogen_hits[organism]['best_identity'], identity)
                    pathogen_hits[organism]['best_evalue'] = min(pathogen_hits[organism]['best_evalue'], evalue)
        
        # Write summary report
        with open(report_path, 'w') as f:
            f.write("PATHOGEN SCREENING SUMMARY\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Total hits found: {total_hits}\n")
            f.write(f"Unique organisms: {len(pathogen_hits)}\n\n")
            
            if pathogen_hits:
                f.write("TOP PATHOGEN HITS:\n")
                f.write("-" * 30 + "\n")
                
                # Sort by number of queries hit
                sorted_hits = sorted(pathogen_hits.items(), 
                                   key=lambda x: len(x[1]['queries']), 
                                   reverse=True)
                
                for organism, data in sorted_hits[:10]:  # Top 10
                    f.write(f"\nOrganism: {organism}\n")
                    f.write(f"  Sequences hit: {len(data['queries'])}\n")
                    f.write(f"  Total hits: {data['count']}\n")
                    f.write(f"  Best identity: {data['best_identity']:.1f}%\n")
                    f.write(f"  Best E-value: {data['best_evalue']:.2e}\n")
            else:
                f.write("No significant pathogen hits found.\n")
        
        print(f"✓ Pathogen summary report created: {report_path}")
        
    except Exception as e:
        print(f"Error generating taxonomy report: {e}")
        # Create basic report
        with open(output_dir / "pathogen_summary.txt", 'w') as f:
            f.write("Error generating pathogen report\n")