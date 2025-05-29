import subprocess
import os
import pandas as pd
from pathlib import Path
from Bio import SeqIO
import numpy as np
from .config import *  # Import all config constants
from .utils import check_dependencies, convert_fastq_to_fasta, compute_sequence_statistics, create_pathogenic_taxids, split_interleaved
from .reporting import generate_final_report
from .taxonomic_analysis import run_kraken, run_bracken, run_fasta_kraken, generate_taxonomy_report
from .pathogen_analysis import *
from .functional_analysis import *
from .visualization import create_visualizations, create_functional_plots, create_sequence_quality_plots, create_pathogen_visualization

def run_analysis(input_file, file_type, output_dir):
    """Main analysis controller"""
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Create pathogenic taxids if missing
    create_pathogenic_taxids()
    
    print(f"Analyzing {input_file} as {file_type}")
    print(f"Output directory: {output_dir}")
    
    if file_type == 'fastq':
        # input_file now is a list: [R1] or [R1, R2]
        analyze_fastq(input_file, output_dir)
    else:
        # still a single-element list
        analyze_fasta(input_file[0], output_dir)

    # Generate final outputs
    generate_final_report(output_dir)
    print(f"\n🎉 Analysis complete! Open {output_dir}/analysis_dashboard.html to explore results")

def analyze_fastq(reads, output_dir):
    """Process FASTQ files"""
    print("\n=== FASTQ Analysis Pipeline ===")
    
    # Taxonomy
    # if we got a single interleaved file, split it first
    if len(reads)==1 and reads[0].endswith(('.fq','.fastq')):
        print("Detected interleaved FASTQ → de-interleaving…")
        reads = split_interleaved(reads[0], output_dir)  # returns [R1, R2]

    print("1. Running taxonomic classification (paired=%s)..." % (len(reads)==2))
    kraken_report = run_kraken(reads, output_dir) 

    print("2. Running Bracken abundance estimation...")
    bracken_report = run_bracken(kraken_report, output_dir)
    
    print("3. Creating taxonomy visualizations...")
    create_visualizations(bracken_report, output_dir)
    
    # Convert FASTQ to FASTA for functional annotation
    print("4. Converting FASTQ to FASTA...")
    # collapse paired or single into one fasta
    fasta_path = convert_fastq_to_fasta(reads if len(reads)>1 else reads[0], output_dir)
    
    # Functional annotation
    print("5. Running gene prediction (Prokka)...")
    prokka_dir = run_prokka(fasta_path, output_dir)
    
    print("6. Running functional annotation (SwissProt)...")
    swissprot_results = run_swissprot_annotation(prokka_dir, output_dir)
    
    print("7. Creating functional plots...")
    create_functional_plots(prokka_dir, swissprot_results, output_dir)

def analyze_fasta(fasta_path, output_dir):
    """Enhanced FASTA analysis pipeline"""
    print("\n=== FASTA Analysis Pipeline ===")
    
    # 1. Sequence statistics
    print("1. Computing sequence statistics...")
    seq_stats = compute_sequence_statistics(fasta_path, output_dir)
    
    # 2. Taxonomic classification
    print("2. Running taxonomic classification...")
    kraken_report = run_fasta_kraken(fasta_path, output_dir)
    bracken_report = run_bracken(kraken_report, output_dir, is_fasta=True)
    
    print("3. Creating taxonomy visualizations...")
    create_visualizations(bracken_report, output_dir)
    
    # 4. Pathogen screening (fixed)
    print("4. Screening for pathogens...")
    blast_report = run_pathogen_scan(fasta_path, output_dir)
    
    # 5. AMR and Virulence scanning (NEW)
    print("5. Scanning for antimicrobial resistance genes...")
    amr_results = run_antimicrobial_resistance_scan(fasta_path, output_dir)
    
    print("6. Scanning for virulence factors...")
    vf_results = run_virulence_factor_scan(fasta_path, output_dir)
    
    # 7. Generate pathogen reports
    print("7. Generating pathogen and resistance reports...")
    generate_taxonomy_report(blast_report, output_dir)
    generate_amr_report(amr_results, output_dir)
    generate_vf_report(vf_results, output_dir)
    
    print("8. Creating pathogen visualization...")
    create_pathogen_visualization(blast_report, output_dir)
    
    # 9. Functional annotation
    print("9. Running gene prediction (Prokka)...")
    prokka_dir = run_prokka(fasta_path, output_dir)
    
    print("10. Running functional annotation (SwissProt)...")
    swissprot_results = run_swissprot_annotation(prokka_dir, output_dir)
    
    print("11. Creating functional plots...")
    create_functional_plots(prokka_dir, swissprot_results, output_dir)
    
    # 12. Enhanced sequence quality plots
    print("12. Creating sequence quality plots...")
    create_sequence_quality_plots(seq_stats, output_dir)
    
    # 13. Assembly quality assessment (NEW)
    print("13. Assessing assembly quality...")
    assess_assembly_quality(seq_stats, prokka_dir, output_dir)

def assess_assembly_quality(seq_stats, prokka_dir, output_dir):
    """Assess assembly quality using sequence statistics and gene prediction results"""
    import json
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    
    quality_metrics = {
        'assembly_stats': {
            'total_sequences': seq_stats['total_sequences'],
            'total_length': seq_stats['total_length'],
            'n50': seq_stats['n50'],
            'mean_length': seq_stats['mean_length'],
            'longest_contig': max(seq_stats['lengths']) if seq_stats['lengths'] else 0
        },
        'quality_assessment': {},
        'gene_prediction_stats': {}
    }
    
    # Assembly quality thresholds
    quality_metrics['quality_assessment'] = {
        'n50_quality': 'Good' if seq_stats['n50'] > 50000 else 'Poor' if seq_stats['n50'] < 10000 else 'Fair',
        'total_length_quality': 'Good' if seq_stats['total_length'] > 2000000 else 'Poor' if seq_stats['total_length'] < 500000 else 'Fair',
        'contig_count_quality': 'Good' if seq_stats['total_sequences'] < 100 else 'Poor' if seq_stats['total_sequences'] > 500 else 'Fair',
        'gc_content_quality': 'Normal' if 30 <= seq_stats['mean_gc_content'] <= 70 else 'Unusual'
    }
    
    # Parse Prokka results if available
    if prokka_dir and prokka_dir.exists():
        gff_file = prokka_dir / "PROKKA_*.gff"
        txt_file = prokka_dir / "PROKKA_*.txt"
        
        # Try to find actual files
        import glob
        gff_files = list(prokka_dir.glob("*.gff"))
        txt_files = list(prokka_dir.glob("*.txt"))
        
        if txt_files:
            try:
                with open(txt_files[0], 'r') as f:
                    content = f.read()
                    
                # Parse Prokka statistics from txt file
                lines = content.split('\n')
                for line in lines:
                    if 'CDS' in line and 'predicted' in line:
                        cds_count = int(line.split()[0])
                        quality_metrics['gene_prediction_stats']['predicted_cds'] = cds_count
                        quality_metrics['gene_prediction_stats']['gene_density'] = cds_count / seq_stats['total_length'] * 1000 if seq_stats['total_length'] > 0 else 0
                    elif 'tRNA' in line and 'predicted' in line:
                        trna_count = int(line.split()[0])
                        quality_metrics['gene_prediction_stats']['predicted_trna'] = trna_count
                    elif 'rRNA' in line and 'predicted' in line:
                        rrna_count = int(line.split()[0])
                        quality_metrics['gene_prediction_stats']['predicted_rrna'] = rrna_count
                        
            except Exception as e:
                print(f"Could not parse Prokka results: {e}")
                quality_metrics['gene_prediction_stats'] = {'error': 'Could not parse Prokka output'}
    
    # Calculate overall quality score
    quality_scores = {
        'Good': 3,
        'Fair': 2, 
        'Normal': 2,
        'Poor': 1,
        'Unusual': 1
    }
    
    total_score = sum(quality_scores.get(score, 1) for score in quality_metrics['quality_assessment'].values())
    max_score = len(quality_metrics['quality_assessment']) * 3
    overall_quality = 'Good' if total_score >= max_score * 0.75 else 'Fair' if total_score >= max_score * 0.5 else 'Poor'
    quality_metrics['overall_quality'] = overall_quality
    
    # Save quality report
    quality_file = output_dir / "assembly_quality_report.json"
    with open(quality_file, 'w') as f:
        json.dump(quality_metrics, f, indent=2)
    
    # Create quality visualization
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=('Assembly Statistics', 'Quality Metrics', 'Contig Length Distribution', 'Quality Summary'),
        specs=[[{"type": "table"}, {"type": "bar"}],
               [{"type": "histogram"}, {"type": "indicator"}]]
    )
    
    # Assembly statistics table
    stats_data = [
        ['Total Contigs', f"{seq_stats['total_sequences']:,}"],
        ['Total Length', f"{seq_stats['total_length']:,} bp"],
        ['N50', f"{seq_stats['n50']:,} bp"],
        ['Mean Length', f"{seq_stats['mean_length']:.0f} bp"],
        ['Longest Contig', f"{max(seq_stats['lengths']):,} bp"],
        ['Mean GC Content', f"{seq_stats['mean_gc_content']:.1f}%"]
    ]
    
    fig.add_trace(
        go.Table(
            header=dict(values=['Metric', 'Value'], fill_color='lightblue'),
            cells=dict(values=list(zip(*stats_data)), fill_color='lavender')
        ),
        row=1, col=1
    )
    
    # Quality metrics bar chart
    quality_names = list(quality_metrics['quality_assessment'].keys())
    quality_values = [quality_scores.get(v, 1) for v in quality_metrics['quality_assessment'].values()]
    quality_colors = ['green' if v == 3 else 'orange' if v == 2 else 'red' for v in quality_values]
    
    fig.add_trace(
        go.Bar(x=quality_names, y=quality_values, marker_color=quality_colors, name='Quality Scores'),
        row=1, col=2
    )
    
    # Contig length distribution
    fig.add_trace(
        go.Histogram(x=seq_stats['lengths'], nbinsx=30, name='Length Distribution'),
        row=2, col=1
    )
    
    # Overall quality indicator
    fig.add_trace(
        go.Indicator(
            mode="gauge+number+delta",
            value=total_score,
            domain={'x': [0, 1], 'y': [0, 1]},
            title={'text': f"Overall Quality: {overall_quality}"},
            gauge={
                'axis': {'range': [None, max_score]},
                'bar': {'color': "darkblue"},
                'steps': [
                    {'range': [0, max_score*0.33], 'color': "lightgray"},
                    {'range': [max_score*0.33, max_score*0.66], 'color': "gray"}
                ],
                'threshold': {
                    'line': {'color': "red", 'width': 4},
                    'thickness': 0.75,
                    'value': max_score*0.75
                }
            }
        ),
        row=2, col=2
    )
    
    fig.update_layout(height=800, showlegend=False, title_text="Assembly Quality Assessment")
    fig.write_html(output_dir / "assembly_quality_report.html")
    
    print(f"✓ Assembly quality assessment completed - Overall quality: {overall_quality}")
    return quality_metrics