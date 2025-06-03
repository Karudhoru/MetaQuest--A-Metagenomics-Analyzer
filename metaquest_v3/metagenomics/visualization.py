import pandas as pd
import subprocess
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
from pathlib import Path
from .config import *

def create_visualizations(bracken_report, output_dir):
    """Create essential taxonomy visualizations - focusing on actionable insights"""
    try:
        # Try reading as Bracken format first
        try:
            df = pd.read_csv(bracken_report, sep='\t')
            if 'fraction_total_reads' in df.columns:
                abundance_col = 'fraction_total_reads'
                name_col = 'name'
            else:
                raise ValueError("Not Bracken format")
        except:
            # Fall back to Kraken format
            df = pd.read_csv(bracken_report, sep='\t', header=None,
                           names=['percentage', 'clade_reads', 'taxon_reads', 'rank', 'taxid', 'name'])
            df['fraction_total_reads'] = df['percentage'] / 100
            abundance_col = 'fraction_total_reads'
            name_col = 'name'
        
        # Filter for significant taxa (>0.1% abundance)
        significant_taxa = df[df[abundance_col] > 0.001].copy()
        
        if len(significant_taxa) > 0:
            # Clean up names
            significant_taxa[name_col] = significant_taxa[name_col].str.strip()
            
            # Single comprehensive taxonomy plot - bar chart (easier to read than pie)
            fig = px.bar(
                significant_taxa.head(10),  # Top 10 is sufficient
                x=abundance_col,
                y=name_col,
                orientation='h',
                title='Top 10 Taxa by Abundance',
                color=abundance_col,
                color_continuous_scale='viridis',
                labels={abundance_col: 'Relative Abundance', name_col: 'Taxon'}
            )
            fig.update_layout(
                height=500,
                yaxis={'categoryorder': 'total ascending'},
                template="plotly_white",
                showlegend=False
            )
            fig.update_traces(
                hovertemplate='<b>%{y}</b><br>Abundance: %{x:.3%}<extra></extra>'
            )
            fig.write_html(output_dir/"taxonomy_overview.html")
        
        # Keep Krona plot as it's uniquely informative for hierarchical data
        create_krona_plot(df, abundance_col, name_col, output_dir)
        
    except Exception as e:
        print(f"Visualization error: {e}")

def create_krona_plot(output_dir):
    """Create Krona hierarchical plot directly from blast_report.txt to ensure consistency"""
    try:
        blast_report_file = output_dir / "blast_report.txt"
        
        if not blast_report_file.exists():
            print("Warning: blast_report.txt not found, cannot create Krona plot")
            return
        
        krona_input = output_dir / "krona_input.txt"
        
        with open(blast_report_file, 'r') as infile, open(krona_input, 'w') as outfile:
            for line in infile:
                parts = line.strip().split('\t')
                if len(parts) >= 6:
                    percentage = float(parts[0])
                    hit_count = int(parts[1])
                    organism = parts[5]
                    
                    # Skip unclassified and very low abundance
                    if organism != "unclassified" and percentage > 0.1:
                        outfile.write(f"{hit_count}\t{organism}\n")
        
        krona_output = output_dir / "taxonomy_krona.html"
        cmd = f"ktImportText {krona_input} -o {krona_output}"
        subprocess.run(cmd, shell=True, check=True)
        print(f"✓ Krona plot created: {krona_output}")
        print(f"  Krona input file: {krona_input}")
        
    except Exception as e:
        print(f"Krona plot error: {e}")

def create_pathogen_visualization(blast_report, output_dir):
    """Create focused pathogen detection summary"""
    try:
        if not Path(blast_report).exists() or Path(blast_report).stat().st_size == 0:
            # Simple message for no pathogens - no need for complex visualization
            with open(output_dir/"pathogen_summary.txt", 'w') as f:
                f.write("No pathogen hits detected in this sample.\n")
            return
            
        df = pd.read_csv(blast_report, sep='\t')
        
        # Simple pathogen count by category - more actionable than sunburst
        if 'Pathogenic' in df.columns:
            pathogen_counts = df['Pathogenic'].value_counts()
            
            fig = px.bar(
                x=pathogen_counts.index,
                y=pathogen_counts.values,
                title='Pathogen Detection Summary',
                labels={'x': 'Pathogen Category', 'y': 'Number of Hits'},
                color=pathogen_counts.values,
                color_continuous_scale='reds'
            )
            fig.update_layout(
                template="plotly_white",
                showlegend=False
            )
            fig.write_html(output_dir/"pathogen_summary.html")
        
        # If there are actual pathogen hits, create a simple table
        if len(df) > 0:
            pathogen_table = df[['staxids', 'pident', 'bitscore']].head(20)
            pathogen_table.to_csv(output_dir/"top_pathogen_hits.csv", index=False)

    except Exception as e:
        print(f"Pathogen visualization error: {e}")

def create_functional_plots(prokka_dir, swissprot_results, output_dir):
    """Create essential functional annotation insights"""
    try:
        # SwissProt annotation results
        if swissprot_results and Path(swissprot_results).exists():
            df = pd.read_csv(swissprot_results, sep='\t', header=None,
                            names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                                  'gapopen', 'qstart', 'qend', 'sstart', 'send', 
                                  'evalue', 'bitscore', 'stitle'])
            
            if len(df) == 0:
                return

            # Focus on annotation quality - most important metric
            # Bin identity scores into quality categories
            df['quality_category'] = pd.cut(df['pident'], 
                                          bins=[0, 30, 50, 70, 90, 100], 
                                          labels=['Poor (<30%)', 'Low (30-50%)', 
                                                'Moderate (50-70%)', 'Good (70-90%)', 
                                                'Excellent (>90%)'])
            
            quality_counts = df['quality_category'].value_counts()
            
            fig = px.pie(
                values=quality_counts.values,
                names=quality_counts.index,
                title='Annotation Quality Distribution',
                color_discrete_sequence=['#ff6b6b', '#feca57', '#48dbfb', '#0abde3', '#00d2d3']
            )
            fig.update_traces(
                textposition='inside', 
                textinfo='percent+label',
                hovertemplate='<b>%{label}</b><br>Count: %{value}<br>Percentage: %{percent}<extra></extra>'
            )
            fig.update_layout(template="plotly_white")
            fig.write_html(output_dir/"annotation_quality.html")
            
            # Summary statistics table instead of complex 3D plot
            summary_stats = {
                'Total Annotations': len(df),
                'Mean Identity %': df['pident'].mean(),
                'High Quality Hits (>70%)': len(df[df['pident'] > 70]),
                'Excellent Hits (>90%)': len(df[df['pident'] > 90])
            }
            
            with open(output_dir/"annotation_summary.txt", 'w') as f:
                f.write("Functional Annotation Summary\n")
                f.write("="*30 + "\n")
                for key, value in summary_stats.items():
                    if isinstance(value, float):
                        f.write(f"{key}: {value:.2f}\n")
                    else:
                        f.write(f"{key}: {value}\n")

    except Exception as e:
        print(f"Functional plot error: {e}")

def generate_analysis_summary(output_dir):
    """Generate a single consolidated summary of key findings"""
    try:
        summary = {
            'files_generated': [],
            'key_findings': []
        }
        
        # Check what files were generated
        viz_files = ['taxonomy_overview.html', 'pathogen_summary.html', 
                    'annotation_quality.html', 'krona.html']
        
        for file in viz_files:
            if (output_dir / file).exists():
                summary['files_generated'].append(file)
        
        # Write simple summary
        with open(output_dir / "analysis_summary.txt", 'w') as f:
            f.write("Microbiome Analysis Summary\n")
            f.write("="*30 + "\n\n")
            f.write("Generated Visualizations:\n")
            for file in summary['files_generated']:
                f.write(f"  - {file}\n")
            f.write(f"\nTotal files: {len(summary['files_generated'])}\n")
            
        print(f"Analysis complete. Generated {len(summary['files_generated'])} visualization files.")
        
    except Exception as e:
        print(f"Summary generation error: {e}")