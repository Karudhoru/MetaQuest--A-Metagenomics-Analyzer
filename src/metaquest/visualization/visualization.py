#!/usr/bin/env python3
"""
MetaQuest Visualization Module
Updated to accommodate new comprehensive reporting system
Handles new data structures and provides robust error handling
"""

import pandas as pd
import subprocess
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from pathlib import Path
import json
from collections import Counter
from ..config import *

class BaseVisualizer:
    # ... (this class remains the same) ...
    def __init__(self, output_dir: Path):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def save_plot(self, fig, filename: str) -> str:
        """Save plotly figure to HTML"""
        filepath = self.output_dir / filename
        fig.write_html(filepath)
        return str(filepath)

class TaxonomicVisualizer(BaseVisualizer):
    """Clean taxonomic classification visualizations"""

    def _prepare_dataframe_from_data(self, data):
        """
        Helper function to process different input data types (BLAST list or Bracken/Kraken file path)
        into a standardized pandas DataFrame.
        """
        if isinstance(data, (str, Path)) and Path(data).exists():
            # Handle Kraken/Bracken report files
            df = pd.read_csv(data, sep='\t')
            if 'fraction_total_reads' in df.columns:
                return df, 'bracken'
            else:
                df = pd.read_csv(data, sep='\t', header=None,
                               names=['percentage', 'clade_reads', 'taxon_reads', 'rank', 'taxid', 'name'])
                df['fraction_total_reads'] = df['percentage'] / 100
                df['new_est_reads'] = df['clade_reads']
                return df, 'kraken'
        elif isinstance(data, list):
            # --- NEW LOGIC TO HANDLE BLAST RESULTS ---
            print("  Processing BLAST results for visualization...")
            all_organisms = []
            for result in data:
                if result and result.get('hits'):
                    for hit in result['hits']:
                        # Get the top hit's organism for abundance counting
                        all_organisms.append(hit.get('organism', 'Unknown'))

            if not all_organisms:
                print("  No identifiable organisms found in BLAST results.")
                return pd.DataFrame(), 'blast'
            
            # Count hits per organism to estimate abundance
            organism_counts = Counter(all_organisms)
            df = pd.DataFrame(organism_counts.items(), columns=['name', 'new_est_reads'])
            df['fraction_total_reads'] = df['new_est_reads'] / df['new_est_reads'].sum()
            return df, 'blast'
        
        # Fallback for other data types like pre-loaded DataFrames
        return pd.DataFrame(data), 'unknown'


    def create_abundance_chart(self, data, title="Taxonomic Abundance Analysis", top_n=15):
        """Create clean abundance bar chart with better insights"""
        try:
            df, data_source = self._prepare_dataframe_from_data(data)

            if df.empty:
                print("No data available to create abundance chart.")
                return None

            # Define standard column names
            abundance_col = 'fraction_total_reads'
            name_col = 'name'
            count_col = 'new_est_reads'
            
            # Filter significant taxa
            significant_taxa = df[df[abundance_col] > 0.0001].copy()
            if significant_taxa.empty:
                print("No significant taxa found for abundance chart.")
                return None
            
            significant_taxa[name_col] = significant_taxa[name_col].str.strip()
            top_taxa = significant_taxa.nlargest(top_n, abundance_col)
            
            # Use appropriate values based on data source
            if data_source == 'blast':
                x_values = top_taxa[count_col]
                x_label = 'Number of BLAST Hits'
                hover_data = None
            else:
                x_values = top_taxa[abundance_col] * 100
                x_label = 'Relative Abundance (%)'
                hover_data = {'Estimated Reads': top_taxa[count_col]}

            fig = px.bar(
                x=x_values,
                y=top_taxa[name_col],
                orientation='h',
                title=f'{title} - Top {top_n} Taxa',
                labels={'x': x_label, 'y': 'Taxon'},
                color=x_values,
                color_continuous_scale='viridis',
                hover_data=hover_data
            )
            
            fig.update_layout(
                height=max(500, top_n * 25),
                yaxis={'categoryorder': 'total ascending'},
                template="plotly_white",
                showlegend=False,
                font=dict(size=11),
                margin=dict(l=250, r=50, t=80, b=50)
            )
            
            return self.save_plot(fig, "taxonomic_abundance_chart.html")
            
        except Exception as e:
            print(f"Error creating abundance chart: {e}")
            return None
    
    def create_krona_plot(self, data, output_filename="taxonomy_krona.html"):
        """Create focused Krona plot for taxonomic data"""
        try:
            df, _ = self._prepare_dataframe_from_data(data)
            
            if df.empty:
                print("No data available to create Krona plot.")
                return None

            # Filter significant taxa
            significant_df = df[df['fraction_total_reads'] > 0.0001].copy()
            if significant_df.empty:
                print("No significant taxa found for Krona plot.")
                return None
            
            krona_input = self.output_dir / "krona_input.txt"
            with open(krona_input, 'w') as f:
                for _, row in significant_df.iterrows():
                    count = int(row['new_est_reads'])
                    taxon_name = row['name'].strip()
                    # Create a simple hierarchy for Krona
                    hierarchy = taxon_name.replace(' ', '\t')
                    f.write(f"{count}\t{hierarchy}\n")
            
            krona_output = self.output_dir / output_filename
            cmd = f"ktImportText {krona_input} -o {krona_output}"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            
            if result.returncode == 0:
                krona_input.unlink()
                print(f"✓ Krona plot created: {krona_output}")
                return str(krona_output)
            else:
                print(f"Krona plot warning: {result.stderr}")
                return None
                
        except Exception as e:
            print(f"Error creating Krona plot: {e}")
            return None


class PathogenVisualizer(BaseVisualizer):
    """Enhanced pathogen detection visualizations for new reporting system"""
    
    def create_risk_detection_chart(self, pathogen_data, title="Pathogen Risk Assessment"):
        """
        Create risk-stratified pathogen detection visualization.
        This function is designed to handle both FASTQ and FASTA+ML report formats.
        """
        try:
            # This function's loading logic is correct and robust.
            # We will copy it to the other function.
            if isinstance(pathogen_data, (str, Path)):
                if not Path(pathogen_data).exists() or Path(pathogen_data).stat().st_size == 0:
                    print("⚠️ Pathogen data file is missing or empty. Skipping visualization.")
                    return None
                with open(pathogen_data, 'r') as f:
                    data = json.load(f)
            elif isinstance(pathogen_data, dict):
                data = pathogen_data
            else:
                print("⚠️ Invalid pathogen data format for visualization.")
                return None

            # Unpack data, which may be nested inside a 'data' key
            if 'data' in data and isinstance(data['data'], dict):
                data = data['data']
            
            detections = {}
            analysis_summary = data.get('analysis_summary', {})
            
            if 'blast_taxonomy_section' in data:
                print("  Visualizing FASTA+ML integrated report...")
                detections = data.get('blast_taxonomy_section', {}).get('detections', {})
            else:
                print("  Visualizing FASTQ traditional report...")
                detections = data.get('pathogen_detections', {})

            if not detections:
                print("⚠️ No pathogen detections found in the report. Skipping chart generation.")
                return None

            # Convert to DataFrame for visualization
            df = pd.DataFrame.from_dict(detections, orient='index').reset_index().rename(columns={'index': 'organism'})
            if df.empty:
                print("⚠️ No pathogen data to visualize.")
                return None

            # Sort by risk level and a relevant metric (blast_hits for FASTA, abundance for FASTQ)
            risk_order = {'HIGH': 3, 'MEDIUM': 2, 'LOW': 1, 'Unknown': 0}
            df['risk_numeric'] = df['risk_level'].map(risk_order)
            sort_metric = 'blast_hits' if 'blast_hits' in df.columns else 'abundance_percentage'
            df_sorted = df.sort_values(['risk_numeric', sort_metric], ascending=[False, False])
            top_pathogens = df_sorted.head(20)

            # Define the primary metric for the y-axis
            if 'blast_hits' in top_pathogens.columns and top_pathogens['blast_hits'].sum() > 0:
                y_values = top_pathogens['blast_hits']
                y_title = "Total BLAST Hits"
            elif 'abundance_percentage' in top_pathogens.columns and top_pathogens['abundance_percentage'].sum() > 0:
                y_values = top_pathogens['abundance_percentage']
                y_title = "Abundance (%)"
            else:
                y_values = pd.Series([1] * len(top_pathogens), index=top_pathogens.index) # Fallback to presence/absence
                y_title = "Detection Count"


            # Create the bar chart
            risk_colors = {'HIGH': '#DC143C', 'MEDIUM': '#FF8C00', 'LOW': '#32CD32', 'Unknown': '#808080'}
            fig = px.bar(
                top_pathogens,
                x='organism',
                y=y_values,
                color='risk_level',
                color_discrete_map=risk_colors,
                labels={'y': y_title, 'organism': 'Pathogen Species'},
                title=f"{title}<br>Overall Risk: {analysis_summary.get('overall_risk_assessment', 'Unknown')}"
            )
            
            fig.update_layout(
                template="plotly_white",
                height=600,
                xaxis={'tickangle': -45, 'categoryorder': 'total descending'},
                margin=dict(b=150) # Increase bottom margin for labels
            )

            return self.save_plot(fig, "pathogen_risk_detection.html")

        except json.JSONDecodeError:
            print(f"Error creating risk detection chart: Invalid JSON format in the report file.")
            return None
        except Exception as e:
            print(f"Error creating risk detection chart: {e}")
            return None

        
    def create_detection_coverage_chart(self, pathogen_data, title="Detection Method Coverage"):
        """Create detection method coverage visualization - UPDATED for robustness"""
        try:
            # --- START: MODIFIED SECTION ---
            # This logic is now identical to the working function above.
            if isinstance(pathogen_data, (str, Path)):
                if not Path(pathogen_data).exists() or Path(pathogen_data).stat().st_size == 0:
                    print("⚠️ Pathogen data file is missing or empty for coverage chart. Skipping.")
                    return None
                with open(pathogen_data, 'r') as f:
                    data = json.load(f)
            elif isinstance(pathogen_data, dict):
                data = pathogen_data
            else:
                return None

            if 'data' in data and isinstance(data['data'], dict):
                data = data['data']
            
            detections = {}
            if 'blast_taxonomy_section' in data:
                # For a FASTA+ML report, all taxonomy comes from one method.
                blast_detections = data.get('blast_taxonomy_section', {}).get('detections', {})
                for organism in blast_detections:
                    detections[organism] = {'detection_methods': ['blast_taxonomy']}
            else:
                # For a FASTQ report, methods are listed per pathogen.
                detections = data.get('pathogen_detections', {})
            # --- END: MODIFIED SECTION ---

            if not detections:
                print("⚠️ No pathogen detections to analyze for coverage chart.")
                return None

            # Count detection methods
            all_methods = []
            for organism, details in detections.items():
                methods = details.get('detection_methods', [])
                all_methods.extend(methods)

            if not all_methods:
                print("⚠️ No detection methods listed in report - skipping coverage chart.")
                return None

            method_counts = Counter(all_methods)

            # Create enhanced method labels
            method_labels = {
                'bracken': 'Taxonomic Classification (Bracken)',
                'blast_taxonomy': 'BLAST Taxonomy Search',
                'taxonomy': 'BLAST Taxonomy Search',
                'sequence': 'Database Sequence Search',
                'diamond': 'DIAMOND Protein Search',
                'ml': 'Machine Learning Prediction'
            }

            methods = list(method_counts.keys())
            counts = list(method_counts.values())
            labels = [method_labels.get(method, method.title()) for method in methods]

            # Create enhanced pie chart
            fig = go.Figure(data=[go.Pie(
                labels=labels,
                values=counts,
                hole=0.3,
                marker_colors=['#4169E1', '#228B22', '#FF6347', '#9932CC', '#20B2AA'][:len(methods)],
                hovertemplate="<b>%{label}</b><br>Pathogens Detected: %{value}<br>Percentage: %{percent}<extra></extra>"
            )])

            fig.update_layout(
                title=f"{title}<br>Multi-method pathogen detection coverage analysis",
                template="plotly_white",
                height=500,
                font=dict(size=11),
                annotations=[dict(text=f"Total<br>Pathogens<br>{sum(counts)}", x=0.5, y=0.5,
                                font_size=14, showarrow=False)]
            )

            return self.save_plot(fig, "detection_method_coverage.html")

        except json.JSONDecodeError:
            print(f"Error creating detection coverage chart: Invalid JSON format in the report file.")
            return None
        except Exception as e:
            print(f"Error creating detection coverage chart: {e}")
            return None


class FunctionalVisualizer(BaseVisualizer):
    """Functional annotation visualizations"""
    
    def create_annotation_quality_chart(self, swissprot_results, title="Annotation Quality Analysis"):
        """Create comprehensive annotation quality visualization"""
        try:
            if not Path(swissprot_results).exists():
                print("⚠️ SwissProt results file does not exist")
                return None
            
            # Check if file is empty
            if Path(swissprot_results).stat().st_size == 0:
                print("⚠️ SwissProt results file is empty")
                return None
            
            df = pd.read_csv(swissprot_results, sep='\t')
            if df.empty:
                print("⚠️ No SwissProt annotation data found")
                return None
            
            # Create subplots for comprehensive analysis
            fig = make_subplots(
                rows=2, cols=2,
                subplot_titles=(
                    'Sequence Identity Distribution',
                    'E-value Distribution (Log Scale)',
                    'Identity vs Bit Score',
                    'Alignment Length Distribution'
                ),
                specs=[[{"type": "histogram"}, {"type": "histogram"}],
                       [{"type": "scatter"}, {"type": "histogram"}]]
            )
            
            # 1. Identity histogram with quality zones
            fig.add_trace(
                go.Histogram(
                    x=df['pident'],
                    nbinsx=25,
                    marker_color='#4169E1',
                    name="Identity",
                    opacity=0.7
                ),
                row=1, col=1
            )
            
            # Add quality zone indicators
            fig.add_vline(x=90, line_dash="dash", line_color="green", 
                         annotation_text="High Quality (90%)", row=1, col=1)
            fig.add_vline(x=75, line_dash="dash", line_color="orange",
                         annotation_text="Medium Quality (75%)", row=1, col=1)
            
            # 2. E-value histogram (log scale)
            fig.add_trace(
                go.Histogram(
                    x=np.log10(df['evalue'] + 1e-300),
                    nbinsx=25,
                    marker_color='#FF6B6B',
                    name="E-value",
                    opacity=0.7
                ),
                row=1, col=2
            )
            
            # 3. Identity vs Bit Score scatter
            fig.add_trace(
                go.Scatter(
                    x=df['pident'],
                    y=df['bitscore'],
                    mode='markers',
                    marker=dict(
                        size=6,
                        color=df['pident'],
                        colorscale='viridis',
                        opacity=0.6,
                        colorbar=dict(title="Identity %", x=1.02)
                    ),
                    name="Quality Correlation",
                    hovertemplate="Identity: %{x:.1f}%<br>Bit Score: %{y:.1f}<extra></extra>"
                ),
                row=2, col=1
            )
            
            # 4. Alignment length histogram
            fig.add_trace(
                go.Histogram(
                    x=df['length'],
                    nbinsx=25,
                    marker_color='#32CD32',
                    name="Alignment Length",
                    opacity=0.7
                ),
                row=2, col=2
            )
            
            # Update layout
            fig.update_layout(
                title_text=f"{title}<br><sub>{len(df)} annotations analyzed | Avg Identity: {df['pident'].mean():.1f}%</sub>",
                template="plotly_white",
                height=700,
                showlegend=False,
                font=dict(size=10)
            )
            
            # Update axes labels
            fig.update_xaxes(title_text="Sequence Identity (%)", row=1, col=1)
            fig.update_xaxes(title_text="Log10(E-value)", row=1, col=2)
            fig.update_xaxes(title_text="Sequence Identity (%)", row=2, col=1)
            fig.update_xaxes(title_text="Alignment Length (aa)", row=2, col=2)
            
            fig.update_yaxes(title_text="Count", row=1, col=1)
            fig.update_yaxes(title_text="Count", row=1, col=2)
            fig.update_yaxes(title_text="Bit Score", row=2, col=1)
            fig.update_yaxes(title_text="Count", row=2, col=2)
            
            return self.save_plot(fig, "annotation_quality_analysis.html")
            
        except Exception as e:
            print(f"Error creating annotation quality chart: {e}")
            return None
    
    def create_protein_length_distribution(self, prokka_results, title="Protein Length Distribution Analysis"):
        """Create protein length distribution visualization"""
        try:
            prokka_path = Path(prokka_results)
            faa_files = list(prokka_path.glob("*.faa"))
            
            if not faa_files:
                print("⚠️ No protein FASTA files found in Prokka results")
                return None
            
            from Bio import SeqIO
            proteins = list(SeqIO.parse(faa_files[0], "fasta"))
            
            if not proteins:
                print("⚠️ No proteins found in FASTA file")
                return None
            
            lengths = [len(seq.seq) for seq in proteins]
            
            # Create distribution plot with quality zones
            fig = go.Figure()
            
            fig.add_trace(go.Histogram(
                x=lengths,
                nbinsx=30,
                marker_color='#4ECDC4',
                opacity=0.7,
                name="Protein Lengths",
                hovertemplate="Length Range: %{x}<br>Count: %{y}<extra></extra>"
            ))
            
            # Add quality zone indicators
            fig.add_vline(x=50, line_dash="dash", line_color="red",
                         annotation_text="Very Short (<50 aa)")
            fig.add_vline(x=100, line_dash="dash", line_color="orange",
                         annotation_text="Short (100 aa)")
            fig.add_vline(x=300, line_dash="dash", line_color="green",
                         annotation_text="Standard (300 aa)")
            
            # Calculate statistics
            mean_length = np.mean(lengths)
            median_length = np.median(lengths)
            
            fig.update_layout(
                title=f"{title}<br><sub>Total Proteins: {len(proteins)} | Mean: {mean_length:.1f} aa | Median: {median_length:.1f} aa</sub>",
                xaxis_title="Protein Length (amino acids)",
                yaxis_title="Number of Proteins",
                template="plotly_white",
                height=500,
                showlegend=False,
                font=dict(size=11)
            )
            
            # Add quality assessment text
            very_short = len([l for l in lengths if l < 50])
            if very_short / len(lengths) > 0.3:
                quality_note = "⚠️ High proportion of very short proteins detected"
                color = "red"
            elif mean_length < 200:
                quality_note = "⚠️ Below average protein lengths observed"
                color = "orange"
            else:
                quality_note = "✓ Normal protein length distribution"
                color = "green"
            
            fig.add_annotation(
                text=quality_note,
                xref="paper", yref="paper",
                x=0.02, y=0.98, xanchor="left", yanchor="top",
                showarrow=False,
                font=dict(size=12, color=color),
                bgcolor="rgba(255,255,255,0.8)",
                bordercolor=color,
                borderwidth=1
            )
            
            return self.save_plot(fig, "protein_length_distribution.html")
            
        except Exception as e:
            print(f"Error creating protein length distribution: {e}")
            return None

# Main visualization functions for pipeline integration - UPDATED
def create_taxonomic_visualizations(output_dir, data, **kwargs):
    """Create all essential taxonomic visualizations"""
    visualizer = TaxonomicVisualizer(Path(output_dir))
    generated_files = {}
    
    print("Creating taxonomic visualizations...")
    
    # Create abundance chart
    abundance_chart = visualizer.create_abundance_chart(data)
    if abundance_chart:
        generated_files['abundance_chart'] = abundance_chart
        print("✓ Taxonomic abundance chart created")
    
    # Create Krona plot
    krona_plot = visualizer.create_krona_plot(data)
    if krona_plot:
        generated_files['krona_plot'] = krona_plot
        print("✓ Krona plot created")
    
    return generated_files

def create_pathogen_visualizations(output_dir, traditional_data=None, **kwargs):
    """Create focused pathogen detection visualizations - UPDATED for new data format"""
    visualizer = PathogenVisualizer(Path(output_dir))
    generated_files = {}
    
    print("Creating pathogen visualizations...")
    
    if traditional_data:
        # Create risk-based detection chart
        risk_chart = visualizer.create_risk_detection_chart(traditional_data)
        if risk_chart:
            generated_files['risk_detection'] = risk_chart
            print("✓ Pathogen risk detection chart created")
        
        # Create detection coverage chart
        coverage_chart = visualizer.create_detection_coverage_chart(traditional_data)
        if coverage_chart:
            generated_files['detection_coverage'] = coverage_chart
            print("✓ Detection method coverage chart created")
    
    return generated_files

def create_functional_visualizations(output_dir, prokka_results=None, swissprot_results=None, **kwargs):
    """Create functional annotation visualizations"""
    visualizer = FunctionalVisualizer(Path(output_dir))
    generated_files = {}
    
    print("Creating functional annotation visualizations...")
    
    # Create annotation quality chart
    if swissprot_results:
        quality_chart = visualizer.create_annotation_quality_chart(swissprot_results)
        if quality_chart:
            generated_files['annotation_quality'] = quality_chart
            print("✓ Annotation quality chart created")
    
    # Create protein length distribution
    if prokka_results:
        length_dist = visualizer.create_protein_length_distribution(prokka_results)
        if length_dist:
            generated_files['protein_length_dist'] = length_dist
            print("✓ Protein length distribution created")
    
    return generated_files

# Legacy functions for backward compatibility (cleaned up)
def create_visualizations(bracken_report, output_dir):
    """Legacy function - now uses clean taxonomic visualizer"""
    return create_taxonomic_visualizations(output_dir, bracken_report)

def create_pathogen_dashboard(output_dir):
    """Legacy function - now uses focused pathogen visualizations"""
    json_report_path = output_dir / "pathogen_detection_report.json"
    if json_report_path.exists():
        return create_pathogen_visualizations(output_dir, traditional_data=json_report_path)
    return {}