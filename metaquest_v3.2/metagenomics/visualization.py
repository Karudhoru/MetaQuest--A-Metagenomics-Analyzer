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
from .config import *

class BaseVisualizer:
    """Base class for all visualizations"""
    
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
    
    def create_abundance_chart(self, data, title="Taxonomic Abundance Analysis", top_n=15):
        """Create clean abundance bar chart with better insights"""
        try:
            # Handle different data formats (Bracken vs Kraken)
            if isinstance(data, (str, Path)):
                df = pd.read_csv(data, sep='\t')
                if 'fraction_total_reads' in df.columns:
                    # Bracken format
                    abundance_col = 'fraction_total_reads'
                    name_col = 'name'
                    count_col = 'new_est_reads'
                else:
                    # Kraken format
                    df = pd.read_csv(data, sep='\t', header=None,
                                   names=['percentage', 'clade_reads', 'taxon_reads', 'rank', 'taxid', 'name'])
                    abundance_col = 'percentage'
                    name_col = 'name'
                    count_col = 'clade_reads'
                    # Convert percentage to fraction for consistency
                    df['fraction_total_reads'] = df['percentage'] / 100
                    abundance_col = 'fraction_total_reads'
            else:
                df = data
                abundance_col = 'fraction_total_reads'
                name_col = 'name'
                count_col = 'new_est_reads'
            
            # Filter significant taxa (>0.01% for better focus)
            significant_taxa = df[
                (df[abundance_col] > 0.0001) & 
                (~df[name_col].str.contains('unclassified', case=False, na=False))
            ].copy()
            
            if len(significant_taxa) == 0:
                print("No significant taxa found for abundance chart")
                return None
            
            # Clean names and get top N
            significant_taxa[name_col] = significant_taxa[name_col].str.strip()
            top_taxa = significant_taxa.nlargest(top_n, abundance_col)
            
            # Convert to percentage for display
            abundance_values = top_taxa[abundance_col] * 100
            
            # Create enhanced bar chart
            fig = px.bar(
                x=abundance_values,
                y=top_taxa[name_col],
                orientation='h',
                title=f'{title} - Top {top_n} Taxa',
                labels={'x': 'Relative Abundance (%)', 'y': 'Taxon'},
                color=abundance_values,
                color_continuous_scale='viridis',
                hover_data={'Estimated Reads': top_taxa[count_col]} if count_col in top_taxa.columns else None
            )
            
            fig.update_layout(
                height=max(500, top_n * 25),
                yaxis={'categoryorder': 'total ascending'},
                template="plotly_white",
                showlegend=False,
                font=dict(size=11),
                margin=dict(l=250, r=50, t=80, b=50)
            )
            
            # Add summary annotation
            total_taxa = len(significant_taxa)
            fig.add_annotation(
                text=f"Showing top {len(top_taxa)} of {total_taxa} detected taxa",
                xref="paper", yref="paper",
                x=0.02, y=0.98, xanchor="left", yanchor="top",
                showarrow=False,
                font=dict(size=10, color="gray")
            )
            
            return self.save_plot(fig, "taxonomic_abundance_chart.html")
            
        except Exception as e:
            print(f"Error creating abundance chart: {e}")
            return None
    
    def create_krona_plot(self, data, output_filename="taxonomy_krona.html"):
        """Create focused Krona plot for taxonomic data"""
        try:
            # Handle different data formats
            if isinstance(data, (str, Path)):
                df = pd.read_csv(data, sep='\t')
                if 'fraction_total_reads' in df.columns:
                    # Bracken format
                    abundance_col = 'fraction_total_reads'
                    name_col = 'name'
                    count_col = 'new_est_reads'
                else:
                    # Kraken format
                    df = pd.read_csv(data, sep='\t', header=None,
                                   names=['percentage', 'clade_reads', 'taxon_reads', 'rank', 'taxid', 'name'])
                    abundance_col = 'percentage'
                    name_col = 'name'
                    count_col = 'clade_reads'
            else:
                df = data
                abundance_col = 'fraction_total_reads'
                name_col = 'name'
                count_col = 'new_est_reads'
            
            # Filter significant taxa (threshold based on format)
            if abundance_col == 'percentage':
                significant_df = df[
                    (df[abundance_col] > 0.01) &
                    (~df[name_col].str.contains('unclassified', case=False, na=False))
                ].copy()
            else:
                significant_df = df[
                    (df[abundance_col] > 0.0001) &
                    (~df[name_col].str.contains('unclassified', case=False, na=False))
                ].copy()
            
            if len(significant_df) == 0:
                print("No significant taxa found for Krona plot")
                return None
            
            # Create Krona input file
            krona_input = self.output_dir / "krona_input.txt"
            with open(krona_input, 'w') as f:
                for _, row in significant_df.iterrows():
                    if abundance_col == 'percentage':
                        count = int(row[count_col]) if row[count_col] > 0 else int(row['taxon_reads'])
                    else:
                        count = int(row.get(count_col, row[abundance_col] * 10000))
                    
                    taxon_name = row[name_col].strip()
                    
                    # Create hierarchical structure for better visualization
                    if ' ' in taxon_name and not taxon_name.startswith('Candidatus'):
                        genus = taxon_name.split(' ')[0]
                        species = taxon_name
                        f.write(f"{count}\t{genus}\t{species}\n")
                    else:
                        f.write(f"{count}\t{taxon_name}\n")
            
            # Generate Krona plot
            krona_output = self.output_dir / output_filename
            cmd = f"ktImportText {krona_input} -o {krona_output}"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            
            if result.returncode == 0:
                krona_input.unlink()  # Clean up
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
        """Create risk-stratified pathogen detection visualization - UPDATED for BLAST+ML format"""
        try:
            # Handle different input formats with enhanced error checking
            if isinstance(pathogen_data, (str, Path)):
                if not Path(pathogen_data).exists():
                    print("⚠️ Pathogen data file does not exist")
                    return None
                with open(pathogen_data, 'r') as f:
                    content = f.read().strip()
                    if not content:
                        print("⚠️ Pathogen data file is empty - no visualizations to create")
                        return None
                    data = json.loads(content)
            elif isinstance(pathogen_data, dict):
                data = pathogen_data
            else:
                return None

            # UPDATED: Handle the new BLAST+ML integrated format
            detections = {}
            analysis_summary = {}
            
            # Check for new BLAST+ML integrated format
            if 'blast_taxonomy_section' in data:
                # New BLAST+ML integrated format
                blast_detections = data.get('blast_taxonomy_section', {}).get('detections', {})
                analysis_summary = data.get('analysis_summary', {})
                
                # Convert BLAST detections to visualization format
                for organism, details in blast_detections.items():
                    detections[organism] = {
                        'risk_level': details.get('risk_level', 'Unknown'),
                        'abundance_percentage': 0,  # BLAST doesn't have abundance
                        'sequence_identity': details.get('blast_identity', 0),
                        'detection_methods': ['blast_taxonomy'],
                        'estimated_reads': 0,
                        'detection_sources': ['BLAST Taxonomic Classification'],
                        'confidence_score': details.get('blast_identity', 0) / 100,
                        'blast_hits': details.get('blast_hits', 0),
                        'blast_evalue': details.get('blast_evalue', 1.0)
                    }
            else:
                # Handle legacy format (metadata wrapper or direct format)
                if 'data' in data and isinstance(data['data'], dict):
                    analysis_summary = data['data'].get('analysis_summary', {})
                    detections = data['data'].get('pathogen_detections', {})
                else:
                    analysis_summary = data.get('analysis_summary', {})
                    detections = data.get('pathogen_detections', data)

            if not detections:
                print("⚠️ No pathogen detections found - skipping risk detection chart")
                return None

            # Convert to DataFrame for visualization
            data_rows = []
            for organism, details in detections.items():
                data_rows.append({
                    'organism': organism,
                    'risk_level': details.get('risk_level', 'Unknown'),
                    'abundance_percentage': details.get('abundance_percentage', 0),
                    'sequence_identity': details.get('sequence_identity', 0) if details.get('sequence_identity') != 'N/A' else 0,
                    'detection_methods': len(details.get('detection_methods', [])),
                    'estimated_reads': details.get('estimated_reads', 0),
                    'detection_sources': ', '.join(details.get('detection_sources', [])),
                    'confidence_score': details.get('confidence_score', 0),
                    'blast_hits': details.get('blast_hits', 0),
                    'blast_evalue': details.get('blast_evalue', 1.0)
                })

            df = pd.DataFrame(data_rows)
            if df.empty:
                print("⚠️ No pathogen data to visualize")
                return None

            # Sort by risk level and sequence identity (for BLAST data)
            risk_order = {'HIGH': 3, 'MEDIUM': 2, 'LOW': 1, 'Unknown': 0}
            df['risk_numeric'] = df['risk_level'].map(risk_order)
            df_sorted = df.sort_values(['risk_numeric', 'sequence_identity'], ascending=[False, False])

            # Take top 20 for clarity
            top_pathogens = df_sorted.head(20)

            # Risk-based color mapping
            risk_colors = {'HIGH': '#DC143C', 'MEDIUM': '#FF8C00', 'LOW': '#32CD32', 'Unknown': '#808080'}
            colors = [risk_colors.get(risk, '#808080') for risk in top_pathogens['risk_level']]

            # Create enhanced bar chart - use BLAST hits instead of abundance for FASTA
            y_values = top_pathogens['blast_hits'] if top_pathogens['blast_hits'].sum() > 0 else top_pathogens['sequence_identity']
            y_title = "BLAST Hits" if top_pathogens['blast_hits'].sum() > 0 else "Sequence Identity (%)"

            fig = go.Figure()
            fig.add_trace(go.Bar(
                x=top_pathogens['organism'],
                y=y_values,
                marker_color=colors,
                name="Pathogen Detection",
                hovertemplate=(
                    "<b>%{x}</b><br>" +
                    "Risk Level: %{customdata[0]}<br>" +
                    "Identity: %{customdata[1]:.1f}%<br>" +
                    "BLAST Hits: %{customdata[2]}<br>" +
                    "E-value: %{customdata[3]:.2e}<br>" +
                    "Sources: %{customdata[4]}<br>" +
                    "Confidence: %{customdata[5]:.3f}<br>" +
                    "<extra></extra>"
                ),
                customdata=list(zip(
                    top_pathogens['risk_level'],
                    top_pathogens['sequence_identity'],
                    top_pathogens['blast_hits'],
                    top_pathogens['blast_evalue'],
                    top_pathogens['detection_sources'],
                    top_pathogens['confidence_score']
                ))
            ))

            # Add risk level indicators
            for i, (risk, organism) in enumerate(zip(top_pathogens['risk_level'], top_pathogens['organism'])):
                if risk == 'HIGH':
                    fig.add_annotation(
                        x=organism, y=y_values.iloc[i],
                        text="⚠️", showarrow=False, yshift=10, font=dict(size=12)
                    )
                elif risk == 'MEDIUM':
                    fig.add_annotation(
                        x=organism, y=y_values.iloc[i],
                        text="🟡", showarrow=False, yshift=10, font=dict(size=10)
                    )

            fig.update_layout(
                title=f"{title}<br>" +
                    f"Overall Risk: {analysis_summary.get('overall_risk_assessment', 'Unknown')} | " +
                    f"Total Detected: {len(df)} organisms",
                xaxis_title="Pathogen Species",
                yaxis_title=y_title,
                template="plotly_white",
                height=600,
                showlegend=False,
                xaxis={'tickangle': 45, 'tickfont': {'size': 10}},
                font=dict(size=11),
                margin=dict(l=60, r=60, t=100, b=120)  # Increase bottom margin for rotated labels
            )

            return self.save_plot(fig, "pathogen_risk_detection.html")

        except Exception as e:
            print(f"Error creating risk detection chart: {e}")
            return None

        
    def create_detection_coverage_chart(self, pathogen_data, title="Detection Method Coverage"):
        """Create detection method coverage visualization - UPDATED for BLAST+ML format"""
        try:
            # Handle different input formats
            if isinstance(pathogen_data, (str, Path)):
                if not Path(pathogen_data).exists():
                    print("⚠️ Pathogen data file does not exist")
                    return None
                with open(pathogen_data, 'r') as f:
                    content = f.read().strip()
                    if not content:
                        print("⚠️ Pathogen data file is empty - no coverage chart to create")
                        return None
                    data = json.loads(content)
            elif isinstance(pathogen_data, dict):
                data = pathogen_data
            else:
                return None

            # UPDATED: Handle the new BLAST+ML integrated format
            detections = {}
            
            # Check for new BLAST+ML integrated format
            if 'blast_taxonomy_section' in data:
                # New BLAST+ML integrated format
                blast_detections = data.get('blast_taxonomy_section', {}).get('detections', {})
                
                # Convert BLAST detections to visualization format
                for organism, details in blast_detections.items():
                    detections[organism] = {
                        'detection_methods': ['blast_taxonomy'],
                        'detection_sources': ['BLAST Taxonomic Classification']
                    }
            else:
                # Handle legacy format
                if 'data' in data and isinstance(data['data'], dict):
                    detections = data['data'].get('pathogen_detections', {})
                else:
                    detections = data.get('pathogen_detections', data)

            if not detections:
                print("⚠️ No pathogen detections found - skipping coverage chart")
                return None

            # Count detection methods
            all_methods = []
            for organism, details in detections.items():
                methods = details.get('detection_methods', [])
                all_methods.extend(methods)

            if not all_methods:
                print("⚠️ No detection methods found")
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
