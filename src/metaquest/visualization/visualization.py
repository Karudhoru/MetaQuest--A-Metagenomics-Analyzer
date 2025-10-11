#!/usr/bin/env python3
"""
MetaQuest Visualization Module v4.0.0
Updated to work with comprehensive reporting system v4.0.0
Provides robust visualizations for taxonomic, pathogen, and functional analysis
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
    """Base class for all visualizers with common functionality"""
    
    def __init__(self, output_dir: Path):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def save_plot(self, fig, filename: str) -> str:
        """Save plotly figure to HTML file"""
        filepath = self.output_dir / filename
        fig.write_html(filepath)
        return str(filepath)
    
    def load_json_report(self, report_path: Path) -> dict:
        """Load and validate JSON report file"""
        try:
            if not report_path.exists():
                print(f"⚠️ Report file not found: {report_path}")
                return {}
            
            if report_path.stat().st_size == 0:
                print(f"⚠️ Report file is empty: {report_path}")
                return {}
            
            with open(report_path, 'r') as f:
                data = json.load(f)
            
            return data
        except json.JSONDecodeError as e:
            print(f"⚠️ Invalid JSON in report file: {e}")
            return {}
        except Exception as e:
            print(f"⚠️ Error loading report: {e}")
            return {}


class TaxonomicVisualizer(BaseVisualizer):
    """Taxonomic classification visualizations aligned with v4.0.0 reports"""

    def _prepare_dataframe_from_data(self, data):
        """
        Process different input data types into standardized DataFrame.
        Supports: Bracken files, Kraken files, BLAST results (list), DataFrames
        """
        # Handle file paths (Bracken/Kraken reports)
        if isinstance(data, (str, Path)) and Path(data).exists():
            file_path = Path(data)
            
            # Try Bracken format first
            try:
                df = pd.read_csv(file_path, sep='\t')
                if 'fraction_total_reads' in df.columns and 'new_est_reads' in df.columns:
                    print("  ✓ Loaded Bracken format data")
                    return df, 'bracken'
            except:
                pass
            
            # Try Kraken format
            try:
                df = pd.read_csv(file_path, sep='\t', header=None,
                               names=['percentage', 'clade_reads', 'taxon_reads', 'rank', 'taxid', 'name'])
                df['fraction_total_reads'] = df['percentage'] / 100
                df['new_est_reads'] = df['clade_reads']
                print("  ✓ Loaded Kraken format data")
                return df, 'kraken'
            except:
                pass
        
        # Handle BLAST results (list of dictionaries)
        elif isinstance(data, list) and len(data) > 0:
            print("  ✓ Processing BLAST results...")
            organisms = []
            
            for result in data:
                if isinstance(result, dict) and 'hits' in result:
                    for hit in result['hits']:
                        organisms.append(hit.get('organism', 'Unknown'))
                elif isinstance(result, dict) and 'organism' in result:
                    organisms.append(result.get('organism', 'Unknown'))
            
            if not organisms:
                print("  ⚠️ No organisms found in BLAST data")
                return pd.DataFrame(), 'blast'
            
            # Count organisms to estimate abundance
            organism_counts = Counter(organisms)
            df = pd.DataFrame(organism_counts.items(), columns=['name', 'new_est_reads'])
            df['fraction_total_reads'] = df['new_est_reads'] / df['new_est_reads'].sum()
            print(f"  ✓ Processed {len(df)} unique organisms from BLAST")
            return df, 'blast'
        
        # Handle pre-loaded DataFrames
        elif isinstance(data, pd.DataFrame):
            return data, 'dataframe'
        
        print("  ⚠️ Unsupported data format")
        return pd.DataFrame(), 'unknown'

    def create_abundance_chart(self, data, title="Taxonomic Abundance Analysis", top_n=15):
        """
        Create professional abundance bar chart aligned with report standards.
        Supports both FASTQ (Bracken) and FASTA (BLAST) pipelines.
        """
        try:
            df, data_source = self._prepare_dataframe_from_data(data)

            if df.empty:
                print("  ⚠️ No data available for abundance chart")
                return None

            # Filter significant taxa (>0.01% for Bracken, any hit for BLAST)
            threshold = 0.0001 if data_source in ['bracken', 'kraken'] else 0
            significant_taxa = df[df['fraction_total_reads'] > threshold].copy()
            
            if significant_taxa.empty:
                print("  ⚠️ No significant taxa found")
                return None
            
            # Clean and sort
            significant_taxa['name'] = significant_taxa['name'].str.strip()
            top_taxa = significant_taxa.nlargest(top_n, 'fraction_total_reads')
            
            # Configure visualization based on data source
            if data_source == 'blast':
                x_values = top_taxa['new_est_reads']
                x_label = 'Number of BLAST Hits'
                hover_template = '<b>%{y}</b><br>Hits: %{x:,}<extra></extra>'
            else:
                x_values = top_taxa['fraction_total_reads'] * 100
                x_label = 'Relative Abundance (%)'
                hover_template = '<b>%{y}</b><br>Abundance: %{x:.3f}%<br>Reads: %{customdata:,}<extra></extra>'

            # Create bar chart with professional styling
            fig = px.bar(
                x=x_values,
                y=top_taxa['name'],
                orientation='h',
                title=f'{title}<br><sub>Top {top_n} Taxa by Abundance</sub>',
                labels={'x': x_label, 'y': 'Taxon'},
                color=x_values,
                color_continuous_scale='viridis',
                custom_data=[top_taxa['new_est_reads']] if data_source != 'blast' else None
            )
            
            fig.update_traces(hovertemplate=hover_template)
            
            fig.update_layout(
                height=max(500, top_n * 30),
                yaxis={'categoryorder': 'total ascending'},
                template="plotly_white",
                showlegend=False,
                font=dict(size=11, family="Arial, sans-serif"),
                margin=dict(l=280, r=80, t=100, b=60),
                coloraxis_colorbar=dict(
                    title=x_label,
                    thickness=15,
                    len=0.7
                )
            )
            
            # Add summary annotation
            total_abundance = top_taxa['fraction_total_reads'].sum() * 100
            fig.add_annotation(
                text=f"These {top_n} taxa represent {total_abundance:.1f}% of classified reads",
                xref="paper", yref="paper",
                x=0.5, y=-0.08,
                showarrow=False,
                font=dict(size=10, color="gray")
            )
            
            filepath = self.save_plot(fig, "taxonomic_abundance_chart.html")
            print(f"  ✓ Abundance chart saved: {filepath}")
            return filepath
            
        except Exception as e:
            print(f"  ✗ Error creating abundance chart: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    def create_diversity_summary(self, data, output_filename="diversity_metrics.html"):
        """Create diversity metrics visualization from taxonomic data"""
        try:
            df, data_source = self._prepare_dataframe_from_data(data)
            
            if df.empty or data_source not in ['bracken', 'kraken']:
                print("  ⚠️ Diversity analysis requires Bracken/Kraken data")
                return None
            
            # Calculate diversity metrics
            abundances = df['fraction_total_reads'].values
            abundances = abundances[abundances > 0]
            
            # Shannon diversity
            shannon = -np.sum(abundances * np.log(abundances + 1e-10))
            
            # Simpson's diversity
            simpson = 1 - np.sum(abundances ** 2)
            
            # Observed richness
            richness = len(df[df['fraction_total_reads'] > 0.0001])
            
            # Pielou's evenness
            evenness = shannon / np.log(richness) if richness > 1 else 0
            
            # Create gauge charts for metrics
            fig = make_subplots(
                rows=2, cols=2,
                specs=[[{"type": "indicator"}, {"type": "indicator"}],
                       [{"type": "indicator"}, {"type": "indicator"}]],
                subplot_titles=("Shannon Diversity", "Simpson's Index", 
                               "Species Richness", "Pielou's Evenness")
            )
            
            # Shannon diversity gauge
            fig.add_trace(go.Indicator(
                mode="gauge+number",
                value=shannon,
                domain={'x': [0, 1], 'y': [0, 1]},
                gauge={
                    'axis': {'range': [0, 5]},
                    'bar': {'color': "#1f77b4"},
                    'steps': [
                        {'range': [0, 2], 'color': "#ffcccc"},
                        {'range': [2, 3.5], 'color': "#ffffcc"},
                        {'range': [3.5, 5], 'color': "#ccffcc"}
                    ],
                    'threshold': {'line': {'color': "red", 'width': 4}, 'thickness': 0.75, 'value': shannon}
                }
            ), row=1, col=1)
            
            # Simpson's index gauge
            fig.add_trace(go.Indicator(
                mode="gauge+number",
                value=simpson,
                domain={'x': [0, 1], 'y': [0, 1]},
                gauge={
                    'axis': {'range': [0, 1]},
                    'bar': {'color': "#2ca02c"},
                    'steps': [
                        {'range': [0, 0.5], 'color': "#ffcccc"},
                        {'range': [0.5, 0.8], 'color': "#ffffcc"},
                        {'range': [0.8, 1], 'color': "#ccffcc"}
                    ]
                }
            ), row=1, col=2)
            
            # Richness indicator
            fig.add_trace(go.Indicator(
                mode="number+delta",
                value=richness,
                title={'text': "Species Count"},
                number={'font': {'size': 50}},
                domain={'x': [0, 1], 'y': [0, 1]}
            ), row=2, col=1)
            
            # Evenness gauge
            fig.add_trace(go.Indicator(
                mode="gauge+number",
                value=evenness,
                domain={'x': [0, 1], 'y': [0, 1]},
                gauge={
                    'axis': {'range': [0, 1]},
                    'bar': {'color': "#ff7f0e"},
                    'steps': [
                        {'range': [0, 0.3], 'color': "#ffcccc"},
                        {'range': [0.3, 0.7], 'color': "#ffffcc"},
                        {'range': [0.7, 1], 'color': "#ccffcc"}
                    ]
                }
            ), row=2, col=2)
            
            fig.update_layout(
                title_text="Alpha Diversity Analysis<br><sub>Ecological metrics of community structure</sub>",
                template="plotly_white",
                height=700,
                font=dict(size=11)
            )
            
            filepath = self.save_plot(fig, output_filename)
            print(f"  ✓ Diversity metrics saved: {filepath}")
            return filepath
            
        except Exception as e:
            print(f"  ✗ Error creating diversity summary: {e}")
            return None
    
    def create_krona_plot(self, data, output_filename="taxonomy_krona.html"):
        """Create Krona hierarchical visualization"""
        try:
            df, data_source = self._prepare_dataframe_from_data(data)
            
            if df.empty:
                print("  ⚠️ No data for Krona plot")
                return None

            # Filter significant taxa
            significant_df = df[df['fraction_total_reads'] > 0.0001].copy()
            if significant_df.empty:
                print("  ⚠️ No significant taxa for Krona")
                return None
            
            # Create Krona input file
            krona_input = self.output_dir / "krona_input.txt"
            with open(krona_input, 'w') as f:
                for _, row in significant_df.iterrows():
                    count = int(row['new_est_reads'])
                    taxon_name = row['name'].strip()
                    # Simple hierarchy for Krona
                    hierarchy = taxon_name.replace(' ', '\t')
                    f.write(f"{count}\t{hierarchy}\n")
            
            krona_output = self.output_dir / output_filename
            cmd = f"ktImportText {krona_input} -o {krona_output}"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            
            if result.returncode == 0:
                krona_input.unlink()  # Clean up temp file
                print(f"  ✓ Krona plot created: {krona_output}")
                return str(krona_output)
            else:
                print(f"  ⚠️ Krona generation warning: {result.stderr}")
                return None
                
        except Exception as e:
            print(f"  ✗ Error creating Krona plot: {e}")
            return None


class PathogenVisualizer(BaseVisualizer):
    """Pathogen detection visualizations aligned with v4.0.0 reporting"""
    
    def create_risk_assessment_chart(self, report_path, title="Pathogen Risk Assessment"):
        """
        Create comprehensive risk assessment visualization.
        Supports both FASTQ and FASTA+ML report formats.
        """
        try:
            # Load report data
            data = self.load_json_report(Path(report_path))
            if not data:
                return None
            
            # Determine report type and extract detections
            summary = data.get('summary', {})
            analysis_type = summary.get('analysis_type', 'unknown')
            
            if 'FASTA' in analysis_type or 'BLAST' in analysis_type:
                detections = data.get('blast_taxonomy', {}).get('organisms', {})
                print("  ✓ Processing FASTA+ML pathogen report")
            else:
                detections = data.get('pathogen_detections', {})
                print("  ✓ Processing FASTQ pathogen report")
            
            if not detections:
                print("  ⚠️ No pathogen detections found")
                return None

            # Convert to DataFrame
            data_rows = []
            for organism, details in detections.items():
                # Handle different data structures
                risk_level = details.get('risk_level', 'Unknown')
                abundance = details.get('abundance_percentage', 0)
                
                # Get identity from various possible keys
                identity = details.get('sequence_identity') or \
                          details.get('max_identity') or \
                          details.get('blast_identity', 0)
                
                hits = details.get('sequence_hits', 0) or \
                      details.get('hits', 0) or \
                      details.get('blast_hits', 0)
                
                data_rows.append({
                    'organism': organism,
                    'risk_level': risk_level,
                    'abundance': abundance,
                    'identity': identity,
                    'hits': hits
                })

            df = pd.DataFrame(data_rows)
            if df.empty:
                print("  ⚠️ No pathogen data to visualize")
                return None

            # Sort by risk and relevant metric
            risk_order = {'CRITICAL': 4, 'HIGH': 3, 'MEDIUM': 2, 'LOW': 1, 'Unknown': 0}
            df['risk_numeric'] = df['risk_level'].map(risk_order)
            df = df.sort_values(['risk_numeric', 'hits'], ascending=[False, False])
            top_pathogens = df.head(20)

            # Create visualization
            risk_colors = {
                'CRITICAL': '#DC143C',
                'HIGH': '#FF6347',
                'MEDIUM': '#FFA500',
                'LOW': '#32CD32',
                'Unknown': '#808080'
            }
            
            # Use appropriate metric for y-axis
            if top_pathogens['hits'].sum() > 0:
                y_values = top_pathogens['hits']
                y_label = "Detection Count"
            elif top_pathogens['abundance'].sum() > 0:
                y_values = top_pathogens['abundance']
                y_label = "Abundance (%)"
            else:
                y_values = pd.Series([1] * len(top_pathogens))
                y_label = "Presence"

            fig = go.Figure()
            
            for risk_level in ['CRITICAL', 'HIGH', 'MEDIUM', 'LOW', 'Unknown']:
                mask = top_pathogens['risk_level'] == risk_level
                if mask.any():
                    subset = top_pathogens[mask]
                    fig.add_trace(go.Bar(
                        name=risk_level,
                        x=subset['organism'],
                        y=y_values[mask],
                        marker_color=risk_colors[risk_level],
                        hovertemplate=(
                            '<b>%{x}</b><br>' +
                            f'Risk: {risk_level}<br>' +
                            'Count: %{y}<br>' +
                            'Identity: %{customdata:.1f}%<extra></extra>'
                        ),
                        customdata=subset['identity']
                    ))
            
            overall_risk = summary.get('overall_risk_assessment', 'Unknown')
            fig.update_layout(
                title=f"{title}<br><sub>Overall Risk Level: <b>{overall_risk}</b></sub>",
                xaxis_title="Pathogen Species",
                yaxis_title=y_label,
                template="plotly_white",
                height=600,
                xaxis={'tickangle': -45},
                margin=dict(b=180, t=100),
                font=dict(size=11),
                legend=dict(
                    title="Risk Level",
                    orientation="v",
                    yanchor="top",
                    y=0.99,
                    xanchor="right",
                    x=0.99
                ),
                barmode='group'
            )

            filepath = self.save_plot(fig, "pathogen_risk_assessment.html")
            print(f"  ✓ Risk assessment chart saved: {filepath}")
            return filepath

        except Exception as e:
            print(f"  ✗ Error creating risk assessment chart: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    def create_who_priority_distribution(self, report_path):
        """Create WHO priority pathogen distribution chart"""
        try:
            data = self.load_json_report(Path(report_path))
            if not data:
                return None
            
            # Extract detections
            summary = data.get('summary', {})
            analysis_type = summary.get('analysis_type', '')
            
            if 'FASTA' in analysis_type or 'BLAST' in analysis_type:
                detections = data.get('blast_taxonomy', {}).get('organisms', {})
            else:
                detections = data.get('pathogen_detections', {})
            
            if not detections:
                return None
            
            # Count by WHO priority
            who_counts = Counter()
            for organism, details in detections.items():
                priority = details.get('who_priority', 'Not Listed')
                who_counts[priority] += 1
            
            if not who_counts:
                return None
            
            # Create pie chart
            labels = list(who_counts.keys())
            values = list(who_counts.values())
            
            colors = {
                'Critical Priority': '#DC143C',
                'High Priority': '#FF6347',
                'Medium Priority': '#FFA500',
                'Not Listed': '#808080'
            }
            color_list = [colors.get(label, '#808080') for label in labels]
            
            fig = go.Figure(data=[go.Pie(
                labels=labels,
                values=values,
                hole=0.4,
                marker_colors=color_list,
                hovertemplate="<b>%{label}</b><br>Count: %{value}<br>Percentage: %{percent}<extra></extra>"
            )])
            
            fig.update_layout(
                title="WHO Priority Pathogen Classification<br><sub>Distribution by WHO 2024 Priority List</sub>",
                template="plotly_white",
                height=500,
                font=dict(size=11),
                annotations=[dict(
                    text=f"Total<br>Pathogens<br>{sum(values)}",
                    x=0.5, y=0.5,
                    font_size=16,
                    showarrow=False
                )]
            )
            
            filepath = self.save_plot(fig, "who_priority_distribution.html")
            print(f"  ✓ WHO priority chart saved: {filepath}")
            return filepath
            
        except Exception as e:
            print(f"  ✗ Error creating WHO priority chart: {e}")
            return None
    
    def create_detection_confidence_chart(self, report_path):
        """Create detection confidence visualization"""
        try:
            data = self.load_json_report(Path(report_path))
            if not data:
                return None
            
            # Extract detections
            summary = data.get('summary', {})
            analysis_type = summary.get('analysis_type', '')
            
            if 'FASTA' in analysis_type:
                detections = data.get('blast_taxonomy', {}).get('organisms', {})
            else:
                detections = data.get('pathogen_detections', {})
            
            if not detections:
                return None
            
            # Prepare data for scatter plot
            organisms = []
            confidences = []
            risk_levels = []
            detection_methods = []
            
            for organism, details in detections.items():
                conf = details.get('confidence_score', 0)
                if conf > 0:  # Only include if confidence is available
                    organisms.append(organism[:30])  # Truncate long names
                    confidences.append(conf)
                    risk_levels.append(details.get('risk_level', 'Unknown'))
                    methods = details.get('detection_methods', [])
                    detection_methods.append(len(methods) if isinstance(methods, list) else 1)
            
            if not confidences:
                print("  ⚠️ No confidence scores available")
                return None
            
            df = pd.DataFrame({
                'organism': organisms,
                'confidence': confidences,
                'risk_level': risk_levels,
                'method_count': detection_methods
            })
            
            # Create scatter plot
            risk_colors = {
                'CRITICAL': '#DC143C',
                'HIGH': '#FF6347',
                'MEDIUM': '#FFA500',
                'LOW': '#32CD32',
                'Unknown': '#808080'
            }
            
            fig = px.scatter(
                df,
                x='confidence',
                y='organism',
                color='risk_level',
                size='method_count',
                color_discrete_map=risk_colors,
                title="Pathogen Detection Confidence Analysis<br><sub>Confidence scores by detection methods</sub>",
                labels={
                    'confidence': 'Confidence Score',
                    'organism': 'Pathogen',
                    'method_count': 'Methods'
                },
                hover_data={'confidence': ':.3f', 'method_count': True}
            )
            
            # Add confidence threshold lines
            fig.add_vline(x=0.7, line_dash="dash", line_color="green",
                         annotation_text="High Confidence", annotation_position="top")
            fig.add_vline(x=0.5, line_dash="dash", line_color="orange",
                         annotation_text="Medium Confidence", annotation_position="top")
            
            fig.update_layout(
                template="plotly_white",
                height=max(400, len(df) * 25),
                xaxis={'range': [0, 1]},
                margin=dict(l=200, r=100, t=100, b=60),
                font=dict(size=10)
            )
            
            filepath = self.save_plot(fig, "detection_confidence.html")
            print(f"  ✓ Confidence chart saved: {filepath}")
            return filepath
            
        except Exception as e:
            print(f"  ✗ Error creating confidence chart: {e}")
            return None


class FunctionalVisualizer(BaseVisualizer):
    """Functional annotation visualizations aligned with v4.0.0 reporting"""
    
    def create_annotation_quality_dashboard(self, swissprot_results, 
                                           title="Annotation Quality Dashboard"):
        """Create comprehensive annotation quality visualization"""
        try:
            swissprot_path = Path(swissprot_results)
            if not swissprot_path.exists() or swissprot_path.stat().st_size == 0:
                print("  ⚠️ SwissProt results file missing or empty")
                return None
            
            df = pd.read_csv(swissprot_path, sep='\t', header=None,
                           names=['query_id', 'subject_id', 'pident', 'length', 
                                  'evalue', 'bitscore', 'stitle'])
            
            if df.empty:
                print("  ⚠️ No annotation data found")
                return None
            
            # Create 2x2 subplot dashboard
            fig = make_subplots(
                rows=2, cols=2,
                subplot_titles=(
                    'Sequence Identity Distribution',
                    'E-value Quality (Log Scale)',
                    'Identity vs Bit Score Correlation',
                    'Alignment Length Distribution'
                ),
                specs=[[{"type": "histogram"}, {"type": "histogram"}],
                       [{"type": "scatter"}, {"type": "histogram"}]]
            )
            
            # 1. Identity histogram with quality zones
            fig.add_trace(
                go.Histogram(
                    x=df['pident'],
                    nbinsx=30,
                    marker_color='#4169E1',
                    name="Identity",
                    opacity=0.75,
                    hovertemplate="Identity: %{x:.1f}%<br>Count: %{y}<extra></extra>"
                ),
                row=1, col=1
            )
            
            # Add quality zone lines
            for threshold, color, label in [(90, 'green', 'Excellent'), 
                                           (80, 'orange', 'Good'), 
                                           (70, 'red', 'Acceptable')]:
                fig.add_vline(x=threshold, line_dash="dash", line_color=color,
                             annotation_text=f"{label} ({threshold}%)",
                             annotation_position="top", row=1, col=1)
            
            # 2. E-value histogram
            log_evalues = np.log10(df['evalue'].replace(0, 1e-300))
            fig.add_trace(
                go.Histogram(
                    x=log_evalues,
                    nbinsx=30,
                    marker_color='#FF6B6B',
                    name="E-value",
                    opacity=0.75,
                    hovertemplate="Log E-value: %{x:.1f}<br>Count: %{y}<extra></extra>"
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
                        size=5,
                        color=df['pident'],
                        colorscale='Viridis',
                        opacity=0.6,
                        showscale=True,
                        colorbar=dict(title="Identity %", x=1.15, len=0.45, y=0.25)
                    ),
                    name="Quality",
                    hovertemplate="Identity: %{x:.1f}%<br>Bit Score: %{y:.1f}<extra></extra>"
                ),
                row=2, col=1
            )
            
            # 4. Alignment length histogram
            fig.add_trace(
                go.Histogram(
                    x=df['length'],
                    nbinsx=30,
                    marker_color='#32CD32',
                    name="Length",
                    opacity=0.75,
                    hovertemplate="Length: %{x} aa<br>Count: %{y}<extra></extra>"
                ),
                row=2, col=2
            )
            
            # Update axes labels
            fig.update_xaxes(title_text="Sequence Identity (%)", row=1, col=1)
            fig.update_xaxes(title_text="Log₁₀(E-value)", row=1, col=2)
            fig.update_xaxes(title_text="Sequence Identity (%)", row=2, col=1)
            fig.update_xaxes(title_text="Alignment Length (aa)", row=2, col=2)
            
            fig.update_yaxes(title_text="Count", row=1, col=1)
            fig.update_yaxes(title_text="Count", row=1, col=2)
            fig.update_yaxes(title_text="Bit Score", row=2, col=1)
            fig.update_yaxes(title_text="Count", row=2, col=2)
            
            # Calculate summary statistics
            avg_identity = df['pident'].mean()
            median_identity = df['pident'].median()
            high_quality = len(df[df['pident'] >= 90])
            
            # Update layout
            fig.update_layout(
                title_text=(
                    f"{title}<br>"
                    f"<sub>{len(df):,} annotations analyzed | "
                    f"Avg Identity: {avg_identity:.1f}% | "
                    f"High Quality (≥90%): {high_quality:,} ({high_quality/len(df)*100:.1f}%)</sub>"
                ),
                template="plotly_white",
                height=800,
                showlegend=False,
                font=dict(size=10)
            )
            
            filepath = self.save_plot(fig, "annotation_quality_dashboard.html")
            print(f"  ✓ Annotation quality dashboard saved: {filepath}")
            return filepath
            
        except Exception as e:
            print(f"  ✗ Error creating annotation quality dashboard: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    def create_protein_length_analysis(self, prokka_results, 
                                      title="Protein Length Distribution"):
        """Create protein length distribution with quality assessment"""
        try:
            prokka_path = Path(prokka_results)
            faa_files = list(prokka_path.glob("*.faa"))
            
            if not faa_files:
                print("  ⚠️ No protein FASTA files found")
                return None
            
            from Bio import SeqIO
            proteins = list(SeqIO.parse(faa_files[0], "fasta"))
            
            if not proteins:
                print("  ⚠️ No proteins in FASTA file")
                return None
            
            lengths = [len(seq.seq) for seq in proteins]
            
            # Calculate statistics
            mean_length = np.mean(lengths)
            median_length = np.median(lengths)
            std_length = np.std(lengths)
            
            # Categorize proteins
            very_short = len([l for l in lengths if l < 50])
            short = len([l for l in lengths if 50 <= l < 100])
            medium = len([l for l in lengths if 100 <= l < 300])
            long = len([l for l in lengths if l >= 300])
            
            # Create main distribution plot
            fig = make_subplots(
                rows=2, cols=1,
                row_heights=[0.7, 0.3],
                subplot_titles=(
                    'Protein Length Distribution',
                    'Length Category Summary'
                ),
                specs=[[{"type": "histogram"}],
                       [{"type": "bar"}]]
            )
            
            # Main histogram
            fig.add_trace(
                go.Histogram(
                    x=lengths,
                    nbinsx=50,
                    marker_color='#4ECDC4',
                    opacity=0.75,
                    name="Lengths",
                    hovertemplate="Length: %{x} aa<br>Count: %{y}<extra></extra>"
                ),
                row=1, col=1
            )
            
            # Add quality zone markers
            zones = [
                (50, 'red', 'Very Short'),
                (100, 'orange', 'Short'),
                (300, 'green', 'Standard'),
                (500, 'blue', 'Long')
            ]
            
            for threshold, color, label in zones:
                fig.add_vline(
                    x=threshold,
                    line_dash="dash",
                    line_color=color,
                    annotation_text=f"{label} ({threshold} aa)",
                    annotation_position="top",
                    row=1, col=1
                )
            
            # Category bar chart
            categories = ['Very Short\n(<50)', 'Short\n(50-99)', 
                         'Medium\n(100-299)', 'Long\n(≥300)']
            counts = [very_short, short, medium, long]
            colors_bar = ['#DC143C', '#FFA500', '#4ECDC4', '#1E88E5']
            
            fig.add_trace(
                go.Bar(
                    x=categories,
                    y=counts,
                    marker_color=colors_bar,
                    text=[f"{c:,}<br>({c/len(lengths)*100:.1f}%)" for c in counts],
                    textposition='auto',
                    hovertemplate="Category: %{x}<br>Count: %{y}<extra></extra>"
                ),
                row=2, col=1
            )
            
            # Quality assessment
            if very_short / len(lengths) > 0.3:
                quality_note = "⚠️ High proportion of very short proteins"
                quality_color = "red"
            elif mean_length < 200:
                quality_note = "⚠️ Below average protein lengths"
                quality_color = "orange"
            else:
                quality_note = "✓ Normal protein length distribution"
                quality_color = "green"
            
            # Update layout
            fig.update_xaxes(title_text="Protein Length (amino acids)", row=1, col=1)
            fig.update_xaxes(title_text="Category", row=2, col=1)
            fig.update_yaxes(title_text="Number of Proteins", row=1, col=1)
            fig.update_yaxes(title_text="Count", row=2, col=1)
            
            fig.update_layout(
                title_text=(
                    f"{title}<br>"
                    f"<sub>Total: {len(proteins):,} proteins | "
                    f"Mean: {mean_length:.1f} ± {std_length:.1f} aa | "
                    f"Median: {median_length:.1f} aa | "
                    f"{quality_note}</sub>"
                ),
                template="plotly_white",
                height=800,
                showlegend=False,
                font=dict(size=11)
            )
            
            # Add quality annotation
            fig.add_annotation(
                text=quality_note,
                xref="paper", yref="paper",
                x=0.98, y=0.98,
                xanchor="right", yanchor="top",
                showarrow=False,
                font=dict(size=12, color=quality_color),
                bgcolor="rgba(255,255,255,0.9)",
                bordercolor=quality_color,
                borderwidth=2,
                borderpad=10
            )
            
            filepath = self.save_plot(fig, "protein_length_analysis.html")
            print(f"  ✓ Protein length analysis saved: {filepath}")
            return filepath
            
        except Exception as e:
            print(f"  ✗ Error creating protein length analysis: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    def create_functional_category_chart(self, prokka_results):
        """Create functional category distribution from GFF annotations"""
        try:
            prokka_path = Path(prokka_results)
            gff_files = list(prokka_path.glob("*.gff"))
            
            if not gff_files:
                print("  ⚠️ No GFF files found")
                return None
            
            # Parse GFF to extract product annotations
            products = []
            with open(gff_files[0], 'r') as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue
                    fields = line.strip().split('\t')
                    if len(fields) >= 9 and fields[2] == 'CDS':
                        attrs = fields[8]
                        if 'product=' in attrs:
                            product = attrs.split('product=')[1].split(';')[0]
                            products.append(product)
            
            if not products:
                print("  ⚠️ No product annotations found")
                return None
            
            # Categorize products (simplified functional categories)
            categories = {
                'Metabolism': ['dehydrogenase', 'synthase', 'metabol', 'kinase'],
                'Transport': ['transporter', 'transport', 'permease', 'channel'],
                'Transcription': ['transcription', 'regulator', 'repressor', 'activator'],
                'Translation': ['ribosomal', 'translation', 'tRNA', 'aminoacyl'],
                'DNA/RNA': ['DNA', 'RNA', 'polymerase', 'helicase', 'topoisomerase'],
                'Cell Structure': ['membrane', 'cell wall', 'flagell', 'pili'],
                'Stress Response': ['stress', 'heat shock', 'chaperone', 'oxidative'],
                'Hypothetical': ['hypothetical', 'unknown', 'uncharacterized']
            }
            
            category_counts = Counter()
            for product in products:
                product_lower = product.lower()
                categorized = False
                for category, keywords in categories.items():
                    if any(keyword in product_lower for keyword in keywords):
                        category_counts[category] += 1
                        categorized = True
                        break
                if not categorized:
                    category_counts['Other'] += 1
            
            # Create bar chart
            sorted_categories = sorted(category_counts.items(), 
                                      key=lambda x: x[1], reverse=True)
            labels = [c[0] for c in sorted_categories]
            values = [c[1] for c in sorted_categories]
            percentages = [v/len(products)*100 for v in values]
            
            fig = go.Figure(data=[go.Bar(
                y=labels,
                x=values,
                orientation='h',
                marker_color='#4ECDC4',
                text=[f"{v:,} ({p:.1f}%)" for v, p in zip(values, percentages)],
                textposition='auto',
                hovertemplate="<b>%{y}</b><br>Proteins: %{x:,}<br>Percentage: %{customdata:.1f}%<extra></extra>",
                customdata=percentages
            )])
            
            fig.update_layout(
                title="Functional Category Distribution<br><sub>Predicted protein functions by category</sub>",
                xaxis_title="Number of Proteins",
                yaxis_title="Functional Category",
                template="plotly_white",
                height=max(400, len(labels) * 40),
                yaxis={'categoryorder': 'total ascending'},
                margin=dict(l=200, r=100, t=100, b=60),
                font=dict(size=11)
            )
            
            filepath = self.save_plot(fig, "functional_categories.html")
            print(f"  ✓ Functional category chart saved: {filepath}")
            return filepath
            
        except Exception as e:
            print(f"  ✗ Error creating functional category chart: {e}")
            return None


# ============================================================================
# Main visualization functions for pipeline integration
# ============================================================================

def create_taxonomic_visualizations(output_dir, data, **kwargs):
    """
    Create all taxonomic visualizations for v4.0.0 reports.
    
    Args:
        output_dir: Output directory for visualizations
        data: Bracken/Kraken file path or BLAST results list
    
    Returns:
        Dictionary of generated file paths
    """
    visualizer = TaxonomicVisualizer(Path(output_dir))
    generated_files = {}
    
    print("\n=== Creating Taxonomic Visualizations ===")
    
    # Abundance chart
    abundance_chart = visualizer.create_abundance_chart(data)
    if abundance_chart:
        generated_files['abundance_chart'] = abundance_chart
    
    # Diversity metrics (only for Bracken/Kraken)
    diversity_chart = visualizer.create_diversity_summary(data)
    if diversity_chart:
        generated_files['diversity_metrics'] = diversity_chart
    
    # Krona plot
    krona_plot = visualizer.create_krona_plot(data)
    if krona_plot:
        generated_files['krona_plot'] = krona_plot
    
    print(f"✓ Generated {len(generated_files)} taxonomic visualizations")
    return generated_files


def create_pathogen_visualizations(output_dir, traditional_data=None, **kwargs):
    """
    Create all pathogen visualizations for v4.0.0 reports.
    
    Args:
        output_dir: Output directory for visualizations
        traditional_data: Path to pathogen detection JSON report
    
    Returns:
        Dictionary of generated file paths
    """
    visualizer = PathogenVisualizer(Path(output_dir))
    generated_files = {}
    
    print("\n=== Creating Pathogen Visualizations ===")
    
    if traditional_data:
        # Risk assessment chart
        risk_chart = visualizer.create_risk_assessment_chart(traditional_data)
        if risk_chart:
            generated_files['risk_assessment'] = risk_chart
        
        # WHO priority distribution
        who_chart = visualizer.create_who_priority_distribution(traditional_data)
        if who_chart:
            generated_files['who_priority'] = who_chart
        
        # Detection confidence
        confidence_chart = visualizer.create_detection_confidence_chart(traditional_data)
        if confidence_chart:
            generated_files['detection_confidence'] = confidence_chart
    
    print(f"✓ Generated {len(generated_files)} pathogen visualizations")
    return generated_files


def create_functional_visualizations(output_dir, prokka_results=None, 
                                    swissprot_results=None, **kwargs):
    """
    Create all functional visualizations for v4.0.0 reports.
    
    Args:
        output_dir: Output directory for visualizations
        prokka_results: Path to Prokka output directory
        swissprot_results: Path to SwissProt annotation file
    
    Returns:
        Dictionary of generated file paths
    """
    visualizer = FunctionalVisualizer(Path(output_dir))
    generated_files = {}
    
    print("\n=== Creating Functional Visualizations ===")
    
    # Annotation quality dashboard
    if swissprot_results:
        quality_dashboard = visualizer.create_annotation_quality_dashboard(swissprot_results)
        if quality_dashboard:
            generated_files['annotation_quality'] = quality_dashboard
    
    # Protein length analysis
    if prokka_results:
        length_analysis = visualizer.create_protein_length_analysis(prokka_results)
        if length_analysis:
            generated_files['protein_length'] = length_analysis
        
        # Functional categories
        category_chart = visualizer.create_functional_category_chart(prokka_results)
        if category_chart:
            generated_files['functional_categories'] = category_chart
    
    print(f"✓ Generated {len(generated_files)} functional visualizations")
    return generated_files


# ============================================================================
# Legacy compatibility functions
# ============================================================================

def create_visualizations(bracken_report, output_dir):
    """Legacy function - redirects to taxonomic visualizations"""
    print("⚠️ Using legacy create_visualizations function")
    return create_taxonomic_visualizations(output_dir, bracken_report)


def create_pathogen_dashboard(output_dir):
    """Legacy function - redirects to pathogen visualizations"""
    print("⚠️ Using legacy create_pathogen_dashboard function")
    json_report_path = Path(output_dir) / "pathogen_detection_report.json"
    if json_report_path.exists():
        return create_pathogen_visualizations(output_dir, traditional_data=str(json_report_path))
    return {}


# ============================================================================
# Utility functions
# ============================================================================

def validate_visualization_data(data, data_type='taxonomic'):
    """
    Validate input data before visualization.
    
    Args:
        data: Input data (file path, list, or DataFrame)
        data_type: Type of data ('taxonomic', 'pathogen', 'functional')
    
    Returns:
        Tuple of (is_valid: bool, message: str)
    """
    if data is None:
        return False, "No data provided"
    
    if isinstance(data, (str, Path)):
        path = Path(data)
        if not path.exists():
            return False, f"File not found: {path}"
        if path.stat().st_size == 0:
            return False, f"File is empty: {path}"
        return True, "Valid file path"
    
    elif isinstance(data, list):
        if len(data) == 0:
            return False, "Empty data list"
        return True, "Valid data list"
    
    elif isinstance(data, pd.DataFrame):
        if data.empty:
            return False, "Empty DataFrame"
        return True, "Valid DataFrame"
    
    else:
        return False, f"Unsupported data type: {type(data)}"


def generate_visualization_summary(output_dir):
    """
    Generate a summary report of all visualizations in output directory.
    
    Args:
        output_dir: Directory containing visualizations
    
    Returns:
        Dictionary with visualization summary
    """
    output_path = Path(output_dir)
    html_files = list(output_path.glob("*.html"))
    
    summary = {
        'total_visualizations': len(html_files),
        'taxonomic': [],
        'pathogen': [],
        'functional': [],
        'other': []
    }
    
    # Categorize visualizations
    for html_file in html_files:
        name = html_file.name
        if 'taxonomic' in name or 'abundance' in name or 'diversity' in name or 'krona' in name:
            summary['taxonomic'].append(name)
        elif 'pathogen' in name or 'risk' in name or 'who' in name or 'confidence' in name:
            summary['pathogen'].append(name)
        elif 'functional' in name or 'annotation' in name or 'protein' in name or 'category' in name:
            summary['functional'].append(name)
        else:
            summary['other'].append(name)
    
    return summary