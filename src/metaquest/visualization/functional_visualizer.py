#!/usr/bin/env python3
"""
MetaQuest Visualization Module v4.0.0
Functional Visualizer Class
"""
from collections import Counter
import json
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from .base_visualizer import BaseVisualizer

class FunctionalVisualizer(BaseVisualizer):
    """Functional annotation visualizations aligned with v4.0.0 reporting"""
    
    def create_annotation_quality_dashboard(self, swissprot_results, 
                                           title="Annotation Quality Dashboard"):
        """Create comprehensive annotation quality visualization"""
        try:
            swissprot_path = Path(swissprot_results)
            if not swissprot_path.exists() or swissprot_path.stat().st_size == 0:
                self.formatter.warning("SwissProt results file missing or empty")
                return None
            
            df = pd.read_csv(swissprot_path, sep='\t', header=None,
                           names=['query_id', 'subject_id', 'pident', 'length', 
                                  'evalue', 'bitscore', 'stitle'])
            
            if df.empty:
                self.formatter.warning("No annotation data found")
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
            self.formatter.success(f"Annotation quality dashboard saved: {Path(filepath).name}")
            return filepath
            
        except Exception as e:
            self.formatter.error(f"Error creating annotation quality dashboard: {e}")
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
                self.formatter.warning("No protein FASTA files found")
                return None
            
            from Bio import SeqIO
            proteins = list(SeqIO.parse(faa_files[0], "fasta"))
            
            if not proteins:
                self.formatter.warning("No proteins in FASTA file")
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
            self.formatter.success(f"Protein length analysis saved: {Path(filepath).name}")
            return filepath
            
        except Exception as e:
            self.formatter.error(f"Error creating protein length analysis: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    def create_functional_category_chart(self, prokka_results):
        """Create functional category distribution from GFF annotations"""
        try:
            prokka_path = Path(prokka_results)
            gff_files = list(prokka_path.glob("*.gff"))
            
            if not gff_files:
                self.formatter.warning("No GFF files found")
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
                self.formatter.warning("No product annotations found")
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
            self.formatter.success(f"Functional category chart saved: {Path(filepath).name}")
            return filepath
            
        except Exception as e:
            self.formatter.error(f"Error creating functional category chart: {e}")
            return None
