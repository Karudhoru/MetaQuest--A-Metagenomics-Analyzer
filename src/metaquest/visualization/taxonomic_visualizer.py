#!/usr/bin/env python3
"""
MetaQuest Visualization Module 
Taxonomic Visualizer Class
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

class TaxonomicVisualizer(BaseVisualizer):
    """Taxonomic classification visualizations aligned with v5.0.0 reports"""

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
                self.formatter.info("Loaded Bracken format data")
                return df, 'kraken'
            except:
                pass
        
        # Handle BLAST results (list of dictionaries)
        elif isinstance(data, list) and len(data) > 0:
            self.formatter.info("  ✓ Processing BLAST results...")
            organisms = []
            
            for result in data:
                if isinstance(result, dict) and 'hits' in result:
                    for hit in result['hits']:
                        organisms.append(hit.get('organism', 'Unknown'))
                elif isinstance(result, dict) and 'organism' in result:
                    organisms.append(result.get('organism', 'Unknown'))
            
            if not organisms:
                self.formatter.warning("No organisms found in BLAST data")
                return pd.DataFrame(), 'blast'
            
            # Count organisms to estimate abundance
            organism_counts = Counter(organisms)
            df = pd.DataFrame(organism_counts.items(), columns=['name', 'new_est_reads'])
            df['fraction_total_reads'] = df['new_est_reads'] / df['new_est_reads'].sum()
            self.formatter.info(f"  ✓ Processed {len(df)} unique organisms from BLAST")
            return df, 'blast'
        
        # Handle pre-loaded DataFrames
        elif isinstance(data, pd.DataFrame):
            return data, 'dataframe'
        
        self.formatter.warning("Unsupported data format for taxonomy viz")
        return pd.DataFrame(), 'unknown'

    def create_abundance_chart(self, data, title="Taxonomic Abundance Analysis", top_n=15):
        """
        Create professional abundance bar chart aligned with report standards.
        Supports both FASTQ (Bracken) and FASTA (BLAST) pipelines.
        """
        try:
            df, data_source = self._prepare_dataframe_from_data(data)

            if df.empty:
                self.formatter.warning("No data available for abundance chart")
                return None

            # Filter significant taxa (>0.01% for Bracken, any hit for BLAST)
            threshold = 0.0001 if data_source in ['bracken', 'kraken'] else 0
            significant_taxa = df[df['fraction_total_reads'] > threshold].copy()
            
            if significant_taxa.empty:
                self.formatter.warning("No significant taxa found for abundance chart")
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
            self.formatter.success(f"Abundance chart saved: {Path(filepath).name}")
            return filepath
            
        except Exception as e:
            self.formatter.error(f"Error creating abundance chart: {e}")
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
