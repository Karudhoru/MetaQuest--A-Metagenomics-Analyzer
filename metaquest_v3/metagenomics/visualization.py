import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.figure_factory as ff
import numpy as np
from pathlib import Path
from .taxonomic_analysis import create_krona_plot
from .reporting import *
from .config import *

def create_visualizations(bracken_report, output_dir):
    """Create interactive taxonomy plots with enhanced visuals"""
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
            
            # Enhanced pie chart
            fig = px.pie(
                significant_taxa.head(15),
                values=abundance_col,
                names=name_col,
                title='Top 15 Taxa by Abundance',
                color_discrete_sequence=TAXONOMIC_COLORS,
                hole=0.3
            )
            fig.update_traces(
                textposition='inside', 
                textinfo='percent+label',
                hovertemplate='<b>%{label}</b><br>Abundance: %{percent}<br>Count: %{value:.6f}<extra></extra>'
            )
            fig.write_html(output_dir/"taxonomy_pie.html")
            
            # Treemap visualization
            fig2 = px.treemap(
                significant_taxa.head(25),
                path=[name_col],
                values=abundance_col,
                title='Taxonomic Abundance (Treemap)',
                color=abundance_col,
                color_continuous_scale='RdYlBu'
            )
            fig2.update_layout(template="plotly_white")
            fig2.write_html(output_dir/"taxonomy_treemap.html")
        
        # Create Krona plot
        create_krona_plot(df, abundance_col, name_col, output_dir)
        
        # Create analysis dashboard
        create_analysis_dashboard(output_dir)

    except Exception as e:
        print(f"Visualization error: {e}")

def create_pathogen_visualization(blast_report, output_dir):
    """Create interactive pathogen summary visual"""
    try:
        if not Path(blast_report).exists() or Path(blast_report).stat().st_size == 0:
            fig = go.Figure()
            fig.add_annotation(
                text='No pathogen hits detected',
                xref="paper", yref="paper",
                x=0.5, y=0.5,
                showarrow=False,
                font=dict(size=20))
            fig.update_layout(title_text="Pathogen Detection Summary")
            fig.write_html(output_dir/"pathogen_summary.html")
            return
            
        df = pd.read_csv(blast_report, sep='\t')
        
        # Sunburst chart for pathogen relationships
        fig = px.sunburst(
            df, 
            path=['Pathogenic', 'staxids'], 
            values='bitscore',
            color='pident',
            color_continuous_scale='RdBu',
            title='Pathogen Hit Hierarchy'
        )
        fig.write_html(output_dir/"pathogen_sunburst.html")

    except Exception as e:
        print(f"Pathogen visualization error: {e}")

def create_functional_plots(prokka_dir, swissprot_results, output_dir):
    """Create interactive functional annotation visualizations"""
    try:
        # SwissProt annotation results
        if swissprot_results and Path(swissprot_results).exists():
            df = pd.read_csv(swissprot_results, sep='\t', header=None,
                            names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                                  'gapopen', 'qstart', 'qend', 'sstart', 'send', 
                                  'evalue', 'bitscore', 'stitle'])
            
            if len(df) == 0:
                return

            # Identity distribution
            fig = px.histogram(
                df, x='pident', nbins=30,
                title='Protein Identity Distribution',
                color_discrete_sequence=['#4ECDC4']
            )
            fig.update_layout(
                xaxis_title="Percent Identity",
                yaxis_title="Number of Hits",
                template="plotly_white"
            )
            fig.write_html(output_dir/"swissprot_identity.html")

            # Interactive 3D scatter plot
            fig2 = px.scatter_3d(
                df.head(500),
                x='pident', y='evalue', z='bitscore',
                color='bitscore',
                title='3D Annotation Overview',
                labels={'pident': 'Identity %', 'evalue': 'E-value', 'bitscore': 'Bit Score'}
            )
            fig2.update_traces(marker_size=3)
            fig2.write_html(output_dir/"3d_annotation.html")

    except Exception as e:
        print(f"Functional plot error: {e}")

def create_swissprot_plots(swissprot_file, output_dir):
    """Create plots from SwissProt annotation results"""
    df = pd.read_csv(swissprot_file, sep='\t', header=None,
                    names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                          'gapopen', 'qstart', 'qend', 'sstart', 'send', 
                          'evalue', 'bitscore', 'stitle'])
    
    if len(df) == 0:
        print("No SwissProt hits found")
        return
    
    # Identity distribution
    plt.figure(figsize=(10, 6))
    plt.hist(df['pident'], bins=30, edgecolor='black', alpha=0.7, color='skyblue')
    plt.xlabel('Percent Identity')
    plt.ylabel('Number of Hits')
    plt.title('SwissProt BLAST Hit Identity Distribution')
    plt.tight_layout()
    plt.savefig(output_dir/"swissprot_identity.png", dpi=300, bbox_inches='tight')
    plt.close()
    print("✓ Created swissprot_identity.png")
    
    # Top functional categories (simplified from protein descriptions)
    df['function'] = df['stitle'].str.extract(r'([A-Z][a-z]+(?:\s+[a-z]+){0,2})')
    top_functions = df['function'].value_counts().head(15)
    
    if len(top_functions) > 0:
        plt.figure(figsize=(12, 8))
        top_functions.plot(kind='barh', color='coral')
        plt.xlabel('Number of Proteins')
        plt.title('Top 15 Protein Functions (SwissProt)')
        plt.tight_layout()
        plt.savefig(output_dir/"swissprot_functions.png", dpi=300, bbox_inches='tight')
        plt.close()
        print("✓ Created swissprot_functions.png")

def create_pathogen_visualization(blast_report, output_dir):
    """Create interactive pathogen summary visual"""
    try:
        if not Path(blast_report).exists() or Path(blast_report).stat().st_size == 0:
            fig = go.Figure()
            fig.add_annotation(
                text='No pathogen hits detected',
                xref="paper", yref="paper",
                x=0.5, y=0.5,
                showarrow=False,
                font=dict(size=20))
            fig.update_layout(title_text="Pathogen Detection Summary")
            fig.write_html(output_dir/"pathogen_summary.html")
            return
            
        df = pd.read_csv(blast_report, sep='\t')
        
        # Sunburst chart for pathogen relationships
        fig = px.sunburst(
            df, 
            path=['Pathogenic', 'staxids'], 
            values='bitscore',
            color='pident',
            color_continuous_scale='RdBu',
            title='Pathogen Hit Hierarchy'
        )
        fig.write_html(output_dir/"pathogen_sunburst.html")

    except Exception as e:
        print(f"Pathogen visualization error: {e}")

def create_sequence_quality_plots(seq_stats, output_dir):
    """Create quality and statistics plots for FASTA sequences"""
    import plotly.express as px
    import plotly.graph_objects as go
    
    # Length distribution
    fig1 = px.histogram(
        x=seq_stats['lengths'],
        nbins=50,
        title='Sequence Length Distribution',
        labels={'x': 'Sequence Length (bp)', 'y': 'Count'},
        color_discrete_sequence=['#FF6B6B']
    )
    fig1.update_layout(template="plotly_white")
    fig1.write_html(output_dir / "length_distribution.html")
    
    # GC content distribution
    fig2 = px.histogram(
        x=seq_stats['gc_contents'],
        nbins=30,
        title='GC Content Distribution',
        labels={'x': 'GC Content (%)', 'y': 'Count'},
        color_discrete_sequence=['#4ECDC4']
    )
    fig2.update_layout(template="plotly_white")
    fig2.write_html(output_dir / "gc_distribution.html")
    
    # Summary statistics
    fig3 = go.Figure()
    
    stats_data = [
        ['Total Sequences', seq_stats['total_sequences']],
        ['Total Length (bp)', f"{seq_stats['total_length']:,}"],
        ['Mean Length (bp)', f"{seq_stats['mean_length']:.0f}"],
        ['N50 (bp)', f"{seq_stats['n50']:,}"],
        ['Mean GC Content (%)', f"{seq_stats['mean_gc_content']:.1f}"],
        ['Mean N Content (%)', f"{seq_stats['mean_n_content']:.2f}"]
    ]
    
    fig3.add_trace(go.Table(
        header=dict(values=['Statistic', 'Value'],
                   fill_color='lightblue',
                   align='left'),
        cells=dict(values=list(zip(*stats_data)),
                  fill_color='lavender',
                  align='left')
    ))
    
    fig3.update_layout(title="Sequence Statistics Summary")
    fig3.write_html(output_dir / "length_summary.html")
    
    print("✓ Created sequence quality plots")