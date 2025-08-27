import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
from pathlib import Path

from .visualization import BaseVisualizer

class ComparativeVisualizer(BaseVisualizer):
    """Generates visualizations for the comparative analysis module."""

    def create_pcoa_plot(self, pcoa_coords: pd.DataFrame, metadata_df: pd.DataFrame, explained_variance: pd.Series):
        """Creates an interactive PCoA plot from beta diversity results."""
        
        # Merge coordinates with metadata to get group information for coloring
        plot_df = pcoa_coords.join(metadata_df.set_index('sample_id'))

        # Get the percentage of variance explained by the first two axes
        pc1_var = explained_variance[0] * 100
        pc2_var = explained_variance[1] * 100

        fig = px.scatter(
            plot_df,
            x='PC1',
            y='PC2',
            color='group', # Color points by the 'group' column from your metadata
            title="Beta Diversity Principal Coordinate Analysis (PCoA)",
            labels={
                "PC1": f"PC1 ({pc1_var:.2f}%)",
                "PC2": f"PC2 ({pc2_var:.2f}%)"
            },
            hover_name=plot_df.index, # Show sample ID on hover
            template="plotly_white"
        )
        
        fig.update_layout(
            height=600,
            width=800,
            legend_title_text='Group'
        )
        fig.update_traces(marker=dict(size=12, line=dict(width=1, color='DarkSlateGrey')))

        self.save_plot(fig, "beta_diversity_pcoa.html")
        print("  ✓ PCoA plot created: beta_diversity_pcoa.html")

    def create_heatmap(self, abundance_df: pd.DataFrame, top_n: int = 50):
        """Creates an interactive heatmap of the top N most abundant species."""
        
        # Calculate relative abundance (proportions) for better color scaling
        relative_abundance = abundance_df.div(abundance_df.sum(axis=1), axis=0)
        
        # Get the top N most abundant species across all samples
        top_species = abundance_df.sum().nlargest(top_n).index
        heatmap_data = relative_abundance[top_species].transpose()

        fig = px.imshow(
            heatmap_data,
            labels=dict(x="Sample ID", y="Species", color="Relative Abundance"),
            x=heatmap_data.columns,
            y=heatmap_data.index,
            aspect="auto",
            color_continuous_scale="Viridis",
            title=f"Heatmap of Top {top_n} Most Abundant Species"
        )
        
        fig.update_layout(
            height=max(600, top_n * 20), # Dynamic height
            width=800
        )

        self.save_plot(fig, "taxonomic_heatmap.html")
        print("  ✓ Heatmap visualization created: taxonomic_heatmap.html")

    def create_volcano_plot(self, diff_abundance_df: pd.DataFrame, p_value_threshold: float = 0.05, fold_change_threshold: float = 2.0):
        """Creates an interactive volcano plot to visualize differential abundance."""
        
        plot_df = diff_abundance_df.copy()
        plot_df['-log10(p-value)'] = -np.log10(plot_df['p_value'])
        plot_df['log2(FoldChange)'] = np.log2(plot_df['fold_change'])
        
        # Determine significance
        plot_df['Significance'] = 'Not Significant'
        plot_df.loc[(plot_df['p_value'] < p_value_threshold) & (plot_df['fold_change'] > fold_change_threshold), 'Significance'] = f"Upregulated (FC > {fold_change_threshold})"
        plot_df.loc[(plot_df['p_value'] < p_value_threshold) & (plot_df['fold_change'] < 1/fold_change_threshold), 'Significance'] = f"Downregulated (FC < {1/fold_change_threshold})"

        fig = px.scatter(
            plot_df,
            x='log2(FoldChange)',
            y='-log10(p-value)',
            color='Significance',
            hover_name='species',
            title='Volcano Plot of Differential Abundance',
            template='plotly_white',
            color_discrete_map={
                'Not Significant': 'grey',
                f"Upregulated (FC > {fold_change_threshold})": 'red',
                f"Downregulated (FC < {1/fold_change_threshold})": 'blue'
            }
        )
        
        fig.add_hline(y=-np.log10(p_value_threshold), line_dash="dash", line_color="black")
        fig.add_vline(x=np.log2(fold_change_threshold), line_dash="dash", line_color="black")
        fig.add_vline(x=-np.log2(fold_change_threshold), line_dash="dash", line_color="black")

        self.save_plot(fig, "volcano_plot.html")
        print("  ✓ Volcano plot created: volcano_plot.html")

    def create_alpha_diversity_boxplot(self, alpha_diversity_df: pd.DataFrame):
        """Creates a box plot comparing alpha diversity between groups."""
        
        fig = px.box(
            alpha_diversity_df,
            x='group',
            y='shannon',
            color='group',
            title='Alpha Diversity (Shannon Index) Comparison',
            labels={"group": "Group", "shannon": "Shannon Diversity Index"},
            template='plotly_white',
            points="all" # Show individual data points
        )
        
        self.save_plot(fig, "alpha_diversity_boxplot.html")
        print("  ✓ Alpha diversity box plot created: alpha_diversity_boxplot.html")