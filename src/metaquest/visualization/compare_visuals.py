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
        
        # FIXED: Handle metadata_df that's already indexed by sample_id
        if 'sample_id' in metadata_df.columns:
            # Original format - has sample_id as column
            plot_df = pcoa_coords.join(metadata_df.set_index('sample_id'))
        else:
            # New format - already indexed by sample_id
            plot_df = pcoa_coords.join(metadata_df)

        # Get the percentage of variance explained by the first two axes
        # FIXED: Handle both Series and ndarray types for explained_variance
        if hasattr(explained_variance, 'iloc'):
            pc1_var = explained_variance.iloc[0] * 100
            pc2_var = explained_variance.iloc[1] * 100
        else:
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
            template="plotly_white",
            color_discrete_sequence=px.colors.qualitative.Set1
        )
        
        fig.update_layout(
            height=600,
            width=800,
            legend_title_text='Group',
            title_x=0.5,  # Center the title
            font=dict(size=12)
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

        # Clean up species names for better readability (keep only species name)
        cleaned_species_names = []
        for species in heatmap_data.index:
            # Extract just the species name (last two parts after splitting by space)
            parts = species.split()
            if len(parts) >= 2:
                cleaned_name = f"{parts[-2]} {parts[-1]}"
            else:
                cleaned_name = species
            cleaned_species_names.append(cleaned_name)
        
        heatmap_data.index = cleaned_species_names

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
            height=max(600, min(top_n * 15, 1200)), # Dynamic height with max limit
            width=900,
            title_x=0.5,
            font=dict(size=10),
            yaxis=dict(tickfont=dict(size=8))  # Smaller font for species names
        )

        self.save_plot(fig, "taxonomic_heatmap.html")
        print("  ✓ Heatmap visualization created: taxonomic_heatmap.html")

    def create_volcano_plot(self, diff_abundance_df: pd.DataFrame, p_value_threshold: float = 0.1, fold_change_threshold: float = 2.0):
        """Creates an interactive volcano plot to visualize differential abundance."""
        
        if diff_abundance_df is None or diff_abundance_df.empty:
            print("  ⚠️ Warning: No differential abundance data available for volcano plot")
            return
        
        plot_df = diff_abundance_df.copy()
        
        # Handle edge cases in p-values
        plot_df['p_value'] = plot_df['p_value'].replace(0, 1e-300)  # Avoid log(0)
        plot_df['-log10(p-value)'] = -np.log10(plot_df['p_value'])
        
        # Handle log2 fold change
        if 'log2_fold_change' in plot_df.columns:
            plot_df['log2(FoldChange)'] = plot_df['log2_fold_change']
        else:
            plot_df['log2(FoldChange)'] = np.log2(plot_df['fold_change'])
        
        # FIXED: Use more lenient thresholds for small sample studies
        # Determine significance categories with multiple levels
        plot_df['Significance'] = 'Not Significant'
        
        # Highly significant (p < 0.01)
        highly_sig_up = (plot_df['p_value'] < 0.01) & (plot_df['fold_change'] > fold_change_threshold)
        highly_sig_down = (plot_df['p_value'] < 0.01) & (plot_df['fold_change'] < 1/fold_change_threshold)
        
        # Significant (p < 0.05)
        sig_up = (plot_df['p_value'] < 0.05) & (plot_df['fold_change'] > fold_change_threshold)
        sig_down = (plot_df['p_value'] < 0.05) & (plot_df['fold_change'] < 1/fold_change_threshold)
        
        # Trend level (p < 0.1) - Important for small sample studies
        trend_up = (plot_df['p_value'] < p_value_threshold) & (plot_df['fold_change'] > fold_change_threshold)
        trend_down = (plot_df['p_value'] < p_value_threshold) & (plot_df['fold_change'] < 1/fold_change_threshold)
        
        plot_df.loc[trend_up, 'Significance'] = f"Trend Up (p < {p_value_threshold})"
        plot_df.loc[trend_down, 'Significance'] = f"Trend Down (p < {p_value_threshold})"
        plot_df.loc[sig_up, 'Significance'] = "Significant Up (p < 0.05)"
        plot_df.loc[sig_down, 'Significance'] = "Significant Down (p < 0.05)"
        plot_df.loc[highly_sig_up, 'Significance'] = "Highly Significant Up (p < 0.01)"
        plot_df.loc[highly_sig_down, 'Significance'] = "Highly Significant Down (p < 0.01)"

        # Create color mapping
        color_map = {
            'Not Significant': 'lightgrey',
            f"Trend Up (p < {p_value_threshold})": 'lightcoral',
            f"Trend Down (p < {p_value_threshold})": 'lightblue',
            'Significant Up (p < 0.05)': 'red',
            'Significant Down (p < 0.05)': 'blue',
            'Highly Significant Up (p < 0.01)': 'darkred',
            'Highly Significant Down (p < 0.01)': 'darkblue'
        }

        # Clean species names for hover text
        plot_df['clean_species'] = plot_df['species'].apply(lambda x: ' '.join(x.split()[-2:]) if len(x.split()) >= 2 else x)

        fig = px.scatter(
            plot_df,
            x='log2(FoldChange)',
            y='-log10(p-value)',
            color='Significance',
            hover_name='clean_species',
            hover_data={
                'p_value': ':.4f',
                'fold_change': ':.2f',
                'Significance': False,
                'log2(FoldChange)': False,
                '-log10(p-value)': False
            },
            title='Volcano Plot of Differential Abundance',
            template='plotly_white',
            color_discrete_map=color_map
        )
        
        # Add significance threshold lines
        fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="black", 
                     annotation_text="p = 0.05", annotation_position="right")
        fig.add_hline(y=-np.log10(p_value_threshold), line_dash="dot", line_color="gray", 
                     annotation_text=f"p = {p_value_threshold}", annotation_position="right")
        fig.add_vline(x=np.log2(fold_change_threshold), line_dash="dash", line_color="black")
        fig.add_vline(x=-np.log2(fold_change_threshold), line_dash="dash", line_color="black")
        
        fig.update_layout(
            height=600,
            width=900,
            title_x=0.5,
            xaxis_title="log₂(Fold Change)",
            yaxis_title="-log₁₀(p-value)",
            legend_title="Significance Level"
        )

        self.save_plot(fig, "volcano_plot.html")
        print("  ✓ Volcano plot created: volcano_plot.html")

    def create_alpha_diversity_boxplot(self, alpha_diversity_df: pd.DataFrame):
        """Creates a comprehensive box plot comparing alpha diversity between groups."""
        
        if alpha_diversity_df is None or alpha_diversity_df.empty:
            print("  ⚠️ Warning: No alpha diversity data available for box plot")
            return
        
        # Available metrics to plot
        metrics = ['shannon', 'simpson', 'chao1', 'sobs']
        available_metrics = [m for m in metrics if m in alpha_diversity_df.columns]
        
        if not available_metrics:
            print("  ⚠️ Warning: No recognized alpha diversity metrics found")
            return
        
        # Create subplot layout
        from plotly.subplots import make_subplots
        
        n_metrics = len(available_metrics)
        cols = 2
        rows = (n_metrics + 1) // 2
        
        metric_titles = {
            'shannon': 'Shannon Diversity Index',
            'simpson': 'Simpson Diversity Index', 
            'chao1': 'Chao1 Species Richness',
            'sobs': 'Observed Species'
        }
        
        fig = make_subplots(
            rows=rows, cols=cols,
            subplot_titles=[metric_titles.get(m, m.title()) for m in available_metrics],
            vertical_spacing=0.12,
            horizontal_spacing=0.1
        )
        
        colors = px.colors.qualitative.Set1
        
        for i, metric in enumerate(available_metrics):
            row = (i // cols) + 1
            col = (i % cols) + 1
            
            # Get unique groups and assign colors
            groups = alpha_diversity_df['group'].unique()
            
            for j, group in enumerate(groups):
                group_data = alpha_diversity_df[alpha_diversity_df['group'] == group]
                
                fig.add_trace(
                    go.Box(
                        y=group_data[metric],
                        name=group,
                        boxpoints='all',
                        jitter=0.3,
                        pointpos=-1.8,
                        marker_color=colors[j % len(colors)],
                        showlegend=(i == 0)  # Only show legend for first subplot
                    ),
                    row=row, col=col
                )
        
        fig.update_layout(
            height=400 * rows,
            width=800,
            title_text="Alpha Diversity Comparison Between Groups",
            title_x=0.5,
            template='plotly_white'
        )
        
        # Update y-axis labels
        for i, metric in enumerate(available_metrics):
            row = (i // cols) + 1
            col = (i % cols) + 1
            fig.update_yaxes(title_text=metric_titles.get(metric, metric.title()), row=row, col=col)
            fig.update_xaxes(title_text="Group", row=row, col=col)

        self.save_plot(fig, "alpha_diversity_boxplot.html")
        print("  ✓ Alpha diversity box plot created: alpha_diversity_boxplot.html")
        
    def create_abundance_barplot(self, abundance_df: pd.DataFrame, metadata_df: pd.DataFrame, top_n: int = 20):
        """Creates a stacked bar plot showing top species abundance by group."""
        
        # Handle metadata indexing
        if 'sample_id' in metadata_df.columns:
            metadata_df = metadata_df.set_index('sample_id')
        
        # Get top N species
        top_species = abundance_df.sum().nlargest(top_n).index
        top_abundance = abundance_df[top_species]
        
        # Calculate relative abundance
        rel_abundance = top_abundance.div(top_abundance.sum(axis=1), axis=0)
        
        # Add group information
        plot_data = rel_abundance.join(metadata_df[['group']])
        
        # Group by condition and calculate mean relative abundance
        grouped_data = plot_data.groupby('group')[top_species].mean()
        
        # Create stacked bar plot
        fig = go.Figure()
        
        colors = px.colors.qualitative.Set3
        
        for i, species in enumerate(top_species):
            # Clean species name
            clean_name = ' '.join(species.split()[-2:]) if len(species.split()) >= 2 else species
            
            fig.add_trace(go.Bar(
                name=clean_name,
                x=grouped_data.index,
                y=grouped_data[species],
                marker_color=colors[i % len(colors)]
            ))
        
        fig.update_layout(
            barmode='stack',
            title=f'Mean Relative Abundance of Top {top_n} Species by Group',
            xaxis_title='Group',
            yaxis_title='Mean Relative Abundance',
            height=600,
            width=900,
            template='plotly_white',
            title_x=0.5,
            legend=dict(
                orientation="v",
                yanchor="top",
                y=1,
                xanchor="left",
                x=1.02
            )
        )
        
        self.save_plot(fig, "abundance_barplot.html")
        print("  ✓ Abundance bar plot created: abundance_barplot.html")