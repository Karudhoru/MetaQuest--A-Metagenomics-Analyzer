import pandas as pd
from pathlib import Path
from typing import List
import numpy as np

from skbio.diversity import beta_diversity, alpha_diversity
from skbio.stats.ordination import pcoa
from scipy.stats import mannwhitneyu
from ..visualization.compare_visuals import ComparativeVisualizer

def run_comparison(input_dirs: List[str], metadata_file: str, output_dir: str):
    """
    Main function to orchestrate the comparative analysis.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    print(f"-> Loading metadata from: {metadata_file}")
    try:
        metadata_df = pd.read_csv(metadata_file, sep='\t')
    except Exception as e:
        print(f"❌ ERROR: Could not parse metadata file. Details: {e}")
        return

    print(f"-> Found {len(metadata_df)} samples in metadata.")
    
    # --- Step 1: Aggregate Taxonomic Data ---
    print("-> Aggregating taxonomic profiles...")
    
    all_samples_data = []
    sample_ids_in_data = []
    for input_dir in input_dirs:
        dir_path = Path(input_dir)
        sample_id = dir_path.name
        bracken_report = dir_path / "bracken_report.tsv"
        if bracken_report.exists():
            try:
                sample_df = pd.read_csv(bracken_report, sep='\t')[['name', 'new_est_reads']]
                sample_df.rename(columns={'new_est_reads': sample_id}, inplace=True)
                sample_df.set_index('name', inplace=True)
                all_samples_data.append(sample_df)
                sample_ids_in_data.append(sample_id)
                print(f"  ✓ Processed: {sample_id}")
            except Exception as e:
                print(f"  ⚠️ Warning: Could not process Bracken report for '{sample_id}'. Skipping. Error: {e}")
        else:
            print(f"  ⚠️ Warning: Bracken report not found for '{sample_id}'. Skipping.")

    if not all_samples_data:
        print("❌ ERROR: No valid Bracken reports found. Cannot proceed.")
        return
        
    abundance_table = pd.concat(all_samples_data, axis=1, join='outer').fillna(0).astype(int)
    abundance_table_transposed = abundance_table.transpose()

    # --- Step 2: Clean and Validate Data ---
    print("\n-> Cleaning and validating abundance table...")

    original_sample_count = abundance_table_transposed.shape[0]
    abundance_table_cleaned = abundance_table_transposed.loc[(abundance_table_transposed.sum(axis=1) > 0)]
    samples_removed = original_sample_count - abundance_table_cleaned.shape[0]
    if samples_removed > 0:
        print(f"  -> Removed {samples_removed} empty samples.")

    original_species_count = abundance_table_cleaned.shape[1]
    abundance_table_cleaned = abundance_table_cleaned.loc[:, (abundance_table_cleaned != 0).any(axis=0)]
    species_removed = original_species_count - abundance_table_cleaned.shape[1]
    if species_removed > 0:
        print(f"  -> Removed {species_removed} species that were absent in all valid samples.")
    
    metadata_df = metadata_df[metadata_df['sample_id'].isin(abundance_table_cleaned.index)]
    
    if len(metadata_df) < 2:
        print("❌ ERROR: Fewer than two valid samples remain after cleaning. Cannot perform comparison.")
        return
        
    abundance_table_cleaned.to_csv(output_path / "taxonomic_abundance_table.tsv", sep='\t')
    print(f"✅ Cleaned abundance table created with {abundance_table_cleaned.shape[1]} species across {abundance_table_cleaned.shape[0]} samples.")

    # --- Step 3: Calculate Alpha Diversity for each sample ---
    print("\n-> Calculating Alpha Diversity...")
    abundance_table_aligned = abundance_table_cleaned.reindex(metadata_df['sample_id'])
    alpha_div_results = []
    for sample_id in abundance_table_aligned.index:
        read_counts = abundance_table_aligned.loc[sample_id].values
        shannon_diversity = alpha_diversity('shannon', [read_counts], validate=False).iloc[0]
        alpha_div_results.append({'sample_id': sample_id, 'shannon': shannon_diversity})
    
    alpha_diversity_df = pd.DataFrame(alpha_div_results).set_index('sample_id')
    alpha_diversity_df = alpha_diversity_df.join(metadata_df.set_index('sample_id'))
    alpha_diversity_df.to_csv(output_path / "alpha_diversity_results.tsv", sep='\t')
    print("  -> Alpha diversity calculations complete.")

    # --- Step 4: Calculate Beta Diversity ---
    print("\n-> Calculating Beta Diversity (Bray-Curtis)...")
    bc_matrix = beta_diversity("braycurtis", abundance_table_aligned.values, abundance_table_aligned.index)
    
    if np.isnan(bc_matrix.data).any():
        print("  -> Warning: NaNs detected in distance matrix, filling with 0.0 to proceed.")
        bc_matrix.data = np.nan_to_num(bc_matrix.data)

    # --- Step 5: Run Statistical Tests for Differential Abundance ---
    print("\n-> Running statistical tests for differential abundance...")
    
    diff_results = []
    groups = metadata_df['group'].unique()
    if len(groups) == 2:
        group1_name, group2_name = groups
        group1_samples = metadata_df[metadata_df['group'] == group1_name]['sample_id'].tolist()
        group2_samples = metadata_df[metadata_df['group'] == group2_name]['sample_id'].tolist()

        for species in abundance_table_aligned.columns:
            group1_values = abundance_table_aligned.loc[group1_samples, species]
            group2_values = abundance_table_aligned.loc[group2_samples, species]
            
            if group1_values.sum() > 0 or group2_values.sum() > 0:
                try:
                    _, p_value = mannwhitneyu(group1_values, group2_values, alternative='two-sided')
                    
                    # Calculate Fold Change
                    mean_g1 = group1_values.mean()
                    mean_g2 = group2_values.mean()
                    fold_change = (mean_g1 + 1e-6) / (mean_g2 + 1e-6) # Add pseudo-count to avoid division by zero
                    
                    enriched_in = group1_name if fold_change > 1 else group2_name
                    
                    diff_results.append({
                        'species': species,
                        'p_value': p_value,
                        'fold_change': fold_change,
                        'enriched_in_group': enriched_in,
                        f'mean_{group1_name}': mean_g1,
                        f'mean_{group2_name}': mean_g2,
                    })
                except ValueError:
                    continue
    
    diff_df = pd.DataFrame(diff_results)
    report_path = output_path / "differential_abundance_report.tsv"
    if not diff_df.empty:
        diff_df = diff_df.sort_values('p_value')
        print(f"  ✓ Found {len(diff_df[diff_df['p_value'] < 0.1])} significantly different species (p < 0.1).")
    else:
        print("  -> No significant differences found.")
    diff_df.to_csv(report_path, sep='\t', index=False)
    print(f"    -> Full report saved to: {report_path}")

    # --- Step 6: Run PCoA ---
    print("\n-> Performing Principal Coordinate Analysis (PCoA)...")
    try:
        pcoa_results = pcoa(bc_matrix)
    except Exception as e:
        print(f"❌ ERROR: PCoA calculation failed. Error: {e}")
        return
    
    # --- Step 7: Generate Visualizations ---
    print("\n-> Generating Visualizations...")
    visualizer = ComparativeVisualizer(output_path)
    visualizer.create_pcoa_plot(pcoa_results.samples, metadata_df, pcoa_results.proportion_explained)
    visualizer.create_heatmap(abundance_table_aligned)
    
    # Add calls for the new plots
    if not alpha_diversity_df.empty:
        visualizer.create_alpha_diversity_boxplot(alpha_diversity_df.reset_index())
    if not diff_df.empty:
        visualizer.create_volcano_plot(diff_df)