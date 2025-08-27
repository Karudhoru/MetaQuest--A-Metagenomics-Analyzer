import pandas as pd
from pathlib import Path
from typing import List, Dict
import numpy as np

# Core scientific libraries
from skbio.diversity import beta_diversity, alpha_diversity
from skbio.stats.ordination import pcoa
from skbio.stats.distance import permanova, anosim
from scipy.stats import mannwhitneyu, kruskal
from statsmodels.stats.multitest import multipletests
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score

from ..visualization.compare_visuals import ComparativeVisualizer

class ComparativeAnalysis:
    """
    A publication-ready comparative analysis pipeline for microbiome data.
    Enhanced with better statistical handling and small sample size considerations.
    """
    def __init__(self, output_dir: str):
        self.output_path = Path(output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        self.abundance_table = None
        self.metadata_df = None

    def run_complete_analysis(self, input_dirs: List[str], metadata_file: str):
        """Main orchestrator for the entire comparative analysis workflow."""
        print("🔬 Starting Enhanced Comparative Microbiome Analysis")
        print("=" * 60)

        if not self.load_and_aggregate_data(input_dirs, metadata_file):
            return

        # --- Diversity and Statistical Analyses ---
        alpha_div_df = self._calculate_alpha_diversity()
        bc_matrix = self._calculate_beta_diversity()
        self._perform_beta_diversity_stats(bc_matrix)
        diff_df = self._run_differential_abundance()
        self._run_ml_biomarker_discovery()

        # --- Visualization ---
        print("\n-> Generating enhanced visualizations...")
        visualizer = ComparativeVisualizer(self.output_path)
        
        if bc_matrix is not None:
            pcoa_results = pcoa(bc_matrix)
            visualizer.create_pcoa_plot(pcoa_results.samples, self.metadata_df, pcoa_results.proportion_explained)
        
        visualizer.create_heatmap(self.abundance_table)
        
        if alpha_div_df is not None and not alpha_div_df.empty:
            visualizer.create_alpha_diversity_boxplot(alpha_div_df)
            
        # Always create volcano plot, even if no significant results
        if diff_df is not None and not diff_df.empty:
            visualizer.create_volcano_plot(diff_df)
        
        # Additional visualization - abundance comparison
        visualizer.create_abundance_barplot(self.abundance_table, self.metadata_df)

        print("\n🎉 Enhanced comparative analysis completed!")
        print(f"📁 Results saved to: {self.output_path}")

    def load_and_aggregate_data(self, input_dirs: List[str], metadata_file: str) -> bool:
        """Load and aggregate taxonomic data from multiple samples."""
        print(f"-> Loading metadata from: {metadata_file}")
        try:
            self.metadata_df = pd.read_csv(metadata_file, sep='\t')
        except Exception as e:
            print(f"❌ ERROR: Could not parse metadata file. Details: {e}")
            return False

        print(f"-> Found {len(self.metadata_df)} samples in metadata.")
        
        # Aggregate taxonomic data
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
                    print(f"  ⚠ Warning: Could not process Bracken report for '{sample_id}'. Error: {e}")
            else:
                print(f"  ⚠ Warning: Bracken report not found for '{sample_id}'. Skipping.")

        if not all_samples_data:
            print("❌ ERROR: No valid Bracken reports found.")
            return False
            
        abundance_table = pd.concat(all_samples_data, axis=1, join='outer').fillna(0).astype(int)
        self.abundance_table = abundance_table.transpose()
        
        return self._clean_and_validate_data()
    
    def _clean_and_validate_data(self) -> bool:
        """Clean and validate the abundance table."""
        print("\n-> Cleaning and validating abundance table...")

        original_sample_count = self.abundance_table.shape[0]
        self.abundance_table = self.abundance_table.loc[(self.abundance_table.sum(axis=1) > 0)]
        samples_removed = original_sample_count - self.abundance_table.shape[0]
        
        if samples_removed > 0:
            print(f"  -> Removed {samples_removed} empty samples.")

        original_species_count = self.abundance_table.shape[1]
        self.abundance_table = self.abundance_table.loc[:, (self.abundance_table != 0).any(axis=0)]
        species_removed = original_species_count - self.abundance_table.shape[1]
        
        if species_removed > 0:
            print(f"  -> Removed {species_removed} species absent in all samples.")
        
        # Filter metadata to match available samples
        self.metadata_df = self.metadata_df[self.metadata_df['sample_id'].isin(self.abundance_table.index)]
        
        if len(self.metadata_df) < 2:
            print("❌ ERROR: Fewer than two valid samples remain after cleaning.")
            return False
            
        # CRITICAL FIX: Ensure proper alignment and indexing
        self.metadata_df = self.metadata_df.set_index('sample_id')
        self.abundance_table = self.abundance_table.reindex(self.metadata_df.index)
        
        # Save cleaned data
        self.abundance_table.to_csv(self.output_path / "taxonomic_abundance_table.tsv", sep='\t')
        print(f"✅ Cleaned abundance table: {self.abundance_table.shape[1]} species × {self.abundance_table.shape[0]} samples")
        
        return True

    def _calculate_alpha_diversity(self) -> pd.DataFrame:
        """Calculate multiple alpha diversity metrics for each sample."""
        print("\n-> Calculating Alpha Diversity...")
        
        metrics = ['shannon', 'simpson', 'chao1', 'sobs']
        alpha_div_results = {}
        for metric in metrics:
            try:
                alpha_div_results[metric] = alpha_diversity(metric, self.abundance_table.values, ids=self.abundance_table.index)
            except Exception as e:
                print(f"  ⚠️ Warning: Could not calculate {metric}: {e}")

        if not alpha_div_results:
            return None
            
        alpha_df = pd.DataFrame(alpha_div_results)
        
        # FIXED: Proper join with metadata
        alpha_df = alpha_df.join(self.metadata_df)
        
        # Perform statistical tests for alpha diversity
        self._perform_alpha_diversity_stats(alpha_df)
        
        alpha_df.to_csv(self.output_path / "alpha_diversity_metrics.tsv", sep='\t')
        print("  ✓ Alpha diversity calculations complete.")
        return alpha_df

    def _perform_alpha_diversity_stats(self, alpha_df):
        """Perform statistical tests on alpha diversity metrics."""
        print("\n-> Running Alpha Diversity Statistical Tests...")
        
        groups = alpha_df['group'].unique()
        if len(groups) != 2:
            print("  ⚠️ Warning: Alpha diversity tests require exactly 2 groups.")
            return
            
        group1_name, group2_name = groups
        results = []
        
        for metric in ['shannon', 'simpson', 'chao1', 'sobs']:
            if metric in alpha_df.columns:
                group1_values = alpha_df[alpha_df['group'] == group1_name][metric].values
                group2_values = alpha_df[alpha_df['group'] == group2_name][metric].values
                
                # Use Mann-Whitney U test (non-parametric)
                try:
                    statistic, p_value = mannwhitneyu(group1_values, group2_values, alternative='two-sided')
                    
                    # Calculate effect size (r = Z / sqrt(N))
                    n1, n2 = len(group1_values), len(group2_values)
                    effect_size = abs(statistic - (n1 * n2 / 2)) / np.sqrt(n1 * n2 * (n1 + n2) / 12)
                    
                    results.append({
                        'metric': metric,
                        'group1_mean': group1_values.mean(),
                        'group2_mean': group2_values.mean(),
                        'statistic': statistic,
                        'p_value': p_value,
                        'effect_size': effect_size,
                        'interpretation': self._interpret_p_value(p_value, effect_size)
                    })
                    
                    print(f"  ✓ {metric}: p = {p_value:.4f} {self._interpret_p_value(p_value, effect_size)}")
                    
                except Exception as e:
                    print(f"  ⚠️ Warning: Could not test {metric}: {e}")
        
        if results:
            results_df = pd.DataFrame(results)
            results_df.to_csv(self.output_path / "alpha_diversity_statistics.tsv", sep='\t', index=False)

    def _interpret_p_value(self, p_value: float, effect_size: float = None) -> str:
        """Interpret p-values with context for small sample sizes."""
        if p_value < 0.001:
            return "(highly significant ***)"
        elif p_value < 0.01:
            return "(significant **)"
        elif p_value < 0.05:
            return "(significant *)"
        elif p_value < 0.1:
            return "(trend/marginally significant †)"
        elif p_value < 0.2:
            return "(weak trend, consider larger n)"
        else:
            return "(not significant)"

    def _calculate_beta_diversity(self):
        """Calculate Bray-Curtis beta diversity."""
        print("\n-> Calculating Beta Diversity (Bray-Curtis)...")
        bc_matrix = beta_diversity("braycurtis", self.abundance_table.values, self.abundance_table.index)
        if np.isnan(bc_matrix.data).any():
            bc_matrix.data = np.nan_to_num(bc_matrix.data)
        return bc_matrix

    def _perform_beta_diversity_stats(self, distance_matrix):
        """Perform PERMANOVA and ANOSIM tests with proper error handling."""
        print("\n-> Running Beta Diversity Statistical Tests...")
        
        try:
            # CRITICAL FIX: Create properly aligned metadata for statistical tests
            # The distance matrix IDs should match the metadata index
            metadata_for_stats = self.metadata_df.copy()
            
            print(f"  -> Distance matrix IDs: {list(distance_matrix.ids)}")
            print(f"  -> Metadata index: {list(metadata_for_stats.index)}")
            
            # Ensure perfect alignment
            common_ids = [id for id in distance_matrix.ids if id in metadata_for_stats.index]
            if len(common_ids) != len(distance_matrix.ids):
                print(f"  ⚠️ Warning: ID mismatch. Common IDs: {len(common_ids)}, Distance matrix: {len(distance_matrix.ids)}")
            
            # PERMANOVA
            try:
                permanova_result = permanova(distance_matrix, metadata_for_stats, column='group', permutations=999)
                p_val = permanova_result['p-value']
                f_stat = permanova_result['test statistic']
                
                print(f"  ✓ PERMANOVA: F = {f_stat:.3f}, p = {p_val:.4f} {self._interpret_p_value(p_val)}")
                
                # Save detailed results
                with open(self.output_path / "permanova_results.txt", 'w') as f:
                    f.write(f"PERMANOVA Results\n")
                    f.write(f"================\n")
                    f.write(f"F-statistic: {f_stat:.6f}\n")
                    f.write(f"p-value: {p_val:.6f}\n")
                    f.write(f"Interpretation: {self._interpret_p_value(p_val)}\n")
                    f.write(f"Permutations: 999\n")
                    f.write(f"Test: Groups differ in microbiome composition\n")
                
            except Exception as e:
                print(f"  ❌ PERMANOVA failed: {e}")

            # ANOSIM
            try:
                anosim_result = anosim(distance_matrix, metadata_for_stats, column='group', permutations=999)
                r_stat = anosim_result['test statistic']
                p_val = anosim_result['p-value']
                
                print(f"  ✓ ANOSIM: R = {r_stat:.3f}, p = {p_val:.4f} {self._interpret_p_value(p_val)}")
                
                # Save detailed results
                with open(self.output_path / "anosim_results.txt", 'w') as f:
                    f.write(f"ANOSIM Results\n")
                    f.write(f"==============\n")
                    f.write(f"R-statistic: {r_stat:.6f}\n")
                    f.write(f"p-value: {p_val:.6f}\n")
                    f.write(f"Interpretation: {self._interpret_p_value(p_val)}\n")
                    f.write(f"Permutations: 999\n")
                    f.write(f"R-value interpretation:\n")
                    f.write(f"  R > 0.75: Well separated\n")
                    f.write(f"  R > 0.5: Overlapping but clearly different\n")
                    f.write(f"  R > 0.25: Barely separable\n")
                    f.write(f"  R ≈ 0: Indistinguishable\n")
                    
            except Exception as e:
                print(f"  ❌ ANOSIM failed: {e}")
                
        except Exception as e:
            print(f"  ❌ Beta diversity statistical tests failed: {e}")

    def _run_differential_abundance(self) -> pd.DataFrame:
        """Run differential abundance testing with multiple significance levels."""
        print("\n-> Running differential abundance analysis...")
        
        groups = self.metadata_df['group'].unique()
        if len(groups) != 2:
            print("  ⚠️ Warning: Differential abundance test requires exactly 2 groups. Skipping.")
            return None
            
        group1_name, group2_name = groups
        group1_samples = self.metadata_df[self.metadata_df['group'] == group1_name].index
        group2_samples = self.metadata_df[self.metadata_df['group'] == group2_name].index

        print(f"  -> Comparing {group1_name} (n={len(group1_samples)}) vs {group2_name} (n={len(group2_samples)})")

        results = []
        for species in self.abundance_table.columns:
            g1_vals = self.abundance_table.loc[group1_samples, species]
            g2_vals = self.abundance_table.loc[group2_samples, species]
            
            # Skip species with no reads in either group
            if g1_vals.sum() == 0 and g2_vals.sum() == 0: 
                continue
            
            try:
                # Mann-Whitney U test
                _, p_value = mannwhitneyu(g1_vals, g2_vals, alternative='two-sided')
                
                mean_g1 = g1_vals.mean()
                mean_g2 = g2_vals.mean()
                fold_change = (mean_g1 + 1e-6) / (mean_g2 + 1e-6)
                log2_fc = np.log2(fold_change)
                
                # Calculate additional metrics
                median_g1 = g1_vals.median()
                median_g2 = g2_vals.median()
                
                results.append({
                    'species': species,
                    'p_value': p_value,
                    'fold_change': fold_change,
                    'log2_fold_change': log2_fc,
                    'mean_group1': mean_g1,
                    'mean_group2': mean_g2,
                    'median_group1': median_g1,
                    'median_group2': median_g2,
                    'group1_name': group1_name,
                    'group2_name': group2_name
                })
            except Exception as e:
                print(f"  ⚠️ Warning: Could not test {species}: {e}")
        
        if not results:
            print("  -> No species were eligible for testing.")
            return None

        diff_df = pd.DataFrame(results)
        
        # Multiple testing correction with different methods
        if len(diff_df) > 0:
            # FDR correction
            try:
                rejected_fdr, p_fdr, _, _ = multipletests(diff_df['p_value'], alpha=0.05, method='fdr_bh')
                diff_df['p_fdr'] = p_fdr
                diff_df['significant_fdr'] = rejected_fdr
            except:
                diff_df['p_fdr'] = diff_df['p_value']
                diff_df['significant_fdr'] = False
            
            # Bonferroni correction
            try:
                rejected_bonf, p_bonf, _, _ = multipletests(diff_df['p_value'], alpha=0.05, method='bonferroni')
                diff_df['p_bonferroni'] = p_bonf
                diff_df['significant_bonferroni'] = rejected_bonf
            except:
                diff_df['p_bonferroni'] = diff_df['p_value']
                diff_df['significant_bonferroni'] = False
        
        # Add significance categories
        diff_df['significance_level'] = diff_df['p_value'].apply(
            lambda p: 'highly_significant' if p < 0.001 
            else 'significant' if p < 0.01
            else 'trend' if p < 0.1
            else 'not_significant'
        )
        
        diff_df = diff_df.sort_values('p_value')
        report_path = self.output_path / "differential_abundance_report.tsv"
        diff_df.to_csv(report_path, sep='\t', index=False)
        
        # Enhanced reporting
        fdr_sig = diff_df['significant_fdr'].sum()
        bonf_sig = diff_df['significant_bonferroni'].sum()
        trend_sig = (diff_df['p_value'] < 0.1).sum()
        nominal_sig = (diff_df['p_value'] < 0.05).sum()
        
        print(f"  ✓ Differential abundance results:")
        print(f"    • FDR significant (q < 0.05): {fdr_sig}")
        print(f"    • Bonferroni significant: {bonf_sig}")
        print(f"    • Nominally significant (p < 0.05): {nominal_sig}")
        print(f"    • Trend level (p < 0.10): {trend_sig}")
        
        if trend_sig > 0:
            print(f"  -> Top candidates (p < 0.10):")
            top_candidates = diff_df[diff_df['p_value'] < 0.10].head(5)
            for _, row in top_candidates.iterrows():
                direction = "↑" if row['fold_change'] > 1 else "↓"
                print(f"    • {row['species']}: p = {row['p_value']:.4f}, FC = {row['fold_change']:.2f} {direction}")
        
        return diff_df

    def _run_ml_biomarker_discovery(self):
        """Use Random Forest to find the most discriminative species."""
        print("\n-> Running Machine Learning for Biomarker Discovery...")
        try:
            X = self.abundance_table
            y = self.metadata_df.loc[X.index, 'group']
            
            # Ensure we have the right data types
            X = X.astype(float)
            
            rf = RandomForestClassifier(n_estimators=100, random_state=42, 
                                      class_weight='balanced')  # Handle class imbalance
            
            # Cross-validation with proper scoring
            cv_scores = cross_val_score(rf, X, y, cv=min(3, len(X)//2), scoring='accuracy')
            
            print(f"  ✓ Random Forest {len(cv_scores)}-fold CV Accuracy: {np.mean(cv_scores):.2f} ± {np.std(cv_scores):.2f}")
            
            # Fit final model
            rf.fit(X, y)
            
            # Feature importance analysis
            importance_df = pd.DataFrame({
                'species': X.columns,
                'importance': rf.feature_importances_
            }).sort_values('importance', ascending=False)
            
            # Add abundance information for top features
            top_features = importance_df.head(20)
            for idx, row in top_features.iterrows():
                species = row['species']
                group_means = X.groupby(y)[species].mean()
                importance_df.loc[idx, 'group1_mean'] = group_means.iloc[0]
                importance_df.loc[idx, 'group2_mean'] = group_means.iloc[1]
                importance_df.loc[idx, 'fold_change'] = (group_means.iloc[0] + 1e-6) / (group_means.iloc[1] + 1e-6)
            
            report_path = self.output_path / "random_forest_feature_importance.tsv"
            importance_df.to_csv(report_path, sep='\t', index=False)
            
            print(f" ✓ Top discriminative species:")
            for i in range(min(2, len(importance_df))):
                print(f"    -> '{importance_df.iloc[i]['species']}' - Importance: {importance_df.iloc[i]['importance']:.4f}")
            print(f"-> Feature importance report saved to: {report_path}")
                
        except Exception as e:
            print(f"  ⚠️ Warning: ML biomarker discovery failed. Error: {e}")

def run_comparison(input_dirs: List[str], metadata_file: str, output_dir: str):
    """Main entry point to run the enhanced comparative analysis."""
    analyzer = ComparativeAnalysis(output_dir)
    analyzer.run_complete_analysis(input_dirs, metadata_file)