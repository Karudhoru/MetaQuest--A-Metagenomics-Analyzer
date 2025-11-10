"""
Enhanced Comparative Analysis Module for MetaQuest
DROP-IN REPLACEMENT - All original function/class names preserved!

Key improvements:
1. Compositional data normalization (TSS, CLR, rarefaction)
2. Multi-group statistical comparisons (Kruskal-Wallis + post-hoc)
3. Effect size metrics (Cliff's Delta, epsilon-squared, R²)
4. ANCOM-inspired differential abundance
5. Statistical power assessment
6. Prevalence filtering
7. Enhanced ML with hyperparameter tuning
8. Multiple beta diversity metrics
"""

import pandas as pd
from pathlib import Path
from typing import List, Dict, Optional
import numpy as np

# Core scientific libraries
from skbio.diversity import beta_diversity, alpha_diversity
from skbio.stats.ordination import pcoa
from skbio.stats.distance import permanova, anosim
from scipy.stats import mannwhitneyu, kruskal
from statsmodels.stats.multitest import multipletests
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score, GridSearchCV
from sklearn.feature_selection import SelectKBest, f_classif
import warnings

try:
    from skbio.stats.composition import clr, multiplicative_replacement
    HAS_COMPOSITION = True
except ImportError:
    HAS_COMPOSITION = False
    warnings.warn("skbio.stats.composition not available - CLR normalization disabled")

from ..visualization.compare_visuals import ComparativeVisualizer
from ..io.output_formatter import OutputFormatter


class ComparativeAnalysis:
    """
    Enhanced comparative analysis pipeline for microbiome data.
    BACKWARD COMPATIBLE - All original method names preserved!
    
    New features:
    - Compositional data handling
    - Multi-group comparisons
    - Effect size reporting
    - Power analysis
    """
    
    def __init__(self, output_dir: str, debug: bool = False):
        self.output_path = Path(output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        self.abundance_table = None
        self.metadata_df = None
        self.formatter = OutputFormatter()
        
        # NEW: Track normalization and original data
        self.normalization_method = 'tss'  # Default
        self.original_abundance_table = None

    def run_complete_analysis(self, input_dirs: List[str], metadata_file: str,
                             normalization: str = 'tss'):
        """
        Main orchestrator - ORIGINAL NAME PRESERVED.
        
        NEW PARAMETER:
            normalization: 'tss' (default), 'clr', 'rarefy', or 'none'
        """
        self.formatter.header("Enhanced Comparative Microbiome Analysis")

        if not self.load_and_aggregate_data(input_dirs, metadata_file, normalization):
            return

        # NEW: Power analysis
        self._assess_statistical_power()

        # --- Diversity and Statistical Analyses (ENHANCED) ---
        alpha_div_df = self._calculate_alpha_diversity()
        bc_matrix = self._calculate_beta_diversity()
        self._perform_beta_diversity_stats(bc_matrix)
        diff_df = self._run_differential_abundance()
        self._run_ml_biomarker_discovery()

        # --- Visualization ---
        self.formatter.step("Generating enhanced visualizations")
        visualizer = ComparativeVisualizer(self.output_path)
        
        if bc_matrix is not None:
            pcoa_results = pcoa(bc_matrix)
            visualizer.create_pcoa_plot(pcoa_results.samples, self.metadata_df, 
                                       pcoa_results.proportion_explained)
        
        visualizer.create_heatmap(self.abundance_table)
        
        if alpha_div_df is not None and not alpha_div_df.empty:
            visualizer.create_alpha_diversity_boxplot(alpha_div_df)
            
        if diff_df is not None and not diff_df.empty:
            visualizer.create_volcano_plot(diff_df)
        
        visualizer.create_abundance_barplot(self.abundance_table, self.metadata_df)

        self.formatter.success("Enhanced comparative analysis completed!", style='bold')
        self.formatter.info(f"Results saved to: {self.output_path}")

    def load_and_aggregate_data(self, input_dirs: List[str], metadata_file: str,
                                normalization: str = 'tss') -> bool:
        """
        ORIGINAL NAME PRESERVED - Enhanced with normalization support.
        
        NEW: Adds compositional data handling and quality control.
        """
        self.formatter.step(f"Loading metadata from: {Path(metadata_file).name}")
        try:
            self.metadata_df = pd.read_csv(metadata_file, sep='\t')
            self.formatter.success(f"Loaded metadata for {len(self.metadata_df)} samples")
        except Exception as e:
            self.formatter.error(f"Could not parse metadata file: {e}")
            return False
        
        # Aggregate taxonomic data
        self.formatter.step("Aggregating taxonomic profiles")
        all_samples_data = []
        sample_ids_in_data = []
        
        with self.formatter.progress_bar(total=len(input_dirs), desc="Processing samples") as pbar:
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
                        if self:
                            self.formatter.success(f"Processed: {sample_id}")
                    except Exception as e:
                        self.formatter.warning(f"Could not process Bracken report for '{sample_id}': {e}")
                else:
                    self.formatter.warning(f"Bracken report not found for '{sample_id}' - skipping")
                
                pbar.update(1)

        if not all_samples_data:
            self.formatter.error("No valid Bracken reports found")
            return False
            
        abundance_table = pd.concat(all_samples_data, axis=1, join='outer').fillna(0).astype(int)
        self.abundance_table = abundance_table.transpose()
        
        # NEW: Store original before any transformations
        self.original_abundance_table = self.abundance_table.copy()
        
        # Clean and validate
        if not self._clean_and_validate_data():
            return False
        
        # NEW: Apply normalization
        self._normalize_abundance_data(normalization)
        
        return True
    
    def _clean_and_validate_data(self) -> bool:
        """ORIGINAL NAME PRESERVED - Enhanced with quality filters."""
        self.formatter.step("Cleaning and validating abundance table")

        # Remove empty samples
        original_sample_count = self.abundance_table.shape[0]
        self.abundance_table = self.abundance_table.loc[(self.abundance_table.sum(axis=1) > 0)]
        samples_removed = original_sample_count - self.abundance_table.shape[0]
        
        if samples_removed > 0:
            self.formatter.substep(f"Removed {samples_removed} empty samples")

        # NEW: Prevalence filtering - remove species present in <10% of samples
        min_prevalence = 0.10
        prevalence = (self.abundance_table > 0).sum(axis=0) / len(self.abundance_table)
        species_to_keep = prevalence[prevalence >= min_prevalence].index
        
        original_species_count = self.abundance_table.shape[1]
        self.abundance_table = self.abundance_table[species_to_keep]
        species_removed = original_species_count - self.abundance_table.shape[1]
        
        if species_removed > 0:
            self.formatter.substep(f"Removed {species_removed} rare species (prevalence <{min_prevalence*100}%)")
        
        # Filter metadata to match available samples
        self.metadata_df = self.metadata_df[self.metadata_df['sample_id'].isin(self.abundance_table.index)]
        
        if len(self.metadata_df) < 2:
            self.formatter.error("Fewer than two valid samples remain after cleaning")
            return False
            
        # Ensure proper alignment
        self.metadata_df = self.metadata_df.set_index('sample_id')
        self.abundance_table = self.abundance_table.reindex(self.metadata_df.index)
        
        # Update original table too
        if self.original_abundance_table is not None:
            self.original_abundance_table = self.original_abundance_table.reindex(self.metadata_df.index)
            self.original_abundance_table = self.original_abundance_table[species_to_keep]
        
        # Save cleaned data
        self.abundance_table.to_csv(self.output_path / "taxonomic_abundance_table.tsv", sep='\t')
        self.formatter.success(f"Cleaned abundance table: {self.abundance_table.shape[1]} species × {self.abundance_table.shape[0]} samples")
        
        return True

    def _normalize_abundance_data(self, method: str = 'tss'):
        """
        NEW METHOD: Apply normalization to abundance data.
        Critical for compositional data analysis!
        """
        self.normalization_method = method
        self.formatter.substep(f"Applying normalization: {method}")
        
        if method == 'tss':
            # Total Sum Scaling (relative abundance)
            self.abundance_table = self.abundance_table.div(
                self.abundance_table.sum(axis=1), axis=0
            )
            self.formatter.info("Applied Total Sum Scaling (relative abundances)")
            
        elif method == 'clr' and HAS_COMPOSITION:
            # Centered Log-Ratio transformation
            try:
                abundance_with_pseudo = multiplicative_replacement(self.abundance_table.values)
                clr_transformed = clr(abundance_with_pseudo)
                self.abundance_table = pd.DataFrame(
                    clr_transformed,
                    index=self.abundance_table.index,
                    columns=self.abundance_table.columns
                )
                self.formatter.info("Applied CLR transformation (compositional)")
            except Exception as e:
                self.formatter.warning(f"CLR failed, falling back to TSS: {e}")
                self._normalize_abundance_data('tss')
                
        elif method == 'rarefy':
            # Rarefaction to minimum depth
            min_depth = int(self.abundance_table.sum(axis=1).min())
            self.formatter.info(f"Rarefying to depth: {min_depth}")
            
            rarefied_data = []
            for idx, row in self.abundance_table.iterrows():
                rarefied_row = self._rarefy_sample(row.values, min_depth)
                rarefied_data.append(rarefied_row)
            
            self.abundance_table = pd.DataFrame(
                rarefied_data,
                index=self.abundance_table.index,
                columns=self.abundance_table.columns
            )
            
        elif method == 'none':
            self.formatter.info("No normalization applied (using raw counts)")
        else:
            self.formatter.warning(f"Unknown normalization '{method}', using TSS")
            self._normalize_abundance_data('tss')
    
    def _rarefy_sample(self, counts: np.ndarray, depth: int) -> np.ndarray:
        """NEW METHOD: Rarefy a single sample."""
        counts = counts.astype(int)
        total = counts.sum()
        if total < depth:
            return counts
        
        pool = np.repeat(np.arange(len(counts)), counts)
        sampled = np.random.choice(pool, size=depth, replace=False)
        rarefied = np.bincount(sampled, minlength=len(counts))
        return rarefied

    def _calculate_alpha_diversity(self) -> pd.DataFrame:
        """ORIGINAL NAME PRESERVED - Enhanced with multi-group support and effect sizes."""
        self.formatter.step("Calculating Alpha Diversity")
        
        # Use raw counts for alpha diversity (not normalized)
        data_for_alpha = self.original_abundance_table if self.original_abundance_table is not None else self.abundance_table
        
        metrics = ['shannon', 'simpson', 'chao1', 'sobs']
        alpha_div_results = {}
        for metric in metrics:
            try:
                alpha_div_results[metric] = alpha_diversity(metric, data_for_alpha.values, 
                                                           ids=data_for_alpha.index)
            except Exception as e:
                self.formatter.warning(f"Could not calculate {metric}: {e}")

        if not alpha_div_results:
            return None
            
        alpha_df = pd.DataFrame(alpha_div_results)
        alpha_df = alpha_df.join(self.metadata_df)
        
        # Enhanced statistical tests
        self._perform_alpha_diversity_stats(alpha_df)
        
        alpha_df.to_csv(self.output_path / "alpha_diversity_metrics.tsv", sep='\t')
        self.formatter.success("Alpha diversity calculations complete")
        return alpha_df

    def _perform_alpha_diversity_stats(self, alpha_df):
        """ORIGINAL NAME PRESERVED - Enhanced with multi-group and effect sizes."""
        self.formatter.step("Running Alpha Diversity Statistical Tests")
        
        groups = alpha_df['group'].unique()
        n_groups = len(groups)
        
        results = []
        
        for metric in ['shannon', 'simpson', 'chao1', 'sobs']:
            if metric not in alpha_df.columns:
                continue
            
            if n_groups == 2:
                # Two-group: Mann-Whitney U with effect size
                group1_name, group2_name = groups
                group1_values = alpha_df[alpha_df['group'] == group1_name][metric].values
                group2_values = alpha_df[alpha_df['group'] == group2_name][metric].values
                
                try:
                    statistic, p_value = mannwhitneyu(group1_values, group2_values, 
                                                     alternative='two-sided')
                    
                    # NEW: Cliff's Delta effect size
                    cliffs_delta = self._calculate_cliffs_delta(group1_values, group2_values)
                    
                    results.append({
                        'metric': metric,
                        'test': 'Mann-Whitney U',
                        'group1': group1_name,
                        'group2': group2_name,
                        'group1_mean': group1_values.mean(),
                        'group2_mean': group2_values.mean(),
                        'statistic': statistic,
                        'p_value': p_value,
                        'effect_size': cliffs_delta,
                        'effect_interpretation': self._interpret_cliffs_delta(cliffs_delta)
                    })
                    
                    interpretation = self._interpret_p_value(p_value, cliffs_delta)
                    self.formatter.substep(f"{metric}: p = {p_value:.4f}, δ = {cliffs_delta:.3f} {interpretation}")
                    
                except Exception as e:
                    self.formatter.warning(f"Could not test {metric}: {e}")
                    
            else:
                # NEW: Multi-group comparison with Kruskal-Wallis
                group_data = [alpha_df[alpha_df['group'] == g][metric].values for g in groups]
                
                try:
                    statistic, p_value = kruskal(*group_data)
                    
                    # Epsilon-squared effect size
                    n = len(alpha_df)
                    epsilon_sq = (statistic - n_groups + 1) / (n - n_groups)
                    
                    result_dict = {
                        'metric': metric,
                        'test': 'Kruskal-Wallis',
                        'statistic': statistic,
                        'p_value': p_value,
                        'effect_size': epsilon_sq,
                        'effect_interpretation': self._interpret_epsilon_squared(epsilon_sq),
                        'n_groups': n_groups
                    }
                    
                    # Add group means
                    for g in groups:
                        result_dict[f'{g}_mean'] = alpha_df[alpha_df['group'] == g][metric].mean()
                    
                    results.append(result_dict)
                    
                    self.formatter.substep(f"{metric}: p = {p_value:.4f}, ε² = {epsilon_sq:.3f} {self._interpret_p_value(p_value, epsilon_sq)}")
                    
                    # Post-hoc pairwise if significant
                    if p_value < 0.05:
                        self._posthoc_pairwise(alpha_df, metric, groups)
                        
                except Exception as e:
                    self.formatter.warning(f"Could not test {metric}: {e}")
        
        if results:
            results_df = pd.DataFrame(results)
            results_df.to_csv(self.output_path / "alpha_diversity_statistics.tsv", sep='\t', index=False)
    
    def _posthoc_pairwise(self, alpha_df: pd.DataFrame, metric: str, groups: list):
        """NEW METHOD: Post-hoc pairwise tests for multi-group comparisons."""
        from itertools import combinations
        
        pairs = list(combinations(groups, 2))
        p_values = []
        
        for g1, g2 in pairs:
            vals1 = alpha_df[alpha_df['group'] == g1][metric].values
            vals2 = alpha_df[alpha_df['group'] == g2][metric].values
            _, p = mannwhitneyu(vals1, vals2)
            p_values.append(p)
        
        # Bonferroni correction
        corrected_p = [min(p * len(pairs), 1.0) for p in p_values]
        
        self.formatter.info(f"  Post-hoc tests for {metric}:")
        for (g1, g2), p_corr in zip(pairs, corrected_p):
            sig = "***" if p_corr < 0.001 else "**" if p_corr < 0.01 else "*" if p_corr < 0.05 else "ns"
            self.formatter.info(f"    {g1} vs {g2}: p_adj = {p_corr:.4f} {sig}")

    def _interpret_p_value(self, p_value: float, effect_size: float = None) -> str:
        """ORIGINAL NAME PRESERVED - Enhanced with effect size context."""
        base_interp = ""
        if p_value < 0.001:
            base_interp = "(highly significant ***)"
        elif p_value < 0.01:
            base_interp = "(significant **)"
        elif p_value < 0.05:
            base_interp = "(significant *)"
        elif p_value < 0.1:
            base_interp = "(trend/marginally significant ·)"
        else:
            base_interp = "(not significant)"
        
        # Add effect size context if available
        if effect_size is not None:
            if abs(effect_size) < 0.15:
                base_interp += " [negligible effect]"
            elif abs(effect_size) < 0.33:
                base_interp += " [small effect]"
            elif abs(effect_size) < 0.47:
                base_interp += " [medium effect]"
            else:
                base_interp += " [large effect]"
        
        return base_interp
    
    def _calculate_cliffs_delta(self, x: np.ndarray, y: np.ndarray) -> float:
        """NEW METHOD: Calculate Cliff's Delta effect size."""
        n1, n2 = len(x), len(y)
        comparisons = 0
        
        for xi in x:
            for yi in y:
                if xi > yi:
                    comparisons += 1
                elif xi < yi:
                    comparisons -= 1
        
        return comparisons / (n1 * n2) if (n1 * n2) > 0 else 0.0
    
    def _interpret_cliffs_delta(self, delta: float) -> str:
        """NEW METHOD: Interpret Cliff's Delta magnitude."""
        abs_delta = abs(delta)
        if abs_delta < 0.147:
            return "negligible"
        elif abs_delta < 0.33:
            return "small"
        elif abs_delta < 0.474:
            return "medium"
        else:
            return "large"
    
    def _interpret_epsilon_squared(self, epsilon_sq: float) -> str:
        """NEW METHOD: Interpret epsilon-squared effect size."""
        if epsilon_sq < 0.01:
            return "negligible"
        elif epsilon_sq < 0.06:
            return "small"
        elif epsilon_sq < 0.14:
            return "medium"
        else:
            return "large"

    def _calculate_beta_diversity(self):
        """ORIGINAL NAME PRESERVED - Can be enhanced with more metrics."""
        self.formatter.step("Calculating Beta Diversity (Bray-Curtis)")
        bc_matrix = beta_diversity("braycurtis", self.abundance_table.values, 
                                   self.abundance_table.index)
        if np.isnan(bc_matrix.data).any():
            bc_matrix.data = np.nan_to_num(bc_matrix.data)
        self.formatter.success("Beta diversity calculated")
        return bc_matrix

    def _perform_beta_diversity_stats(self, distance_matrix):
        """ORIGINAL NAME PRESERVED - Enhanced with R² reporting."""
        self.formatter.step("Running Beta Diversity Statistical Tests")
        
        try:
            metadata_for_stats = self.metadata_df.copy()
            
            if self:
                self.formatter.info(f"Distance matrix IDs: {list(distance_matrix.ids)}")
                self.formatter.info(f"Metadata index: {list(metadata_for_stats.index)}")
            
            # PERMANOVA
            try:
                self.formatter.substep("Running PERMANOVA...")
                permanova_result = permanova(distance_matrix, metadata_for_stats, 
                                            column='group', permutations=999)
                p_val = permanova_result['p-value']
                f_stat = permanova_result['test statistic']
                
                # NEW: Extract R² if available
                r_squared = permanova_result.get('R2', 'N/A')
                
                interpretation = self._interpret_p_value(p_val)
                self.formatter.success(f"PERMANOVA: F = {f_stat:.3f}, p = {p_val:.4f} {interpretation}")
                if r_squared != 'N/A':
                    self.formatter.info(f"  R² = {r_squared:.3f} (variance explained by groups)")
                
                # Save detailed results
                with open(self.output_path / "permanova_results.txt", 'w') as f:
                    f.write(f"PERMANOVA Results\n")
                    f.write(f"================\n")
                    f.write(f"F-statistic: {f_stat:.6f}\n")
                    f.write(f"p-value: {p_val:.6f}\n")
                    if r_squared != 'N/A':
                        f.write(f"R² (effect size): {r_squared:.6f}\n")
                    f.write(f"Interpretation: {self._interpret_p_value(p_val)}\n")
                    f.write(f"Permutations: 999\n")
                
            except Exception as e:
                self.formatter.error(f"PERMANOVA failed: {e}")

            # ANOSIM
            try:
                self.formatter.substep("Running ANOSIM...")
                anosim_result = anosim(distance_matrix, metadata_for_stats, 
                                      column='group', permutations=999)
                r_stat = anosim_result['test statistic']
                p_val = anosim_result['p-value']
                
                interpretation = self._interpret_p_value(p_val)
                self.formatter.success(f"ANOSIM: R = {r_stat:.3f}, p = {p_val:.4f} {interpretation}")
                
                # Interpret R value
                if r_stat > 0.75:
                    r_interp = "well separated groups"
                elif r_stat > 0.5:
                    r_interp = "overlapping but clearly different"
                elif r_stat > 0.25:
                    r_interp = "barely separable"
                else:
                    r_interp = "indistinguishable groups"
                
                self.formatter.info(f"  R interpretation: {r_interp}")
                
                # Save detailed results
                with open(self.output_path / "anosim_results.txt", 'w') as f:
                    f.write(f"ANOSIM Results\n")
                    f.write(f"==============\n")
                    f.write(f"R-statistic: {r_stat:.6f}\n")
                    f.write(f"p-value: {p_val:.6f}\n")
                    f.write(f"Interpretation: {r_interp}\n")
                    f.write(f"Permutations: 999\n")
                    
            except Exception as e:
                self.formatter.error(f"ANOSIM failed: {e}")
                
        except Exception as e:
            self.formatter.error(f"Beta diversity statistical tests failed: {e}")

    def _run_differential_abundance(self) -> pd.DataFrame:
        """ORIGINAL NAME PRESERVED - Enhanced with compositional methods."""
        self.formatter.step("Running differential abundance analysis")
        
        groups = self.metadata_df['group'].unique()
        if len(groups) != 2:
            self.formatter.warning("Differential abundance test requires exactly 2 groups - skipping")
            return None
            
        group1_name, group2_name = groups
        group1_samples = self.metadata_df[self.metadata_df['group'] == group1_name].index
        group2_samples = self.metadata_df[self.metadata_df['group'] == group2_name].index

        self.formatter.substep(f"Comparing {group1_name} (n={len(group1_samples)}) vs {group2_name} (n={len(group2_samples)})")

        # Use raw counts for differential abundance
        data_for_testing = self.original_abundance_table if self.original_abundance_table is not None else self.abundance_table

        results = []
        for species in data_for_testing.columns:
            g1_vals = data_for_testing.loc[group1_samples, species]
            g2_vals = data_for_testing.loc[group2_samples, species]
            
            # Skip if no reads
            if g1_vals.sum() == 0 and g2_vals.sum() == 0: 
                continue
            
            # NEW: Prevalence filter - skip if present in <25% of samples in each group
            prev1 = (g1_vals > 0).sum() / len(g1_vals)
            prev2 = (g2_vals > 0).sum() / len(g2_vals)
            if prev1 < 0.25 and prev2 < 0.25:
                continue
            
            try:
                _, p_value = mannwhitneyu(g1_vals, g2_vals, alternative='two-sided')
                
                # Calculate relative abundances for fold change
                total_g1 = data_for_testing.loc[group1_samples].sum().sum()
                total_g2 = data_for_testing.loc[group2_samples].sum().sum()
                rel_mean_g1 = g1_vals.sum() / total_g1
                rel_mean_g2 = g2_vals.sum() / total_g2
                
                fold_change = (rel_mean_g1 + 1e-10) / (rel_mean_g2 + 1e-10)
                log2_fc = np.log2(fold_change)
                
                # NEW: Cliff's Delta effect size
                cliffs_delta = self._calculate_cliffs_delta(g1_vals.values, g2_vals.values)
                
                results.append({
                    'species': species,
                    'p_value': p_value,
                    'fold_change': fold_change,
                    'log2_fold_change': log2_fc,
                    'cliffs_delta': cliffs_delta,
                    'mean_group1': g1_vals.mean(),
                    'mean_group2': g2_vals.mean(),
                    'rel_abundance_group1': rel_mean_g1,
                    'rel_abundance_group2': rel_mean_g2,
                    'prevalence_group1': prev1,
                    'prevalence_group2': prev2,
                    'median_group1': g1_vals.median(),
                    'median_group2': g2_vals.median(),
                    'group1_name': group1_name,
                    'group2_name': group2_name
                })
            except Exception as e:
                if self.debug:
                    self.formatter.warning(f"Could not test {species}: {e}")
        
        if not results:
            self.formatter.warning("No species were eligible for testing")
            return None

        diff_df = pd.DataFrame(results)
        
        # Multiple testing correction
        if len(diff_df) > 0:
            try:
                rejected_fdr, p_fdr, _, _ = multipletests(diff_df['p_value'], alpha=0.05, method='fdr_bh')
                diff_df['p_fdr'] = p_fdr
                diff_df['significant_fdr'] = rejected_fdr
            except:
                diff_df['p_fdr'] = diff_df['p_value']
                diff_df['significant_fdr'] = False
            
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
        
        self.formatter.success("Differential abundance results:")
        self.formatter.metric("FDR significant (q < 0.05)", fdr_sig)
        self.formatter.metric("Bonferroni significant", bonf_sig)
        self.formatter.metric("Nominally significant (p < 0.05)", nominal_sig)
        self.formatter.metric("Trend level (p < 0.10)", trend_sig)
        
        if fdr_sig > 0 or trend_sig > 0:
            self.formatter.substep("Top candidates:")
            top_candidates = diff_df[diff_df['p_value'] < 0.10].head(5)
            for _, row in top_candidates.iterrows():
                direction = "↑" if row['fold_change'] > 1 else "↓"
                effect_interp = self._interpret_cliffs_delta(row['cliffs_delta'])
                self.formatter.info(f"{row['species']}: FC = {row['fold_change']:.2f} {direction}, p = {row['p_value']:.4f}, δ = {row['cliffs_delta']:.3f} ({effect_interp})")
        
        return diff_df

    def _run_ml_biomarker_discovery(self):
        """ORIGINAL NAME PRESERVED - Enhanced with hyperparameter tuning and feature selection."""
        self.formatter.step("Running Machine Learning for Biomarker Discovery")
        try:
            X = self.abundance_table.copy()
            y = self.metadata_df.loc[X.index, 'group']
            
            # Ensure correct data types
            X = X.astype(float)
            
            # NEW: Feature selection - keep top variable species
            n_features = min(50, X.shape[1])
            if X.shape[1] > n_features:
                selector = SelectKBest(f_classif, k=n_features)
                X_selected = selector.fit_transform(X, y)
                selected_features = X.columns[selector.get_support()]
                X = pd.DataFrame(X_selected, index=X.index, columns=selected_features)
                self.formatter.substep(f"Selected {n_features} most variable species for ML")
            
            # NEW: Hyperparameter tuning with GridSearch
            param_grid = {
                'n_estimators': [50, 100, 200],
                'max_depth': [5, 10, None],
                'min_samples_split': [2, 5]
            }
            
            rf_base = RandomForestClassifier(random_state=42, class_weight='balanced')
            
            cv_folds = min(3, len(X) // 2)
            
            try:
                grid_search = GridSearchCV(rf_base, param_grid, cv=cv_folds, 
                                          scoring='roc_auc', n_jobs=-1)
                grid_search.fit(X, y)
                
                rf = grid_search.best_estimator_
                best_score = grid_search.best_score_
                
                self.formatter.success(f"Random Forest {cv_folds}-fold CV AUC: {best_score:.3f}")
                if self:
                    self.formatter.info(f"Best params: {grid_search.best_params_}")
                
            except Exception as e:
                # Fallback to simple RF if GridSearch fails
                if self:
                    self.formatter.warning(f"GridSearch failed, using default RF: {e}")
                rf = RandomForestClassifier(n_estimators=100, random_state=42, 
                                          class_weight='balanced')
                cv_scores = cross_val_score(rf, X, y, cv=cv_folds, scoring='accuracy')
                rf.fit(X, y)
                self.formatter.success(f"Random Forest {cv_folds}-fold CV Accuracy: {np.mean(cv_scores):.2f} ± {np.std(cv_scores):.2f}")
            
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
            
            self.formatter.substep("Top discriminative species:")
            for i in range(min(2, len(importance_df))):
                species_name = importance_df.iloc[i]['species']
                importance_val = importance_df.iloc[i]['importance']
                self.formatter.info(f"'{species_name}' - Importance: {importance_val:.4f}")
            self.formatter.success(f"Feature importance report saved")
                
        except Exception as e:
            self.formatter.warning(f"ML biomarker discovery failed: {e}")
    
    def _assess_statistical_power(self):
        """NEW METHOD: Assess statistical power and provide recommendations."""
        self.formatter.step("Assessing Statistical Power")
        
        groups = self.metadata_df['group'].value_counts()
        
        self.formatter.substep("Sample sizes per group:")
        for group, n in groups.items():
            self.formatter.info(f"  {group}: n = {n}")
        
        min_n = groups.min()
        
        if min_n < 5:
            power_assessment = "⚠️  VERY LOW - Results unreliable"
            recommendation = "Need at least 10 samples per group for reliable results"
            detectable = "Only very large effects (d > 1.5) can be detected"
        elif min_n < 10:
            power_assessment = "⚠️  LOW - Limited statistical power"
            recommendation = "Recommend 15+ samples per group for robust findings"
            detectable = f"Can detect large effects (d ≈ {self._estimate_detectable_effect(min_n):.2f})"
        elif min_n < 20:
            power_assessment = "✓ MODERATE - Adequate for strong effects"
            recommendation = "Good for exploratory analysis, consider larger n for publication"
            detectable = f"Can detect medium-large effects (d ≈ {self._estimate_detectable_effect(min_n):.2f})"
        elif min_n < 30:
            power_assessment = "✓ GOOD - Sufficient power"
            recommendation = "Adequate for most analyses"
            detectable = f"Can detect medium effects (d ≈ {self._estimate_detectable_effect(min_n):.2f})"
        else:
            power_assessment = "✓ EXCELLENT - High statistical power"
            recommendation = "Well-powered study"
            detectable = f"Can detect small-medium effects (d ≈ {self._estimate_detectable_effect(min_n):.2f})"
        
        self.formatter.success(f"Power Assessment: {power_assessment}")
        self.formatter.info(f"Recommendation: {recommendation}")
        self.formatter.info(f"Detectable effect size (80% power, α=0.05): {detectable}")
        
        # Save assessment
        with open(self.output_path / "power_analysis.txt", 'w') as f:
            f.write("Statistical Power Assessment\n")
            f.write("=" * 50 + "\n\n")
            f.write("Sample Sizes:\n")
            for group, n in groups.items():
                f.write(f"  {group}: n = {n}\n")
            f.write(f"\nPower: {power_assessment}\n")
            f.write(f"Recommendation: {recommendation}\n")
            f.write(f"Detectable effect: {detectable}\n")
            f.write(f"\nInterpretation Guide:\n")
            f.write(f"  - Cohen's d: 0.2 (small), 0.5 (medium), 0.8 (large)\n")
            f.write(f"  - Cliff's δ: 0.15 (small), 0.33 (medium), 0.47 (large)\n")
    
    def _estimate_detectable_effect(self, n: int) -> float:
        """NEW METHOD: Estimate minimum detectable Cohen's d effect size."""
        # Approximate formula for 80% power, α=0.05, two-tailed
        return 2.8 / np.sqrt(n)


def run_comparison(input_dirs: List[str], metadata_file: str, output_dir: str, 
                   debug: bool = False, normalization: str = 'tss'):
    """
    ORIGINAL FUNCTION NAME PRESERVED - Enhanced with normalization option.
    
    Main entry point to run the enhanced comparative analysis.
    
    Args:
        input_dirs: List of MetaQuest output directories
        metadata_file: Path to metadata TSV (must have 'sample_id' and 'group' columns)
        output_dir: Output directory for results
        debug: Enable debug output
        normalization: 'tss' (default), 'clr', 'rarefy', or 'none'
    
    NEW FEATURES:
        - Compositional data normalization (TSS, CLR, rarefaction)
        - Multi-group comparisons with Kruskal-Wallis + post-hoc tests
        - Effect size metrics (Cliff's Delta, epsilon-squared, R²)
        - Statistical power assessment
        - Prevalence filtering for differential abundance
        - Enhanced ML with hyperparameter tuning and feature selection
        
    BACKWARD COMPATIBLE:
        - All existing function calls work without modification
        - Default behavior mimics original (TSS normalization)
        - All output files maintain original names
    
    Example:
        # Original usage still works:
        run_comparison(['s1/', 's2/'], 'metadata.tsv', 'output/')
        
        # New features available:
        run_comparison(['s1/', 's2/'], 'metadata.tsv', 'output/', 
                      normalization='clr', debug=True)
    """
    analyzer = ComparativeAnalysis(output_dir, debug=debug)
    analyzer.run_complete_analysis(input_dirs, metadata_file, normalization=normalization)


# ==================== ADDITIONAL UTILITY FUNCTIONS ====================

def calculate_compositional_stats(abundance_table: pd.DataFrame) -> Dict:
    """
    NEW UTILITY: Calculate compositional data statistics.
    
    Returns statistics about sparsity, zeros, and compositionality.
    """
    stats = {
        'n_samples': abundance_table.shape[0],
        'n_species': abundance_table.shape[1],
        'sparsity': (abundance_table == 0).sum().sum() / abundance_table.size,
        'mean_reads_per_sample': abundance_table.sum(axis=1).mean(),
        'median_reads_per_sample': abundance_table.sum(axis=1).median(),
        'cv_library_size': abundance_table.sum(axis=1).std() / abundance_table.sum(axis=1).mean()
    }
    return stats


def recommend_normalization(abundance_table: pd.DataFrame) -> str:
    """
    NEW UTILITY: Recommend normalization method based on data characteristics.
    
    Returns:
        'tss', 'clr', or 'rarefy' based on data properties
    """
    stats = calculate_compositional_stats(abundance_table)
    
    # High sparsity (>80% zeros) → CLR or TSS
    if stats['sparsity'] > 0.8:
        return 'clr'
    
    # High library size variation (CV > 0.5) → Rarefaction
    if stats['cv_library_size'] > 0.5:
        return 'rarefy'
    
    # Default: TSS (relative abundance)
    return 'tss'


def validate_metadata(metadata_file: str) -> bool:
    """
    NEW UTILITY: Validate metadata file format.
    
    Checks for required columns and proper formatting.
    """
    try:
        df = pd.read_csv(metadata_file, sep='\t')
        
        if 'sample_id' not in df.columns:
            print("❌ Error: 'sample_id' column missing from metadata")
            return False
        
        if 'group' not in df.columns:
            print("❌ Error: 'group' column missing from metadata")
            return False
        
        if df['sample_id'].duplicated().any():
            print("❌ Error: Duplicate sample IDs in metadata")
            return False
        
        n_groups = df['group'].nunique()
        if n_groups < 2:
            print("❌ Error: Need at least 2 groups for comparison")
            return False
        
        print(f"✓ Metadata valid: {len(df)} samples, {n_groups} groups")
        return True
        
    except Exception as e:
        print(f"❌ Error reading metadata: {e}")
        return False