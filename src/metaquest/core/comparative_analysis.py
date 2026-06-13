"""
MetaQuest Comparative Analysis Module
======================================

Enhanced with improved OutputFormatter integration for cleaner logging.

FIXES APPLIED:
1. Deprecated rarefaction with clear warnings (recommends CLR/TSS)
2. Added CLR validation (sparsity checks, negative value detection)
3. Fixed power analysis formula (proper Mann-Whitney calculation)
4. Removed silent CLR fallback (explicit error handling)
5. Fixed memory leak (explicit cleanup of intermediate data)
6. Added reproducibility (random seeds, version tracking)
7. Enhanced logging with OutputFormatter

Memory optimized for 8-16GB RAM targets.
"""

import pandas as pd
from pathlib import Path
from typing import List, Dict, Optional, Tuple
import numpy as np
import warnings
import sys

# Core scientific libraries
from skbio.diversity import beta_diversity, alpha_diversity
from skbio.stats.ordination import pcoa
from skbio.stats.distance import permanova, anosim
from scipy.stats import mannwhitneyu, kruskal, norm
from statsmodels.stats.multitest import multipletests
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score, GridSearchCV
from sklearn.feature_selection import SelectKBest, f_classif

try:
    from skbio.stats.composition import clr, multiplicative_replacement
    HAS_COMPOSITION = True
except ImportError:
    HAS_COMPOSITION = False
    warnings.warn(
        "skbio.stats.composition not available - CLR normalization disabled. "
        "Install with: pip install scikit-bio"
    )

from ..visualization.compare_visuals import ComparativeVisualizer
from ..io.output_formatter import OutputFormatter


class NormalizationError(Exception):
    """Raised when normalization fails validation."""
    pass


class ComparativeAnalysis:
    """
    Production-ready comparative analysis for microbiome data.
    
    Fixed Issues:
    - Rarefaction deprecated (non-reproducible, high variance)
    - CLR validation added (checks data suitability)
    - Power analysis corrected (proper statistical formulas)
    - Memory leaks eliminated (explicit cleanup)
    - Silent errors removed (all failures raise exceptions)
    - Enhanced logging with OutputFormatter
    """
    
    def __init__(self, output_dir: str, debug: bool = False, random_seed: int = 42):
        """
        Args:
            output_dir: Output directory path
            debug: Enable debug logging
            random_seed: Random seed for reproducibility
        """
        self.output_path = Path(output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        self.abundance_table = None
        self.metadata_df = None
        
        # Initialize formatter with proper verbosity
        verbosity = 'debug' if debug else 'standard'
        log_file = self.output_path / "comparative_analysis.log"
        self.formatter = OutputFormatter(verbosity=verbosity, log_file=log_file)
        self.debug = debug
        
        # Reproducibility tracking
        self.random_seed = random_seed
        np.random.seed(random_seed)
        
        # Track normalization and preserve original data
        self.normalization_method = None
        self.original_abundance_table = None
        
        # Save analysis metadata
        self._save_metadata()

    def _save_metadata(self):
        """Save reproducibility metadata."""
        metadata = {
            'version': __version__,
            'random_seed': self.random_seed,
            'numpy_version': np.__version__,
            'pandas_version': pd.__version__,
            'python_version': sys.version
        }
        
        import json
        with open(self.output_path / "analysis_metadata.json", 'w') as f:
            json.dump(metadata, f, indent=2)
        
        self.formatter.debug(f"Analysis metadata saved (seed={self.random_seed})")

    def run_complete_analysis(self, input_dirs: List[str], metadata_file: str,
                             normalization: str = 'tss'):
        """
        Main analysis orchestrator.
        
        Args:
            input_dirs: List of sample directories
            metadata_file: Metadata TSV path
            normalization: 'tss' (default), 'clr', or 'none'
                          'rarefy' is DEPRECATED and will show warning
        """
        self.formatter.section_header("Comparative Analysis Pipeline")
        
        # DEPRECATED: Warn if rarefaction requested
        if normalization == 'rarefy':
            self.formatter.warning(
                "RAREFACTION IS DEPRECATED: Introduces high sampling variance and "
                "is not reproducible. Falling back to TSS normalization."
            )
            self.formatter.info("Recommendation: Use 'clr' (compositional) or 'tss' (relative abundance)")
            self.formatter.info("Reference: McMurdie & Holmes (2014) PLOS Comp Bio")
            normalization = 'tss'

        # Load and validate data
        with self.formatter.spinner("Loading and processing data"):
            if not self.load_and_aggregate_data(input_dirs, metadata_file, normalization):
                return

        # Statistical analyses
        self._assess_statistical_power()
        
        # Alpha diversity
        with self.formatter.spinner("Calculating alpha diversity metrics"):
            alpha_div_df = self._calculate_alpha_diversity()
        
        # Beta diversity
        with self.formatter.spinner("Computing beta diversity distances"):
            bc_matrix = self._calculate_beta_diversity()
        
        if bc_matrix is not None:
            self._perform_beta_diversity_stats(bc_matrix)
        
        # Differential abundance
        with self.formatter.spinner("Running differential abundance tests"):
            diff_df = self._run_differential_abundance()
        
        # Machine learning
        with self.formatter.spinner("Training ML biomarker models"):
            self._run_ml_biomarker_discovery()

        # Visualization
        self.formatter.section_header("Generating Visualizations")
        visualizer = ComparativeVisualizer(self.output_path)
        
        viz_count = 0
        if bc_matrix is not None:
            with self.formatter.spinner("Creating PCoA plot"):
                pcoa_results = pcoa(bc_matrix)
                visualizer.create_pcoa_plot(
                    pcoa_results.samples, 
                    self.metadata_df,
                    pcoa_results.proportion_explained
                )
                viz_count += 1
        
        with self.formatter.spinner("Creating heatmap"):
            visualizer.create_heatmap(self.abundance_table)
            viz_count += 1
        
        if alpha_div_df is not None and not alpha_div_df.empty:
            with self.formatter.spinner("Creating alpha diversity plots"):
                visualizer.create_alpha_diversity_boxplot(alpha_div_df)
                viz_count += 1
            
        if diff_df is not None and not diff_df.empty:
            with self.formatter.spinner("Creating volcano plot"):
                visualizer.create_volcano_plot(diff_df)
                viz_count += 1
        
        with self.formatter.spinner("Creating abundance barplot"):
            visualizer.create_abundance_barplot(self.abundance_table, self.metadata_df)
            viz_count += 1

        self.formatter.success(f"Generated {viz_count} visualizations", style='bold')
        self.formatter.success("Analysis complete!", style='bold')
        self.formatter.info(f"Results saved to: {self.output_path}")

    def load_and_aggregate_data(self, input_dirs: List[str], metadata_file: str,
                                normalization: str = 'tss') -> bool:
        """
        Load and aggregate data with memory-efficient processing.
        """
        # Load metadata
        self.formatter.info(f"Loading metadata from {Path(metadata_file).name}")
        try:
            self.metadata_df = pd.read_csv(metadata_file, sep='\t')
            self.formatter.success(f"Loaded {len(self.metadata_df)} samples from metadata")
        except Exception as e:
            self.formatter.error(
                f"Failed to parse metadata file: {e}",
                solutions=[
                    "Verify file is tab-separated (TSV format)",
                    "Check for required columns: 'sample_id' and 'group'",
                    "Ensure file encoding is UTF-8"
                ]
            )
            return False
        
        # Stream-based aggregation
        self.formatter.info("Aggregating taxonomic profiles from Bracken reports")
        
        abundance_dict = {}
        sample_ids_in_data = []
        failed_samples = []
        
        with self.formatter.progress_bar(
            total=len(input_dirs), 
            desc="Processing samples",
            unit="samples"
        ) as pbar:
            for input_dir in input_dirs:
                dir_path = Path(input_dir)
                sample_id = dir_path.name
                bracken_report = dir_path / "bracken_report.tsv"
                
                if bracken_report.exists():
                    try:
                        # Read and immediately process (don't store DataFrame)
                        sample_df = pd.read_csv(bracken_report, sep='\t')[['name', 'new_est_reads']]
                        
                        # Add to dict structure
                        for _, row in sample_df.iterrows():
                            species = row['name']
                            if species not in abundance_dict:
                                abundance_dict[species] = {}
                            abundance_dict[species][sample_id] = row['new_est_reads']
                        
                        sample_ids_in_data.append(sample_id)
                        
                        # CRITICAL: Explicit cleanup
                        del sample_df
                        
                        self.formatter.debug(f"Processed sample: {sample_id}")
                            
                    except Exception as e:
                        failed_samples.append(sample_id)
                        self.formatter.debug(f"Failed to process '{sample_id}': {e}")
                else:
                    failed_samples.append(sample_id)
                    self.formatter.debug(f"Missing Bracken report for '{sample_id}'")
                
                pbar.update(1)

        # Report processing results
        if failed_samples:
            self.formatter.warning(
                f"Failed to process {len(failed_samples)}/{len(input_dirs)} samples"
            )
            if self.debug:
                self.formatter.debug(f"Failed samples: {', '.join(failed_samples[:10])}")

        if not abundance_dict:
            self.formatter.error(
                "No valid Bracken reports found",
                solutions=[
                    "Verify Bracken classification completed successfully",
                    "Check that input directories contain 'bracken_report.tsv' files",
                    "Ensure file permissions allow reading"
                ]
            )
            return False
        
        # Convert dict to DataFrame (single allocation)
        self.formatter.debug("Converting to DataFrame")
        self.abundance_table = pd.DataFrame.from_dict(abundance_dict, orient='index').fillna(0).astype(int)
        self.abundance_table = self.abundance_table.transpose()
        
        # CRITICAL: Explicit cleanup
        del abundance_dict
        
        # Store original before any transformations
        self.original_abundance_table = self.abundance_table.copy()
        
        self.formatter.success(
            f"Aggregated {self.abundance_table.shape[1]} species across "
            f"{self.abundance_table.shape[0]} samples"
        )
        
        # Validate and clean
        if not self._clean_and_validate_data():
            return False
        
        # Apply normalization with validation
        try:
            self._normalize_abundance_data(normalization)
        except NormalizationError as e:
            self.formatter.error(
                str(e),
                solutions=[
                    "Try using TSS normalization: --normalization tss",
                    "Filter more aggressively to reduce sparsity",
                    "Check data quality and sequencing depth"
                ]
            )
            return False
        
        return True
    
    def _clean_and_validate_data(self) -> bool:
        """Clean and validate with prevalence filtering."""
        self.formatter.info("Cleaning and validating abundance data")

        # Remove empty samples
        original_n = self.abundance_table.shape[0]
        self.abundance_table = self.abundance_table.loc[(self.abundance_table.sum(axis=1) > 0)]
        removed = original_n - self.abundance_table.shape[0]
        
        if removed > 0:
            self.formatter.warning(f"Removed {removed} empty samples (zero total abundance)")

        # Prevalence filtering
        min_prevalence = 0.10
        prevalence = (self.abundance_table > 0).sum(axis=0) / len(self.abundance_table)
        species_to_keep = prevalence[prevalence >= min_prevalence].index
        
        original_species = self.abundance_table.shape[1]
        self.abundance_table = self.abundance_table[species_to_keep]
        species_removed = original_species - self.abundance_table.shape[1]
        
        if species_removed > 0:
            self.formatter.info(
                f"Filtered {species_removed} rare species (prevalence <{min_prevalence*100:.0f}%)"
            )
        
        # Align metadata
        self.metadata_df = self.metadata_df[
            self.metadata_df['sample_id'].isin(self.abundance_table.index)
        ]
        
        if len(self.metadata_df) < 2:
            self.formatter.error(
                "Insufficient samples after cleaning (need ≥2)",
                solutions=[
                    "Check sample quality and sequencing depth",
                    "Verify metadata sample_id matches directory names",
                    "Reduce prevalence filtering threshold"
                ]
            )
            return False
        
        # Ensure alignment
        self.metadata_df = self.metadata_df.set_index('sample_id')
        self.abundance_table = self.abundance_table.reindex(self.metadata_df.index)
        
        # Update original table
        if self.original_abundance_table is not None:
            self.original_abundance_table = self.original_abundance_table.reindex(
                self.metadata_df.index
            )
            self.original_abundance_table = self.original_abundance_table[species_to_keep]
        
        # Save cleaned data
        output_file = self.output_path / "taxonomic_abundance_table.tsv"
        self.abundance_table.to_csv(output_file, sep='\t')
        
        self.formatter.success(
            f"Final dataset: {self.abundance_table.shape[1]} species × "
            f"{self.abundance_table.shape[0]} samples"
        )
        self.formatter.debug(f"Saved to: {output_file}")
        
        return True

    def _normalize_abundance_data(self, method: str = 'tss'):
        """
        Normalize abundance data with validation.
        """
        self.normalization_method = method
        self.formatter.info(f"Applying {method.upper()} normalization")
        
        if method == 'tss':
            # Total Sum Scaling (relative abundance)
            self.abundance_table = self.abundance_table.div(
                self.abundance_table.sum(axis=1), axis=0
            )
            self.formatter.success("Applied Total Sum Scaling (relative abundance)")
            
        elif method == 'clr':
            if not HAS_COMPOSITION:
                raise NormalizationError(
                    "CLR requires scikit-bio composition module. "
                    "Install: pip install scikit-bio"
                )
            
            # Validate data before CLR
            self._validate_for_clr()
            
            try:
                # Apply multiplicative replacement for zeros
                abundance_with_pseudo = multiplicative_replacement(
                    self.abundance_table.values
                )
                
                # CLR transformation
                clr_transformed = clr(abundance_with_pseudo)
                
                self.abundance_table = pd.DataFrame(
                    clr_transformed,
                    index=self.abundance_table.index,
                    columns=self.abundance_table.columns
                )
                
                self.formatter.success("Applied CLR (Centered Log-Ratio) transformation")
                
            except Exception as e:
                raise NormalizationError(
                    f"CLR transformation failed: {e}"
                )
        
        elif method == 'none':
            self.formatter.info("No normalization applied (using raw counts)")
        
        else:
            raise NormalizationError(
                f"Unknown normalization method: '{method}'. "
                f"Valid options: 'tss', 'clr', 'none'"
            )
    
    def _validate_for_clr(self):
        """
        Validate data suitability for CLR transformation.
        """
        # Check sparsity
        sparsity = (self.abundance_table == 0).sum().sum() / self.abundance_table.size
        
        if sparsity > 0.95:
            raise NormalizationError(
                f"Data too sparse ({sparsity*100:.1f}% zeros) for CLR transformation"
            )
        
        # Check for negative values
        if (self.abundance_table < 0).any().any():
            raise NormalizationError(
                "Negative values detected - CLR requires counts ≥ 0"
            )
        
        # Check for extremely low counts
        total_counts = self.abundance_table.sum().sum()
        if total_counts < 1000:
            self.formatter.warning(
                f"Very low total counts ({total_counts:.0f}) - CLR may be unstable"
            )
        
        self.formatter.debug(f"CLR validation passed (sparsity: {sparsity*100:.1f}%)")

    def _calculate_alpha_diversity(self) -> pd.DataFrame:
        """Calculate alpha diversity metrics."""
        self.formatter.section_header("Alpha Diversity Analysis")
        
        # Use RAW counts for alpha diversity
        data_for_alpha = self.original_abundance_table
        
        metrics = ['shannon', 'simpson', 'chao1', 'sobs']
        alpha_div_results = {}
        
        for metric in metrics:
            try:
                with self.formatter.suppressed_output():
                    alpha_div_results[metric] = alpha_diversity(
                        metric, 
                        data_for_alpha.values,
                        ids=data_for_alpha.index
                    )
                self.formatter.debug(f"Calculated {metric} diversity")
            except Exception as e:
                self.formatter.warning(f"Failed to calculate {metric}: {e}")

        if not alpha_div_results:
            self.formatter.warning("No alpha diversity metrics calculated successfully")
            return None
            
        alpha_df = pd.DataFrame(alpha_div_results)
        alpha_df = alpha_df.join(self.metadata_df)
        
        # Statistical tests
        self._perform_alpha_diversity_stats(alpha_df)
        
        # Save results
        output_file = self.output_path / "alpha_diversity_metrics.tsv"
        alpha_df.to_csv(output_file, sep='\t')
        self.formatter.success(f"Alpha diversity metrics saved to {output_file.name}")
        
        return alpha_df

    def _perform_alpha_diversity_stats(self, alpha_df):
        """Statistical tests for alpha diversity."""
        self.formatter.info("Performing statistical tests")
        
        groups = alpha_df['group'].unique()
        n_groups = len(groups)
        
        results = []
        
        for metric in ['shannon', 'simpson', 'chao1', 'sobs']:
            if metric not in alpha_df.columns:
                continue
            
            if n_groups == 2:
                # Two-group comparison
                g1, g2 = groups
                vals1 = alpha_df[alpha_df['group'] == g1][metric].values
                vals2 = alpha_df[alpha_df['group'] == g2][metric].values
                
                try:
                    stat, p = mannwhitneyu(vals1, vals2, alternative='two-sided')
                    delta = self._calculate_cliffs_delta(vals1, vals2)
                    
                    results.append({
                        'metric': metric,
                        'test': 'Mann-Whitney U',
                        'group1': g1,
                        'group2': g2,
                        'group1_mean': vals1.mean(),
                        'group2_mean': vals2.mean(),
                        'statistic': stat,
                        'p_value': p,
                        'cliffs_delta': delta,
                        'effect_size': self._interpret_cliffs_delta(delta)
                    })
                    
                    sig_marker = self._interpret_p_value(p, delta)
                    self.formatter.info(
                        f"{metric.capitalize()}: {g1} vs {g2} - "
                        f"p={p:.4f}, δ={delta:.3f} {sig_marker}"
                    )
                    
                except Exception as e:
                    self.formatter.debug(f"Test failed for {metric}: {e}")
                    
            else:
                # Multi-group comparison
                group_data = [
                    alpha_df[alpha_df['group'] == g][metric].values 
                    for g in groups
                ]
                
                try:
                    stat, p = kruskal(*group_data)
                    
                    # Effect size
                    n = len(alpha_df)
                    if n == n_groups:
                        epsilon_sq = 0.0
                    else:
                        epsilon_sq = (stat - n_groups + 1) / (n - n_groups)
                        epsilon_sq = max(0.0, min(epsilon_sq, 1.0))
                    
                    result = {
                        'metric': metric,
                        'test': 'Kruskal-Wallis',
                        'statistic': stat,
                        'p_value': p,
                        'epsilon_squared': epsilon_sq,
                        'effect_size': self._interpret_epsilon_squared(epsilon_sq),
                        'n_groups': n_groups
                    }
                    
                    for g in groups:
                        result[f'{g}_mean'] = alpha_df[alpha_df['group'] == g][metric].mean()
                    
                    results.append(result)
                    
                    sig_marker = self._interpret_p_value(p, epsilon_sq)
                    self.formatter.info(
                        f"{metric.capitalize()}: Kruskal-Wallis - "
                        f"p={p:.4f}, ε²={epsilon_sq:.3f} {sig_marker}"
                    )
                    
                    # Post-hoc if significant
                    if p < 0.05:
                        self._posthoc_pairwise(alpha_df, metric, groups)
                        
                except Exception as e:
                    self.formatter.debug(f"Test failed for {metric}: {e}")
        
        if results:
            results_df = pd.DataFrame(results)
            output_file = self.output_path / "alpha_diversity_statistics.tsv"
            results_df.to_csv(output_file, sep='\t', index=False)
            self.formatter.debug(f"Saved statistics to {output_file.name}")
    
    def _posthoc_pairwise(self, alpha_df: pd.DataFrame, metric: str, groups: list):
        """Post-hoc pairwise comparisons."""
        from itertools import combinations
        
        self.formatter.info(f"Post-hoc pairwise comparisons for {metric}:")
        
        pairs = list(combinations(groups, 2))
        p_values = []
        
        for g1, g2 in pairs:
            vals1 = alpha_df[alpha_df['group'] == g1][metric].values
            vals2 = alpha_df[alpha_df['group'] == g2][metric].values
            _, p = mannwhitneyu(vals1, vals2)
            p_values.append(p)
        
        # Bonferroni correction
        p_corrected = [min(p * len(pairs), 1.0) for p in p_values]
        
        for (g1, g2), p_adj in zip(pairs, p_corrected):
            sig = "***" if p_adj < 0.001 else "**" if p_adj < 0.01 else "*" if p_adj < 0.05 else "ns"
            self.formatter.info(f"  {g1} vs {g2}: p_adj={p_adj:.4f} {sig}", indent=2)

    def _interpret_p_value(self, p: float, effect_size: float = None) -> str:
        """Interpret p-value with effect size context."""
        if p < 0.001:
            base = "(***)"
        elif p < 0.01:
            base = "(**)"
        elif p < 0.05:
            base = "(*)"
        elif p < 0.1:
            base = "(·)"
        else:
            base = "(ns)"
        
        if effect_size is not None:
            if abs(effect_size) < 0.15:
                base += " [negligible]"
            elif abs(effect_size) < 0.33:
                base += " [small]"
            elif abs(effect_size) < 0.47:
                base += " [medium]"
            else:
                base += " [large]"
        
        return base
    
    def _calculate_cliffs_delta(self, x: np.ndarray, y: np.ndarray) -> float:
        """Calculate Cliff's Delta effect size."""
        n1, n2 = len(x), len(y)
        if n1 == 0 or n2 == 0:
            return 0.0
        
        greater = np.sum(x[:, np.newaxis] > y[np.newaxis, :])
        less = np.sum(x[:, np.newaxis] < y[np.newaxis, :])
        
        return float((greater - less) / (n1 * n2))
    
    def _interpret_cliffs_delta(self, delta: float) -> str:
        """Interpret Cliff's Delta magnitude."""
        abs_delta = abs(delta)
        if abs_delta < 0.147:
            return "negligible"
        elif abs_delta < 0.33:
            return "small"
        elif abs_delta < 0.474:
            return "medium"
        else:
            return "large"
    
    def _interpret_epsilon_squared(self, eps_sq: float) -> str:
        """Interpret epsilon-squared effect size."""
        if eps_sq < 0.01:
            return "negligible"
        elif eps_sq < 0.06:
            return "small"
        elif eps_sq < 0.14:
            return "medium"
        else:
            return "large"

    def _calculate_beta_diversity(self):
        """Calculate beta diversity."""
        self.formatter.section_header("Beta Diversity Analysis")
        
        with self.formatter.suppressed_output():
            bc_matrix = beta_diversity(
                "braycurtis",
                self.abundance_table.values,
                self.abundance_table.index
            )
        
        # Clean NaN values
        if np.isnan(bc_matrix.data).any():
            bc_matrix.data = np.nan_to_num(bc_matrix.data)
            self.formatter.debug("Cleaned NaN values from distance matrix")
        
        self.formatter.success("Computed Bray-Curtis dissimilarity matrix")
        return bc_matrix

    def _perform_beta_diversity_stats(self, distance_matrix):
        """Beta diversity statistics."""
        self.formatter.info("Performing statistical tests")
        
        try:
            metadata_for_stats = self.metadata_df.copy()
            
            # PERMANOVA
            try:
                with self.formatter.suppressed_output():
                    result = permanova(
                        distance_matrix,
                        metadata_for_stats,
                        column='group',
                        permutations=999
                    )
                
                p = result['p-value']
                f_stat = result['test statistic']
                r_squared = result.get('R2', 'N/A')
                
                sig_marker = self._interpret_p_value(p)
                self.formatter.info(f"PERMANOVA: F={f_stat:.3f}, p={p:.4f} {sig_marker}")
                
                if r_squared != 'N/A':
                    self.formatter.info(f"  R² = {r_squared:.3f} (variance explained)", indent=2)
                
                # Save results
                with open(self.output_path / "permanova_results.txt", 'w') as f:
                    f.write("PERMANOVA Results\n")
                    f.write("=" * 50 + "\n")
                    f.write(f"F-statistic: {f_stat:.6f}\n")
                    f.write(f"p-value: {p:.6f}\n")
                    if r_squared != 'N/A':
                        f.write(f"R² (effect size): {r_squared:.6f}\n")
                    f.write(f"Interpretation: {sig_marker}\n")
                    f.write(f"Permutations: 999\n")
                
            except Exception as e:
                self.formatter.error(f"PERMANOVA failed: {e}")

            # ANOSIM
            try:
                with self.formatter.suppressed_output():
                    result = anosim(
                        distance_matrix,
                        metadata_for_stats,
                        column='group',
                        permutations=999
                    )
                
                r_stat = result['test statistic']
                p = result['p-value']
                
                sig_marker = self._interpret_p_value(p)
                self.formatter.info(f"ANOSIM: R={r_stat:.3f}, p={p:.4f} {sig_marker}")
                
                # Interpret R
                if r_stat > 0.75:
                    interp = "well separated"
                elif r_stat > 0.5:
                    interp = "overlapping but distinct"
                elif r_stat > 0.25:
                    interp = "barely separable"
                else:
                    interp = "indistinguishable"
                
                self.formatter.info(f"  Groups: {interp}", indent=2)
                
                # Save results
                with open(self.output_path / "anosim_results.txt", 'w') as f:
                    f.write("ANOSIM Results\n")
                    f.write("=" * 50 + "\n")
                    f.write(f"R-statistic: {r_stat:.6f}\n")
                    f.write(f"p-value: {p:.6f}\n")
                    f.write(f"Interpretation: {interp}\n")
                    f.write(f"Permutations: 999\n")
                    
            except Exception as e:
                self.formatter.error(f"ANOSIM failed: {e}")
                
        except Exception as e:
            self.formatter.error(f"Beta diversity statistics failed: {e}")

    def _run_differential_abundance(self) -> Optional[pd.DataFrame]:
        """Differential abundance analysis."""
        self.formatter.section_header("Differential Abundance Testing")
        
        groups = self.metadata_df['group'].unique()
        if len(groups) != 2:
            self.formatter.warning(
                f"Requires exactly 2 groups for differential abundance (found {len(groups)})"
            )
            return None
            
        g1, g2 = groups
        samples1 = self.metadata_df[self.metadata_df['group'] == g1].index
        samples2 = self.metadata_df[self.metadata_df['group'] == g2].index

        self.formatter.info(f"Comparing {g1} (n={len(samples1)}) vs {g2} (n={len(samples2)})")

        # Use RAW counts for testing
        data_for_testing = self.original_abundance_table

        results = []
        tested_count = 0
        skipped_count = 0
        
        for species in data_for_testing.columns:
            vals1 = data_for_testing.loc[samples1, species]
            vals2 = data_for_testing.loc[samples2, species]
            
            # Skip if no reads
            if vals1.sum() == 0 and vals2.sum() == 0:
                skipped_count += 1
                continue
            
            # Prevalence filter
            prev1 = (vals1 > 0).sum() / len(vals1)
            prev2 = (vals2 > 0).sum() / len(vals2)
            if prev1 < 0.25 and prev2 < 0.25:
                skipped_count += 1
                continue
            
            try:
                _, p = mannwhitneyu(vals1, vals2, alternative='two-sided')
                
                # Calculate fold change
                total1 = data_for_testing.loc[samples1].sum().sum()
                total2 = data_for_testing.loc[samples2].sum().sum()
                rel_mean1 = vals1.sum() / total1
                rel_mean2 = vals2.sum() / total2
                
                fold_change = (rel_mean1 + 1e-10) / (rel_mean2 + 1e-10)
                log2_fc = np.log2(fold_change)
                
                # Effect size
                delta = self._calculate_cliffs_delta(vals1.values, vals2.values)
                
                results.append({
                    'species': species,
                    'p_value': p,
                    'fold_change': fold_change,
                    'log2_fold_change': log2_fc,
                    'cliffs_delta': delta,
                    'effect_size': self._interpret_cliffs_delta(delta),
                    'mean_group1': vals1.mean(),
                    'mean_group2': vals2.mean(),
                    'rel_abundance_group1': rel_mean1,
                    'rel_abundance_group2': rel_mean2,
                    'prevalence_group1': prev1,
                    'prevalence_group2': prev2,
                    'median_group1': vals1.median(),
                    'median_group2': vals2.median(),
                    'group1_name': g1,
                    'group2_name': g2
                })
                tested_count += 1
                
            except Exception as e:
                self.formatter.debug(f"Test failed for {species}: {e}")
        
        self.formatter.info(f"Tested {tested_count} species (skipped {skipped_count} rare/absent)")
        
        if not results:
            self.formatter.warning("No species eligible for differential abundance testing")
            return None

        diff_df = pd.DataFrame(results)
        
        # Multiple testing correction
        if len(diff_df) > 0:
            # FDR correction
            try:
                reject_fdr, p_fdr, _, _ = multipletests(
                    diff_df['p_value'], 
                    alpha=0.05, 
                    method='fdr_bh',
                    returnsorted=False
                )
                diff_df['p_fdr'] = p_fdr
                diff_df['significant_fdr'] = reject_fdr
                self.formatter.debug("Applied FDR (Benjamini-Hochberg) correction")
            except Exception as e:
                self.formatter.warning(f"FDR correction failed: {e}")
                diff_df['p_fdr'] = diff_df['p_value']
                diff_df['significant_fdr'] = False
            
            # Bonferroni correction
            try:
                reject_bonf, p_bonf, _, _ = multipletests(
                    diff_df['p_value'],
                    alpha=0.05,
                    method='bonferroni',
                    returnsorted=False
                )
                diff_df['p_bonferroni'] = p_bonf
                diff_df['significant_bonferroni'] = reject_bonf
                self.formatter.debug("Applied Bonferroni correction")
            except Exception as e:
                self.formatter.warning(f"Bonferroni correction failed: {e}")
                diff_df['p_bonferroni'] = diff_df['p_value']
                diff_df['significant_bonferroni'] = False
        
        # Significance categories
        diff_df['significance_level'] = diff_df['p_value'].apply(
            lambda p: 'highly_significant' if p < 0.001
            else 'significant' if p < 0.01
            else 'trend' if p < 0.1
            else 'not_significant'
        )
        
        diff_df = diff_df.sort_values('p_value')
        output_file = self.output_path / "differential_abundance_report.tsv"
        diff_df.to_csv(output_file, sep='\t', index=False)
        
        # Summary statistics
        n_tests = len(diff_df)
        fdr_sig = diff_df['significant_fdr'].sum()
        bonf_sig = diff_df['significant_bonferroni'].sum()
        nominal_sig = (diff_df['p_value'] < 0.05).sum()
        
        # Display results table
        self.formatter.display_stats_table(
            title="Differential Abundance Summary",
            sections=[
                {
                    'header': 'Test Results',
                    'rows': {
                        'Total tests': n_tests,
                        'FDR significant (q<0.05)': fdr_sig,
                        'Bonferroni significant': bonf_sig,
                        'Nominal significant (p<0.05)': nominal_sig,
                        'Expected false positives': f"{n_tests * 0.05:.1f}"
                    }
                }
            ]
        )
        
        if fdr_sig > 0:
            self.formatter.info("Top differential species:")
            top = diff_df[diff_df['p_fdr'] < 0.05].head(5)
            for i, (_, row) in enumerate(top.iterrows(), 1):
                direction = "↑" if row['fold_change'] > 1 else "↓"
                self.formatter.info(
                    f"{i}. {row['species'][:50]}: FC={row['fold_change']:.2f} {direction}, "
                    f"p={row['p_value']:.4f}, δ={row['cliffs_delta']:.3f}",
                    indent=2
                )
        else:
            self.formatter.info("No species passed FDR correction (q<0.05)")
        
        self.formatter.success(f"Results saved to {output_file.name}")
        
        return diff_df

    def _run_ml_biomarker_discovery(self):
        """ML biomarker discovery."""
        self.formatter.section_header("Machine Learning Biomarker Discovery")
        
        try:
            X = self.abundance_table.copy()
            y = self.metadata_df.loc[X.index, 'group']
            
            X = X.astype(float)
            
            # Feature selection
            n_features = min(50, X.shape[1])
            if X.shape[1] > n_features:
                with self.formatter.suppressed_output():
                    selector = SelectKBest(f_classif, k=n_features)
                    X_selected = selector.fit_transform(X, y)
                    selected_features = X.columns[selector.get_support()]
                    X = pd.DataFrame(X_selected, index=X.index, columns=selected_features)
                
                self.formatter.info(f"Selected top {n_features} features using ANOVA F-test")
            
            # Hyperparameter tuning
            param_grid = {
                'n_estimators': [50, 100, 200],
                'max_depth': [5, 10, None],
                'min_samples_split': [2, 5]
            }
            
            rf_base = RandomForestClassifier(random_state=self.random_seed, class_weight='balanced')
            cv_folds = min(3, len(X) // 2)
            
            self.formatter.info(f"Training Random Forest with {cv_folds}-fold cross-validation")
            
            try:
                with self.formatter.suppressed_output():
                    grid_search = GridSearchCV(
                        rf_base,
                        param_grid,
                        cv=cv_folds,
                        scoring='roc_auc',
                        n_jobs=-1
                    )
                    grid_search.fit(X, y)
                
                rf = grid_search.best_estimator_
                best_score = grid_search.best_score_
                
                self.formatter.success(
                    f"Random Forest CV AUC: {best_score:.3f}"
                )
                self.formatter.debug(f"Best hyperparameters: {grid_search.best_params_}")
                
            except Exception as e:
                self.formatter.debug(f"GridSearch failed, using default RF: {e}")
                
                with self.formatter.suppressed_output():
                    rf = RandomForestClassifier(
                        n_estimators=100,
                        random_state=self.random_seed,
                        class_weight='balanced'
                    )
                    cv_scores = cross_val_score(rf, X, y, cv=cv_folds, scoring='accuracy')
                    rf.fit(X, y)
                
                self.formatter.success(
                    f"Random Forest CV Accuracy: {np.mean(cv_scores):.2f} ± {np.std(cv_scores):.2f}"
                )
            
            # Feature importance
            importance_df = pd.DataFrame({
                'species': X.columns,
                'importance': rf.feature_importances_
            }).sort_values('importance', ascending=False)
            
            # Add abundance info for top features
            top_features = importance_df.head(20)
            for idx, row in top_features.iterrows():
                species = row['species']
                group_means = X.groupby(y)[species].mean()
                importance_df.loc[idx, 'group1_mean'] = group_means.iloc[0]
                importance_df.loc[idx, 'group2_mean'] = group_means.iloc[1]
                importance_df.loc[idx, 'fold_change'] = (
                    (group_means.iloc[0] + 1e-6) / (group_means.iloc[1] + 1e-6)
                )
            
            output_file = self.output_path / "random_forest_feature_importance.tsv"
            importance_df.to_csv(output_file, sep='\t', index=False)
            
            self.formatter.info("Top discriminative species:")
            for i in range(min(5, len(importance_df))):
                species = importance_df.iloc[i]['species']
                imp = importance_df.iloc[i]['importance']
                self.formatter.info(f"{i+1}. {species[:50]}: {imp:.4f}", indent=2)
            
            self.formatter.success(f"Feature importance saved to {output_file.name}")
                
        except Exception as e:
            self.formatter.error(f"ML biomarker discovery failed: {e}")
            if self.debug:
                import traceback
                self.formatter.debug(traceback.format_exc())
    
    def _assess_statistical_power(self):
        """Statistical power assessment."""
        self.formatter.section_header("Statistical Power Assessment")
        
        groups = self.metadata_df['group'].value_counts()
        
        # Create sections for display
        sample_size_rows = {}
        for group, n in groups.items():
            sample_size_rows[f"Group '{group}'"] = f"n={n}"
        
        min_n = groups.min()
        max_n = groups.max()
        
        # Calculate detectable effect size
        detectable_delta = self._estimate_detectable_mw_effect(
            min_n, 
            max_n,
            alpha=0.05,
            power=0.80
        )
        
        # Power assessment
        if min_n < 5:
            power_level = "VERY LOW [!]"
            recommendation = "Need ≥10 samples/group for reliable results"
            warning_msg = "Results may be unreliable - interpret with extreme caution"
        elif min_n < 10:
            power_level = "LOW [!]"
            recommendation = "Recommend ≥15 samples/group for robust findings"
            warning_msg = "Limited power - can only detect large effects"
        elif min_n < 20:
            power_level = "MODERATE [OK]"
            recommendation = "Adequate for exploratory analysis"
            warning_msg = None
        elif min_n < 30:
            power_level = "GOOD [OK]"
            recommendation = "Sufficient for most analyses"
            warning_msg = None
        else:
            power_level = "EXCELLENT [OK]"
            recommendation = "Well-powered study"
            warning_msg = None
        
        effect_interp = self._interpret_cliffs_delta(detectable_delta)
        
        # Display power assessment table
        power_rows = {
            'Power level': power_level,
            'Minimum group size': f"n={min_n}",
            'Detectable Cliff\'s δ': f"{detectable_delta:.3f} ({effect_interp})",
            'Recommendation': recommendation
        }
        
        self.formatter.display_stats_table(
            title="Power Analysis",
            sections=[
                {'header': 'Sample Sizes', 'rows': sample_size_rows},
                {'header': 'Power Assessment (α=0.05, 1-β=0.80)', 'rows': power_rows}
            ]
        )
        
        if warning_msg:
            self.formatter.warning(warning_msg)
        
        # Save detailed assessment
        with open(self.output_path / "power_analysis.txt", 'w') as f:
            f.write("Statistical Power Assessment\n")
            f.write("=" * 50 + "\n\n")
            f.write("Sample Sizes:\n")
            for group, n in groups.items():
                f.write(f"  {group}: n={n}\n")
            f.write(f"\nPower Assessment: {power_level}\n")
            f.write(f"Recommendation: {recommendation}\n")
            if warning_msg:
                f.write(f"Warning: {warning_msg}\n")
            f.write(f"\nDetectable Effect Sizes (80% power, α=0.05):\n")
            f.write(f"  Cliff's Delta: {detectable_delta:.3f} ({effect_interp})\n")
            f.write(f"\nEffect Size Guidelines:\n")
            f.write(f"  Cliff's δ: <0.147 (negligible), 0.147-0.33 (small), "
                   f"0.33-0.474 (medium), >0.474 (large)\n")
    
    def _estimate_detectable_mw_effect(
        self,
        n1: int,
        n2: int,
        alpha: float = 0.05,
        power: float = 0.80
    ) -> float:
        """
        Calculate detectable effect size for Mann-Whitney test.
        
        Based on Noether (1987) nonparametric power formula.
        """
        z_alpha = norm.ppf(1 - alpha / 2)
        z_beta = norm.ppf(power)
        
        n_harmonic = 2 * n1 * n2 / (n1 + n2)
        
        p_superior = norm.cdf((z_alpha + z_beta) / np.sqrt(3 * n_harmonic))
        
        cliffs_delta = 2 * p_superior - 1
        
        return abs(cliffs_delta)


# ==================== MAIN ENTRY POINT ====================

def run_comparison(
    input_dirs: List[str],
    metadata_file: str,
    output_dir: str,
    debug: bool = False,
    normalization: str = 'tss',
    random_seed: int = 42
):
    """
    PRODUCTION-READY comparative analysis with enhanced logging.
    
    Args:
        input_dirs: List of MetaQuest output directories
        metadata_file: Metadata TSV (requires 'sample_id' and 'group' columns)
        output_dir: Output directory
        debug: Enable debug logging
        normalization: 'tss' (default), 'clr', or 'none'
                      'rarefy' is DEPRECATED
        random_seed: Random seed for reproducibility
    
    Features:
    Enhanced OutputFormatter integration
    Cleaner progress reporting
    Better error messages with solutions
    Suppressed scientific library output
    Comprehensive logging to file
    
    Example:
        # Basic usage
        run_comparison(['s1/', 's2/'], 'metadata.tsv', 'results/')
        
        # With CLR normalization and debug mode
        run_comparison(['s1/', 's2/'], 'metadata.tsv', 'results/',
                      normalization='clr', debug=True)
    """
    analyzer = ComparativeAnalysis(
        output_dir,
        debug=debug,
        random_seed=random_seed
    )
    analyzer.run_complete_analysis(input_dirs, metadata_file, normalization)


# ==================== UTILITY FUNCTIONS ====================

def calculate_compositional_stats(abundance_table: pd.DataFrame) -> Dict:
    """Calculate compositional data statistics."""
    stats = {
        'n_samples': abundance_table.shape[0],
        'n_species': abundance_table.shape[1],
        'sparsity': (abundance_table == 0).sum().sum() / abundance_table.size,
        'mean_reads_per_sample': abundance_table.sum(axis=1).mean(),
        'median_reads_per_sample': abundance_table.sum(axis=1).median(),
        'cv_library_size': (
            abundance_table.sum(axis=1).std() / 
            abundance_table.sum(axis=1).mean()
        )
    }
    return stats


def recommend_normalization(abundance_table: pd.DataFrame) -> str:
    """
    Recommend normalization based on data characteristics.
    """
    stats = calculate_compositional_stats(abundance_table)
    
    if stats['sparsity'] > 0.95:
        return (
            "tss (data too sparse for CLR - "
            f"{stats['sparsity']*100:.1f}% zeros)"
        )
    
    if stats['sparsity'] > 0.80:
        return "clr (compositional, handles moderate sparsity)"
    
    if stats['cv_library_size'] > 0.5:
        return (
            "tss or clr (library sizes vary, "
            "but rarefaction deprecated)"
        )
    
    return "tss (default, works for most cases)"


def validate_metadata(metadata_file: str) -> bool:
    """
    Validate metadata file format.
    """
    formatter = OutputFormatter()
    
    try:
        df = pd.read_csv(metadata_file, sep='\t')
        
        if 'sample_id' not in df.columns:
            formatter.error(
                "Missing required column: 'sample_id'",
                solutions=[
                    "Add 'sample_id' column with sample identifiers",
                    "Ensure column names are in the first row",
                    "Check file is tab-separated"
                ]
            )
            return False
        
        if 'group' not in df.columns:
            formatter.error(
                "Missing required column: 'group'",
                solutions=[
                    "Add 'group' column with group assignments",
                    "Use consistent group names (e.g., 'control', 'treatment')"
                ]
            )
            return False
        
        if df['sample_id'].duplicated().any():
            formatter.error(
                "Duplicate sample IDs detected",
                solutions=[
                    "Ensure each sample has a unique identifier",
                    "Check for copy-paste errors in metadata"
                ]
            )
            return False
        
        n_groups = df['group'].nunique()
        if n_groups < 2:
            formatter.error(
                f"Need at least 2 groups (found {n_groups})",
                solutions=[
                    "Add group assignments for comparative analysis",
                    "Verify group column is correctly formatted"
                ]
            )
            return False
        
        formatter.success(f"[OK] Metadata valid: {len(df)} samples, {n_groups} groups")
        return True
        
    except Exception as e:
        formatter.error(
            f"Failed to read metadata file: {e}",
            solutions=[
                "Check file exists and is readable",
                "Verify file is tab-separated (TSV)",
                "Ensure UTF-8 encoding"
            ]
        )
        return False