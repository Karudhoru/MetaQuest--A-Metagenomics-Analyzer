"""
Comparative Analysis Subpackage
================================
Split from the monolithic comparative_analysis.py into focused modules:
- normalization: TSS, CLR, validation
- diversity: Alpha/beta diversity, ordination
- statistics: Differential abundance, effect sizes, power analysis
- biomarkers: ML-based biomarker discovery (Random Forest)
- runner: Orchestrator that composes all modules
"""

from .runner import ComparativeAnalysis, run_comparison

__all__ = ["ComparativeAnalysis", "run_comparison"]
