"""
Comparative Analysis Runner
============================
Re-exports from the original module for backward compatibility.
A full decomposition into normalization/diversity/statistics/biomarkers
modules is planned as a follow-up refactor.
"""

# For now, import directly from the monolithic module.
# This will be decomposed further once the pipeline architecture stabilizes.
from metaquest.core.comparative_analysis import (
    ComparativeAnalysis,
    run_comparison,
    calculate_compositional_stats,
    recommend_normalization,
    validate_metadata,
)

__all__ = [
    "ComparativeAnalysis",
    "run_comparison",
    "calculate_compositional_stats",
    "recommend_normalization",
    "validate_metadata",
]
