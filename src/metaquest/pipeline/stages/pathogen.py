"""
Pathogen Detection Stage — DIAMOND database scan + ML prediction.
"""

from __future__ import annotations

import logging
from pathlib import Path

from metaquest.exceptions import PathogenDetectionError
from metaquest.pipeline.context import PipelineContext, PathogenResult

logger = logging.getLogger(__name__)


def run_pathogen_stage(ctx: PipelineContext) -> PipelineContext:
    """Run pathogen database scan and ML-based prediction."""
    from metaquest.core.pathogen_analysis import run_pathogen_scan
    from metaquest.ml.pathogen_predictor import run_ml_pathogen_prediction
    from metaquest.io.utils import should_use_ml_prediction
    from metaquest.reporting.risk_scoring import calculate_all_risks

    if ctx.annotation is None or ctx.classification is None:
        raise PathogenDetectionError("Classification and annotation must complete first")

    config = ctx.config.pathogen_detection
    ml_config = ctx.config.ml

    protein_file = ctx.annotation.protein_file
    bracken_file = ctx.classification.bracken_file
    prokka_dir = ctx.annotation.prokka_dir
    output_dir = ctx.output_dir

    pathogen_hits_file = None
    ml_predictions_file = None
    ml_summary = None

    # Database scan
    if protein_file and protein_file.exists():
        logger.info(
            "Scanning pathogen database (identity>=%.0f%%, evalue<=%s)",
            config.min_identity, config.evalue_threshold,
        )
        try:
            pathogen_hits_file = run_pathogen_scan(
                protein_file, output_dir,
                bracken_results=bracken_file,
                taxonomy_results=None,
            )
        except Exception as e:
            logger.warning("Pathogen database scan failed: %s", e)

    # ML prediction
    if prokka_dir and ml_config.enabled:
        from metaquest.io.output_formatter import get_formatter
        fmt = get_formatter()

        ml_suitable = should_use_ml_prediction(prokka_dir, fmt)
        if ml_suitable:
            logger.info("Running ML pathogen prediction")
            try:
                ml_results, ml_summary = run_ml_pathogen_prediction(prokka_dir, output_dir)
                ml_predictions_file = output_dir / "ml_pathogen_predictions.json"
                if not ml_predictions_file.exists():
                    ml_predictions_file = None
            except Exception as e:
                logger.warning("ML prediction failed: %s", e)

    # Comprehensive risk assessment
    logger.info("Calculating integrated risk scores")
    functional_file = ctx.annotation.functional_annotations

    risk_data = calculate_all_risks(
        bracken_file=bracken_file,
        functional_file=functional_file,
        pathogen_hits_file=pathogen_hits_file,
        ml_predictions_file=ml_predictions_file,
    )

    ctx.pathogen = PathogenResult(
        pathogen_hits_file=pathogen_hits_file,
        ml_predictions_file=ml_predictions_file,
        ml_summary=ml_summary,
        risk_data=risk_data,
    )

    risk_level = risk_data.get("integrated", {}).get("risk_level", "Unknown")
    risk_score = risk_data.get("integrated", {}).get("final_score", 0)
    logger.info("Risk assessment: %s (score=%.0f/100)", risk_level, risk_score)
    return ctx
