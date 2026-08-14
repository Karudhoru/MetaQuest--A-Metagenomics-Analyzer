"""Stable reporting stage for descriptive research outputs."""

from __future__ import annotations

import logging

from metaquest.exceptions import ReportingError
from metaquest.pipeline.context import PipelineContext


logger = logging.getLogger(__name__)


def run_reporting_stage(ctx: PipelineContext) -> PipelineContext:
    """Generate descriptive reports and non-risk visualizations."""
    from metaquest.reporting.stable_reporter import generate_stable_reports
    from metaquest.visualization.main_visualizer import (
        create_functional_visualizations,
        create_taxonomic_visualizations,
    )

    try:
        generate_stable_reports(ctx)
    except Exception as exc:
        raise ReportingError(f"Stable report generation failed: {exc}") from exc

    if ctx.classification:
        try:
            create_taxonomic_visualizations(
                ctx.output_dir,
                ctx.classification.bracken_file,
            )
        except Exception as exc:
            logger.warning("Taxonomic visualizations failed: %s", exc)

    if ctx.annotation and ctx.annotation.functional_annotations:
        try:
            create_functional_visualizations(
                ctx.output_dir,
                ctx.annotation.gene_prediction_dir,
                ctx.annotation.functional_annotations,
            )
        except Exception as exc:
            logger.warning("Functional visualizations failed: %s", exc)

    logger.info("Reporting complete")
    return ctx
