"""Stable reporting stage for descriptive research outputs."""

from __future__ import annotations

import logging

from metaquest.exceptions import ReportingError
from metaquest.pipeline.context import PipelineContext


logger = logging.getLogger(__name__)


def run_reporting_stage(ctx: PipelineContext) -> PipelineContext:
    """Generate descriptive reports from validated pipeline outputs."""
    from metaquest.reporting.stable_reporter import generate_stable_reports

    try:
        generate_stable_reports(ctx)
    except Exception as exc:
        raise ReportingError(f"Stable report generation failed: {exc}") from exc

    logger.info("Reporting complete")
    return ctx
