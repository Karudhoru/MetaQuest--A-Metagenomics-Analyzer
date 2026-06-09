"""
Reporting Stage — generates text, JSON, and HTML reports + visualizations.
"""

from __future__ import annotations

import logging
from pathlib import Path

from metaquest.exceptions import ReportingError
from metaquest.pipeline.context import PipelineContext

logger = logging.getLogger(__name__)


def run_reporting_stage(ctx: PipelineContext) -> PipelineContext:
    """Generate all reports and visualizations."""
    from metaquest.reporting.main_reporter import MainReporter
    from metaquest.reporting.risk_scoring import RiskScorer
    from metaquest.visualization.main_visualizer import (
        create_taxonomic_visualizations,
        create_pathogen_visualizations,
        create_functional_visualizations,
    )
    from metaquest.visualization.dashboard import create_dashboard

    output_dir = ctx.output_dir
    report_config = ctx.config.reporting

    try:
        reporter = MainReporter(output_dir, view_mode="both")
        risk_scorer = RiskScorer()

        # Taxonomic report
        if ctx.classification:
            _generate_taxonomy_report(ctx, reporter, risk_scorer)
            create_taxonomic_visualizations(output_dir, ctx.classification.bracken_file)

        # Functional report
        if ctx.annotation and ctx.annotation.functional_annotations:
            _generate_functional_report(ctx, reporter, risk_scorer)
            create_functional_visualizations(
                output_dir,
                ctx.annotation.prokka_dir,
                ctx.annotation.functional_annotations,
            )

        # Pathogen report
        if ctx.pathogen and ctx.pathogen.risk_data:
            _generate_pathogen_report(ctx, reporter)
            create_pathogen_visualizations(output_dir, risk_data=ctx.pathogen.risk_data)

        # Comprehensive report
        _generate_comprehensive_report(ctx, reporter)

        # Dashboard
        create_dashboard(analysis_type="fastq", output_dir=output_dir)

    except Exception as e:
        raise ReportingError(f"Report generation failed: {e}") from e

    logger.info("Reporting complete")
    return ctx


def _generate_taxonomy_report(ctx: PipelineContext, reporter, risk_scorer) -> None:
    from metaquest.io.data_loaders import load_bracken_report

    bracken_df = load_bracken_report(ctx.classification.bracken_file)
    taxonomic_risk = risk_scorer.calculate_taxonomic_risk(bracken_df)

    risk_data = {
        "taxonomic": taxonomic_risk,
        "functional": {"score": 0, "details": {}},
        "ml": {"score": 0, "details": []},
        "integrated": {
            "final_score": taxonomic_risk["score"],
            "risk_level": risk_scorer.get_risk_level(taxonomic_risk["score"]),
            "tier_scores": {"taxonomic": taxonomic_risk["score"], "functional": 0, "ml": 0},
        },
    }

    report = reporter.generate_taxonomy_report(
        bracken_file=ctx.classification.bracken_file,
        risk_data=risk_data,
    )
    with open(ctx.output_dir / "01_taxonomic_report.txt", "w") as f:
        f.write(report)


def _generate_functional_report(ctx: PipelineContext, reporter, risk_scorer) -> None:
    from metaquest.io.data_loaders import load_annotation_file, load_pathogen_hits
    import pandas as pd

    functional_file = ctx.annotation.functional_annotations
    functional_df = load_annotation_file(functional_file)
    pathogen_df = pd.DataFrame()

    if ctx.pathogen and ctx.pathogen.pathogen_hits_file:
        try:
            pathogen_df = load_pathogen_hits(ctx.pathogen.pathogen_hits_file)
        except Exception:
            pass

    total_cds = len(functional_df["query_id"].unique()) if not functional_df.empty and "query_id" in functional_df.columns else None

    taxonomic_risk = risk_scorer.calculate_taxonomic_risk(
        load_bracken_report_safe(ctx.classification.bracken_file)
    ) if ctx.classification else {"score": 0}

    functional_risk = risk_scorer.calculate_functional_risk(functional_df, pathogen_df, total_cds=total_cds)

    risk_data = {
        "taxonomic": taxonomic_risk,
        "functional": functional_risk,
        "ml": {"score": 0, "details": []},
        "integrated": {
            "final_score": taxonomic_risk["score"] * 0.4 + functional_risk["score"] * 0.3,
            "risk_level": risk_scorer.get_risk_level(
                taxonomic_risk["score"] * 0.4 + functional_risk["score"] * 0.3
            ),
            "tier_scores": {
                "taxonomic": taxonomic_risk["score"],
                "functional": functional_risk["score"],
                "ml": 0,
            },
        },
    }

    sample_info_file = ctx.annotation.prokka_dir / "sample.txt"
    report = reporter.generate_functional_report(
        sample_info_file=sample_info_file,
        functional_annotations_file=functional_file,
        risk_data=risk_data,
    )
    with open(ctx.output_dir / "02_functional_report.txt", "w") as f:
        f.write(report)


def _generate_pathogen_report(ctx: PipelineContext, reporter) -> None:
    ml_predictions_file = ctx.pathogen.ml_predictions_file
    pathogen_hits_file = ctx.pathogen.pathogen_hits_file or ctx.output_dir / "pathogen_results.tsv"

    report = reporter.generate_pathogen_report(
        risk_data=ctx.pathogen.risk_data,
        pathogen_hits_file=pathogen_hits_file,
        ml_predictions_file=ml_predictions_file,
    )
    with open(ctx.output_dir / "03_pathogen_risk_report.txt", "w") as f:
        f.write(report)


def _generate_comprehensive_report(ctx: PipelineContext, reporter) -> None:
    bracken_file = ctx.classification.bracken_file if ctx.classification else None
    sample_info_file = ctx.annotation.prokka_dir / "sample.txt" if ctx.annotation else None
    functional_file = ctx.annotation.functional_annotations if ctx.annotation else None
    pathogen_hits_file = ctx.pathogen.pathogen_hits_file if ctx.pathogen else None
    ml_predictions_file = ctx.pathogen.ml_predictions_file if ctx.pathogen else None

    try:
        reporter.generate_report(
            bracken_file=bracken_file,
            sample_info_file=sample_info_file,
            functional_annotations_file=functional_file,
            pathogen_hits_file=pathogen_hits_file,
            ml_predictions_file=ml_predictions_file,
        )

        if bracken_file and ctx.pathogen and ctx.pathogen.risk_data:
            reporter.export_tables(
                bracken_file=bracken_file,
                annotation_file=functional_file,
                risk_data=ctx.pathogen.risk_data,
            )
    except Exception as e:
        logger.warning("Comprehensive report generation failed: %s", e)


def load_bracken_report_safe(bracken_file):
    """Load bracken report, return empty DataFrame on failure."""
    import pandas as pd
    try:
        from metaquest.io.data_loaders import load_bracken_report
        return load_bracken_report(bracken_file)
    except Exception:
        return pd.DataFrame()
