"""Annotation stage — metagenomic gene prediction with Pyrodigal."""

from __future__ import annotations

import logging
from pathlib import Path

from metaquest.exceptions import AnnotationError
from metaquest.pipeline.context import AnnotationResult, PipelineContext

logger = logging.getLogger(__name__)


def run_annotation_stage(ctx: PipelineContext) -> PipelineContext:
    """Predict protein-coding genes from assembled metagenomic contigs."""
    from metaquest.core.functional_analysis import run_gene_prediction

    if ctx.assembly is None:
        raise AnnotationError("Assembly stage must complete before gene prediction")

    config = ctx.config.annotation
    logger.info("Running Pyrodigal in metagenomic mode")
    try:
        prediction_dir, gene_count = run_gene_prediction(
            ctx.assembly.contigs_fasta,
            ctx.output_dir,
            min_contig_length=config.min_contig_length,
        )
    except Exception as exc:
        raise AnnotationError(
            f"Pyrodigal gene prediction failed: {exc}",
            cause=exc,
        ) from exc

    prediction_dir = Path(prediction_dir)
    ctx.annotation = AnnotationResult(
        gene_prediction_dir=prediction_dir,
        protein_file=prediction_dir / "genes.faa",
        gene_count=gene_count,
    )
    logger.info("Gene prediction complete: %d genes", gene_count)
    return ctx
