"""Annotation stage — Pyrodigal gene prediction and DIAMOND search."""

from __future__ import annotations

import logging
from pathlib import Path

from metaquest.exceptions import AnnotationError
from metaquest.pipeline.context import PipelineContext, AnnotationResult

logger = logging.getLogger(__name__)


def run_annotation_stage(ctx: PipelineContext) -> PipelineContext:
    """Run gene calling followed by provisional functional annotation."""
    from metaquest.core.functional_analysis import (
        run_gene_prediction,
        run_functional_annotation,
    )

    if ctx.assembly is None:
        raise AnnotationError("Assembly stage must complete before annotation")

    config = ctx.config.annotation
    contigs_fasta = ctx.assembly.contigs_fasta

    logger.info("Running Pyrodigal in metagenomic mode")
    try:
        prediction_dir, gene_count = run_gene_prediction(
            contigs_fasta,
            ctx.output_dir,
            min_contig_length=config.min_contig_length,
        )
    except Exception as e:
        raise AnnotationError(f"Pyrodigal gene prediction failed: {e}", cause=e) from e

    prediction_dir = Path(prediction_dir)
    protein_file = prediction_dir / "genes.faa"

    if gene_count == 0:
        logger.warning("No proteins predicted (contigs may be too short)")
        ctx.annotation = AnnotationResult(
            gene_prediction_dir=prediction_dir,
            protein_file=protein_file,
        )
        return ctx

    logger.info("Pyrodigal predicted %d genes", gene_count)

    # Run DIAMOND functional annotation
    logger.info("Running DIAMOND annotation (threads=%d)", config.threads)
    try:
        annotation_file = run_functional_annotation(
            prediction_dir, ctx.output_dir,
            threads=config.threads,
            evalue=config.evalue,
        )
    except Exception as e:
        raise AnnotationError(f"DIAMOND annotation failed: {e}", cause=e) from e

    annotated_count = 0
    if annotation_file and Path(annotation_file).exists():
        with open(annotation_file) as f:
            annotated_count = len(set(line.split("\t")[0] for line in f if line.strip()))

    ctx.annotation = AnnotationResult(
        gene_prediction_dir=prediction_dir,
        protein_file=protein_file,
        functional_annotations=Path(annotation_file) if annotation_file else None,
        gene_count=gene_count,
        annotated_count=annotated_count,
    )

    logger.info("Annotation complete: %d/%d genes annotated", annotated_count, gene_count)
    return ctx
