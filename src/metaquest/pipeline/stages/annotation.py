"""
Annotation Stage — Prokka gene prediction + DIAMOND functional annotation.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path

from Bio import SeqIO

from metaquest.exceptions import AnnotationError
from metaquest.pipeline.context import PipelineContext, AnnotationResult

logger = logging.getLogger(__name__)


def run_annotation_stage(ctx: PipelineContext) -> PipelineContext:
    """Run Prokka for gene calling and DIAMOND for functional annotation."""
    from metaquest.core.functional_analysis import (
        run_prokka,
        run_functional_annotation,
    )

    if ctx.assembly is None:
        raise AnnotationError("Assembly stage must complete before annotation")

    config = ctx.config.annotation
    contigs_fasta = ctx.assembly.contigs_fasta

    # Run Prokka
    logger.info("Running Prokka (threads=%d, mode=%s)", config.threads, config.mode)
    try:
        prokka_dir = run_prokka(
            contigs_fasta,
            ctx.output_dir,
            threads=config.threads,
            mode=config.mode,
        )
    except Exception as e:
        raise AnnotationError(f"Prokka gene prediction failed: {e}", cause=e) from e

    prokka_dir = Path(prokka_dir)
    protein_file = prokka_dir / "sample.faa"

    if not protein_file.exists() or os.path.getsize(protein_file) == 0:
        logger.warning("No proteins predicted (contigs may be too short)")
        ctx.annotation = AnnotationResult(prokka_dir=prokka_dir)
        return ctx

    gene_count = sum(1 for _ in SeqIO.parse(protein_file, "fasta"))
    logger.info("Prokka predicted %d genes", gene_count)

    # Run DIAMOND functional annotation
    logger.info("Running DIAMOND annotation (threads=%d)", config.threads)
    try:
        annotation_file = run_functional_annotation(
            prokka_dir, ctx.output_dir,
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
        prokka_dir=prokka_dir,
        protein_file=protein_file,
        functional_annotations=Path(annotation_file) if annotation_file else None,
        gene_count=gene_count,
        annotated_count=annotated_count,
    )

    logger.info("Annotation complete: %d/%d genes annotated", annotated_count, gene_count)
    return ctx
