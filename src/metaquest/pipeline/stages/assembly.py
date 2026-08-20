"""
Assembly Stage — metagenomic assembly with MEGAHIT.
"""

from __future__ import annotations

import gc
import logging
from pathlib import Path

import numpy as np
from Bio import SeqIO

from metaquest.exceptions import AssemblyError
from metaquest.pipeline.context import PipelineContext, AssemblyResult

logger = logging.getLogger(__name__)


def run_assembly_stage(ctx: PipelineContext) -> PipelineContext:
    """Run metagenomic assembly and compute statistics."""
    from metaquest.core.assembly import assemble_reads_to_fasta
    from metaquest.io.output_formatter import get_formatter

    config = ctx.config.assembly
    fmt = get_formatter()

    reads = ctx.input_files
    if isinstance(reads, (str, Path)):
        reads = [reads]

    logger.info("Assembling with %s (threads=%d)", config.assembler, config.threads)

    existing_contigs = ctx.output_dir / "megahit_assembly" / "final.contigs.fa"
    if ctx.resume and existing_contigs.is_file() and existing_contigs.stat().st_size:
        contigs_fasta = existing_contigs
        fmt.info("Reusing completed MEGAHIT assembly")
    else:
        try:
            contigs_fasta = assemble_reads_to_fasta(
                reads, ctx.output_dir, fmt,
                threads=config.threads,
            )
        except Exception as e:
            raise AssemblyError(f"Assembly failed: {e}", cause=e) from e

    # Compute stats (streaming, memory-efficient)
    stats = _compute_assembly_stats(contigs_fasta)

    ctx.assembly = AssemblyResult(
        contigs_fasta=contigs_fasta,
        total_contigs=stats["total_contigs"],
        total_bases=stats["total_bases"],
        n50=stats["n50"],
        max_length=stats["max_length"],
        mean_length=stats["mean_length"],
    )

    logger.info(
        "Assembly complete: %d contigs, N50=%d, total=%d bp",
        stats["total_contigs"], stats["n50"], stats["total_bases"],
    )
    return ctx


def _compute_assembly_stats(contigs_fasta: Path) -> dict:
    """Calculate assembly statistics without loading all sequences into memory."""
    lengths = []
    for record in SeqIO.parse(contigs_fasta, "fasta"):
        lengths.append(len(record.seq))

    if not lengths:
        return {"total_contigs": 0, "total_bases": 0, "n50": 0, "max_length": 0, "mean_length": 0.0}

    total_contigs = len(lengths)
    total_bases = sum(lengths)
    max_length = max(lengths)
    mean_length = total_bases / total_contigs

    # N50
    sorted_lengths = sorted(lengths, reverse=True)
    cumsum = np.cumsum(sorted_lengths)
    threshold = total_bases / 2.0
    n50 = 0
    for i, cs in enumerate(cumsum):
        if cs >= threshold:
            n50 = sorted_lengths[i]
            break

    del lengths, sorted_lengths, cumsum
    gc.collect()

    return {
        "total_contigs": total_contigs,
        "total_bases": total_bases,
        "n50": int(n50),
        "max_length": max_length,
        "mean_length": mean_length,
    }
