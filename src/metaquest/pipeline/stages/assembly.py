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

    reads = ctx.analysis_input_files or ctx.input_files
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
        l50=stats["l50"],
        n90=stats["n90"],
        l90=stats["l90"],
        gc_percent=stats["gc_percent"],
        contigs_ge_1000=stats["contigs_ge_1000"],
        contigs_ge_5000=stats["contigs_ge_5000"],
        contigs_ge_10000=stats["contigs_ge_10000"],
        lengths=stats["lengths"],
    )

    logger.info(
        "Assembly complete: %d contigs, N50=%d, total=%d bp",
        stats["total_contigs"], stats["n50"], stats["total_bases"],
    )
    return ctx


def _compute_assembly_stats(contigs_fasta: Path) -> dict:
    """Calculate assembly statistics without loading all sequences into memory."""
    lengths = []
    gc_bases = 0
    for record in SeqIO.parse(contigs_fasta, "fasta"):
        lengths.append(len(record.seq))
        sequence = str(record.seq).upper()
        gc_bases += sequence.count("G") + sequence.count("C")

    if not lengths:
        return {"total_contigs": 0, "total_bases": 0, "n50": 0, "l50": 0,
                "n90": 0, "l90": 0, "max_length": 0, "mean_length": 0.0,
                "gc_percent": 0.0, "contigs_ge_1000": 0, "contigs_ge_5000": 0,
                "contigs_ge_10000": 0, "lengths": []}

    total_contigs = len(lengths)
    total_bases = sum(lengths)
    max_length = max(lengths)
    mean_length = total_bases / total_contigs

    # N50
    sorted_lengths = sorted(lengths, reverse=True)
    cumsum = np.cumsum(sorted_lengths)
    threshold = total_bases / 2.0
    n50 = l50 = n90 = l90 = 0
    for i, cs in enumerate(cumsum):
        if not n50 and cs >= threshold:
            n50 = sorted_lengths[i]
            l50 = i + 1
        if not n90 and cs >= total_bases * 0.9:
            n90 = sorted_lengths[i]
            l90 = i + 1
            break

    retained_lengths = list(lengths)
    del lengths, sorted_lengths, cumsum
    gc.collect()

    return {
        "total_contigs": total_contigs,
        "total_bases": total_bases,
        "n50": int(n50),
        "l50": l50,
        "n90": int(n90),
        "l90": l90,
        "max_length": max_length,
        "mean_length": mean_length,
        "gc_percent": gc_bases / total_bases * 100,
        "contigs_ge_1000": sum(length >= 1000 for length in retained_lengths),
        "contigs_ge_5000": sum(length >= 5000 for length in retained_lengths),
        "contigs_ge_10000": sum(length >= 10000 for length in retained_lengths),
        "lengths": retained_lengths,
    }
