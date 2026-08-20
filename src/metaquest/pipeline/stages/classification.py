"""
Classification Stage — Kraken2 + Bracken taxonomic profiling.
"""

from __future__ import annotations

import logging
from pathlib import Path

from metaquest.exceptions import ClassificationError
from metaquest.pipeline.context import PipelineContext, ClassificationResult

logger = logging.getLogger(__name__)


def run_classification_stage(ctx: PipelineContext) -> PipelineContext:
    """
    Run taxonomic classification with Kraken2 and abundance estimation with Bracken.
    """
    from metaquest.core.taxonomic_analysis import (
        infer_read_length,
        parse_kraken_read_counts,
        run_bracken,
        run_kraken,
        select_bracken_read_length,
    )
    from metaquest.io.data_loaders import load_bracken_report
    from metaquest.io.utils import split_interleaved
    from metaquest.io.output_formatter import get_formatter

    config = ctx.config.classification
    output_dir = ctx.output_dir
    reads = ctx.analysis_input_files or ctx.input_files

    # Handle interleaved reads
    if ctx.read_mode == "interleaved" and not ctx.analysis_input_files:
        fmt = get_formatter()
        existing_split = [
            output_dir / "split_R1.fastq.gz",
            output_dir / "split_R2.fastq.gz",
        ]
        if ctx.resume and all(path.is_file() for path in existing_split):
            reads = existing_split
            fmt.info("Reusing synchronized interleaved read split")
        else:
            reads = split_interleaved(str(reads[0]), output_dir, fmt)
        if isinstance(reads, (str, Path)):
            reads = [Path(reads)]

    # Ensure reads are list of paths
    if isinstance(reads, (str, Path)):
        reads = [reads]
    reads = [Path(r) if not isinstance(r, Path) else r for r in reads]

    observed_read_length = infer_read_length(reads)
    bracken_read_length = select_bracken_read_length(
        ctx.config.databases.kraken_db, observed_read_length
    )
    if bracken_read_length != observed_read_length:
        get_formatter().warning(
            f"No Bracken model for {observed_read_length} bp reads; "
            f"using the closest installed model ({bracken_read_length} bp)"
        )
    logger.info(
        "Running Bracken (level=%s, observed=%d, model=%d, threshold=%d)",
        config.taxonomic_level,
        observed_read_length,
        bracken_read_length,
        config.bracken_threshold,
    )
    existing_kraken = output_dir / "kraken_report.txt"
    existing_bracken = output_dir / "bracken_report.tsv"
    if ctx.resume and existing_kraken.is_file() and existing_bracken.is_file():
        kraken_report = existing_kraken
        bracken_file = existing_bracken
        get_formatter().info("Reusing completed Kraken2 and Bracken results")
    else:
        logger.info("Running Kraken2 (threads=%d)", config.threads)
        try:
            kraken_report = run_kraken(
                reads,
                output_dir,
                threads=config.threads,
                min_hit_groups=config.min_hit_groups,
                memory_mapping=ctx.low_memory,
            )
        except Exception as e:
            raise ClassificationError(
                f"Kraken2 classification failed: {e}",
                cause=e,
            ) from e

        try:
            bracken_file = run_bracken(
                kraken_report,
                output_dir,
                read_length=bracken_read_length,
                level=config.taxonomic_level,
                threshold=config.bracken_threshold,
            )
        except Exception as e:
            raise ClassificationError(
                f"Bracken abundance estimation failed: {e}",
                cause=e,
            ) from e

    # Collect stats
    bracken_df = load_bracken_report(bracken_file)
    species_count = len(bracken_df)
    species_assigned_reads = (
        int(bracken_df["new_est_reads"].sum())
        if "new_est_reads" in bracken_df.columns
        else 0
    )
    total_reads, kraken_classified_reads, unclassified_reads = parse_kraken_read_counts(
        kraken_report
    )

    ctx.classification = ClassificationResult(
        kraken_report=Path(kraken_report),
        bracken_file=Path(bracken_file),
        species_count=species_count,
        species_assigned_reads=species_assigned_reads,
        total_reads=total_reads,
        kraken_classified_reads=kraken_classified_reads,
        unclassified_reads=unclassified_reads,
        observed_read_length=observed_read_length,
        bracken_read_length=bracken_read_length,
    )

    logger.info(
        "Classification complete: %d species, %d/%d Kraken-classified reads",
        species_count,
        kraken_classified_reads,
        total_reads,
    )
    return ctx
