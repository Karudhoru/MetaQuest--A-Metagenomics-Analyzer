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
    from metaquest.core.taxonomic_analysis import run_kraken, run_bracken
    from metaquest.io.data_loaders import load_bracken_report
    from metaquest.io.utils import split_interleaved

    config = ctx.config.classification
    output_dir = ctx.output_dir
    reads = ctx.input_files

    # Handle interleaved reads
    if ctx.read_mode == "interleaved":
        reads = split_interleaved(str(reads[0]), output_dir)
        if isinstance(reads, (str, Path)):
            reads = [Path(reads)]

    # Ensure reads are list of paths
    if isinstance(reads, (str, Path)):
        reads = [reads]
    reads = [Path(r) if not isinstance(r, Path) else r for r in reads]

    # Run Kraken2
    logger.info("Running Kraken2 (threads=%d)", config.threads)
    try:
        kraken_report = run_kraken(
            reads, output_dir,
            threads=config.threads,
            min_hit_groups=config.min_hit_groups,
        )
    except Exception as e:
        raise ClassificationError(
            f"Kraken2 classification failed: {e}",
            cause=e,
        ) from e

    # Run Bracken
    logger.info(
        "Running Bracken (level=%s, read_len=%d, threshold=%d)",
        config.taxonomic_level, config.read_length, config.bracken_threshold,
    )
    try:
        bracken_file = run_bracken(
            kraken_report, output_dir,
            read_length=config.read_length,
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
    classified_reads = int(bracken_df["new_est_reads"].sum()) if "new_est_reads" in bracken_df.columns else 0

    ctx.classification = ClassificationResult(
        kraken_report=Path(kraken_report),
        bracken_file=Path(bracken_file),
        species_count=species_count,
        classified_reads=classified_reads,
    )

    logger.info("Classification complete: %d species, %d classified reads", species_count, classified_reads)
    return ctx
