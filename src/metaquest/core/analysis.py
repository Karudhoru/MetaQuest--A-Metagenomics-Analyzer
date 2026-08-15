"""
MetaQuest Core Analysis Module — FASTQ-only v6.0.0
===================================================
Delegates to pipeline.runner for actual orchestration.
This module provides backward-compatible entry points.
"""

from pathlib import Path
from dataclasses import replace

from ..settings import get_config, load_config
from ..pipeline.runner import build_default_pipeline
from ..pipeline.context import PipelineContext


def run_analysis(input_file, output_dir, cli_args=None):
    """
    Main analysis entry point — backward compatible.

    Args:
        input_file: FASTQ file path(s).
        output_dir: Output directory path.
        cli_args: Command-line arguments namespace.
    """
    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)

    # Load config (from cli_args.config if available, else defaults)
    config_path = getattr(cli_args, "config", None)
    db_dir = getattr(cli_args, "db_dir", None)
    if config_path or db_dir:
        config = load_config(
            Path(config_path) if config_path else None,
            db_dir=Path(db_dir) if db_dir else None,
        )
    else:
        config = get_config()

    # Preserve the existing CLI while ensuring command-line values actually
    # reach the immutable configuration used by pipeline stages.
    annotation_threads = getattr(cli_args, "annotation_threads", None)
    if annotation_threads is not None:
        config = replace(
            config,
            annotation=replace(config.annotation, threads=annotation_threads),
        )

    # Determine read mode
    read_mode = "paired"
    if cli_args:
        if getattr(cli_args, "single", False):
            read_mode = "single"
        elif getattr(cli_args, "interleaved", False):
            read_mode = "interleaved"

    # Normalize input files
    if isinstance(input_file, str):
        input_files = [Path(input_file)]
    elif isinstance(input_file, (list, tuple)):
        input_files = [Path(f) for f in input_file]
    else:
        input_files = [Path(input_file)]

    skip_annotation = getattr(cli_args, "skip_annotation", False)

    # Build and run pipeline
    ctx = PipelineContext(
        config=config,
        input_files=input_files,
        output_dir=output_dir_path,
        read_mode=read_mode,
        skip_annotation=skip_annotation,
        low_memory=getattr(cli_args, "low_memory", False),
    )

    pipeline = build_default_pipeline(config, skip_annotation=skip_annotation)

    ctx = pipeline.run(ctx)
