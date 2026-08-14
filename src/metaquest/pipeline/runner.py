"""
Pipeline Runner — orchestrates stages in sequence.
"""

from __future__ import annotations

import gc
import json
import logging
import time
from datetime import datetime
from pathlib import Path
from typing import Any

from metaquest.exceptions import MetaQuestError, PipelineStageError
from metaquest.settings import MetaQuestConfig
from metaquest.pipeline.context import PipelineContext

logger = logging.getLogger(__name__)


class PipelineRunner:
    """
    Runs the MetaQuest analysis pipeline as a sequence of stages.

    Each stage receives the PipelineContext and returns it (possibly mutated).
    Stages are functions: (ctx: PipelineContext) -> PipelineContext
    """

    def __init__(self, config: MetaQuestConfig):
        self.config = config
        self._stages: list[tuple[str, Any]] = []

    def add_stage(self, name: str, stage_fn):
        """Register a stage function."""
        self._stages.append((name, stage_fn))
        return self

    def run(self, ctx: PipelineContext) -> PipelineContext:
        """Execute all registered stages in order."""
        from metaquest.io.output_formatter import get_formatter

        total = len(self._stages)
        formatter = get_formatter()
        logger.info("Pipeline starting with %d stage(s)", total)
        formatter.section_header("PIPELINE PROGRESS")

        for i, (name, stage_fn) in enumerate(self._stages, 1):
            logger.info("[%d/%d] %s", i, total, name)
            formatter.info(f"[{i}/{total}] {name}")
            started = time.monotonic()
            try:
                with formatter.spinner(name):
                    ctx = stage_fn(ctx)
                ctx.completed_stages.append(name)
                elapsed = time.monotonic() - started
                formatter.success(
                    f"{name} completed in {formatter._format_time(elapsed)}"
                )
            except PipelineStageError:
                raise
            except MetaQuestError:
                raise
            except Exception as e:
                raise PipelineStageError(name, str(e), cause=e) from e
            finally:
                gc.collect()

        self._save_metadata(ctx)
        logger.info("Pipeline complete: %d stage(s) finished", len(ctx.completed_stages))
        return ctx

    def _save_metadata(self, ctx: PipelineContext) -> None:
        """Save pipeline metadata for reproducibility."""
        from metaquest.settings import get_config

        metadata = {
            "timestamp": datetime.now().isoformat(),
            "metaquest_version": _get_version(),
            "analysis_type": "fastq",
            "completed_stages": ctx.completed_stages,
            "config_summary": {
                "assembler": ctx.config.assembly.assembler,
                "threads": ctx.config.assembly.threads,
            },
            "input_files": [str(f) for f in ctx.input_files],
        }
        metadata.update(ctx.metadata)

        metadata_file = ctx.output_dir / "analysis_metadata.json"
        with open(metadata_file, "w") as f:
            json.dump(metadata, f, indent=2)


def _get_version() -> str:
    try:
        from metaquest.config import __version__
        return __version__
    except ImportError:
        return "unknown"


def build_default_pipeline(config: MetaQuestConfig, skip_annotation: bool = False) -> PipelineRunner:
    """
    Build the standard MetaQuest analysis pipeline.

    Stages:
      1. Taxonomic classification (Kraken2 + Bracken)
      2. Metagenomic assembly (MEGAHIT) [full workflow only]
      3. Gene prediction (Pyrodigal)
      4. Stable reporting

    Custom pathogen detection, ML, HMM, ESM, island detection, and risk
    scoring are intentionally excluded until independently validated.
    """
    from metaquest.pipeline.stages.classification import run_classification_stage
    from metaquest.pipeline.stages.reporting import run_reporting_stage

    runner = PipelineRunner(config)
    runner.add_stage("Taxonomic Classification", run_classification_stage)

    if not skip_annotation:
        from metaquest.pipeline.stages.assembly import run_assembly_stage
        from metaquest.pipeline.stages.annotation import run_annotation_stage

        runner.add_stage("Metagenomic Assembly", run_assembly_stage)
        runner.add_stage("Gene Prediction", run_annotation_stage)

    runner.add_stage("Reporting", run_reporting_stage)
    return runner
