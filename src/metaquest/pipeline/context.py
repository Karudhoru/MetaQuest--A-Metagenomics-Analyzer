"""
Pipeline Context — carries state between stages.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from metaquest.settings import MetaQuestConfig


@dataclass
class AssemblyResult:
    contigs_fasta: Path
    total_contigs: int = 0
    total_bases: int = 0
    n50: int = 0
    max_length: int = 0
    mean_length: float = 0.0


@dataclass
class ClassificationResult:
    kraken_report: Path
    bracken_file: Path
    species_count: int = 0
    classified_reads: int = 0


@dataclass
class AnnotationResult:
    gene_prediction_dir: Path
    protein_file: Path | None = None
    functional_annotations: Path | None = None
    gene_count: int = 0
    annotated_count: int = 0


@dataclass
class PipelineContext:
    """Accumulates results as stages execute sequentially."""

    config: MetaQuestConfig
    input_files: list[Path]
    output_dir: Path
    read_mode: str = "paired"  # paired | single | interleaved
    skip_annotation: bool = False

    # Filled by stages
    classification: ClassificationResult | None = None
    assembly: AssemblyResult | None = None
    annotation: AnnotationResult | None = None

    # Metadata for reproducibility
    metadata: dict[str, Any] = field(default_factory=dict)
    completed_stages: list[str] = field(default_factory=list)
