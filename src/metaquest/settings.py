"""
MetaQuest Settings
==================
YAML-based configuration system. Single source of truth for all pipeline parameters.

Resolution order for database path:
  1. METAQUEST_DB_DIR environment variable
  2. databases.base_dir in config YAML
  3. ./databases (relative to working directory)
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from metaquest.exceptions import ConfigError


# ---------------------------------------------------------------------------
# Config dataclasses
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class DatabasesConfig:
    base_dir: Path
    kraken_db: Path
    pathogen_markers: Path
    swissprot_cog: Path
    kraken_db_version: str = "Standard-16 (Built: 2024-01-10)"


@dataclass(frozen=True)
class AssemblyConfig:
    assembler: str = "megahit"
    threads: int = 8
    min_contig_length: int = 500
    memory_limit_gb: int = 8
    kmer_min: int = 21
    kmer_max: int = 141
    kmer_step: int = 12


@dataclass(frozen=True)
class ClassificationConfig:
    threads: int = 8
    read_length: int = 150
    taxonomic_level: str = "S"
    min_hit_groups: int = 2
    bracken_threshold: int = 10


@dataclass(frozen=True)
class AnnotationConfig:
    tool: str = "prokka"
    threads: int = 8
    evalue: float = 1e-6
    mode: str = "metagenome"
    kill_tbl2asn: bool = False
    tbl2asn_timeout: int = 30
    rfam: bool = False
    min_contig_length: int = 200


@dataclass(frozen=True)
class PathogenDetectionConfig:
    min_fragment_length: int = 80
    min_identity: float = 80.0
    min_query_coverage: float = 0.7
    min_subject_coverage: float = 0.7
    evalue_threshold: float = 1e-5
    diamond_threads: int = 4
    diamond_sensitivity: str = "more-sensitive"
    confidence_high: int = 70
    confidence_medium: int = 45
    # Novel detection system
    use_hmm: bool = False
    use_esm: bool = False
    use_island_detection: bool = False
    use_bayesian_risk: bool = False
    hmm_evalue: float = 1e-5
    hmm_threads: int = 4
    esm_model: str = "esm2_t12_35M_UR50D"
    esm_batch_size: int = 32
    esm_min_delta_score: float = 0.0
    island_window_size: int = 20000
    island_step_size: int = 5000
    island_min_vf_genes: int = 3


@dataclass(frozen=True)
class MLConfig:
    enabled: bool = False
    confidence_threshold: float = 0.7
    batch_size: int = 100
    feature_cache: bool = True
    random_seed: int = 42
    model_artifacts_dir: Path | None = None
    uncertainty_threshold: float = 0.3


@dataclass(frozen=True)
class ComparativeConfig:
    normalization: str = "tss"
    min_prevalence: float = 0.10
    random_seed: int = 42
    permutations: int = 999
    alpha_level: float = 0.05
    fdr_method: str = "fdr_bh"


@dataclass(frozen=True)
class MemoryConfig:
    max_sequences_in_memory: int = 10000
    chunk_size: int = 1000
    streaming_enabled: bool = True


@dataclass(frozen=True)
class QCConfig:
    min_read_quality: int = 20
    min_read_length: int = 50
    min_classification_rate: float = 0.5
    min_species_level_rate: float = 0.3
    max_unclassified_rate: float = 0.5
    min_contig_coverage: float = 2.0


@dataclass(frozen=True)
class ReportingConfig:
    formats: list[str] = field(default_factory=lambda: ["txt", "json", "html"])
    detail_level: str = "full"
    include_plots: bool = True
    include_tables: bool = True
    confidence_warnings: bool = True
    max_table_rows: int = 100
    decimal_precision: int = 2


@dataclass(frozen=True)
class LoggingConfig:
    level: str = "INFO"
    log_to_file: bool = True
    log_to_console: bool = True
    format: str = "%(asctime)s [%(levelname)s] %(message)s"


@dataclass(frozen=True)
class MetaQuestConfig:
    databases: DatabasesConfig
    assembly: AssemblyConfig = field(default_factory=AssemblyConfig)
    classification: ClassificationConfig = field(default_factory=ClassificationConfig)
    annotation: AnnotationConfig = field(default_factory=AnnotationConfig)
    pathogen_detection: PathogenDetectionConfig = field(default_factory=PathogenDetectionConfig)
    ml: MLConfig = field(default_factory=MLConfig)
    comparative: ComparativeConfig = field(default_factory=ComparativeConfig)
    memory: MemoryConfig = field(default_factory=MemoryConfig)
    qc: QCConfig = field(default_factory=QCConfig)
    reporting: ReportingConfig = field(default_factory=ReportingConfig)
    logging: LoggingConfig = field(default_factory=LoggingConfig)


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

_DEFAULT_CONFIG_PATH = Path(__file__).parent / "metaquest_default.yaml"

_cached_config: MetaQuestConfig | None = None


def _resolve_db_dir(raw: dict) -> Path:
    """Resolve database directory from env > yaml > default."""
    env_val = os.environ.get("METAQUEST_DB_DIR")
    if env_val:
        return Path(env_val)
    yaml_val = raw.get("databases", {}).get("base_dir")
    if yaml_val:
        return Path(yaml_val)
    return Path("./databases")


def _build_databases_config(raw: dict, base_dir: Path) -> DatabasesConfig:
    db = raw.get("databases", {})
    kraken = Path(db.get("kraken_db") or str(base_dir))
    pathogen = base_dir / db.get("pathogen_markers", "metaquest_pathogen_markers.dmnd")
    swissprot = base_dir / db.get("swissprot_cog", "SwissProt_COG_db.dmnd")
    version = db.get("kraken_db_version", "Standard-16 (Built: 2024-01-10)")
    return DatabasesConfig(
        base_dir=base_dir,
        kraken_db=kraken,
        pathogen_markers=pathogen,
        swissprot_cog=swissprot,
        kraken_db_version=version,
    )


def _build_section(cls, raw: dict, key: str):
    """Build a frozen dataclass from a YAML section, ignoring unknown keys."""
    section = raw.get(key, {}) or {}
    valid_fields = {f.name for f in cls.__dataclass_fields__.values()}
    filtered = {k: v for k, v in section.items() if k in valid_fields}
    # Convert Path fields
    for fname, fobj in cls.__dataclass_fields__.items():
        if fname in filtered and filtered[fname] is not None:
            if fobj.type in ("Path", "Path | None"):
                filtered[fname] = Path(filtered[fname])
    return cls(**filtered)


def load_config(config_path: Path | None = None) -> MetaQuestConfig:
    """
    Load configuration from YAML file.

    Resolution:
      1. User-provided config_path
      2. Default config bundled with package

    Values in user config override defaults (shallow merge per section).
    """
    global _cached_config

    # Load defaults
    if not _DEFAULT_CONFIG_PATH.exists():
        raise ConfigError(f"Default config missing: {_DEFAULT_CONFIG_PATH}")

    with open(_DEFAULT_CONFIG_PATH) as f:
        defaults = yaml.safe_load(f) or {}

    # Merge user config on top
    raw = dict(defaults)
    if config_path:
        if not config_path.exists():
            raise ConfigError(f"Config file not found: {config_path}")
        with open(config_path) as f:
            user = yaml.safe_load(f) or {}
        for key, val in user.items():
            if isinstance(val, dict) and key in raw and isinstance(raw[key], dict):
                raw[key] = {**raw[key], **val}
            else:
                raw[key] = val

    db_dir = _resolve_db_dir(raw)
    databases = _build_databases_config(raw, db_dir)

    # Resolve ML model artifacts dir
    ml_section = raw.get("ml", {}) or {}
    if not ml_section.get("model_artifacts_dir"):
        ml_section["model_artifacts_dir"] = str(
            Path(__file__).parent / "ml" / "model_artifacts"
        )
    raw["ml"] = ml_section

    config = MetaQuestConfig(
        databases=databases,
        assembly=_build_section(AssemblyConfig, raw, "assembly"),
        classification=_build_section(ClassificationConfig, raw, "classification"),
        annotation=_build_section(AnnotationConfig, raw, "annotation"),
        pathogen_detection=_build_section(PathogenDetectionConfig, raw, "pathogen_detection"),
        ml=_build_section(MLConfig, raw, "ml"),
        comparative=_build_section(ComparativeConfig, raw, "comparative"),
        memory=_build_section(MemoryConfig, raw, "memory"),
        qc=_build_section(QCConfig, raw, "qc"),
        reporting=_build_section(ReportingConfig, raw, "reporting"),
        logging=_build_section(LoggingConfig, raw, "logging"),
    )

    _cached_config = config
    return config


def get_config() -> MetaQuestConfig:
    """Return cached config or load defaults."""
    global _cached_config
    if _cached_config is None:
        _cached_config = load_config()
    return _cached_config


def reset_config() -> None:
    """Clear cached config (for testing)."""
    global _cached_config
    _cached_config = None


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def validate_config(config: MetaQuestConfig | None = None) -> tuple[bool, list[str]]:
    """Validate configuration, return (is_valid, errors)."""
    cfg = config or get_config()
    errors: list[str] = []

    if cfg.memory.max_sequences_in_memory < 1000:
        errors.append("memory.max_sequences_in_memory must be >= 1000")
    if cfg.assembly.memory_limit_gb < 4:
        errors.append("assembly.memory_limit_gb should be >= 4")
    if not (0.0 <= cfg.ml.confidence_threshold <= 1.0):
        errors.append("ml.confidence_threshold must be between 0 and 1")
    if cfg.assembly.threads < 1:
        errors.append("assembly.threads must be >= 1")
    if cfg.pathogen_detection.min_identity < 0 or cfg.pathogen_detection.min_identity > 100:
        errors.append("pathogen_detection.min_identity must be 0-100")
    if cfg.comparative.normalization not in ("tss", "clr", "none"):
        errors.append("comparative.normalization must be tss, clr, or none")
    for fmt in cfg.reporting.formats:
        if fmt not in ("txt", "json", "html", "csv"):
            errors.append(f"reporting.formats: unknown format '{fmt}'")

    return len(errors) == 0, errors
