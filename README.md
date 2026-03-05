# MetaQuest

**A Comprehensive Metagenomics Analysis Pipeline**

---

## Overview

MetaQuest is an integrated bioinformatics pipeline for metagenomic data analysis. It combines state-of-the-art tools and databases to provide researchers with a streamlined workflow for understanding microbial community composition, identifying potential pathogens, and assessing antimicrobial resistance profiles from sequencing data.

**Key Capabilities:**

- **Taxonomic Classification** — High-accuracy taxonomic profiling via Kraken2/Bracken.
- **Pathogen Detection** — Fragment-aware confidence scoring with 100% specificity (v4.1 DB).
- **Functional Annotation** — Dual-database system: COG + SwissProt via DIAMOND/Prokka.
- **Machine Learning** — ML-enhanced pathogen prediction and biomarker discovery.
- **Comparative Analysis** — Compositional normalization (CLR/TSS), multi-group stats, effect sizes, power analysis.
- **Validation Engine** — Multi-evidence gene-to-species linkage with Bracken confidence scoring.
- **Professional Reporting** — Publication-ready text reports with scientific formatting.
- **Interactive Visualizations** — Modular, publication-ready dashboard and plots.

---

## Table of Contents

- [Development Status](#development-status)
- [What's New in v5.0.0](#whats-new-in-v500)
- [Features](#features)
- [System Architecture](#system-architecture)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Documentation](#documentation)
- [Contributing](#contributing)
- [Support](#support)

---

## Development Status

**Current Version:** 5.0.0 (March 2026)

### Completed Features

| Component | Status | Description |
|-----------|--------|-------------|
| **File Validation & Quality Control** | ✅ Complete | Comprehensive validation for both FASTQ and FASTA inputs with detailed statistics. |
| **Taxonomic Classification** | ✅ Complete | FASTQ and FASTA taxonomic profiling via Kraken2/Bracken. |
| **Pathogen Detection** | ✅ Complete | Fragment-aware confidence scoring with 100% specificity (housekeeping genes filtered). |
| **Validation Engine** | ✅ Complete | Vectorized gene-to-species linkage, Bracken confidence scoring, multi-evidence validation. |
| **Machine Learning Integration** | ✅ Complete | ML-based pathogen prediction with feature extraction and pre-trained model artifacts. |
| **Comparative Analysis** | ✅ Complete | v2.0: CLR/TSS normalization, power analysis, Kruskal-Wallis, Cliff's Delta, R². |
| **Beta Diversity** | ✅ Complete | PERMANOVA (with R²) and ANOSIM with PCoA visualization. |
| **Alpha Diversity** | ✅ Complete | Multi-group Kruskal-Wallis with post-hoc tests and Cliff's Delta effect size. |
| **Differential Abundance** | ✅ Complete | Prevalence filtering and Cliff's Delta effect size. |
| **ML Biomarker Discovery** | ✅ Complete | GridSearchCV hyperparameter tuning + SelectKBest feature selection. |
| **Functional Annotation** | ✅ Complete | COG + SwissProt dual-database with methodology-accurate reporting. |
| **Gene Prediction** | ✅ Complete | Prokka with contig filtering, tbl2asn timeout controls, configurable threading. |
| **Professional Logging** | ✅ Complete | Dual-mode (standard/debug) via OutputFormatter. |
| **Reporting System** | ✅ Complete | Publication-ready taxonomic, functional, and pathogen risk reports. |

### In Development

| Component | Status | Target |
|-----------|--------|--------|
| **Pathogen DB Benchmarking** | In Progress | Q2 2026 |
| **Metabolic Pathway Reconstruction** | Planning | Q3 2026 |
| **Advanced ML Models (v4.2)** | Planning | Q2 2026 |

### Release History

See [CHANGELOG.md](CHANGELOG.md) for the full version history.

| Version | Date | Summary |
|---------|------|---------|
| **v5.0.0** | Mar 2026 | Validation engine, fragment-aware pathogen scoring, publication-ready reports, scientific audit fixes, code cleanup. |
| **v4.1.0** | Nov 2025 | Comparative Analysis v2.0, Pathogen DB v4.1 (100% sensitivity & specificity), effect sizes. |
| **v4.0.0** | Oct 2025 | COG+SwissProt annotation, ML model v2.1, professional logging system. |
| **v3.6.0** | Oct 2025 | SPAdes assembly, modular visualization system. |

---

## What's New in v5.0.0

### 1. Validation Engine

A new `validation_engine.py` provides multi-evidence pathogen validation:

- **Vectorized gene-to-species linkage** via pandas merge (5–10× faster than previous loop implementation).
- **Bracken confidence scoring** — distinguishes direct Kraken hits from inferred reads to assign HIGH/MODERATE/LOW species confidence. Risk scores reduced 30% for low-confidence dominant taxa.
- **Exact species matching** — no genus-level mismatches.

### 2. Fragment-Aware Pathogen Detection

Replaced simple identity/coverage filter with a 7-factor confidence model:

- Fragment completeness (sigmoid curve), sequence identity × coverage, organism abundance (log₁₀ scale), contig sequencing depth, critical motif detection, Prokka cross-validation, signal peptide detection.
- Tiers: **HIGH** (≥70), **MEDIUM** (45–69), **LOW** (<45).

### 3. Scientific Report Fixes

- **Risk scoring**: z-score baseline (gut microbiome-specific) replaced with sample-agnostic absolute functional density scoring.
- **Taxonomic reporter**: classification rate now from Kraken2 `unclassified` row, not Bracken estimates.
- **Functional reporter**: annotation counts properly distinguish unique annotated proteins from total DIAMOND alignments.
- **All reporters**: version unified to 5.0.0.

### 4. New I/O Modules

- `io/data_loaders.py` — centralized, memory-efficient data loading.
- `io/text_parsers.py` — shared parsing utilities.
- `config/pathogen_config.json` — externalized configuration for easier DB updates.

---

## Features

### Core Analysis Capabilities

**File Validation & QC**
- FASTQ: quality score analysis, contamination detection, MD5 checksums.
- FASTA: sequence type detection, composition analysis, N50 metrics.

**Taxonomic Classification**
- Kraken2/Bracken classification with confidence-weighted species assignment.
- True classification rate from Kraken2 unclassified read counts.

**Pathogen Detection**
- v4.1 database: housekeeping genes filtered, 100% sensitivity and specificity.
- 7-factor fragment confidence scoring with tiered output (HIGH/MEDIUM/LOW).
- Bracken confidence integration to flag low-evidence detections.

**Functional Annotation**
- Prokka gene prediction with configurable contig filtering and timeout management.
- DIAMOND against COG + SwissProt databases.
- COG category analysis, IS family classification, mobile genetic element tracking.

**Comparative Analysis (v2.0)**
- Normalization: TSS (relative abundance), CLR (compositional). Rarefaction deprecated.
- Alpha diversity: Shannon, Simpson, Chao1 with Kruskal-Wallis and post-hoc tests.
- Beta diversity: Bray-Curtis, PERMANOVA (R²), ANOSIM.
- Differential abundance with prevalence filtering and Cliff's Delta.
- ML biomarker discovery with GridSearchCV + SelectKBest.
- Statistical power assessment with detectable effect size estimates.

**Reporting**
- Taxonomic report: classification metrics, BSL lookup, confidence-annotated species tables.
- Functional report: annotation coverage, COG categories, methodology notes.
- Pathogen risk report: three-tier scoring (taxonomic + functional + ML), sample-agnostic interpretation.

**Visualization**
- Interactive HTML dashboard.
- Taxonomic, functional, pathogen, and comparative visualizations.
- PCoA, volcano plots, heatmaps, alpha diversity boxplots.

---

## System Architecture

```
MetaQuest/
├── README.md
├── CHANGELOG.md
├── setup.py
├── requirements.txt
├── environment/
│   └── environment.yml
├── scripts/
│   └── setup_databases.sh
├── docs/
│   ├── installation.md
│   └── annotation.md
└── src/metaquest/
    ├── cli.py                              # CLI entry point and argument parsing
    ├── config.py                           # Central configuration and constants
    │
    ├── config/                             # Externalized configuration
    │   └── pathogen_config.json            # Pathogen DB config and critical motifs
    │
    ├── core/                               # Core analysis modules
    │   ├── analysis.py                     # Main pipeline orchestration
    │   ├── taxonomic_analysis.py           # Kraken2/Bracken classification
    │   ├── pathogen_analysis.py            # Fragment-aware pathogen detection (v2)
    │   ├── functional_analysis.py          # Prokka + DIAMOND annotation
    │   ├── assembly.py                     # MEGAHIT/SPAdes assembly wrapper
    │   └── comparative_analysis.py         # Multi-sample statistical comparison (v2.0)
    │
    ├── io/                                 # Input/Output modules
    │   ├── file_validator.py               # FASTQ/FASTA validation and QC
    │   ├── output_formatter.py             # Professional logging and output system
    │   ├── data_loaders.py                 # Centralized data loading utilities
    │   ├── text_parsers.py                 # Shared text parsing utilities
    │   └── utils.py                        # General I/O utilities
    │
    ├── ml/                                 # Machine Learning
    │   ├── feature_extractor.py            # Feature extraction for ML models
    │   ├── pathogen_predictor.py           # ML pathogen classification
    │   └── model_artifacts/               # Pre-trained model files
    │
    ├── reporting/                          # Report generation
    │   ├── main_reporter.py               # Report orchestration
    │   ├── base_reporter.py               # Base class with formatting utilities
    │   ├── validation_engine.py           # Multi-evidence validation and confidence
    │   ├── taxonomy_reporter.py           # Taxonomic reports
    │   ├── functional_reporter.py         # Functional annotation reports
    │   ├── pathogen_risk_reporter.py      # Pathogen risk reports
    │   └── risk_scoring.py                # Risk score calculation engine
    │
    └── visualization/                     # Visualization system
        ├── main_visualizer.py             # Visualization coordinator
        ├── base_visualizer.py             # Base class
        ├── dashboard.py                   # Interactive HTML dashboard
        ├── taxonomic_visualizer.py        # Taxonomic plots
        ├── functional_visualizer.py       # Functional annotation plots
        ├── pathogenic_visualizer.py       # Pathogen detection plots
        └── compare_visuals.py             # Comparative analysis plots
```

---

## Installation

MetaQuest requires a Linux/macOS environment with conda.

**For detailed instructions, see [docs/installation.md](docs/installation.md)**

### System Requirements

- Linux or macOS
- Conda package manager
- Minimum 16 GB RAM (32 GB recommended)
- ~100 GB disk space for databases

### Quick Installation

```bash
git clone https://github.com/Karudhoru/MetaQuest.git
cd MetaQuest

conda env create -f environment/environment.yml
conda activate metaquest

pip install -e .

# Download databases (COG + SwissProt + Kraken2)
./scripts/setup_databases.sh --swissprot --kraken

# Build pathogen database (v4.1)
python scripts/custom_pathogen_db.py
./scripts/setup_databases.sh --pathogen

# Verify
metaquest check
```

---

## Quick Start

### Validate Input Files

```bash
metaquest validate fastq --single sample.fastq.gz
metaquest validate fasta genome.fasta
```

### Run Analysis

```bash
# FASTQ analysis
metaquest analyze fastq --single sample.fastq.gz -o results/

# FASTA analysis
metaquest analyze fasta genome.fasta -o results/ -s 100

# Custom contig filtering
metaquest analyze fastq --single reads.fq --min-contig-length 500 -o results/

# Skip annotation (taxonomy only)
metaquest analyze fastq --single sample.fastq.gz --skip-annotation -o results/

# Debug mode
metaquest --debug analyze fastq --single reads.fq -o results/
```

### Comparative Analysis

```bash
# metadata.tsv format:
# sample_id    group
# sample1      Healthy
# sample2      Diseased

# TSS normalization (default)
metaquest compare -i sample1/ sample2/ sample3/ -m metadata.tsv -o comparison/

# CLR normalization (compositional)
metaquest compare -i sample*/ -m metadata.tsv -o comparison/ --normalization clr
```

### Output Files

| Category | Files |
|----------|-------|
| **Reports** | `taxonomic_report.txt`, `functional_report.txt`, `pathogen_risk_report.txt` |
| **Validation** | `pathogen_detections_validated.json`, `pathogen_results.tsv` |
| **Taxonomy** | `bracken_report.tsv`, `kraken2_report.txt` |
| **Annotation** | `functional_annotations.tsv`, `prokka_annotation/` |
| **Comparative** | `alpha_diversity_metrics.tsv`, `permanova_results.txt`, PCoA/volcano plots |
| **Visualization** | `dashboard.html`, `*.png` plots |

---

## Documentation

- **[Installation Guide](docs/installation.md)**
- **[Annotation Guide](docs/annotation.md)**
- **[Changelog](CHANGELOG.md)**
- API Documentation — *coming soon*

---

## Contributing

Contributions are welcome from the scientific community.

### Priority Areas

| Area | Impact |
|------|--------|
| Pathogen DB benchmarking (formal publication) | High |
| Advanced ML models beyond Random Forest | High |
| Additional databases (KEGG, Pfam) | High |
| Clinical validation on real patient samples | Critical |

### How to Contribute

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/your-feature`
3. Commit with descriptive messages
4. Submit a pull request with a detailed description

---

## Support

- **Bug Reports**: GitHub Issues with reproduction steps
- **Feature Requests**: GitHub Discussions
- **Email**: metaquest-support@example.org
- **Troubleshooting**: Use `--debug` flag; check `metaquest.log` in output directory

---

## Citation

*Citation information will be provided upon publication.*

## License

*License information will be added upon publication.*

---

**Last Updated:** March 2026 — v5.0.0