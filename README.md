# MetaQuest

**A Comprehensive Metagenomics Analysis Pipeline**

-----

## Overview

MetaQuest is an integrated bioinformatics pipeline that addresses the complex challenges of metagenomic data analysis. By combining state-of-the-art tools and databases, MetaQuest provides researchers with a streamlined workflow for understanding microbial community composition, identifying potential pathogens, and assessing antimicrobial resistance profiles from sequencing data.

**Key Capabilities:**

  - **Taxonomic Classification**: High-accuracy taxonomic profiling.
  - **Pathogen Detection**: **Perfected 100% specificity** by filtering housekeeping genes.
  - **Machine learning-enhanced pathogen prediction** (Also in FASTQ).
  - **Advanced Comparative Analysis (v2.0)**: Compositional data normalization (CLR, rarefaction), multi-group stats (Kruskal-Wallis), effect sizes (Cliff's Delta, R²), and statistical power analysis.
  - **Comprehensive File Validation** and quality control.
  - **Interactive Visualizations** and professional reporting via a new **Output Formatter**.
  - **Enhanced Functional Annotation** with COG and SwissProt databases.
  - **Professional Logging System** with debug mode.

-----

## Table of Contents

  - [Development Status](https://www.google.com/search?q=%23development-status)
  - [What's New in Version 4.1.0](https://www.google.com/search?q=%23whats-new-in-version-410)
  - [Features](https://www.google.com/search?q=%23features)
  - [System Architecture](https://www.google.com/search?q=%23system-architecture)
  - [Benchmark Performance](https://www.google.com/search?q=%23benchmark-performance)
  - [Installation](https://www.google.com/search?q=%23installation)
  - [Quick Start](https://www.google.com/search?q=%23quick-start)
  - [Documentation](https://www.google.com/search?q=%23documentation)
  - [Contributing](https://www.google.com/search?q=%23contributing)
  - [Support](https://www.google.com/search?q=%23support)

-----

## Development Status

**Current Version:** 4.1.0 (November 2025)

MetaQuest continues active development with major improvements to statistical analysis, database accuracy, and user experience.

### Completed Features

| Component | Status | Description |
|-----------|---------|-------------|
| **File Validation & Quality Control** | ✅ Complete | Comprehensive validation with detailed statistics for both FASTQ and FASTA inputs. |
| **Taxonomic Classification** | ✅ Complete | Both FASTQ and FASTA taxonomic profiling are accurate and fully functional. |
| **Pathogen Detection** | ✅ Complete | **v4.1**: Comprehensive detection with **100% specificity** (housekeeping genes filtered). |
| **Machine Learning Integration** | ✅ Complete | Advanced ML-based pathogen prediction with feature extraction and model artifacts. |
| **Comparative Analysis** | ✅ Complete | **v2.0**: Advanced statistical comparison with normalization, power analysis, multi-group stats, and effect sizes. |
| **Beta Diversity Analysis** | ✅ Complete | PERMANOVA (with **R²**) and ANOSIM statistical tests with PCoA visualization. |
| **Alpha Diversity Analysis** | ✅ Complete | **Multi-group stats (Kruskal-Wallis)**, post-hoc tests, and **effect sizes (Cliff's Delta)**. |
| **Differential Abundance Testing** | ✅ Complete | Enhanced with **prevalence filtering** and **Cliff's Delta effect size**. |
| **Machine Learning Biomarker Discovery** | ✅ Complete | Enhanced with **hyperparameter tuning (GridSearchCV)** and **feature selection (SelectKBest)**. |
| **Enhanced Functional Annotation** | ✅ Complete | Dual-database annotation system combining COG and SwissProt with detailed reporting. |
| **Advanced Gene Prediction** | ✅ Complete | Prokka annotation with customizable contig filtering and timeout controls. |
| **Professional Logging System** | ✅ Complete | Dual-mode output (standard/debug) via **OutputFormatter**. |
| **Enhanced Reporting System** | ✅ Complete | Comprehensive text-based reports with clinical and research perspectives. |

### In Development

| Component | Status | Target Completion |
|-----------|---------|-------------------|
| **Pathogen DB Benchmarking** | In Progress | Q1 2026 |
| **Metabolic Pathway Reconstruction** | Planning | Q3 2026 |
| **Advanced ML Models (v4.2)** | Planning | Q2 2026 |

### Recent Releases

| Version | Date | Key Features |
|---------|------|-------------|
| **v4.1.0** | Nov 2025 | **Advanced Comparative Analysis (v2.0)**: Complete overhaul with compositional data normalization (CLR, rarefaction), multi-group stats (Kruskal-Wallis), effect sizes (Cliff's Delta, R²), statistical power analysis, and enhanced ML biomarker discovery (hyperparameter tuning). **Pathogen DB v4.1**: Perfected 100% specificity by filtering housekeeping genes. **Code Refinement**: Standardized output formatting and bug fixes. |
| **v4.0.0** | Oct 2025 | Enhanced COG+SwissProt annotation with comprehensive reporting, new and updated ML model v2.1, advanced gene prediction controls, professional logging system with debug mode. |
| **v3.6.0** | Oct 2025 | SPAdes metagenomic assembly, enhanced COG+SwissProt annotation, modular visualization system with specialized plotters. |
| **v3.5.0** | Aug 2025 | Enhanced statistical testing, advanced ML biomarker discovery, comprehensive comparative visualizations. |
| **v3.3.0** | Aug 2025 | Comparative analysis pipeline, statistical testing, beta diversity visualization, alpha diversity visualization. |

-----

## What's New in Version 4.1.0

### Major Improvements

**1. Advanced Comparative Analysis (v2.0)**

  - **Compositional Data Handling**: Full support for **TSS** (relative abundance), **CLR** (Centered Log-Ratio), and **rarefaction** normalization.
  - **Statistical Power Assessment**: Automatically analyzes your cohort size and provides recommendations and estimates of detectable effect sizes.
  - **Multi-Group Statistics**: Upgraded alpha diversity to use **Kruskal-Wallis** for \>2 groups, with automated pairwise post-hoc tests.
  - **Effect Size Metrics**: Integrated robust effect size calculations (**Cliff's Delta** for differential abundance, **epsilon-squared** for alpha diversity, **R²** for PERMANOVA) for measuring the magnitude of findings.
  - **Enhanced ML Biomarkers**: The biomarker discovery pipeline now uses **SelectKBest** for feature selection and **GridSearchCV** for hyperparameter tuning, resulting in more accurate and reliable models.
  - **Prevalence Filtering**: Smarter differential abundance testing that filters out extremely rare species for more meaningful results.

**2. Pathogen Database v4.1 (Critical Update)**

  - **Perfected Specificity**: The database build process now actively filters common housekeeping genes (*rpoB*, *gyrA*, etc.), eliminating major sources of false positives and achieving **100% specificity**.
  - **Benchmark-Proven**: Achieves **100% sensitivity** and **100% specificity** on our critical marker test suite (see Benchmark section).
  - **Comprehensive Coverage**: Retains 100% coverage for all critical markers, including **mecA (MRSA)**, **mcr-1 (Colistin)**, and all major carbapenemases (KPC, NDM, OXA).

**3. Standardized Output & Logging**

  - The professional `OutputFormatter` is now fully integrated into the **reporting** and **visualization** modules, ensuring all console output is consistent, clean, and professional.
  - Better progress tracking for comparative analysis and visualization steps.

**4. Code Refinement & Bug Fixes**

  - General code cleanup and performance optimizations.
  - Fixed bugs in command-line argument parsing (e.g., removed non-existent `ident` arg).

-----

## Features

### Core Analysis Capabilities

**Advanced File Validation**

  - FASTQ quality control with quality score analysis and contamination detection.
  - FASTA format validation with sequence type detection and composition analysis.
  - Comprehensive statistics including MD5 checksums and N50 metrics.

**Enhanced Functional Annotation**

  - Dual-database annotation system combining COG and SwissProt.
  - COG database: Comprehensive functional categories and orthologous groups.
  - SwissProt database: High-quality, manually curated protein annotations.
  - Increased annotation coverage and functional depth.
  - Mobile genetic element analysis with IS family classification.

**Advanced Gene Prediction Controls**

  - Prokka annotation with customizable parameters.
  - Contig filtering options (enable/disable, custom length thresholds).
  - tbl2asn timeout management with auto-kill functionality.
  - Configurable threading for optimal performance.

**Professional Logging System**

  - **Standard Mode**: Clean, user-friendly progress output powered by the **OutputFormatter**.
  - **Debug Mode** (`--debug`): Comprehensive diagnostic output.
  - Full command-line tool invocations.

**Machine Learning Integration**

  - Automated feature extraction for ML models.
  - Pre-trained pathogen classification models.
  - Random Forest biomarker discovery with cross-validation.

**Pathogen Screening (v4.1 DB)**

  - **100% Specificity**: New database filters housekeeping genes to prevent false positives.
  - **100% Sensitivity**: Correctly identifies all critical markers like `mecA` and `mcr-1`.
  - Clinical risk assessment with comprehensive risk stratification.

**Advanced Statistical Analysis (v2.0)**

  - **Statistical Power Assessment**: Provides recommendations based on sample size.
  - **Compositional Data Normalization**: Supports **TSS**, **CLR**, and **rarefaction**.
  - **Alpha Diversity**: Multi-group comparisons (Kruskal-Wallis) with post-hoc tests and **Cliff's Delta** effect size.
  - **Beta Diversity**: PERMANOVA (with **R² effect size**) and ANOSIM statistical tests.
  - **Differential Abundance**: Enhanced with prevalence filtering and **Cliff's Delta** effect size.
  - **Machine Learning Biomarkers**: Upgraded with **hyperparameter tuning** (GridSearchCV) and **feature selection** (SelectKBest).

**Enhanced Reporting System**

  - **Powered by OutputFormatter** for consistent, professional console output.
  - **Taxonomic Reports**: Clinical and research perspectives with diversity metrics.
  - **Functional Reports**: Detailed COG category analysis, mobile element tracking.
  - **Pathogen Risk Reports**: Three-tier assessment with integrated scoring.

**Modular Visualization System**

  - **Powered by OutputFormatter** for standardized logging during plot generation.
  - Completely redesigned with modular architecture.
  - Modern, publication-ready visualizations.

### Analysis Workflows

**FASTQ Analysis (Clinical Focus)**

  - Quality control and preprocessing.
  - Rapid Kraken2/Bracken classification.
  - **High-specificity** pathogen screening with clinical recommendations.
  - Enhanced functional annotation with COG+SwissProt.
  - Detailed text reports for all analysis components.

**FASTA Analysis (Research Focus)**

  - High-accuracy BLAST classification.
  - Machine learning pathogen prediction.
  - Comprehensive taxonomic profiling.
  - Detailed functional annotation with dual databases.
  - Research-oriented detailed reports.

**Comparative Analysis (v2.0)**

  - Statistical comparison across sample groups with **power analysis**.
  - Beta diversity analysis (PCoA) with **R² effect size**.
  - Alpha diversity comparison with **multi-group stats** and **effect sizes**.
  - Differential abundance testing with volcano plots.
  - Interactive taxonomic heatmaps with enhanced styling.
  - Machine learning biomarker discovery with **hyperparameter tuning**.
  - Publication-ready figure generation.

-----

## System Architecture

```
src/metaquest/
├── cli.py                          # Enhanced CLI with annotation controls
├── config.py                       # Central configuration management
│
├── core/                           # Core analysis modules
│   ├── analysis.py                 # Main analysis orchestration
│   ├── taxonomic_analysis.py       # Taxonomic classification
│   ├── pathogen_analysis.py        # Pathogen detection
│   ├── comparative_analysis.py     # v2.0: Advanced comparison (normalization, power, multi-group stats)
│   └── functional_analysis.py      # Enhanced functional annotation (COG+SwissProt)
│
├── io/                             # Input/Output and validation
│   ├── file_validator.py           # File validation & QC
│   ├── output_formatter.py         # Professional logging & output system
│   └── utils.py                    # I/O utilities
│
├── ml/                             # Machine Learning components
│   ├── feature_extractor.py        # Feature extraction for ML
│   ├── pathogen_predictor.py       # ML pathogen prediction
│   └── model_artifacts/            # Pre-trained models
│
├── reporting/                      # Enhanced reporting engine (uses OutputFormatter)
│   ├── main_reporter.py            # Report orchestrator
│   ├── base_reporter.py            # Base class with formatting utilities
│   ├── taxonomy_reporter.py        # Enhanced taxonomic reports
│   ├── pathogen_risk_reporter.py   # Pathogen risk reports
│   ├── functional_reporter.py      # Expanded functional reports
│   └── risk_scoring.py             # Risk calculation engine
│
└── visualization/                  # Modular visualization system (uses OutputFormatter)
    ├── base_visualizer.py          # Base class for all visualizers
    ├── taxonomic_visualizer.py     # Taxonomic visualization module
    ├── pathogenic_visualizer.py    # Pathogen detection plots
    ├── functional_visualizer.py    # Functional annotation visualizations
    ├── compare_visuals.py          # Multi-sample comparison plots
    ├── dashboard.py                # Interactive dashboard orchestrator
    └── main_visualizer.py          # Main visualization coordinator
```

-----

-----

## Benchmark Performance

MetaQuest has been rigorously benchmarked against established tools using standardized datasets. The results demonstrate superior performance across all analysis categories, especially with the new **Pathogen DB v4.1**.

### Key Performance Metrics

#### Taxonomic Classification Excellence

  - **3.1× lower error rate** than MetaPhlAn4 (1.69% vs 5.23% MAE).
  - **100% species detection rate** for all expected organisms.
  - **Perfect consistency**: All predictions within ±3.5% of expected values.
  - **Statistically indistinguishable** from gold standard (p = 0.920).

#### Functional Annotation Leadership

  - **98.9% annotation coverage** using COG database.
  - **99.3% gene detection accuracy** compared to NCBI reference.
  - **Near-reference quality** with only 1.1% annotation gap.
  - **18-29% higher coverage** than single-database approaches.

#### Pathogen Detection: Perfected Accuracy & Specificity (v4.1 DB)

The new database (v4.1) was rebuilt to filter housekeeping genes and focus on curated markers. This fixed all previous sensitivity and specificity issues.

| Metric | Old Database | **New Database (v4.1)** | Improvement |
|--------|--------------|-------------------------|-------------|
| **Sensitivity** | 33.3% | **100.0%** | **+200%** |
| **Specificity** | 33.3% | **100.0%** | **+200%** |
| **Coverage (Colistin Res.)** | 0.0% | **100.0%** | **+100%** |
| **Coverage (ESKAPE)** | 87.5% | **100.0%** | **+12.5%** |
| **Database Size** | 70,752 seqs | **21,152 seqs** | **-70.1%** |
| **Performance** | 1.5 q/s | **3.9 q/s** | **2.6× faster** |

### Benchmarking Standards

MetaQuest was validated using industry-standard datasets:

  - **ZymoBIOMICS Microbial Community Standard** (taxonomic classification).
  - **E. coli K-12 MG1655 Reference Genome** (functional annotation).
  - **Critical Marker Test Suite** (mecA, mcr-1, KPC-2, etc.) for pathogen detection.

### Detailed Results

For comprehensive benchmarking methodology, detailed performance metrics, and statistical analyses, see the complete **[Benchmarking Report](https://www.google.com/search?q=benchmark/benchmark.md)**.

-----

## Installation

MetaQuest requires a Linux/macOS environment with conda package manager.

**For detailed installation instructions, see [installation.md](https://www.google.com/search?q=installation.md)**

### System Requirements

  - Linux/macOS operating system
  - Conda package manager
  - Minimum 16GB RAM (32GB recommended)
  - 100GB available disk space for databases

### Quick Installation

```bash
# Clone repository
git clone https://github.com/Karudhoru/MetaQuest.git
cd MetaQuest

# Create conda environment
conda env create -f environment.yml
conda activate metaquest

# Install MetaQuest
pip install -e .

# Download annotation databases (COG + SwissProt)
./scripts/setup_databases.sh --swissprot --kraken

# Download and build the new pathogen database (v4.1)
# (This script downloads sequences and builds the DB)
python scripts/custom_pathogen_db.py 
# (This script formats the DB with DIAMOND)
./scripts/setup_databases.sh --pathogen

# Verify installation
metaquest check
```

-----

## Quick Start

### 1\. Verify Installation

```bash
metaquest check
```

### 2\. Validate Input Files

```bash
# FASTQ validation
metaquest validate fastq --single sample.fastq.gz

# FASTA validation  
metaquest validate fasta genome.fasta
```

### 3\. Run Analysis

#### Basic Analysis

```bash
# FASTQ analysis (uses default annotation settings)
metaquest analyze fastq --single sample.fastq.gz -o results/

# FASTA analysis with enhanced annotation
metaquest analyze fasta genome.fasta -o results/ -s 100
```

#### Advanced Annotation Controls

```bash
# Custom contig filtering threshold
metaquest analyze fastq --single reads.fq --min-contig-length 500 -o results/

# Disable contig filtering (annotate all contigs)
metaquest analyze fastq --single reads.fq --no-filter-contigs -o results/

# Skip annotation for faster taxonomic-only analysis
metaquest analyze fastq --single sample.fastq.gz --skip-annotation -o fast_results/
```

#### Debug Mode for Troubleshooting

```bash
# Run with full diagnostic output
metaquest --debug analyze fastq --single reads.fq -o debug_results/
```

### 4\. Compare Multiple Samples (v2.0)

```bash
# Create metadata file (metadata.tsv)
# sample_id    group
# sample1      Healthy
# sample2      Diseased
# sample3      Healthy
# sample4      Diseased

# Run comparative analysis (default: TSS normalization)
metaquest compare -i sample1_results/ sample2_results/ sample3_results/ sample4_results/ -m metadata.tsv -o comparison/

# Run with CLR normalization for compositional data
metaquest compare -i sample*_results/ -m metadata.tsv -o comparison_clr/ --normalization clr
```

### 5\. View Results

The analysis generates comprehensive outputs:

**Enhanced Text Reports**

  - **Taxonomic Report**: Clinical summary, researcher view, diversity metrics.
  - **Functional Report**: COG categories, mobile elements, annotation quality.
  - **Pathogen Risk Report**: Three-tier assessment with clinical interpretation.

**Functional Annotation** (COG + SwissProt)

  - Expanded functional category coverage.
  - High-confidence protein annotations.
  - Detailed COG category distributions.
  - Mobile genetic element analysis.

**Interactive Visualizations**

  - Modern, publication-ready dashboards.
  - Enhanced color schemes and layouts.

-----

## Documentation

### Core Documentation

  - **[Usage Guide](usage.md)** - Comprehensive usage examples and command reference
  - **[Installation Guide](https://www.google.com/search?q=installation.md)** - Detailed installation instructions
  - **[Annotation Guide](https://www.google.com/search?q=annotation.md)** - COG and SwissProt database information
  - **API Documentation** - Coming soon

### Analysis Outputs

**FASTQ Results**

  - Interactive dashboard with modern styling.
  - Pathogen risk assessment with clinical recommendations.
  - Taxonomic classification reports (text + JSON).
  - Enhanced functional annotation reports (COG + SwissProt).
  - Mobile genetic element analysis.
  - Publication-ready visualizations.

**FASTA Results**

  - BLAST+ML integrated pathogen reports.
  - High-accuracy taxonomic classification.
  - Machine learning predictions.
  - Comprehensive functional annotation with dual databases.
  - Detailed text reports for all components.

**Comparative Results (v2.0)**

  - **Statistical power assessment report**.
  - Alpha diversity box plots with **multi-group stats (Kruskal-Wallis)**.
  - Beta diversity PCoA plots with **PERMANOVA (R²) / ANOSIM results**.
  - Differential abundance volcano plots with **Cliff's Delta effect size**.
  - Interactive taxonomic heatmaps with modern color schemes.
  - **Enhanced ML feature importance reports** (from tuned models).
  - Publication-ready figure exports.

-----

## Contributing

We welcome contributions from the scientific community.

### Priority Areas

| Area | Description | Impact |
|------|-------------|---------|
| **Pathogen DB Benchmarking** | Publish comprehensive benchmarks for the new v4.1 DB | High |
| **Advanced ML Models (v4.2)** | Implement new models beyond Random Forest | High |
| **Database Expansion** | Add additional functional databases (KEGG, Pfam) | High |
| **Clinical Validation** | Real-world testing of pathogen detection | Critical |
| **Statistical Methods** | Additional statistical tests and corrections | High |
| **Visualization Features** | New plot types and interactive elements | Medium |
| **Documentation** | User guides and clinical interpretation | Medium |

### How to Contribute

1.  Fork the repository
2.  Create a feature branch (`git checkout -b feature/new-analysis`)
3.  Make your changes with appropriate tests
4.  Submit a pull request with detailed description
5.  Participate in code review process

-----

## Support

### Documentation Resources

  - **Installation Guide**: [installation.md](https://www.google.com/search?q=installation.md)
  - **Usage Examples**: [usage.md](usage.md)
  - **Annotation Guide**: [annotation.md](https://www.google.com/search?q=annotation.md)
  - **GitHub Wiki**: Additional tutorials and examples

### Getting Help

  - **Bug Reports**: Submit via GitHub issues with detailed reproduction steps
  - **Feature Requests**: Use GitHub discussions for new feature proposals
  - **Email Support**: metaquest-support@example.org

### Troubleshooting

  - Use `--debug` flag for detailed diagnostic output.
  - Check `metaquest.log` in output directory for error details.
  - Run `metaquest check` to verify system dependencies.
  - See [usage.md](usage.md) for common issues and solutions.

-----

## Citation

*Citation information will be provided upon publication*

-----

## License

*License information will be added upon publication*

-----

**Last Updated**: November 2025 - Version 4.1.0 with Advanced Comparative Analysis v2.0 and Pathogen DB v4.1

-----

*MetaQuest is actively developed with regular updates and improvements. Check our GitHub repository for the latest releases and features.*