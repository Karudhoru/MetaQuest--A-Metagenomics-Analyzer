# MetaQuest

**A Comprehensive Metagenomics Analysis Pipeline**

---

## Overview

MetaQuest is an integrated bioinformatics pipeline that addresses the complex challenges of metagenomic data analysis. By combining state-of-the-art tools and databases, MetaQuest provides researchers with a streamlined workflow for understanding microbial community composition, identifying potential pathogens, and assessing antimicrobial resistance profiles from sequencing data.

**Key Capabilities:**
- Taxonomic classification and pathogen detection
- Machine learning-enhanced pathogen prediction (Also in FASTQ)
- Advanced comparative analysis with statistical testing
- Comprehensive file validation and quality control
- Interactive visualizations and reporting
- Enhanced functional annotation with COG and SwissProt databases
- Professional logging system with debug mode for troubleshooting

---

## Table of Contents

- [Development Status](#development-status)
- [Features](#features)
- [System Architecture](#system-architecture)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Documentation](#documentation)
- [Contributing](#contributing)
- [Support](#support)

---

## Development Status

**Current Version:** 4.0.0 (October 2025)

MetaQuest continues active development with major improvements to annotation, reporting, and user experience.

### Completed Features

| Component | Status | Description |
|-----------|---------|-------------|
| **File Validation & Quality Control** | ✅ Complete | Comprehensive validation with detailed statistics for both FASTQ and FASTA inputs |
| **Taxonomic Classification** | ✅ Complete | Both FASTQ and FASTA taxonomic profiling are accurate and fully functional |
| **Pathogen Detection** | ✅ Complete | Comprehensive pathogen detection workflows with clinical risk assessment |
| **Machine Learning Integration** | ✅ Complete | Advanced ML-based pathogen prediction with feature extraction and model artifacts |
| **Comparative Analysis** | ✅ Complete | Statistical comparison across multiple samples with differential abundance testing and visualization |
| **Beta Diversity Analysis** | ✅ Complete | PERMANOVA and ANOSIM statistical tests with PCoA visualization |
| **Alpha Diversity Analysis** | ✅ Complete | Shannon, Simpson, Chao1, and Observed Species metrics with statistical testing |
| **Differential Abundance Testing** | ✅ Complete | Mann-Whitney U tests with multiple testing correction (FDR, Bonferroni) |
| **Machine Learning Biomarker Discovery** | ✅ Complete | Random Forest classification with feature importance analysis |
| **Enhanced Functional Annotation** | ✅ Complete | Dual-database annotation system combining COG and SwissProt with detailed reporting |
| **Advanced Gene Prediction** | ✅ Complete | Prokka annotation with customizable contig filtering and timeout controls |
| **Professional Logging System** | ✅ Complete | Dual-mode output (standard/debug) with formatted progress tracking |
| **Enhanced Reporting System** | ✅ Complete | Comprehensive text-based reports with clinical and research perspectives |

### In Development

| Component | Status | Target Completion |
|-----------|---------|-------------------|
| **Metabolic Pathway Reconstruction** | Planning | Q3 2026 |

### Recent Releases

| Version | Date | Key Features |
|---------|------|-------------|
| **v4.0.0** | Oct 2025 | Enhanced COG+SwissProt annotation with comprehensive reporting, new and updated ML model v2.1, advanced gene prediction controls, professional logging system with debug mode |
| **v3.6.0** | Oct 2025 | SPAdes metagenomic assembly, enhanced COG+SwissProt annotation, modular visualization system with specialized plotters |
| **v3.5.0** | Aug 2025 | Enhanced statistical testing, advanced ML biomarker discovery, comprehensive comparative visualizations |
| **v3.3.0** | Aug 2025 | Comparative analysis pipeline, statistical testing, beta diversity visualization, alpha diversity visualization |
| **v3.2.2** | Aug 2025 | Professional CLI revamp, advanced reporting with scikit-bio integration |
| **v3.2.1** | Jul 2025 | Complete ML pipeline, advanced file validation system |
| **v3.2.0** | Jun 2025 | Pathogen detection & clinical risk assessment |

---

## Features

### Core Analysis Capabilities

**Advanced File Validation**
- FASTQ quality control with quality score analysis and contamination detection
- FASTA format validation with sequence type detection and composition analysis
- Comprehensive statistics including MD5 checksums and N50 metrics

**Enhanced Functional Annotation**
- Dual-database annotation system combining COG and SwissProt
- COG database: Comprehensive functional categories and orthologous groups
- SwissProt database: High-quality, manually curated protein annotations
- Increased annotation coverage and functional depth
- Detailed functional category distribution and pathway insights
- Mobile genetic element analysis with IS family classification

**Advanced Gene Prediction Controls**
- Prokka annotation with customizable parameters
- Contig filtering options (enable/disable, custom length thresholds)
- tbl2asn timeout management with auto-kill functionality
- Configurable threading for optimal performance
- Enhanced error handling and recovery

**Professional Logging System**
- **Standard Mode**: Clean, user-friendly progress output
  - Progress indicators with spinners
  - Formatted section headers
  - Success/error/warning notifications
  - Time-formatted completion messages
- **Debug Mode** (`--debug`): Comprehensive diagnostic output
  - Full command-line tool invocations
  - Complete tool output streams
  - Detailed error traces
  - Performance metrics

**Machine Learning Integration**
- Automated feature extraction for ML models
- Pre-trained pathogen classification models
- Complete model artifacts with scalers and feature selectors
- Random Forest biomarker discovery with cross-validation

**Pathogen Screening**
- FASTQ processing with multi-source screening and Bracken integration
- FASTA processing with integrated BLAST taxonomy and ML predictions
- Clinical risk assessment with comprehensive risk stratification
- Three-tier pathogenicity assessment system

**Advanced Statistical Analysis**
- Alpha diversity metrics (Shannon, Simpson, Chao1, Observed Species)
- Beta diversity analysis with Bray-Curtis dissimilarity
- PERMANOVA and ANOSIM statistical tests
- Differential abundance testing with multiple correction methods
- Machine learning-based biomarker identification

**Enhanced Reporting System**
- **Taxonomic Reports**: Clinical and research perspectives with diversity metrics
- **Functional Reports**: Detailed COG category analysis, mobile element tracking, annotation quality assessment
- **Pathogen Risk Reports**: Three-tier risk assessment with integrated scoring and clinical interpretation
- Professional formatting with emojis and visual indicators
- Structured sections for different user roles (clinician vs researcher)

**Modular Visualization System**
- Completely redesigned with modular architecture
- Specialized visualization modules for different analysis types
- Modern, publication-ready visualizations
- Enhanced color schemes and data presentation
- Dynamic plots with improved interactivity
- Comprehensive multi-panel layouts
- Professional styling and responsive design

### Analysis Workflows

**FASTQ Analysis (Clinical Focus)**
- Quality control and preprocessing
- Rapid Kraken2/Bracken classification
- Pathogen screening with clinical recommendations
- Enhanced functional annotation with COG+SwissProt
- Detailed text reports for all analysis components

**FASTA Analysis (Research Focus)**
- High-accuracy BLAST classification
- Machine learning pathogen prediction
- Comprehensive taxonomic profiling
- Detailed functional annotation with dual databases
- Research-oriented detailed reports

**Comparative Analysis**
- Statistical comparison across sample groups
- Beta diversity analysis with PCoA visualization
- Alpha diversity comparison with box plots
- Differential abundance testing with volcano plots
- Interactive taxonomic heatmaps with enhanced styling
- Abundance bar plots by group with modern aesthetics
- Machine learning biomarker discovery
- Publication-ready figure generation

---

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
│   ├── comparative_analysis.py     # Multi-sample comparison
│   └── functional_analysis.py      # Enhanced functional annotation (COG+SwissProt)
│
├── io/                             # Input/Output and validation
│   ├── file_validator.py           # File validation & QC
│   ├── output_formatter.py         # Professional logging system (NEW)
│   └── utils.py                    # I/O utilities
│
├── ml/                             # Machine Learning components
│   ├── feature_extractor.py        # Feature extraction for ML
│   ├── pathogen_predictor.py       # ML pathogen prediction
│   └── model_artifacts/            # Pre-trained models
│
├── reporting/                      # Enhanced reporting engine
│   ├── main_reporter.py            # Report orchestrator
│   ├── base_reporter.py            # Base class with formatting utilities
│   ├── taxonomy_reporter.py        # Enhanced taxonomic reports (UPDATED)
│   ├── pathogen_risk_reporter.py   # Pathogen risk reports (UPDATED)
│   ├── functional_reporter.py      # Expanded functional reports (UPDATED)
│   └── risk_scoring.py             # Risk calculation engine
│
└── visualization/                  # Modular visualization system
    ├── base_visualizer.py          # Base class for all visualizers
    ├── taxonomic_visualizer.py     # Taxonomic visualization module
    ├── pathogenic_visualizer.py    # Pathogen detection plots
    ├── functional_visualizer.py    # Functional annotation visualizations
    ├── compare_visuals.py          # Multi-sample comparison plots
    ├── dashboard.py                # Interactive dashboard orchestrator
    └── main_visualizer.py          # Main visualization coordinator
```

---

## Installation

MetaQuest requires a Linux/macOS environment with conda package manager.

**For detailed installation instructions, see [installation.md](installation.md)**

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
./scripts/setup_databases.sh

# Download pathogen database (Custom Database)
python scripts/custom_pathogen_db.py

# Verify installation
metaquest check
```

---

## Quick Start

### 1. Verify Installation
```bash
metaquest check
```

### 2. Validate Input Files
```bash
# FASTQ validation
metaquest validate fastq --single sample.fastq.gz

# FASTA validation  
metaquest validate fasta genome.fasta
```

### 3. Run Analysis

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

# Extended tbl2asn timeout for large datasets
metaquest analyze fasta assembly.fasta --tbl2asn-timeout 600 -o results/

# Use more threads for faster annotation
metaquest analyze fastq --paired R1.fq R2.fq --annotation-threads 16 -o results/

# Disable tbl2asn auto-kill (let it run indefinitely)
metaquest analyze fasta assembly.fasta --no-kill-tbl2asn -o results/

# Skip annotation for faster taxonomic-only analysis
metaquest analyze fastq --single sample.fastq.gz --skip-annotation -o fast_results/
```

#### Debug Mode for Troubleshooting
```bash
# Run with full diagnostic output
metaquest --debug analyze fastq --single reads.fq -o debug_results/

# Debug mode shows:
# - All tool command invocations
# - Complete tool output
# - Detailed error messages
# - Performance timing
```

### 4. Compare Multiple Samples
```bash
# Create metadata file (metadata.tsv)
# sample_id    group
# sample1      Healthy
# sample2      Diseased

# Run comparative analysis with enhanced visualizations
metaquest compare -i sample1_results/ sample2_results/ -m metadata.tsv -o comparison/
```

### 5. View Results

The analysis generates comprehensive outputs:

**Enhanced Text Reports**
- **Taxonomic Report**: Clinical summary, researcher view, diversity metrics
- **Functional Report**: COG categories, mobile elements, annotation quality
- **Pathogen Risk Report**: Three-tier assessment with clinical interpretation

**Functional Annotation** (COG + SwissProt)
- Expanded functional category coverage
- High-confidence protein annotations
- Pathway and process enrichment
- Detailed COG category distributions
- Mobile genetic element analysis
- IS family classification

**Interactive Visualizations**
- Modern, publication-ready dashboards
- Enhanced color schemes and layouts
- Improved data presentation
- Responsive design for all screen sizes

---

## Documentation

### Core Documentation
- **[Usage Guide](usage.md)** - Comprehensive usage examples and command reference
- **[Installation Guide](installation.md)** - Detailed installation instructions
- **[Annotation Guide](annotation.md)** - COG and SwissProt database information
- **API Documentation** - Coming soon

### Analysis Outputs

**FASTQ Results**
- Interactive dashboard with modern styling
- Pathogen risk assessment with clinical recommendations
- Taxonomic classification reports (text + JSON)
- Enhanced functional annotation reports (COG + SwissProt)
- Mobile genetic element analysis
- Publication-ready visualizations

**FASTA Results**  
- BLAST+ML integrated pathogen reports
- High-accuracy taxonomic classification
- Machine learning predictions
- Comprehensive functional annotation with dual databases
- Organism comparison data
- Detailed text reports for all components

**Comparative Results**
- Alpha diversity box plots with statistical tests
- Beta diversity PCoA plots with PERMANOVA/ANOSIM results
- Differential abundance volcano plots with enhanced styling
- Interactive taxonomic heatmaps with modern color schemes
- Abundance bar plots by group with improved aesthetics
- Machine learning feature importance reports
- Statistical comparison summaries
- Publication-ready figure exports

---

## What's New in Version 3.6.0

### Major Improvements

**1. Enhanced Annotation System**
- **Dual-database approach**: Combined COG and SwissProt annotation for maximum coverage
- **COG Database**: Comprehensive functional categories with detailed pathway analysis
- **SwissProt Database**: High-quality, manually curated protein annotations
- **Significantly increased annotation coverage** compared to previous versions
- **Mobile genetic element tracking**: IS family classification and transposase analysis

**2. Advanced Gene Prediction Controls**
- **Contig filtering**: Enable/disable filtering with custom length thresholds
  - Default: Filter contigs <1000bp before annotation
  - Customize: `--min-contig-length 500` for different thresholds
  - Disable: `--no-filter-contigs` to annotate all sequences
- **tbl2asn management**: Automatic timeout and termination controls
  - Default 300s timeout with auto-kill for stuck processes
  - Customize: `--tbl2asn-timeout 600` for longer runs
  - Disable: `--no-kill-tbl2asn` for unlimited runtime
- **Threading optimization**: `--annotation-threads` for parallel processing

**3. Professional Logging System**
- **Standard Mode**: Clean, informative progress tracking
  - Formatted section headers with visual separators
  - Progress spinners for long-running operations
  - Color-coded status messages (success/error/warning)
  - Time-formatted completion summaries
- **Debug Mode** (`--debug` flag): Complete diagnostic output
  - Full command-line invocations for all tools
  - Complete stdout/stderr from external programs
  - Detailed error traces and stack information
  - Performance metrics and timing data
- **Structured logging**: Consistent format across all operations
- **Log file support**: Automatic logging to `metaquest.log` in output directory

**4. Enhanced Reporting System**
- **Comprehensive text reports** for all analysis components:
  - **Taxonomic Reports**: Clinical and research perspectives, diversity metrics
  - **Functional Reports**: COG category analysis, mobile element tracking, quality scores
  - **Pathogen Risk Reports**: Three-tier assessment with clinical interpretation
- **Professional formatting**: Emojis, visual indicators, structured sections
- **Role-based views**: Tailored information for clinicians vs researchers
- **Quality metrics**: Annotation coverage, identity scores, functional diversity

**5. Improved Error Handling**
- Better error messages with actionable solutions
- Graceful handling of annotation failures
- Recovery mechanisms for common issues
- Detailed troubleshooting guidance in debug mode

---

## Contributing

We welcome contributions from the scientific community.

### Priority Areas

| Area | Description | Impact |
|------|-------------|---------|
| **Annotation Optimization** | Improve COG/SwissProt coverage and accuracy | High |
| **Database Expansion** | Add additional functional databases (KEGG, Pfam) | High |
| **Machine Learning Enhancement** | Improve ML models with additional training data | High |
| **Clinical Validation** | Real-world testing of pathogen detection | Critical |
| **Statistical Methods** | Additional statistical tests and corrections | High |
| **Visualization Features** | New plot types and interactive elements | Medium |
| **Documentation** | User guides and clinical interpretation | Medium |

### How to Contribute

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/new-analysis`)
3. Make your changes with appropriate tests
4. Submit a pull request with detailed description
5. Participate in code review process

### Development Guidelines

- Follow existing code style and structure
- Include comprehensive tests for new features
- Update documentation for user-facing changes
- Ensure compatibility with existing workflows
- Test with diverse metagenomic datasets

---

## Support

### Documentation Resources
- **Installation Guide**: [installation.md](installation.md)
- **Usage Examples**: [usage.md](usage.md)
- **Annotation Guide**: [annotation.md](annotation.md)
- **GitHub Wiki**: Additional tutorials and examples

### Getting Help
- **Bug Reports**: Submit via GitHub issues with detailed reproduction steps
- **Feature Requests**: Use GitHub discussions for new feature proposals
- **Community Support**: Join our discussion forum for user questions
- **Email Support**: metaquest-support@example.org

### Troubleshooting
- Use `--debug` flag for detailed diagnostic output
- Check `metaquest.log` in output directory for error details
- Run `metaquest check` to verify system dependencies
- See [usage.md](usage.md) for common issues and solutions

### Development Team
**MetaQuest Development Team**  
*Advancing metagenomics through integrated computational solutions*

---

## Citation

*Citation information will be provided upon publication*

---

## License

*License information will be added upon publication*

---

**Last Updated**: October 2025 - Version 3.6.0 with enhanced annotation, advanced gene prediction controls, and professional logging system

---


*MetaQuest is actively developed with regular updates and improvements. Check our GitHub repository for the latest releases and features.*


