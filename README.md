# MetaQuest

**A Comprehensive Metagenomics Analysis Pipeline**

---

## Overview

MetaQuest is an integrated bioinformatics pipeline that addresses the complex challenges of metagenomic data analysis. By combining state-of-the-art tools and databases, MetaQuest provides researchers with a streamlined workflow for understanding microbial community composition, identifying potential pathogens, and assessing antimicrobial resistance profiles from sequencing data.

**Key Capabilities:**
- Taxonomic classification and pathogen detection
- Machine learning-enhanced pathogen prediction
- Advanced comparative analysis with statistical testing
- Comprehensive file validation and quality control
- Interactive visualizations and reporting

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

**Current Version:** 3.5.0 (August 2025)

MetaQuest is under active development with core functionality operational and major improvements completed across all processing capabilities.

### Completed Features

| Component | Status | Description |
|-----------|---------|-------------|
| **File Validation & Quality Control** | Complete | Comprehensive validation with detailed statistics for both FASTQ and FASTA inputs |
| **Taxonomic Classification** | Complete | Both FASTQ and FASTA taxonomic profiling are accurate and fully functional |
| **Pathogen Detection** | Complete | Comprehensive pathogen detection workflows with clinical risk assessment |
| **Machine Learning Integration** | Complete | Advanced ML-based pathogen prediction with feature extraction and model artifacts |
| **Comparative Analysis** | Complete | Statistical comparison across multiple samples with differential abundance testing and visualization |
| **Beta Diversity Analysis** | Complete | PERMANOVA and ANOSIM statistical tests with PCoA visualization |
| **Alpha Diversity Analysis** | Complete | Shannon, Simpson, Chao1, and Observed Species metrics with statistical testing |
| **Differential Abundance Testing** | Complete | Mann-Whitney U tests with multiple testing correction (FDR, Bonferroni) |
| **Machine Learning Biomarker Discovery** | Complete | Random Forest classification with feature importance analysis |

### In Development

| Component | Status | Target Completion |
|-----------|---------|-------------------|
| **Virulence Factor Analysis** | In Progress | Q3 2025 |
| **AMR Analysis** | In Progress | Q4 2025 |

### Recent Releases

| Version | Date | Key Features |
|---------|------|-------------|
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

**Machine Learning Integration**
- Automated feature extraction for ML models
- Pre-trained pathogen classification models
- Complete model artifacts with scalers and feature selectors
- Random Forest biomarker discovery with cross-validation

**Pathogen Screening**
- FASTQ processing with multi-source screening and Bracken integration
- FASTA processing with integrated BLAST taxonomy and ML predictions
- Clinical risk assessment with comprehensive risk stratification

**Advanced Statistical Analysis**
- Alpha diversity metrics (Shannon, Simpson, Chao1, Observed Species)
- Beta diversity analysis with Bray-Curtis dissimilarity
- PERMANOVA and ANOSIM statistical tests
- Differential abundance testing with multiple correction methods
- Machine learning-based biomarker identification

**Analysis & Reporting**
- Species-level taxonomic profiling with diversity metrics
- Interactive HTML dashboards with dynamic visualizations
- Statistical analysis and quality assessment
- Functional annotation with gene prediction
- Multi-sample comparative analysis with publication-ready outputs

### Analysis Workflows

**FASTQ Analysis (Clinical Focus)**
- Rapid Kraken2/Bracken classification
- Pathogen screening with clinical recommendations
- Quality control and contamination detection

**FASTA Analysis (Research Focus)**
- High-accuracy BLAST classification
- Machine learning pathogen prediction
- Comprehensive taxonomic profiling

**Comparative Analysis**
- Statistical comparison across sample groups
- Beta diversity analysis with PCoA visualization
- Alpha diversity comparison with box plots
- Differential abundance testing with volcano plots
- Interactive heatmaps and abundance bar plots
- Machine learning biomarker discovery

---

## System Architecture

```
src/metaquest/
├── cli.py                          # Command-line interface
├── config.py                       # Configuration management
├── core/                           # Core analysis modules
│   ├── analysis.py                 # Main analysis orchestration
│   ├── taxonomic_analysis.py       # Taxonomic classification
│   ├── pathogen_analysis.py        # Pathogen detection
│   ├── comparative_analysis.py     # Multi-sample comparison
│   └── functional_analysis.py      # Functional annotation
├── io/                             # Input/Output handling
│   ├── file_validator.py           # File validation & QC
│   └── utils.py                    # I/O utilities
├── ml/                             # Machine Learning components
│   ├── feature_extractor.py        # Feature extraction for ML
│   ├── pathogen_predictor.py       # ML pathogen prediction
│   └── model_artifacts/            # Pre-trained models
├── reporting/                      # Report generation
│   ├── reporting.py               # Comprehensive reporting
│   └── dashboard.py               # Interactive dashboards
└── visualization/                  # Data visualization
    ├── visualization.py           # Interactive visualizations
    └── compare_visuals.py         # Comparative analysis plots
```

---

## Installation

MetaQuest requires a Linux/macOS environment with conda package manager.

**For detailed installation instructions, see [installation.md](installation.md)**

### System Requirements

- Linux/macOS operating system
- Conda package manager
- Minimum 8GB RAM (16GB recommended)
- 50GB available disk space for databases

### Quick Installation

```bash
# Clone repository
git clone https://github.com/your-org/metaquest.git
cd metaquest

# Create conda environment
conda env create -f environment.yml
conda activate metaquest

# Install MetaQuest
pip install -e .

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
```bash
# FASTQ analysis (clinical focus)
metaquest analyze fastq --single sample.fastq.gz -o results/

# FASTA analysis (research focus)
metaquest analyze fasta genome.fasta -o results/ -s 100

# Skip annotation for faster analysis
metaquest analyze fastq --single sample.fastq.gz --skip-annotation -o fast_results/
```

### 4. Compare Multiple Samples
```bash
# Create metadata file (metadata.tsv)
# sample_id    group
# sample1      Healthy
# sample2      Diseased

# Run comparative analysis
metaquest compare -i sample1_results/ sample2_results/ -m metadata.tsv -o comparison/
```

### 5. Statistical Testing Results

The comparative analysis provides comprehensive statistical testing:

**Alpha Diversity Testing**
- Shannon diversity: significant differences (p < 0.01)
- Simpson diversity: significant differences (p < 0.01)  
- Chao1 richness: significant differences (p < 0.05)
- Observed species: significant differences (p < 0.05)

**Beta Diversity Testing**
- PERMANOVA: F = 5.418, p = 0.004 (significant)
- ANOSIM: R = 0.496, p = 0.018 (overlapping but clearly different)

**Differential Abundance**
- 150 nominally significant species (p < 0.05)
- 241 trend-level significant species (p < 0.10)
- Multiple testing correction (FDR and Bonferroni)

**Machine Learning Biomarker Discovery**
- Random Forest cross-validation accuracy: 83% ± 12%
- Feature importance ranking for discriminative species

---

## Documentation

### Core Documentation
- **[Usage Guide](usage.md)** - Comprehensive usage examples and command reference
- **[Installation Guide](installation.md)** - Detailed installation instructions
- **API Documentation** - Coming soon

### Analysis Outputs

**FASTQ Results**
- Interactive dashboard with pathogen risk assessment
- Taxonomic classification reports
- Pathogen detection summaries
- Functional annotation results

**FASTA Results**  
- BLAST+ML integrated pathogen reports
- High-accuracy taxonomic classification
- Machine learning predictions
- Organism comparison data

**Comparative Results**
- Alpha diversity box plots with statistical tests
- Beta diversity PCoA plots with PERMANOVA/ANOSIM results
- Differential abundance volcano plots
- Interactive taxonomic heatmaps
- Abundance bar plots by group
- Machine learning feature importance reports
- Statistical comparison summaries

---

## Contributing

We welcome contributions from the scientific community.

### Priority Areas

| Area | Description | Impact |
|------|-------------|---------|
| **Machine Learning Enhancement** | Improve ML models with additional training data | High |
| **Clinical Validation** | Real-world testing of pathogen detection | Critical |
| **Statistical Methods** | Additional statistical tests and corrections | High |
| **Virulence Factor Analysis** | Optimize VF detection workflows | High |
| **AMR Enhancement** | Improve antimicrobial resistance detection | High |
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

---

## Support

### Documentation Resources
- **Installation Guide**: [installation.md](installation.md)
- **Usage Examples**: [usage.md](usage.md)
- **GitHub Wiki**: Additional tutorials and examples

### Getting Help
- **Bug Reports**: Submit via GitHub issues with detailed reproduction steps
- **Feature Requests**: Use GitHub discussions for new feature proposals
- **Community Support**: Join our discussion forum for user questions

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

**Last Updated**: August 2025 - Version 3.5.0 with enhanced statistical testing and ML biomarker discovery

---

*MetaQuest is actively developed with regular updates and improvements. Check our GitHub repository for the latest releases and features.*