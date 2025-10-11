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
- De novo metagenomic assembly with SPAdes
- Enhanced functional annotation with combined COG and SwissProt databases

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

**Current Version:** 3.6.0 (October 2025)

MetaQuest continues active development with major improvements to assembly, annotation, and modular visualization systems.

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
| **Metagenomic Assembly** | Complete | SPAdes-based de novo assembly optimized for metagenomic data |
| **Enhanced Functional Annotation** | Complete | Dual-database annotation system combining COG and SwissProt |
| **Advanced Visualization System** | Complete | Completely revamped reporting, visualization, and interactive dashboards |

### In Development

| Component | Status | Target Completion |
|-----------|---------|-------------------|
| **Virulence Factor Analysis** | In Progress | Q1 2026 |
| **AMR Analysis** | In Progress | Q2 2026 |
| **Metabolic Pathway Reconstruction** | Planning | Q3 2026 |

### Recent Releases

| Version | Date | Key Features |
|---------|------|-------------|
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

**De Novo Metagenomic Assembly**
- SPAdes metagenomic assembler integration for high-quality contig generation
- Optimized assembly parameters for complex microbial communities
- Quality metrics including N50, L50, and assembly statistics
- Support for both single-end and paired-end sequencing data

**Enhanced Functional Annotation**
- Dual-database annotation system combining COG and SwissProt
- COG database: Comprehensive functional categories and orthologous groups
- SwissProt database: High-quality, manually curated protein annotations
- Increased annotation coverage and functional depth
- Detailed functional category distribution and pathway insights

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

**Revamped Modular Visualization System**
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
- SPAdes metagenomic assembly for improved accuracy
- Rapid Kraken2/Bracken classification
- Pathogen screening with clinical recommendations
- Enhanced functional annotation with COG+SwissProt

**FASTA Analysis (Research Focus)**
- High-accuracy BLAST classification
- Machine learning pathogen prediction
- Comprehensive taxonomic profiling
- Detailed functional annotation with dual databases

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
├── cli.py                          # Command-line interface logic
├── config.py                       # Central configuration management
│
├── core/                           # Core analysis modules (the "engine")
│   ├── analysis.py                 # Main analysis orchestration
│   ├── assembly.py                 # SPAdes metagenomic assembly
│   ├── taxonomic_analysis.py       # Taxonomic classification
│   ├── pathogen_analysis.py        # Pathogen detection
│   ├── comparative_analysis.py     # Multi-sample comparison
│   └── functional_analysis.py      # Enhanced functional annotation (COG+SwissProt)
│
├── io/                             # Input/Output and validation
│   ├── file_validator.py           # File validation & QC
│   └── utils.py                    # I/O utilities
│
├── ml/                             # Machine Learning components
│   ├── feature_extractor.py        # Feature extraction for ML
│   ├── pathogen_predictor.py       # ML pathogen prediction
│   └── model_artifacts/            # Pre-trained models
│
├── reporting/                      # Modular reporting engine (REVAMPED)
│   ├── main_reporter.py            # Orchestrator with enhanced layouts
│   ├── base_reporter.py            # Base class with modern styling
│   ├── taxonomic_reporter.py       # Enhanced taxonomic section
│   ├── pathogen_reporter.py        # Improved pathogen section
│   ├── functional_reporter.py      # Expanded functional section (COG+SwissProt)
│   └── assembly_reporter.py        # Assembly statistics and quality metrics
│
└── visualization/                  # Modular visualization system (REVAMPED)
    ├── base_visualizer.py          # Base class for all visualizers
    ├── taxonomic_visualizer.py     # Taxonomic visualization module
    ├── pathogenic_visualizer.py    # Pathogen detection plots
    ├── functional_visualizer.py    # Functional annotation visualizations
    ├── compare_visuals.py          # Multi-sample comparison plots
    └── dashboard.py                # Modern interactive dashboard orchestrator

```

---

## Installation

MetaQuest requires a Linux/macOS environment with conda package manager.

**For detailed installation instructions, see [installation.md](installation.md)**

### System Requirements

- Linux/macOS operating system
- Conda package manager
- Minimum 16GB RAM (32GB recommended for assembly)
- 100GB available disk space for databases and assembly

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

# Download annotation databases (COG + SwissProt)
metaquest download-databases

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
# FASTQ analysis with assembly (recommended)
metaquest analyze fastq --single sample.fastq.gz -o results/ --assemble

# Paired-end FASTQ with assembly
metaquest analyze fastq --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz -o results/ --assemble

# FASTA analysis with enhanced annotation
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

# Run comparative analysis with enhanced visualizations
metaquest compare -i sample1_results/ sample2_results/ -m metadata.tsv -o comparison/
```

### 5. View Results

The analysis generates comprehensive outputs:

**Assembly Metrics** (when using --assemble)
- N50, L50, and total assembly size
- Contig length distribution
- Assembly quality assessment
- Improved taxonomic classification from assembled contigs

**Functional Annotation** (COG + SwissProt)
- Expanded functional category coverage
- High-confidence protein annotations
- Pathway and process enrichment
- Detailed COG category distributions

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
- **[Assembly Guide](assembly.md)** - SPAdes assembly parameters and optimization
- **[Annotation Guide](annotation.md)** - COG and SwissProt database information
- **API Documentation** - Coming soon

### Analysis Outputs

**FASTQ Results**
- Interactive dashboard with modern styling
- Assembly statistics and quality metrics
- Pathogen risk assessment with clinical recommendations
- Taxonomic classification reports
- Enhanced functional annotation (COG + SwissProt)
- Publication-ready visualizations

**FASTA Results**  
- BLAST+ML integrated pathogen reports
- High-accuracy taxonomic classification
- Machine learning predictions
- Comprehensive functional annotation
- Organism comparison data

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

## What's New in Version 3.6

### Major Improvements

**1. SPAdes Metagenomic Assembly**
- Integrated SPAdes assembler for de novo contig generation
- Optimized parameters for complex metagenomic samples
- Improved taxonomic classification accuracy from assembled sequences
- Comprehensive assembly quality metrics

**2. Enhanced Annotation System**
- Dual-database approach combining COG and SwissProt
- Significantly increased annotation coverage
- Higher-quality functional predictions
- Detailed functional category analysis

**3. Modular Visualization Architecture**
- Redesigned with modular architecture matching reporting system
- Specialized plotter modules for each analysis type
- Base plotter class with shared functionality
- Modern, publication-ready output
- Enhanced color schemes and styling
- Improved plot layouts and interactivity

**4. Improved Reporting**
- Comprehensive assembly statistics section
- Expanded functional annotation reports
- Enhanced data presentation
- Better organization and readability
- Seamless integration with modular visualization system

---

## Contributing

We welcome contributions from the scientific community.

### Priority Areas

| Area | Description | Impact |
|------|-------------|---------|
| **Assembly Optimization** | Fine-tune SPAdes parameters for specific sample types | High |
| **Database Expansion** | Add additional functional databases | High |
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
- **Assembly Guide**: [assembly.md](assembly.md)
- **Annotation Guide**: [annotation.md](annotation.md)
- **GitHub Wiki**: Additional tutorials and examples

### Getting Help
- **Bug Reports**: Submit via GitHub issues with detailed reproduction steps
- **Feature Requests**: Use GitHub discussions for new feature proposals
- **Community Support**: Join our discussion forum for user questions
- **Email Support**: metaquest-support@example.org

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

**Last Updated**: October 2025 - Version 3.6.0 with SPAdes assembly, enhanced annotation, and modular visualization system

---

*MetaQuest is actively developed with regular updates and improvements. Check our GitHub repository for the latest releases and features.*