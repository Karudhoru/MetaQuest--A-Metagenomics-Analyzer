# MetaQuest

> **A comprehensive metagenomics analysis pipeline for taxonomic classification, pathogen detection, and functional annotation of microbial communities.**

---

## 📋 Table of Contents

- [Overview](#overview)
- [Development Status](#development-status)
- [Architecture & Features](#architecture--features)
- [Installation](#installation)
- [Usage](#usage)
- [Output Structure](#output-structure)
- [Contributing](#contributing)
- [Support & Contact](#support--contact)

---

## 🔬 Overview

MetaQuest is an integrated bioinformatics pipeline that addresses the complex challenges of metagenomic data analysis. By combining state-of-the-art tools and databases, MetaQuest provides researchers with a streamlined workflow for understanding microbial community composition, identifying potential pathogens, and assessing antimicrobial resistance profiles from sequencing data.

---

## ⚠️ Development Status

**MetaQuest is currently under active development** with core functionality operational and major improvements completed across all processing capabilities.

### ✅ Completed Components

| Component | Status | Description |
|-----------|---------|-------------|
| **File Validation & Quality Control** | ✅ Complete | Comprehensive validation with detailed statistics for both FASTQ and FASTA inputs |
| **Taxonomic Classification** | ✅ Complete | Both FASTQ and FASTA taxonomic profiling are accurate and fully functional |
| **Pathogenicity Assessment** | ✅ Complete | Comprehensive pathogen detection workflows with clinical risk assessment |
| **Machine Learning Integration** | ✅ Complete | Advanced ML-based pathogen prediction with feature extraction and model artifacts |

### 🔄 In Development

| Component | Status | Target Completion |
|-----------|---------|-------------------|
| **Virulence Factor Analysis** | 🔄 In Progress | Q3 2025 |
| **AMR Analysis** | 🔄 In Progress | Q4 2025 |
| **scikit-bio Integration** | 🔄 Planned | Q4 2025 |

### 🎉 Major Achievements (August 2025)

- ✅ **Advanced File Validation**: Comprehensive quality control with contamination detection
- ✅ **Machine Learning Pipeline**: Complete ML pathogen prediction system
- ✅ **Dual-Method Pathogen Detection**: Optimized workflows for FASTQ and FASTA inputs
- ✅ **Clinical Risk Assessment**: Comprehensive risk stratification with emergency protocols
- ✅ **Analysis-Specific Dashboards**: Separate optimized interfaces with themed visualizations
- ✅ **Modular Architecture**: Clean separation of concerns with dedicated modules
- ✅ **Enhanced Visualization**: Interactive dashboards with comprehensive analysis reports

### 📈 Recent Version History

| Version | Release | Key Features |
|---------|---------|-------------|
| **v3.2.2** | Aug 2025 | Professional CLI revamp, advanced reporting with scikit-bio integration |
| **v3.2.1** | Aug 2025 | Complete ML pipeline, advanced file validation system |
| **v3.2.0** | Aug 2025 | Pathogen detection & clinical risk assessment |

---

## 🏗️ Architecture & Features

### System Architecture

```
src/metaquest/
├── cli.py                    # Command-line interface
├── config.py                 # Configuration management
├── core/                     # Core analysis modules
│   ├── analysis.py          # Main analysis orchestration
│   ├── taxonomic_analysis.py # Taxonomic classification
│   ├── pathogen_analysis.py  # Pathogen detection
│   └── functional_analysis.py # Functional annotation
├── io/                       # Input/Output handling
│   ├── file_validator.py     # File validation & QC
│   └── utils.py             # I/O utilities
├── ml/                       # Machine Learning components
│   ├── feature_extractor.py  # Feature extraction for ML
│   ├── pathogen_predictor.py # ML pathogen prediction
│   └── model_artifacts/      # Pre-trained models
├── reporting/                # Report generation
│   ├── reporting.py         # Comprehensive reporting
│   └── dashboard.py         # Interactive dashboards
└── visualization/            # Data visualization
    └── visualization.py     # Interactive visualizations
```

### 🔑 Key Features

#### 🔍 Advanced File Validation
- **FASTQ Quality Control**: Quality score analysis, overrepresented sequence detection, adapter contamination screening
- **FASTA Format Validation**: Sequence type detection, composition analysis, duplicate ID checking
- **Comprehensive Statistics**: MD5 checksums, N50/length metrics, GC content analysis

#### 🤖 Machine Learning Integration
- **Feature Extraction**: Automated extraction of sequence features for ML models
- **Pathogen Prediction**: Pre-trained models for pathogen classification
- **Model Artifacts**: Complete set of trained models, scalers, and feature selectors

#### 🦠 Advanced Pathogen Screening
- **FASTQ Processing**: Multi-source traditional screening with Bracken integration
- **FASTA Processing**: Integrated BLAST taxonomy + ML pathogen predictions
- **Clinical Risk Assessment**: Comprehensive risk stratification with clinical recommendations

#### 📊 Analysis & Reporting
- **Taxonomic Profiling**: Species-level identification with scikit-bio diversity metrics
- **Interactive Dashboards**: Analysis-specific HTML reports with dynamic visualizations
- **Quality Assessment**: Statistical analysis and quality metrics
- **Functional Annotation**: Gene prediction and functional characterization

#### 🔬 Specialized Analysis (In Development)
- **Antimicrobial Resistance (AMR)**: Resistance gene detection
- **Virulence Factor Assessment**: Virulence factor identification

---

## 🚀 Installation

MetaQuest requires a Linux/macOS environment with conda package manager.

**📖 For detailed installation instructions, see [installation.md](installation.md)**

---

## 💻 Usage

MetaQuest supports optimized workflows for both FASTQ and FASTA inputs with advanced validation and ML integration.

### Command Structure
```bash
# Analyze samples
metaquest analyze [options]

# Validate files
metaquest validate [options]

# System check
metaquest check
```

**📖 For comprehensive usage examples and validation options, see [usage.md](usage.md)**

---

## 📁 Output Structure

MetaQuest generates comprehensive outputs organized by analysis type:

### FASTQ Analysis Outputs
```
results/
├── analysis_dashboard.html              # Interactive green-themed dashboard
├── pathogen_summary.txt                 # Definitive pathogen report
├── taxonomic_classification_report.txt  # Kraken2/Bracken results
├── quality_stats.txt                    # File validation statistics
├── functional_annotation/               # Gene prediction results
└── ... (additional reports and visualizations)
```

### FASTA Analysis Outputs
```
results/
├── analysis_dashboard.html              # Interactive blue-themed dashboard
├── blast_ml_pathogen_summary.txt        # BLAST+ML pathogen report
├── ml_pathogen_predictions.csv          # Detailed ML predictions
├── taxonomic_classification_report.txt  # BLAST taxonomy results
├── quality_stats.txt                    # File validation statistics
├── functional_annotation/               # Gene prediction results
└── ... (additional reports and visualizations)
```

---

## 🤝 Contributing

We welcome contributions from the scientific community! 

### 🎯 Priority Areas

| Area | Description | Impact |
|------|-------------|---------|
| **Machine Learning Enhancement** | Improve ML models with additional training data | High |
| **Advanced Validation** | Enhance file validation with quality metrics | Medium |
| **Virulence Factor Analysis** | Optimize VF detection workflows | High |
| **AMR Enhancement** | Improve antimicrobial resistance detection | High |
| **Clinical Validation** | Real-world testing of pathogen detection | Critical |
| **Documentation** | User guides and clinical interpretation | Medium |

### 📋 How to Contribute
1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request
5. Participate in code review

---

## 📞 Support & Contact

### 📚 Documentation
- **Installation Guide**: [installation.md](installation.md)
- **Usage Examples**: [usage.md](usage.md)

### 🐛 Issues & Support
- **Bug Reports**: Submit via GitHub issues
- **Feature Requests**: Use GitHub discussions
- **Community Forum**: Join our discussion forum

### 👥 Development Team
**MetaQuest Development Team**  
*Advancing metagenomics through integrated computational solutions*

---

## 📄 License & Citation

*Documentation and citation information will be added upon publication*

---

**Last Updated**: August 2025 - Machine learning integration with comprehensive pathogen prediction capabilities

---

*MetaQuest is actively developed with regular updates and improvements. Check our GitHub repository for the latest releases and features.*
