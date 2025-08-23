# MetaQuest

A comprehensive metagenomics analysis pipeline for taxonomic classification, pathogen detection, and functional annotation of microbial communities.

## Overview

MetaQuest is an integrated bioinformatics pipeline that addresses the complex challenges of metagenomic data analysis. By combining state-of-the-art tools and databases, MetaQuest provides researchers with a streamlined workflow for understanding microbial community composition, identifying potential pathogens, and assessing antimicrobial resistance profiles from sequencing data.

## ⚠️ Development Status

**MetaQuest is currently under active development.** The core functionality is operational with major improvements completed across all processing capabilities:

### Current Processing Status
- **File Validation & Quality Control**: ✅ **COMPLETED** - Comprehensive validation with detailed statistics for both FASTQ and FASTA inputs
- **Taxonomic Classification**: ✅ **COMPLETED** - Both FASTQ and FASTA taxonomic profiling are now accurate and fully functional
- **Pathogenicity Assessment**: ✅ **COMPLETED** - Comprehensive pathogen detection workflows for both FASTQ and FASTA inputs with clinical risk assessment
- **Machine Learning Integration**: ✅ **COMPLETED** - Advanced ML-based pathogen prediction with feature extraction and model artifacts
- **Virulence Factor Analysis**: 🔄 **IN DEVELOPMENT** - Virulence factor identification workflows require optimization and integration
- **AMR Analysis**: 🔄 **IN DEVELOPMENT** - Antimicrobial resistance detection needs enhancement and clinical integration

**Major Achievements (August 2025)**: 
- ✅ **Advanced File Validation**: Comprehensive quality control with detailed statistics, format validation, quality thresholds, adapter contamination and overrepresented sequences.
- ✅ **Machine Learning Pipeline**: Complete ML pathogen prediction system with feature extraction, model training, and inference
- ✅ **Dual-Method Pathogen Detection**: FASTQ uses multi-source traditional screening, FASTA uses integrated BLAST+ML predictions
- ✅ **Clinical Risk Assessment**: Comprehensive risk stratification with emergency protocols and clinical recommendations
- ✅ **Analysis-Specific Dashboards**: Separate optimized interfaces for FASTQ (green theme) and FASTA (blue theme) workflows
- ✅ **Modular Architecture**: Clean separation of concerns with dedicated modules for IO, ML, visualization, and reporting
- ✅ **Enhanced Visualization**: Fixed pathogen risk charts and detection coverage analysis with interactive dashboards

### Development Timeline
- **Updated CLI**: ✅ **COMPLETED** (v3.2.2) - *New major revamp*
- **File Validation System**: ✅ **COMPLETED** (v3.2.1) - *New major feature*
- **Machine Learning Integration**: ✅ **COMPLETED** (v3.2.1) - *Major milestone achieved*
- **Pathogen Detection & Risk Assessment**: ✅ **COMPLETED** (v3.2.0) - *Major milestone achieved*
- **Enhanced Virulence Factor Detection**: Target completion Q3 2025
- **Integrating "scikit-bio" Package**: Target completion Q4 2025
- **Comprehensive AMR Analysis**: Target completion Q4 2025

**Current Status**: Both FASTQ and FASTA inputs now provide comprehensive analysis with reliable taxonomic classification, advanced pathogen detection, and ML-powered predictions. The pipeline successfully identifies high-risk pathogens (e.g., Salmonella enterica, E. coli, Klebsiella pneumoniae) with proper clinical risk assessment and emergency protocol recommendations. Advanced file validation ensures data quality before analysis begins.

## Architecture & Modules

MetaQuest follows a modular architecture with clear separation of concerns:

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
│   |── reporting.py         # Comprehensive reporting
|   ├── dashboard.py         # Comprehensive Dashboard
└── visualization/            # Data visualization
    └── visualization.py     # Interactive dashboards


### Key Features
**Advanced File Validation:**

**FASTQ Quality Control:** Quality score analysis, overrepresented sequence detection, and adapter contamination screening (Fully operational)

**FASTA Format Validation:** Sequence type detection, composition analysis, duplicate ID checking (Fully operational)

**Comprehensive Statistics:** MD5 checksums, N50/length metrics, GC content analysis (Fully operational)

#### Machine Learning Integration:

**Feature Extraction:** Automated extraction of sequence features for ML models (Fully operational)

**Pathogen Prediction:** Pre-trained models for pathogen classification (Fully operational)

**Model Artifacts:** Complete set of trained models, scalers, and feature selectors (Fully operational)

**Taxonomic Profiling:** Species-level identification and abundance estimation with robust scikit-bio diversity metrics (Fully supported for both FASTQ and FASTA)

#### Advanced Pathogen Screening:

**FASTQ:** Multi-source traditional screening with Bracken integration (Fully operational)

**FASTA:** Integrated BLAST taxonomy + ML pathogen predictions (Fully operational)

**Clinical Risk Assessment:** Comprehensive risk stratification with clinical recommendations (Fully operational)

**Antimicrobial Resistance (AMR) Analysis:** Basic resistance gene detection (In Development)

**Virulence Factor Assessment:** Basic virulence factor identification (In Development)

**Functional Annotation:** Gene prediction and functional characterization (Fully supported for both formats)

**Interactive Dashboards:** Sleek, professional, and analysis-specific HTML reports with dynamic visualizations (Fully operational)

**Quality Assessment:** Statistical analysis and quality metrics (Fully operational)

### Recent Updates

#### CLI Revamp & Reporting Polish (v3.2.2) - August 2025
##### ✅ Professional Command-Line Interface
**New Subcommand Structure:** The CLI is now organized into clear, intuitive commands: analyze, validate, and check.

**Streamlined Arguments:** Input for FASTQ is handled with explicit flags (--single, --paired, --interleaved) for clarity.

**Convenience Features:** Added short-form arguments (-o, -s, etc.) and a professional header for all commands.

##### ✅ Advanced Reporting & Insights
**scikit-bio Integration:** Upgraded alpha diversity reporting to include robust metrics like Shannon, Simpson, and Chao1 for more rigorous ecological analysis.

**Sequence Hit Summaries:** Added a new, dedicated reporter that transforms raw DIAMOND output (e.g., from AMR or pathogen screens) into human-readable summary tables and insights.

**Sleek Dashboards:** Completely redesigned the HTML dashboard with a modern, professional, and animated interface with distinct themes for FASTQ and FASTA results.

#### Machine Learning Integration (v3.2.1) - August 2025
##### ✅ Complete ML Pipeline
**Feature Extraction Module:** Automated extraction of sequence-based features for pathogen prediction.

**Pathogen Prediction Engine:** Pre-trained machine learning models for accurate pathogen classification.

**Model Artifacts:** Complete set of trained models including best_model.pkl, scaler.pkl, feature_selector.pkl, and feature name files.

#### Advanced File Validation System (v3.2.1) - August 2025
##### ✅ Comprehensive Quality Control
**Pre-Analysis Validation:** Thorough file format and quality checking before analysis begins.

**Contamination & Duplication Checks:** Added robust detection for adapter contamination and overrepresented sequences, with actionable recommendations.

**Intelligent Quality Thresholds:** Implemented customizable criteria for all key validation metrics via CLI flags.

## Installation

MetaQuest requires a Linux/macOS environment with conda package manager.

**📖 For detailed installation instructions, see [installation.md](installation.md)**

## Usage

MetaQuest now supports optimized workflows for both FASTQ and FASTA inputs with advanced validation and ML integration.

**📖 For comprehensive usage examples and validation options, see [usage.md](usage.md)**

## Output Structure

MetaQuest generates comprehensive outputs organized by analysis type:

### FASTQ Analysis Outputs
results/
├── analysis_dashboard.html          # Interactive green-themed dashboard
├── pathogen_summary.txt            # Definitive pathogen report
├── taxonomic_classification_report.txt     # Kraken2/Bracken results summary
├── quality_stats.txt               # File validation statistics
├── functional_annotation/          # Gene prediction results
└── ... (other reports and visualization files)


### FASTA Analysis Outputs
results/
├── analysis_dashboard.html          # Interactive blue-themed dashboard
├── blast_ml_pathogen_summary.txt   # Definitive BLAST+ML pathogen report
├── ml_pathogen_predictions.csv     # Detailed ML predictions
├── taxonomic_classification_report.txt     # BLAST taxonomy results summary
├── quality_stats.txt               # File validation statistics
├── functional_annotation/          # Gene prediction results
└── ... (other reports and visualization files)


## Contributing

We welcome contributions from the scientific community. Current priority areas include:

- **Machine Learning Enhancement**: Improve ML models with additional training data and feature engineering
- **Advanced Validation**: Help enhance file validation with additional quality metrics and contamination detection
- **Virulence Factor Analysis**: Help optimize VF detection workflows and clinical integration
- **AMR Enhancement**: Improve antimicrobial resistance detection and clinical mapping
- **Clinical Validation**: Real-world testing of pathogen detection accuracy and risk assessment
- **Model Training**: Contribute to expanding ML training datasets and model optimization
- **Documentation**: User guides and clinical interpretation guidelines

## Contact and Support

- **Documentation**: [installation.md](installation.md) | [usage.md](usage.md)
- **Issues**: Report bugs and feature requests via GitHub issues
- **Community**: Join our discussion forum for questions and collaboration

---

**MetaQuest Development Team** *Advancing metagenomics through integrated computational solutions*

---

*Last updated: August 2025 - Machine learning integration with comprehensive pathogen prediction capabilities*
