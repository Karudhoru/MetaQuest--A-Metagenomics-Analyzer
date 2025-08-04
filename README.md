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
- ✅ **Advanced File Validation**: Comprehensive quality control with detailed statistics, format validation, and quality thresholds
- ✅ **Machine Learning Pipeline**: Complete ML pathogen prediction system with feature extraction, model training, and inference
- ✅ **Dual-Method Pathogen Detection**: FASTQ uses multi-source traditional screening, FASTA uses integrated BLAST+ML predictions
- ✅ **Clinical Risk Assessment**: Comprehensive risk stratification with emergency protocols and clinical recommendations
- ✅ **Analysis-Specific Dashboards**: Separate optimized interfaces for FASTQ (green theme) and FASTA (blue theme) workflows
- ✅ **Modular Architecture**: Clean separation of concerns with dedicated modules for IO, ML, visualization, and reporting
- ✅ **Enhanced Visualization**: Fixed pathogen risk charts and detection coverage analysis with interactive dashboards

### Development Timeline
- **File Validation System**: ✅ **COMPLETED** (v3.2.1) - *New major feature*
- **Machine Learning Integration**: ✅ **COMPLETED** (v3.2.1) - *Major milestone achieved*
- **Pathogen Detection & Risk Assessment**: ✅ **COMPLETED** (v3.2.0) - *Major milestone achieved*
- **Taxonomic Classification**: ✅ **COMPLETED** (v3.1.2)
- **Enhanced Virulence Factor Detection**: Target completion Q3 2025
- **Comprehensive AMR Analysis**: Target completion Q4 2025
- **Integrated Clinical Decision Support**: Target completion Q1 2026

**Current Status**: Both FASTQ and FASTA inputs now provide comprehensive analysis with reliable taxonomic classification, advanced pathogen detection, and ML-powered predictions. The pipeline successfully identifies high-risk pathogens (e.g., Salmonella enterica, E. coli, Klebsiella pneumoniae) with proper clinical risk assessment and emergency protocol recommendations. Advanced file validation ensures data quality before analysis begins.

## Architecture & Modules

MetaQuest follows a modular architecture with clear separation of concerns:

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
│   └── reporting.py         # Comprehensive reporting
└── visualization/            # Data visualization
    └── visualization.py     # Interactive dashboards
```

### Key Features

- **Advanced File Validation**: 
  - **FASTQ Quality Control**: Quality score analysis, encoding detection, adapter contamination screening *(Fully operational)*
  - **FASTA Format Validation**: Sequence type detection, composition analysis, duplicate ID checking *(Fully operational)*
  - **Comprehensive Statistics**: MD5 checksums, N50/length metrics, GC content analysis *(Fully operational)*

- **Machine Learning Integration**:
  - **Feature Extraction**: Automated extraction of sequence features for ML models *(Fully operational)*
  - **Pathogen Prediction**: Pre-trained models for pathogen classification *(Fully operational)*
  - **Model Artifacts**: Complete set of trained models, scalers, and feature selectors *(Fully operational)*

- **Taxonomic Profiling**: Species-level identification and abundance estimation *(Fully supported for both FASTQ and FASTA)*

- **Advanced Pathogen Screening**: 
  - **FASTQ**: Multi-source traditional screening with Bracken integration *(Fully operational)*
  - **FASTA**: Integrated BLAST taxonomy + ML pathogen predictions *(Fully operational)*

- **Clinical Risk Assessment**: Comprehensive risk stratification with clinical recommendations *(Fully operational)*

- **Antimicrobial Resistance (AMR) Analysis**: Basic resistance gene detection *(Under enhancement)*

- **Virulence Factor Assessment**: Basic pathogenicity determinant identification *(Under development)*

- **Functional Annotation**: Gene prediction and functional characterization *(Fully supported for both formats)*

- **Interactive Dashboards**: Analysis-specific HTML reports with dynamic visualizations *(Fully operational)*

- **Quality Assessment**: Statistical analysis and quality metrics *(Fully operational)*

## Recent Updates

### Machine Learning Integration (v3.2.1) - August 2025

#### ✅ **Complete ML Pipeline**
- **Feature Extraction Module**: Automated extraction of sequence-based features for pathogen prediction
- **Pathogen Prediction Engine**: Pre-trained machine learning models for accurate pathogen classification
- **Model Artifacts**: Complete set of trained models including:
  - `best_model.pkl`: Primary classification model
  - `scaler.pkl`: Feature scaling transformation
  - `feature_selector.pkl`: Feature selection model
  - `feature_names.pkl`: Selected feature names
  - `all_feature_names.pkl`: Complete feature vocabulary

#### ✅ **Enhanced Analysis Workflows**
- **FASTA ML Integration**: Seamless integration of ML predictions with BLAST-based taxonomy
- **Feature Engineering**: Advanced sequence feature extraction including composition, k-mer, and structural features
- **Model Inference**: Real-time pathogen prediction during analysis pipeline
- **Confidence Scoring**: Prediction confidence scores for clinical decision support

### Advanced File Validation System (v3.2.1) - August 2025

#### ✅ **Comprehensive Quality Control**
- **Pre-Analysis Validation**: Thorough file format and quality checking before analysis begins
- **FASTQ Analysis**: Quality encoding detection, Q20/Q30 statistics, adapter contamination screening, sequence length distribution
- **FASTA Analysis**: Sequence type detection (protein/nucleotide), composition analysis, duplicate ID detection, GC content variance
- **Statistical Reporting**: N50 calculations, coefficient of variation, MD5 checksums, file size analysis

#### ✅ **Intelligent Quality Thresholds**
- **Customizable Criteria**: User-defined minimum quality scores and sequence counts
- **Analysis-Specific Validation**: Different validation rules for FASTQ vs FASTA inputs
- **Quality Warnings**: Clear identification of potential issues with actionable recommendations
- **Validation-Only Mode**: Standalone file analysis without running full pipeline

#### ✅ **Enhanced User Experience**
- **Beautiful Terminal Output**: Color-coded statistics with progress indicators and clear pass/fail status
- **Detailed Diagnostics**: Comprehensive file statistics including sequence diversity, quality metrics, and composition analysis
- **Skip Validation Option**: Advanced users can bypass validation when needed
- **Integration with CLI**: Seamless validation integrated into main analysis workflow

### Major Pathogen Analysis Completion (v3.2.0) - July 2025

#### ✅ **Comprehensive Pathogen Detection System**
- **Dual-Method Architecture**: Separate optimized workflows for FASTQ and FASTA inputs
- **FASTQ Pathogen Analysis**: Multi-source traditional screening using Bracken taxonomy, custom databases, AMR/VF scanning
- **FASTA Pathogen Analysis**: Integrated BLAST taxonomy + Machine Learning predictions with separated reporting sections
- **Clinical Risk Assessment**: Emergency protocol recommendations for CRITICAL/HIGH risk samples

#### ✅ **Enhanced Reporting & Visualization**
- **Streamlined Reports**: Single definitive pathogen report per analysis type (pathogen_summary.txt for FASTQ, blast_ml_pathogen_summary.txt for FASTA)
- **Analysis-Specific Dashboards**: Green-themed FASTQ dashboard vs. Blue-themed FASTA dashboard with optimized layouts
- **Fixed Visualizations**: Proper pathogen risk detection charts and method coverage analysis
- **Redundancy Elimination**: Removed all duplicate and overlapping report files

#### ✅ **Production-Ready Pipeline**
- **Error-Free Execution**: Resolved all analysis failures and missing method issues
- **Proper File Detection**: Fixed dashboard generation with correct analysis type identification
- **Clean Output Structure**: Users receive exactly the reports they need without confusion
- **Validated Results**: Successfully detects dangerous pathogens with proper CRITICAL risk classification

## Installation

MetaQuest requires a Linux/macOS environment with conda package manager. The installation process includes:

1. Environment setup with conda
2. Database downloads (~8GB total)
3. Tool integration and testing
4. ML model setup

**📖 For detailed installation instructions, see [Installation Guide](docs/installation.md)**

## Usage

### Comprehensive Input Support

MetaQuest now supports optimized workflows for both FASTQ and FASTA inputs with advanced validation and ML integration:

#### **FASTQ Analysis** *(Recommended for clinical samples)*
- **Pre-Analysis Validation**: Comprehensive quality control with detailed statistics
- **Taxonomic Classification**: Kraken2/Bracken rapid classification
- **Pathogen Detection**: Multi-source traditional screening with clinical risk assessment
- **Output**: Green-themed dashboard with pathogen_summary.txt as the definitive report

#### **FASTA Analysis** *(Optimized for assembled genomes)*
- **Format Validation**: Sequence type detection and composition analysis
- **Taxonomic Classification**: BLAST-based high-accuracy classification
- **ML Pathogen Prediction**: Advanced machine learning-based pathogen classification
- **Pathogen Detection**: Integrated BLAST taxonomy + ML predictions with separated reporting
- **Output**: Blue-themed dashboard with blast_ml_pathogen_summary.txt as the definitive report

#### **Validation Features**
- **Quality Thresholds**: Customizable minimum quality and sequence count requirements
- **Statistical Analysis**: Comprehensive file statistics including N50, GC content, and quality distributions
- **Format Compliance**: Strict validation of FASTQ/FASTA format requirements
- **Pre-flight Checks**: Validation-only mode for quick file assessment

**📖 For comprehensive usage examples and validation options, see [Usage Guide](docs/usage.md)**

## Output Structure

MetaQuest generates comprehensive outputs organized by analysis type:

### FASTQ Analysis Outputs
```
results/
├── analysis_dashboard.html          # Interactive green-themed dashboard
├── pathogen_summary.txt            # Definitive pathogen report
├── taxonomic_classification.tsv     # Kraken2/Bracken results
├── quality_stats.txt               # File validation statistics
├── functional_annotation/          # Gene prediction results
└── visualization/                   # Charts and plots
```

### FASTA Analysis Outputs
```
results/
├── analysis_dashboard.html          # Interactive blue-themed dashboard
├── blast_ml_pathogen_summary.txt   # Definitive BLAST+ML pathogen report
├── ml_pathogen_predictions.csv     # Detailed ML predictions
├── taxonomic_classification.tsv     # BLAST taxonomy results
├── quality_stats.txt               # File validation statistics
├── functional_annotation/          # Gene prediction results
└── visualization/                   # Charts and plots
```

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

- **Documentation**: [Installation Guide](docs/installation.md) | [Usage Guide](docs/usage.md)
- **Issues**: Report bugs and feature requests via GitHub issues
- **Community**: Join our discussion forum for questions and collaboration

---

**MetaQuest Development Team**  
*Advancing metagenomics through integrated computational solutions*

---
*Last updated: August 2025 - Machine learning integration with comprehensive pathogen prediction capabilities*