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
- **Virulence Factor Analysis**: 🔄 **IN DEVELOPMENT** - Virulence factor identification workflows require optimization and integration
- **AMR Analysis**: 🔄 **IN DEVELOPMENT** - Antimicrobial resistance detection needs enhancement and clinical integration

**Major Achievements (August 2025)**: 
- ✅ **Advanced File Validation**: Comprehensive quality control with detailed statistics, format validation, and quality thresholds
- ✅ **Dual-Method Pathogen Detection**: FASTQ uses multi-source traditional screening, FASTA uses integrated BLAST+ML predictions
- ✅ **Clinical Risk Assessment**: Comprehensive risk stratification with emergency protocols and clinical recommendations
- ✅ **Analysis-Specific Dashboards**: Separate optimized interfaces for FASTQ (green theme) and FASTA (blue theme) workflows
- ✅ **Redundancy Elimination**: Streamlined reporting with single definitive pathogen reports per analysis type
- ✅ **Enhanced Visualization**: Fixed pathogen risk charts and detection coverage analysis

### Development Timeline
- **File Validation System**: ✅ **COMPLETED** (v3.2.1) - *New major feature*
- **Taxonomic Classification**: ✅ **COMPLETED** (v3.1.2)
- **Pathogen Detection & Risk Assessment**: ✅ **COMPLETED** (v3.2.0) - *Major milestone achieved*
- **Enhanced Virulence Factor Detection**: Target completion Q3 2025
- **Comprehensive AMR Analysis**: Target completion Q4 2025
- **Integrated Clinical Decision Support**: Target completion Q1 2026

**Current Status**: Both FASTQ and FASTA inputs now provide comprehensive analysis with reliable taxonomic classification and advanced pathogen detection. The pipeline successfully identifies high-risk pathogens (e.g., Salmonella enterica, E. coli, Klebsiella pneumoniae) with proper clinical risk assessment and emergency protocol recommendations. Advanced file validation ensures data quality before analysis begins.

### Key Features

- **Advanced File Validation**: 
  - **FASTQ Quality Control**: Quality score analysis, encoding detection, adapter contamination screening *(Fully operational)*
  - **FASTA Format Validation**: Sequence type detection, composition analysis, duplicate ID checking *(Fully operational)*
  - **Comprehensive Statistics**: MD5 checksums, N50/length metrics, GC content analysis *(Fully operational)*
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

### FASTA Processing Improvements (v3.1.2)

#### Enhanced Taxonomic Classification
- **Complete FASTA Support**: Taxonomic profiling for FASTA files is now fully functional and accurate
- **Improved k-mer Classification**: Optimized Kraken2 integration for assembled sequences
- **Enhanced Visualization**: Accurate taxonomic visualizations (Krona, pie charts, treemaps)
- **Quality Metrics**: Comprehensive assembly quality assessment and reporting

### CLI & Pipeline Enhancements (v3.1.1)

#### Enhanced FASTQ Input Support
- **Mutually Exclusive Input Groups**: Added support for single-end, paired-end, and interleaved FASTQ inputs
- **Automatic Interleaved Processing**: Seamless detection and splitting of interleaved FASTQ files
- **Improved Validation**: Enhanced file existence checking and format validation
- **Flexible Command Interface**: Multiple input options with clear usage patterns

## Validation Features

### FASTQ Validation Capabilities

- **Quality Score Analysis**: Automatic detection of quality encoding (Phred+33/Phred+64)
- **Statistical Metrics**: Mean, median, N50 length calculations with quality score distributions
- **Quality Thresholds**: Q20/Q30 base percentage analysis with customizable thresholds
- **Contamination Detection**: Adapter sequence screening in first 1000 reads
- **Sequence Composition**: GC content analysis and ambiguous base detection
- **File Integrity**: MD5 checksum calculation and compression detection

### FASTA Validation Capabilities

- **Sequence Type Detection**: Automatic classification as protein or nucleotide sequences
- **Composition Analysis**: GC content statistics with per-sequence variance analysis
- **Quality Assessment**: Duplicate ID detection, gap analysis, and N-base counting
- **Length Distribution**: Comprehensive length statistics with coefficient of variation
- **Category Classification**: Intelligent sequence categorization (genome, contigs, genes, etc.)
- **Format Validation**: Strict FASTA format compliance checking

## Installation

MetaQuest requires a Linux/macOS environment with conda package manager. The installation process includes:

1. Environment setup with conda
2. Database downloads (~8GB total)
3. Tool integration and testing

**📖 For detailed installation instructions, see [Installation Guide](docs/installation.md)**

## Usage

### Comprehensive Input Support

MetaQuest now supports optimized workflows for both FASTQ and FASTA inputs with advanced validation:

#### **FASTQ Analysis** *(Recommended for clinical samples)*
- **Pre-Analysis Validation**: Comprehensive quality control with detailed statistics
- **Taxonomic Classification**: Kraken2/Bracken rapid classification
- **Pathogen Detection**: Multi-source traditional screening with clinical risk assessment
- **Output**: Green-themed dashboard with pathogen_summary.txt as the definitive report

#### **FASTA Analysis** *(Optimized for assembled genomes)*
- **Format Validation**: Sequence type detection and composition analysis
- **Taxonomic Classification**: BLAST-based high-accuracy classification
- **Pathogen Detection**: Integrated BLAST taxonomy + ML predictions with separated reporting
- **Output**: Blue-themed dashboard with blast_ml_pathogen_summary.txt as the definitive report

#### **Validation Features**
- **Quality Thresholds**: Customizable minimum quality and sequence count requirements
- **Statistical Analysis**: Comprehensive file statistics including N50, GC content, and quality distributions
- **Format Compliance**: Strict validation of FASTQ/FASTA format requirements
- **Pre-flight Checks**: Validation-only mode for quick file assessment

**📖 For comprehensive usage examples and validation options, see [Usage Guide](docs/usage.md)**

## Contributing

We welcome contributions from the scientific community. Current priority areas include:

- **Advanced Validation**: Help enhance file validation with additional quality metrics and contamination detection
- **Virulence Factor Analysis**: Help optimize VF detection workflows and clinical integration
- **AMR Enhancement**: Improve antimicrobial resistance detection and clinical mapping
- **Clinical Validation**: Real-world testing of pathogen detection accuracy and risk assessment
- **Machine Learning**: Expand ML capabilities for enhanced pathogen prediction
- **Documentation**: User guides and clinical interpretation guidelines

## Contact and Support

- **Documentation**: [Installation Guide](docs/installation.md) | [Usage Guide](docs/usage.md)
- **Issues**: Report bugs and feature requests via GitHub issues
- **Community**: Join our discussion forum for questions and collaboration

---

**MetaQuest Development Team**  
*Advancing metagenomics through integrated computational solutions*

---
*Last updated: August 2025 - Advanced file validation system with comprehensive quality control*