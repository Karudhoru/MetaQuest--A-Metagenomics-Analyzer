# MetaQuest

A comprehensive metagenomics analysis pipeline for taxonomic classification, pathogen detection, and functional annotation of microbial communities.

## Overview

MetaQuest is an integrated bioinformatics pipeline that addresses the complex challenges of metagenomic data analysis. By combining state-of-the-art tools and databases, MetaQuest provides researchers with a streamlined workflow for understanding microbial community composition, identifying potential pathogens, and assessing antimicrobial resistance profiles from sequencing data.

## ⚠️ Development Status

**MetaQuest is currently under active development.** The core functionality is operational with major improvements completed across all processing capabilities:

### Current Processing Status
- **Taxonomic Classification**: ✅ **COMPLETED** - Both FASTQ and FASTA taxonomic profiling are now accurate and fully functional
- **Pathogenicity Assessment**: ✅ **COMPLETED** - Comprehensive pathogen detection workflows for both FASTQ and FASTA inputs with clinical risk assessment
- **Virulence Factor Analysis**: 🔄 **IN DEVELOPMENT** - Virulence factor identification workflows require optimization and integration
- **AMR Analysis**: 🔄 **IN DEVELOPMENT** - Antimicrobial resistance detection needs enhancement and clinical integration

**Major Achievements (July 2025)**: 
- ✅ **Dual-Method Pathogen Detection**: FASTQ uses multi-source traditional screening, FASTA uses integrated BLAST+ML predictions
- ✅ **Clinical Risk Assessment**: Comprehensive risk stratification with emergency protocols and clinical recommendations
- ✅ **Analysis-Specific Dashboards**: Separate optimized interfaces for FASTQ (green theme) and FASTA (blue theme) workflows
- ✅ **Redundancy Elimination**: Streamlined reporting with single definitive pathogen reports per analysis type
- ✅ **Enhanced Visualization**: Fixed pathogen risk charts and detection coverage analysis

### Development Timeline
- **Taxonomic Classification**: ✅ **COMPLETED** (v3.1.2)
- **Pathogen Detection & Risk Assessment**: ✅ **COMPLETED** (v3.2.0) - *Major milestone achieved*
- **Enhanced Virulence Factor Detection**: Target completion Q3 2025
- **Comprehensive AMR Analysis**: Target completion Q4 2025
- **Integrated Clinical Decision Support**: Target completion Q1 2026

**Current Status**: Both FASTQ and FASTA inputs now provide comprehensive analysis with reliable taxonomic classification and advanced pathogen detection. The pipeline successfully identifies high-risk pathogens (e.g., Salmonella enterica, E. coli, Klebsiella pneumoniae) with proper clinical risk assessment and emergency protocol recommendations.

### Key Features

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

## Future Directions

### Immediate Priorities (v3.3-3.4)

#### **Virulence Factor Analysis Enhancement** *(High Priority)*
- **Advanced VF Detection**: Optimize virulence factor identification workflows for both FASTQ and FASTA
- **Clinical Integration**: Map virulence factors to clinical severity and treatment protocols
- **Interactive VF Visualization**: Enhanced charts showing virulence profiles with risk assessment
- **Pathogenesis Pathway Analysis**: Integrate VF findings with pathogen detection for comprehensive threat assessment

#### **AMR Analysis Optimization** *(High Priority)*
- **Enhanced Resistance Detection**: Improve antimicrobial resistance gene identification accuracy
- **Clinical Breakpoint Integration**: Map resistance genes to clinical susceptibility testing
- **Treatment Recommendation Engine**: Suggest appropriate antimicrobial therapy based on resistance profiles
- **Resistance Mechanism Classification**: Detailed categorization by resistance mechanisms and drug classes

#### **Integrated Clinical Decision Support**
- **Unified Risk Assessment**: Combine pathogen, virulence, and resistance data for comprehensive threat evaluation
- **Emergency Protocol Enhancement**: Expand clinical recommendations based on integrated findings
- **Multi-sample Comparison**: Comparative analysis dashboards for outbreak investigation
- **Real-time Alerting**: Automated high-risk pathogen detection notifications

### Medium-term Developments (v4.0+)

#### **Machine Learning Integration**
- **Enhanced ML Pathogen Prediction**: Expand ML capabilities beyond current FASTA integration
- **Resistance Phenotype Prediction**: ML models for clinical resistance prediction from genotype
- **Outbreak Detection**: Automated surveillance and clustering for epidemiological investigation
- **Treatment Outcome Prediction**: Integration of resistance/virulence data for therapy success prediction

#### **Advanced Clinical Features**
- **Multi-omics Integration**: Combine metagenomic data with host genomics and clinical parameters
- **Longitudinal Analysis**: Track pathogen evolution and resistance development over time
- **Population Surveillance**: Large-scale monitoring for public health applications
- **Point-of-Care Integration**: Rapid diagnostic platform compatibility

## Installation

MetaQuest requires a Linux/macOS environment with conda package manager. The installation process includes:

1. Environment setup with conda
2. Database downloads (~8GB total)
3. Tool integration and testing

**📖 For detailed installation instructions, see [Installation Guide](docs/installation.md)**

## Usage

### Comprehensive Input Support

MetaQuest now supports optimized workflows for both FASTQ and FASTA inputs:

#### **FASTQ Analysis** *(Recommended for clinical samples)*
- **Taxonomic Classification**: Kraken2/Bracken rapid classification
- **Pathogen Detection**: Multi-source traditional screening with clinical risk assessment
- **Output**: Green-themed dashboard with pathogen_summary.txt as the definitive report

#### **FASTA Analysis** *(Optimized for assembled genomes)*
- **Taxonomic Classification**: BLAST-based high-accuracy classification
- **Pathogen Detection**: Integrated BLAST taxonomy + ML predictions with separated reporting
- **Output**: Blue-themed dashboard with blast_ml_pathogen_summary.txt as the definitive report

**📖 For comprehensive usage examples and options, see [Usage Guide](docs/usage.md)**

## Contributing

We welcome contributions from the scientific community. Current priority areas include:

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
*Last updated: July 2025 - Major pathogen analysis completion with comprehensive risk assessment and clinical decision support*
