# MetaQuest

A comprehensive metagenomics analysis pipeline for taxonomic classification, pathogen detection, and functional annotation of microbial communities.

## Overview

MetaQuest is an integrated bioinformatics pipeline that addresses the complex challenges of metagenomic data analysis. By combining state-of-the-art tools and databases, MetaQuest provides researchers with a streamlined workflow for understanding microbial community composition, identifying potential pathogens, and assessing antimicrobial resistance profiles from sequencing data.

## ⚠️ Development Status

**MetaQuest is currently under active development.** The core functionality is operational with significant improvements in FASTA processing capabilities:

### Current FASTA Processing Status
- **Taxonomic Classification**: ✅ **COMPLETED** - FASTA taxonomic profiling is now accurate and fully functional
- **Pathogenicity Assessment**: 🔄 **IN DEVELOPMENT** - Pathogen detection workflows for FASTA inputs require optimization
- **Virulence Factor Analysis**: 🔄 **IN DEVELOPMENT** - Virulence factor identification from FASTA files needs refinement
- **AMR Analysis**: 🔄 **IN DEVELOPMENT** - Antimicrobial resistance detection for FASTA inputs under active development

**Current Status**: FASTA files now provide reliable taxonomic classification. Pathogen screening, AMR analysis, and virulence factor detection are functional but undergoing optimization for improved accuracy.

### Development Timeline
- **FASTA Taxonomic Classification**: ✅ **COMPLETED** (v3.1.2)
- **Enhanced Pathogen Detection**: Target completion Q3 2025
- **AMR & Virulence Optimization**: Target completion Q4 2025
- **Comprehensive Validation**: Ongoing testing across file formats

**Recommendation**: Both FASTQ and FASTA inputs are now viable, with FASTA providing excellent taxonomic results. For pathogen detection and AMR analysis, FASTQ inputs remain more reliable until optimization is complete.

### Key Features

- **Taxonomic Profiling**: Species-level identification and abundance estimation *(Fully supported for both FASTQ and FASTA)*
- **Pathogen Screening**: Detection of bacterial, viral, and fungal pathogens *(FASTQ recommended, FASTA under optimization)*
- **Antimicrobial Resistance (AMR) Analysis**: Comprehensive resistance gene detection *(FASTQ recommended, FASTA under optimization)*
- **Virulence Factor Assessment**: Identification of pathogenicity determinants *(FASTQ recommended, FASTA under optimization)*
- **Functional Annotation**: Gene prediction and functional characterization *(Fully supported for both formats)*
- **Quality Assessment**: Statistical analysis and quality metrics
- **Interactive Visualization**: Rich HTML reports with dynamic plots

## Theoretical Background

### Metagenomics Workflow

MetaQuest implements a comprehensive metagenomics analysis workflow based on established computational biology principles:

#### 1. **Sequence Quality Assessment**
- Statistical analysis of read length distribution
- GC content profiling for contamination detection
- Base quality scoring and filtering
- Assembly quality metrics (N50, coverage, contiguity)

#### 2. **Taxonomic Classification**
The pipeline employs a k-mer based classification approach:
- **K-mer Analysis**: Uses Kraken2's exact k-mer matching against a curated database
- **Lowest Common Ancestor (LCA) Algorithm**: Resolves taxonomic assignments when k-mers match multiple taxa
- **Abundance Estimation**: Bracken's Bayesian approach redistributes reads to species level
- **Confidence Scoring**: Statistical confidence assessment for each taxonomic assignment

#### 3. **Pathogen Detection Strategy**
Multi-layered approach for pathogen identification:
- **Primary Screening**: Taxonomic classification against known pathogen databases
- **Sequence Similarity**: DIAMOND/BLAST alignment for validation
- **Risk Assessment**: Integration of pathogenicity databases (CARD, VFDB)
- **Clinical Relevance**: Priority scoring based on medical significance

#### 4. **Antimicrobial Resistance Analysis**
Comprehensive AMR profiling methodology:
- **Gene Homology**: Sequence similarity search against CARD database
- **Resistance Mechanisms**: Classification by mechanism (efflux, modification, target alteration)
- **Drug Classes**: Mapping to specific antimicrobial categories
- **Clinical Context**: Integration with resistance breakpoints and clinical guidelines

#### 5. **Functional Annotation Framework**
Multi-step functional characterization:
- **Gene Prediction**: Open reading frame identification using Prokka
- **Homology Search**: Sequence similarity against SwissProt database
- **Functional Classification**: GO terms, KEGG pathways, and enzyme classification
- **Metabolic Reconstruction**: Pathway completeness assessment

## Tools and Technologies

### Core Bioinformatics Tools

#### **Kraken2** - Taxonomic Classification Engine
- **Algorithm**: Exact k-mer matching with minimizer-based indexing
- **Database**: Pre-built taxonomic database with bacterial, archaeal, viral genomes
- **Performance**: Ultra-fast classification (~1M reads/minute)
- **Accuracy**: High precision taxonomic assignment with confidence scoring

#### **Bracken** - Abundance Estimation
- **Method**: Bayesian re-estimation of taxonomic abundances
- **Input**: Kraken2 classification results
- **Output**: Species-level abundance profiles with statistical confidence
- **Application**: Corrects for database composition bias

#### **DIAMOND** - High-Performance Sequence Aligner
- **Algorithm**: Double Index Alignment of Next-generation sequencing Data
- **Speed**: 100-20,000x faster than BLASTX
- **Sensitivity**: Comparable to BLAST with optimized scoring matrices
- **Use Case**: Pathogen detection and functional annotation

#### **Prokka** - Rapid Prokaryotic Genome Annotation
- **Pipeline**: Integrates multiple annotation tools (Prodigal, HMMER, BLAST)
- **Features**: Gene prediction, tRNA/rRNA detection, functional annotation
- **Output**: Standard formats (GFF3, GenBank, FASTA)
- **Speed**: Complete bacterial genome annotation in <10 minutes

#### **Seqkit** - Fast Sequence Processing Toolkit
- **Purpose**: High-performance FASTQ/FASTA file manipulation
- **Features**: Format conversion, quality filtering, sequence splitting
- **Application**: Interleaved FASTQ splitting and preprocessing
- **Speed**: Optimized for large-scale sequence processing

### Specialized Databases

#### **CARD** - Comprehensive Antibiotic Resistance Database
- **Content**: >6,000 reference sequences across all major AMR gene families
- **Organization**: Hierarchical classification by mechanism and drug class
- **Curation**: Expert-curated with regular updates
- **Integration**: Direct mapping to clinical resistance phenotypes

#### **VFDB** - Virulence Factor Database
- **Scope**: Comprehensive collection of bacterial virulence factors
- **Classification**: Organized by pathogenesis mechanism
- **Validation**: Experimentally verified virulence associations
- **Coverage**: >2,500 virulence factor genes from >70 bacterial genera

#### **MiniKraken2** - Optimized Taxonomic Database
- **Size**: ~8GB compressed database for efficient processing
- **Content**: Representative genomes from bacteria, archaea, viruses
- **Performance**: Balanced between speed and taxonomic resolution
- **Maintenance**: Regular updates with new reference genomes

### Visualization and Reporting

#### **Krona** - Interactive Taxonomic Visualization
- **Format**: HTML5-based hierarchical pie charts
- **Interactivity**: Zoom, filter, and explore taxonomic distributions
- **Integration**: Direct import from Kraken2 results
- **Customization**: Multiple visualization modes and color schemes

#### **Plotly** - Interactive Scientific Plotting
- **Technology**: JavaScript-based interactive plotting library
- **Features**: Zoom, pan, hover tooltips, data export
- **Chart Types**: Bar plots, heatmaps, scatter plots, treemaps
- **Output**: Self-contained HTML files for easy sharing

## Scientific Applications

### Clinical Metagenomics
- **Infectious Disease Diagnosis**: Rapid pathogen identification from clinical samples
- **AMR Surveillance**: Population-level resistance monitoring
- **Outbreak Investigation**: Source tracking and transmission analysis
- **Personalized Medicine**: Treatment selection based on resistance profiles

### Environmental Microbiology
- **Microbiome Studies**: Community structure and function analysis
- **Contamination Assessment**: Pathogen detection in water/food samples
- **Biodegradation Studies**: Functional potential assessment
- **Ecosystem Monitoring**: Microbial community health indicators

### Agricultural Applications
- **Plant Microbiome**: Beneficial and pathogenic microorganism detection
- **Soil Health**: Microbial diversity and functional assessment
- **Food Safety**: Pathogen screening in agricultural products
- **Livestock Health**: Gut microbiome and pathogen monitoring

## Installation

MetaQuest requires a Linux/macOS environment with conda package manager. The installation process includes:

1. Environment setup with conda
2. Database downloads (~8GB total)
3. Tool integration and testing

**📖 For detailed installation instructions, see [Installation Guide](docs/installation.md)**

## Usage

### Enhanced FASTQ Input Support

MetaQuest now supports multiple FASTQ input modes for improved flexibility:

**📖 For detailed installation instructions, see [Usage Guide](docs/usage.md)**


### Command-line Options

#### Required Arguments
- `--type` (`-t`): Input file type (`fastq` or `fasta`)
- `--output` (`-o`): Output directory path

#### FASTQ Input Options (Mutually Exclusive)
- `--reads` (`-r`): Single-end FASTQ file
- `--reads1` (`-1`) + `--reads2` (`-2`): Paired-end FASTQ files
- `--interleaved` (`-i`): Interleaved paired-end FASTQ file

#### FASTA Input Options
- Input file specified as positional argument

### Input Validation and Processing

MetaQuest automatically:
- Validates file existence and format
- Detects interleaved FASTQ files and splits them using seqkit
- Configures analysis parameters based on input type
- Provides informative error messages for invalid inputs

**📖 For comprehensive usage examples and options, see [Usage Guide](docs/usage.md)**

## Recent Updates

### FASTA Processing Improvements (v3.1.2)

#### Enhanced Taxonomic Classification
- **Complete FASTA Support**: Taxonomic profiling for FASTA files is now fully functional and accurate
- **Improved k-mer Classification**: Optimized Kraken2 integration for assembled sequences
- **Enhanced Visualization**: Accurate taxonomic visualizations (Krona, pie charts, treemaps)
- **Quality Metrics**: Comprehensive assembly quality assessment and reporting

#### Ongoing Developments
- **Pathogen Detection**: Active optimization of pathogen screening workflows for FASTA inputs
- **AMR Analysis**: Enhanced antimicrobial resistance detection algorithms under development
- **Virulence Factors**: Improved virulence factor identification accuracy in progress
- **Performance Tuning**: Continued optimization for large FASTA datasets

### CLI & Pipeline Enhancements (v3.1.1)

#### Enhanced FASTQ Input Support
- **Mutually Exclusive Input Groups**: Added support for single-end, paired-end, and interleaved FASTQ inputs
- **Automatic Interleaved Processing**: Seamless detection and splitting of interleaved FASTQ files
- **Improved Validation**: Enhanced file existence checking and format validation
- **Flexible Command Interface**: Multiple input options with clear usage patterns

#### Technical Improvements
- **Seqkit Integration**: Added seqkit for efficient interleaved FASTQ splitting
- **Enhanced Logging**: Improved progress reporting and paired-end detection
- **Streamlined Processing**: Optimized workflow for different input formats
- **Error Handling**: Better error messages and validation feedback

#### Backend Updates
- **CLI Module**: Restructured argument parsing with mutually exclusive groups
- **Analysis Pipeline**: Enhanced to handle list-based input file processing
- **Taxonomic Classification**: Improved Kraken2 integration for paired-end reads
- **Utility Functions**: New helper functions for format conversion and file splitting

## Future Directions

### Short-term Enhancements (v3.2-3.3)

#### **FASTA Processing Optimization** *(High Priority)*
- **Pathogen Detection**: Complete optimization of assembly-aware pathogen screening
- **AMR Analysis**: Enhanced antimicrobial resistance detection for longer sequences
- **Virulence Analysis**: Improved virulence factor detection accuracy
- **Clinical Validation**: Comprehensive benchmarking against known datasets

#### **Performance Optimization**
- **GPU Acceleration**: CUDA-enabled alignment for DIAMOND searches
- **Memory Management**: Streaming algorithms for large dataset processing
- **Parallel Processing**: Multi-sample batch processing with job queuing
- **Database Compression**: Advanced indexing for reduced memory footprint

#### **Extended Analysis Capabilities**
- **Plasmid Detection**: Identification and characterization of mobile genetic elements
- **Prophage Analysis**: Integrated viral sequence detection and annotation
- **Metabolic Profiling**: KEGG pathway reconstruction and completeness scoring
- **Strain-level Resolution**: Sub-species identification using marker genes

#### **Enhanced Visualization**
- **3D Network Graphs**: Interactive microbial interaction networks
- **Time-series Analysis**: Longitudinal microbiome dynamics visualization
- **Comparative Analysis**: Multi-sample comparison dashboards
- **Publication-ready Figures**: Export to vector formats (SVG, PDF)

### Medium-term Developments (v4.0-4.5)

#### **Machine Learning Integration**
- **Pathogenicity Prediction**: ML models for novel pathogen risk assessment
- **Resistance Prediction**: Phenotype prediction from genotype data
- **Community Classification**: Microbiome state classification and clustering
- **Anomaly Detection**: Automated identification of unusual microbial signatures

#### **Advanced Genomics Features**
- **Pangenome Analysis**: Core and accessory genome characterization
- **Horizontal Gene Transfer**: Detection and visualization of HGT events
- **Evolutionary Analysis**: Phylogenetic reconstruction and molecular evolution
- **Functional Redundancy**: Assessment of metabolic pathway robustness

#### **Clinical Decision Support**
- **Treatment Recommendations**: Evidence-based therapy suggestions
- **Resistance Prediction**: Clinical breakpoint integration
- **Risk Stratification**: Patient-specific pathogen risk assessment
- **Outbreak Detection**: Automated surveillance and alerting systems

### Long-term Vision (v5.0+)

#### **Real-time Analysis Platform**
- **Streaming Analysis**: Live analysis of MinION/GridION sequencing data
- **Cloud Integration**: Scalable cloud-based processing infrastructure
- **Mobile Interface**: Tablet/smartphone access for field applications
- **API Development**: RESTful API for integration with LIMS systems

#### **Multi-omics Integration**
- **Proteomics**: Integration with metaproteomic analysis pipelines
- **Metabolomics**: Correlation with metabolic profiling data
- **Transcriptomics**: RNA-seq integration for activity assessment
- **Host-Microbe Interactions**: Integrated host genomics analysis

#### **Artificial Intelligence**
- **Natural Language Processing**: Automated literature mining for pathogen information
- **Computer Vision**: Integration of microscopy and imaging data
- **Predictive Modeling**: Population-level outbreak prediction
- **Knowledge Graphs**: Semantic integration of multi-source biological data

#### **Global Health Applications**
- **Surveillance Networks**: International pathogen monitoring platforms
- **Rapid Diagnostics**: Point-of-care sequencing integration
- **Antimicrobial Stewardship**: Population-level resistance management
- **One Health Integration**: Environmental-clinical-veterinary data fusion

## Contributing

We welcome contributions from the scientific community. Areas of particular interest include:

- **Algorithm Development**: New methods for metagenomic analysis
- **Database Curation**: Expansion and validation of reference databases
- **Visualization Tools**: Novel approaches for data presentation
- **Clinical Validation**: Real-world testing and benchmarking studies
- **Documentation**: User guides, tutorials, and best practices
- **FASTA Processing**: Help optimize pathogen detection and AMR analysis workflows

## Reporting Issues

If you encounter problems, please report them with:
- Input file format and size
- Error messages or unexpected output
- System specifications
- MetaQuest version information

**Note**: For FASTA files, taxonomic classification is now fully reliable. If you experience issues with pathogen detection, AMR analysis, or virulence factor identification, please note that these features are under active optimization.

## Contact and Support

- **Documentation**: [Installation Guide](docs/installation.md) | [Usage Guide](docs/usage.md)
- **Issues**: Report bugs and feature requests via GitHub issues
- **Community**: Join our discussion forum for questions and collaboration

---

**MetaQuest Development Team**  
*Advancing metagenomics through integrated computational solutions*

---
*Last updated: June 2025 - Development version with enhanced FASTA taxonomic classification and ongoing pathogen detection optimization*