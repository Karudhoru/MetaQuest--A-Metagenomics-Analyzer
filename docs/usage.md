# MetaQuest Usage Guide

This guide covers how to use MetaQuest for metagenomics analysis after installation is complete.

## Table of Contents
- [Quick Start](#quick-start)
- [Command Overview](#command-overview)
- [Usage Examples](#usage-examples)
- [Command Line Options](#command-line-options)
- [Input File Formats](#input-file-formats)
- [Output Files](#output-files)
- [Advanced Usage](#advanced-usage)
- [Performance Optimization](#performance-optimization)
- [Troubleshooting](#troubleshooting)

## Quick Start

### Basic Analysis

Activate your environment and run a complete metagenomics analysis:

```bash
# Activate environment
conda activate metagenomics_app

# Run complete analysis
metaquest examples/example.fastq.gz -t fastq -o results/ # for single-end FASTQ files

metaquest examples/sequence.fasta -t fasta -o results/ # for FASTA files
```

### View Results

```bash
# Open the main HTML report in your browser
firefox results/final_report.html

# View taxonomic visualizations
firefox results/taxonomy_krona.html
firefox results/taxonomy_pie.html
firefox results/taxonomy_treemap.html

# Check pathogen findings
less results/pathogen_summary.txt
less results/amr_hits.txt
less results/virulence_hits.txt

# Review functional annotations
less results/swissprot_annotation.tsv
firefox results/prokka_annotation/sample.gbk  # In genome browser
```

## Command Overview

MetaQuest provides several analysis modes:

```bash
# Show all available commands
metaquest --help

# Complete analysis (recommended for most users)
metaquest example.fastq.gz -t fastq -o results
metaquest examples/sequence.fasta -t fasta -o results
```

## Usage Examples

### Example 1: Single-end FASTQ Analysis

Analyze a single FASTQ file with default settings:

```bash
# Single-end sequencing data
metaquest sample_R1.fastq -t fastq -o results/
metaquest sample_R1.fastq.gz -t fastq -o results/

# Alternative syntax
metaquest --type fastq --reads sample_R1.fastq -o results/
metaquest --type fastq -r sample_R1.fastq.gz -o results/
```

### Example 2: Paired-end FASTQ Analysis

Analyze paired-end FASTQ files:

```bash
# Paired-end sequencing data
metaquest --type fastq --reads1 sample_R1.fastq --reads2 sample_R2.fastq -o results/
metaquest --type fastq -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -o results/

# With additional options
metaquest \
    --type fastq \
    --reads1 /path/to/sample_R1.fastq.gz \
    --reads2 /path/to/sample_R2.fastq.gz \
    --output /path/to/results/
```

### Example 3: Interleaved FASTQ Analysis

Analyze interleaved paired-end data (both reads in single file):

```bash
# Interleaved paired-end data
metaquest --type fastq --interleaved sample_interleaved.fastq -o results/
metaquest --type fastq -i sample_interleaved.fastq.gz -o results/
```

### Example 4: FASTA Analysis

Analyze pre-assembled contigs or sequences:

```bash
# FASTA input (limited functionality - see development status below)
metaquest --type fasta sample.fasta -o results/
metaquest sample.fasta -t fasta -o results/
```

### Example 5: Functional Analysis of Assembled Contigs

Analyze pre-assembled contigs for functional content:

```bash
metaquest functional \
    /path/to/assembled_contigs.fasta \
    -t fasta \
    -o functional_results/
```

## Command Line Options

### Global Options

Available for all commands:

- `-h, --help`: Show help information
- `-t, --type`: Specify the file type {fasta/fastq} (required)
- `-o, --output`: Output directory (required)
- `--threads`: Number of threads to use (optional, default: auto-detect)

### FASTQ-specific Options

For FASTQ input files:

- `-r, --reads`: Single-end FASTQ file
- `-1, --reads1`: First paired-end FASTQ file (R1)
- `-2, --reads2`: Second paired-end FASTQ file (R2)
- `-i, --interleaved`: Interleaved paired-end FASTQ file

## Input File Formats

### Supported Input Formats

MetaQuest accepts the following input formats:

1. **FASTQ files** (recommended for raw sequencing data)
   - Single-end: `sample.fastq`, `sample.fastq.gz`
   - Paired-end: `sample_R1.fastq`, `sample_R2.fastq` (separate files)
   - Interleaved: `sample_interleaved.fastq` (both reads in one file)
 
2. **FASTA files** (for assembled sequences)
   - Contigs: `contigs.fasta`, `contigs.fasta.gz`
   - Gene sequences: `genes.fna`

### Quality Requirements

- **Minimum read length**: 35 bp (configurable)
- **Minimum base quality**: Q20 (configurable)
- **Supported encodings**: Illumina 1.8+ (Phred+33)

### File Naming Conventions

For paired-end data, common naming conventions are supported:
- `sample_R1.fastq` and `sample_R2.fastq`
- `sample_1.fastq` and `sample_2.fastq`
- `sample.1.fastq` and `sample.2.fastq`

## Output Files

MetaQuest generates a comprehensive set of output files:

### Directory Structure

#### FASTQ Input Results
```
results/
├── 3d_annotation.html                    # 3D visualization of annotations
├── amr_classes_distribution.html         # AMR class distribution plots
├── amr_hits.txt                         # Antimicrobial resistance hits
├── amr_identity_distribution.html        # AMR identity distribution plots
├── analysis_dashboard.html              # Analysis dashboard with key metrics
├── antimicrobial_resistance_report.json # AMR analysis results (JSON)
├── assembly_quality_report.html         # Assembly quality metrics (HTML)
├── assembly_quality_report.json         # Assembly quality metrics (JSON)
├── bracken_report.tsv                   # Bracken abundance estimation (TSV)
├── bracken_report.txt                   # Bracken abundance estimation (TXT)
├── converted.fasta                      # Converted/processed FASTA
├── fasta_kraken_classified.txt          # Kraken2 classified sequences
├── fasta_kraken_report.txt              # Kraken2 classification report
├── final_report.html                    # Main comprehensive report
├── gc_distribution.html                 # GC content distribution plots
├── kraken_classified.txt                # Kraken2 classified reads
├── kraken_report.txt                    # Kraken2 report for reads
├── krona_input.txt                      # Krona visualization input
├── length_distribution.html             # Sequence length distribution
├── length_summary.html                  # Length statistics summary
├── pathogen_blast_results.txt           # Pathogen BLAST search results
├── pathogen_summary.txt                 # Pathogen identification summary
├── prokka_annotation/                   # Prokka gene annotation results
│   ├── sample.err                       # Error log
│   ├── sample.faa                       # Protein sequences (FASTA)
│   ├── sample.ffn                       # Gene sequences (FASTA)
│   ├── sample.fna                       # Nucleotide sequences
│   ├── sample.fsa                       # Contig sequences
│   ├── sample.gbk                       # GenBank format
│   ├── sample.gff                       # Gene feature format
│   ├── sample.log                       # Annotation log
│   ├── sample.sqn                       # Sequin format
│   ├── sample.tbl                       # Feature table
│   ├── sample.tsv                       # Tab-separated annotations
│   └── sample.txt                       # Text summary
├── sequence_statistics.json             # Sequence statistics (JSON)
├── swissprot_annotation.tsv             # SwissProt functional annotations
├── swissprot_identity.html              # SwissProt identity distribution
├── taxonomy_krona.html                  # Krona taxonomic visualization
├── taxonomy_pie.html                    # Taxonomic pie chart
├── taxonomy_treemap.html                # Taxonomic treemap visualization
└── virulence_hits.txt                   # Virulence factor hits
```

#### FASTA Input Results
```
results/
├── 3d_annotation.html                    # 3D visualization of annotations
├── amr_classes_distribution.html         # AMR class distribution plots
├── amr_hits.txt                         # Antimicrobial resistance hits
├── amr_identity_distribution.html        # AMR identity distribution plots
├── antimicrobial_resistance_report.json # AMR analysis results (JSON)
├── assembly_quality_report.html         # Assembly quality metrics (HTML)
├── assembly_quality_report.json         # Assembly quality metrics (JSON)
├── fasta_kraken_classified.txt          # Kraken2 classified sequences
├── fasta_kraken_report.txt              # Kraken2 classification report
├── final_report.html                    # Main comprehensive report
├── gc_distribution.html                 # GC content distribution plots
├── krona_input.txt                      # Krona visualization input
├── length_distribution.html             # Sequence length distribution
├── length_summary.html                  # Length statistics summary
├── pathogen_blast_results.txt           # Pathogen BLAST search results
├── pathogen_summary.txt                 # Pathogen identification summary
├── prokka_annotation/                   # Prokka gene annotation results
│   ├── sample.err                       # Error log
│   ├── sample.faa                       # Protein sequences (FASTA)
│   ├── sample.ffn                       # Gene sequences (FASTA)
│   ├── sample.fna                       # Nucleotide sequences
│   ├── sample.fsa                       # Contig sequences
│   ├── sample.gbk                       # GenBank format
│   ├── sample.gff                       # Gene feature format
│   ├── sample.log                       # Annotation log
│   ├── sample.sqn                       # Sequin format
│   ├── sample.tbl                       # Feature table
│   ├── sample.tsv                       # Tab-separated annotations
│   └── sample.txt                       # Text summary
├── sequence_statistics.json             # Sequence statistics (JSON)
├── swissprot_annotation.tsv             # SwissProt functional annotations
├── swissprot_identity.html              # SwissProt identity distribution
├── taxonomy_krona.html                  # Krona taxonomic visualization
├── taxonomy_pie.html                    # Taxonomic pie chart
├── taxonomy_treemap.html                # Taxonomic treemap visualization
└── virulence_hits.txt                   # Virulence factor hits
```

### Key Output Files

#### Main Reports
- **final_report.html**: Comprehensive interactive HTML report with all results
- **analysis_dashboard.html**: Analysis dashboard with key metrics (FASTQ only)

#### Taxonomic Results
- **kraken_report.txt** / **fasta_kraken_report.txt**: Kraken2 classification reports
- **kraken_classified.txt** / **fasta_kraken_classified.txt**: Classified sequences
- **bracken_report.tsv** / **bracken_report.txt**: Bracken abundance estimation (FASTQ only)
- **taxonomy_krona.html**: Interactive Krona taxonomic visualization
- **taxonomy_pie.html**: Taxonomic composition pie chart
- **taxonomy_treemap.html**: Taxonomic treemap visualization

#### Pathogen & Resistance Analysis
- **pathogen_summary.txt**: Summary of pathogenic organism findings
- **pathogen_blast_results.txt**: Detailed pathogen BLAST search results
- **amr_hits.txt**: Antimicrobial resistance gene hits
- **virulence_hits.txt**: Virulence factor identifications
- **antimicrobial_resistance_report.json**: Structured AMR analysis results
- **amr_classes_distribution.html**: AMR class distribution visualization
- **amr_identity_distribution.html**: AMR identity distribution plots

#### Functional Annotation
- **prokka_annotation/**: Complete Prokka gene annotation results
  - **sample.tsv**: Tab-separated gene annotations
  - **sample.gff**: Gene feature format file
  - **sample.faa**: Protein sequences
  - **sample.ffn**: Gene sequences
  - **sample.gbk**: GenBank format annotation
- **swissprot_annotation.tsv**: SwissProt functional annotations
- **swissprot_identity.html**: SwissProt identity distribution
- **3d_annotation.html**: 3D visualization of annotations

#### Quality & Statistics
- **assembly_quality_report.html** / **assembly_quality_report.json**: Assembly quality metrics
- **sequence_statistics.json**: Comprehensive sequence statistics
- **gc_distribution.html**: GC content distribution analysis
- **length_distribution.html**: Sequence length distribution
- **length_summary.html**: Length statistics summary

### Output Format Details

#### Taxonomic Classification Results
```txt
# kraken_report.txt format
 15.23  1234    1234    S    562      Escherichia coli
  8.45   678     678    S    1280     Staphylococcus aureus
```

#### Pathogen Summary Format
```txt
# pathogen_summary.txt
Organism: Salmonella enterica
Confidence: High
Number of hits: 234
Risk assessment: Pathogenic
Classification: Bacterial pathogen
```

#### AMR Hits Format
```txt
# amr_hits.txt
Gene: blaTEM-1
Class: Beta-lactamase
Mechanism: Hydrolysis
Identity: 98.5%
Coverage: 100%
Organism: E. coli
```

#### Virulence Factors Format
```txt
# virulence_hits.txt
Factor: Shiga toxin
Type: Cytotoxin
Organism: E. coli O157:H7
Identity: 97.2%
Function: Cell damage
```


## Development Status and Future Directions

### Current Limitations

- **Paired-end reads**: Full analysis pipeline is implemented but may need optimization
- **Pathogenicity reports**: The reports for pathogenicity and virulence are still not refined and may fetch inaccurate results
- **Taxonomic Classification**: For FASTA analysis, the taxonomic report is inaccurate and is being improved
- **Integration with other tools**: Future versions will include integration with FastQC and Trimmomatic

### Planned Improvements

- Enhanced pathogen detection accuracy
- Improved taxonomic classification for FASTA inputs
- Better quality control integration
- Performance optimizations for large datasets
- Support for additional file formats

## Troubleshooting

### Common Issues

1. **Memory errors**: Reduce thread count or process smaller file chunks
2. **Incomplete results**: Check input file format and ensure sufficient disk space
3. **Slow performance**: Increase thread count and ensure adequate RAM

### Getting Help

For issues not covered in this guide, please check the project repository or contact the development team.

This completes the MetaQuest usage guide. For installation instructions, see the [Installation Guide](installation.md).