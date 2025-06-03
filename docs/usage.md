# MetaQuest Usage Guide

This guide covers how to use MetaQuest for metagenomics analysis after installation is complete.

## Table of Contents
- [Quick Start](#quick-start)
- [Command Overview](#command-overview)
- [Usage Examples](#usage-examples)
- [Command Line Options](#command-line-options)
- [Input File Formats](#input-file-formats)
- [Output Files](#output-files)

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
├── bracken_report.tsv                   # Bracken abundance estimation (TSV)
├── bracken_report.txt                   # Bracken abundance estimation (TXT)
├── converted.fasta                      # Converted/processed FASTA
├── kraken_classified.txt                # Kraken2 classified reads
├── kraken_report.txt                    # Kraken2 report for reads
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
├── swissprot_annotation.tsv             # SwissProt functional annotations
└── taxonomy_overview.html               # Krona taxonomic visualization
```

#### FASTA Input Results
```
results/
├── amr_hits.txt                         # AMR hits (under development)
├── annotation_quality.html             # Annotation quality metrics
├── annotation_summary.txt               # Annotation summary
├── blast_cache/                         # BLAST cache directory
│   └── blast_cache.json                 # Cached BLAST results
├── blast_report.txt                     # BLAST analysis report
├── blast_taxonomy_results.json          # BLAST taxonomy results (JSON)
├── blast_taxonomy_summary.txt           # BLAST taxonomy summary
├── krona_input.txt                      # Krona visualization input
├── organism_comparison_data.csv          # Organism comparison (CSV)
├── organism_comparison_data.json         # Organism comparison (JSON)
├── pathogen_blast_results.txt           # Pathogen BLAST results (under development)
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
├── swissprot_annotation.tsv             # SwissProt functional annotations
├── taxonomy_krona.html                  # Krona taxonomic visualization
└── virulence_hits.txt                   # Virulence hits (under development)
```

### Key Output Files

#### Main Reports
- **taxonomy_krona.html**: Interactive Krona taxonomic visualization (FASTA)
- **taxonomy_overview.html**: Krona taxonomic visualization (FASTQ)

#### Taxonomic Results
- **blast_report.txt**: BLAST-based taxonomic classification report (FASTA)
- **blast_taxonomy_summary.txt**: Detailed BLAST taxonomic analysis summary (FASTA)
- **blast_taxonomy_results.json**: BLAST taxonomy results in JSON format (FASTA)
- **kraken_report.txt**: Kraken2 classification report (FASTQ)
- **kraken_classified.txt**: Kraken2 classified sequences (FASTQ)
- **bracken_report.tsv**: Bracken abundance estimation in TSV format (FASTQ)
- **bracken_report.txt**: Bracken abundance estimation in text format (FASTQ)
- **organism_comparison_data.csv**: Organism comparison data in CSV format (FASTA)
- **organism_comparison_data.json**: Organism comparison data in JSON format (FASTA)

#### Pathogen & Resistance Analysis (Under Development)
- **pathogen_blast_results.txt**: Pathogen BLAST search results (FASTA)
- **amr_hits.txt**: Antimicrobial resistance gene hits (FASTA)
- **virulence_hits.txt**: Virulence factor identifications (FASTA)

#### Functional Annotation
- **prokka_annotation/**: Complete Prokka gene annotation results (Available for both FASTA and FASTQ inputs)
  - **sample.txt**: Annotation summary with organism info and statistics
  - **sample.tsv**: Tab-separated gene annotations with detailed functional information
  - **sample.gff**: Gene feature format file
  - **sample.faa**: Protein sequences
  - **sample.ffn**: Gene sequences
  - **sample.gbk**: GenBank format annotation
- **swissprot_annotation.tsv**: SwissProt functional annotations
- **annotation_quality.html**: Annotation quality metrics (FASTA)
- **annotation_summary.txt**: Annotation summary report (FASTA)

#### Quality & Statistics
- **converted.fasta**: Converted/processed FASTA sequences (FASTQ)
- **krona_input.txt**: Input file for Krona visualization (FASTA)

#### Cache and Support Files
- **blast_cache/**: BLAST results caching directory (FASTA)
  - **blast_cache.json**: Cached BLAST results for faster re-analysis

### Output Format Details

#### BLAST Taxonomic Classification (FASTA Input)
```txt
# blast_report.txt format
0.00	0	0	U	0	unclassified
50.00	10	1	S	0	Salmonella enterica
30.00	6	1	S	0	Escherichia coli
10.00	2	1	S	0	Klebsiella pneumoniae
10.00	2	1	S	0	Cloning vector
```

#### BLAST Taxonomy Summary (FASTA Input)
```txt
# blast_taxonomy_summary.txt format
BLAST TAXONOMIC CLASSIFICATION SUMMARY
==================================================

Total sequences analyzed: 1
Sequences with hits: 1
Total BLAST hits: 20
Unique organisms identified: 4

TOP ORGANISMS BY TOTAL HITS:
----------------------------------------
Organism                       Total Hits Sequences  Avg Hits/Seq
----------------------------------------
Salmonella enterica            10         1          10.0        
Escherichia coli               6          1          6.0         
Klebsiella pneumoniae          2          1          2.0         
Cloning vector                 2          1          2.0
```

#### Organism Comparison Data (FASTA Input)
```csv
# organism_comparison_data.csv format
organism,total_hits,sequences_with_hits,avg_hits_per_sequence
Salmonella enterica,10,1,10.0
Escherichia coli,6,1,6.0
Klebsiella pneumoniae,2,1,2.0
Cloning vector,2,1,2.0
```

#### Bracken Abundance Report (FASTQ Input)
```tsv
# bracken_report.tsv format
name	taxonomy_id	taxonomy_lvl	kraken_assigned_reads	added_reads	new_est_reads	fraction_total_reads
Marinilactibacillus sp. 15R	1911586	S	120	4406	4526	0.24906
Paucilactobacillus nenjiangensis	1296540	S	41	573	614	0.03379
Amylolactobacillus amylophilus	1603	S	37	4840	4877	0.26838
Aerococcus urinae	1376	S	48	345	393	0.02163
Suicoccus acidiformans	2036206	S	10	283	293	0.01612
Tetragenococcus osmophilus	526944	S	11	596	607	0.03340
Peribacillus psychrosaccharolyticus	1407	S	69	308	377	0.02075
Peribacillus butanolivorans	421767	S	17	65	82	0.00451
Bacillus thuringiensis	1428	S	12	70	82	0.00451
Bacillus velezensis	492670	S	10	464	474	0.02608
Jeotgalicoccus saudimassiliensis	1461582	S	11	25	36	0.00198
Finegoldia magna	1260	S	3884	639	4523	0.24890
Anaerococcus mediterraneensis	1870984	S	28	8	36	0.00198
Gudongella oleilytica	1582259	S	707	18	725	0.03990
Mycoplasma sp. (ex Biomphalaria glabrata)	1749074	S	499	21	520	0.02862
```

#### Prokka Annotation Results (Both FASTA and FASTQ Input)
```tsv
# sample.tsv format (Prokka functional annotations)
locus_tag	ftype	length_bp	gene	EC_number	COG	product
KIEHJLIC_00001	CDS	555	rdmC	3.1.1.95		Aclacinomycin methylesterase RdmC
KIEHJLIC_00002	CDS	285				hypothetical protein
KIEHJLIC_00003	CDS	2130	mobA			Mobilization protein A
KIEHJLIC_00004	CDS	213				hypothetical protein
KIEHJLIC_00005	CDS	207				hypothetical protein
KIEHJLIC_00006	CDS	840	repA			Regulatory protein RepA
KIEHJLIC_00007	CDS	852				hypothetical protein
```

#### Prokka Annotation Summary
```txt
# sample.txt format
organism: Genus species strain 
contigs: 19643
bases: 8479628
CDS: 44150
```

### Getting Help

For issues not covered in this guide, please check the project repository or contact the development team.

This completes the MetaQuest usage guide. For installation instructions, see the [Installation Guide](installation.md).