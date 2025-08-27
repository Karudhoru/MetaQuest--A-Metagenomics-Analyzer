# MetaQuest Usage Guide

**Version 3.3.0** | Metagenomics Analysis Pipeline

---

## 📋 Table of Contents

- [Overview](#overview)
- [Quick Start](#quick-start)
- [Command Reference](#command-reference)
- [Individual Commands](#individual-commands)
  - [analyze](#analyze-command)
  - [validate](#validate-command)
  - [compare](#compare-command)
  - [check](#check-command)
- [Output Structure](#output-structure)
- [Examples & Use Cases](#examples--use-cases)
- [Troubleshooting](#troubleshooting)

---

## 🔬 Overview

MetaQuest provides **analysis-specific workflows** optimized for different input types and use cases:

| Analysis Type | Optimized For | Key Features |
|---------------|---------------|-------------|
| **FASTQ Analysis** | Clinical samples | Rapid Kraken2/Bracken classification, pathogen screening |
| **FASTA Analysis** | Research applications | High-accuracy BLAST classification, ML pathogen prediction |
| **Comparative Analysis** | Multi-sample studies | Statistical comparison, differential abundance, beta diversity |
| **File Validation** | Quality control | Comprehensive validation with detailed statistics |

---

## ⚡ Quick Start

### 1. System Check (Recommended)
```bash
metaquest check
```

### 2. File Validation (Recommended)
```bash
# Single FASTQ file
metaquest validate fastq --single sample.fastq.gz

# FASTA file
metaquest validate fasta genome.fasta
```

### 3. Run Analysis
```bash
# FASTQ analysis
metaquest analyze fastq --single sample.fastq.gz -o results/

# FASTA analysis
metaquest analyze fasta genome.fasta -o results/
```

### 4. Compare Multiple Samples
```bash
metaquest compare -i sample1_results/ sample2_results/ sample3_results/ -m metadata.tsv -o comparison/
```

---

## 📖 Command Reference

### Main Commands

| Command | Purpose | Input Types |
|---------|---------|-------------|
| `metaquest analyze` | Full analysis pipeline | FASTQ, FASTA |
| `metaquest validate` | File validation only | FASTQ, FASTA |
| `metaquest compare` | Multi-sample comparison | MetaQuest results |
| `metaquest check` | System dependencies | None |
| `metaquest --version` | Version information | None |
| `metaquest --help` | General help | None |

### Global Options

```bash
-v, --version    Display version information
-h, --help      Show help message
```

---

## 🔧 Individual Commands

## `analyze` Command

**Purpose**: Run the complete MetaQuest analysis pipeline on a single sample.

### Command Syntax
```bash
metaquest analyze <type> <input_files> [options]
```

### Subcommands

#### FASTA Analysis
```bash
metaquest analyze fasta <input_file> [options]
```

**Arguments:**
- `<input_file>` *(required)*: Path to input FASTA file

**Options:**
- `-o, --output <dir>`: Output directory (default: `results`)
- `-s, --blast-sample-size <int>`: Number of sequences to BLAST (default: `50`)
- `--skip-validation`: Skip pre-analysis file validation
- `--skip-annotation`: Skip functional annotation for faster analysis

**Example:**
```bash
metaquest analyze fasta genome.fasta -o fasta_results/ -s 100
```

**Help Output:**
```bash
$ metaquest analyze fasta --help
usage: metaquest analyze fasta [-h] [-o OUTPUT] [--skip-validation] [--skip-annotation] 
                              [-s BLAST_SAMPLE_SIZE] input_file

Analyze a single FASTA file.

positional arguments:
  input_file            Path to the input FASTA file.

options:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output directory name (default: results).
  --skip-validation     Skip input file validation.
  --skip-annotation     Skip functional and pathogen annotation steps for a faster 
                        taxonomic-only analysis.
  -s BLAST_SAMPLE_SIZE, --blast-sample-size BLAST_SAMPLE_SIZE
                        Number of sequences to BLAST (default: 50).
```

#### FASTQ Analysis
```bash
metaquest analyze fastq <mode> [options]
```

**Modes** (choose one):
- `--single <file>`: Single-end reads
- `--paired <r1> <r2>`: Paired-end reads  
- `--interleaved <file>`: Interleaved reads

**Options:**
- `-o, --output <dir>`: Output directory (default: `results`)
- `--skip-validation`: Skip pre-analysis file validation
- `--skip-annotation`: Skip functional annotation for faster analysis

**FASTQ Validation Options:**
- `-q, --min-quality <int>`: Minimum mean quality score (default: `20`)
- `-n, --min-sequences <int>`: Minimum sequence count (default: `100`)
- `--overrep-threshold <float>`: Overrepresentation threshold % (default: `0.1`)

**Examples:**
```bash
# Single-end analysis
metaquest analyze fastq --single reads.fastq.gz -o single_results/

# Paired-end analysis with custom quality threshold
metaquest analyze fastq --paired R1.fastq.gz R2.fastq.gz -q 25 -o paired_results/

# Interleaved analysis skipping validation
metaquest analyze fastq --interleaved interleaved.fastq.gz --skip-validation -o fast_results/
```

**Help Output:**
```bash
$ metaquest analyze fastq --help
usage: metaquest analyze fastq [-h] [-o OUTPUT] [--skip-validation] [--skip-annotation] 
                              [-q MIN_QUALITY] [-n MIN_SEQUENCES] 
                              [--overrep-threshold OVERREP_THRESHOLD] 
                              (--single READS.fastq | --paired R1.fastq R2.fastq | --interleaved INTERLEAVED.fastq)

Analyze FASTQ files.

options:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output directory name (default: results).
  --skip-validation     Skip input file validation.
  --skip-annotation     Skip functional and pathogen annotation steps for a faster 
                        taxonomic-only analysis.
  --single READS.fastq  Single-end FASTQ file.
  --paired R1.fastq R2.fastq
                        Paired-end FASTQ files.
  --interleaved INTERLEAVED.fastq
                        Interleaved paired-end FASTQ file.

FASTQ Validation Options:
  -q MIN_QUALITY, --min-quality MIN_QUALITY
                        Minimum mean quality score (default: 20).
  -n MIN_SEQUENCES, --min-sequences MIN_SEQUENCES
                        Minimum number of sequences (default: 100).
  --overrep-threshold OVERREP_THRESHOLD
                        Percentage threshold to flag a sequence as overrepresented 
                        (default: 0.1).
```

---

## `validate` Command

**Purpose**: Validate input files without running the full analysis pipeline.

### Command Syntax
```bash
metaquest validate <type> <input_files> [options]
```

### Subcommands

#### FASTA Validation
```bash
metaquest validate fasta <input_file>
```

**Arguments:**
- `<input_file>` *(required)*: Path to input FASTA file

**Example:**
```bash
metaquest validate fasta genome.fasta
```

#### FASTQ Validation
```bash
metaquest validate fastq <mode> [options]
```

**Modes** (choose one):
- `--single <file>`: Single-end reads
- `--paired <r1> <r2>`: Paired-end reads
- `--interleaved <file>`: Interleaved reads

**Options:**
- `-q, --min-quality <int>`: Minimum mean quality score (default: `20`)
- `-n, --min-sequences <int>`: Minimum sequence count (default: `100`)
- `--overrep-threshold <float>`: Overrepresentation threshold % (default: `0.1`)

**Examples:**
```bash
# Basic validation
metaquest validate fastq --single reads.fastq.gz

# Custom thresholds
metaquest validate fastq --paired R1.fastq R2.fastq -q 30 --overrep-threshold 0.2

# Interleaved validation
metaquest validate fastq --interleaved reads.fastq.gz -n 500
```

**Help Output:**
```bash
$ metaquest validate fastq --help
usage: metaquest validate fastq [-h] [-q MIN_QUALITY] [-n MIN_SEQUENCES] 
                                [--overrep-threshold OVERREP_THRESHOLD] 
                                (--single READS.fastq | --paired R1.fastq R2.fastq | --interleaved INTERLEAVED.fastq)

Validate FASTQ files.

options:
  -h, --help            show this help message and exit
  --single READS.fastq  Single-end FASTQ file.
  --paired R1.fastq R2.fastq
                        Paired-end FASTQ files (R1 and R2).
  --interleaved INTERLEAVED.fastq
                        Interleaved paired-end FASTQ file.

FASTQ Validation Options:
  -q MIN_QUALITY, --min-quality MIN_QUALITY
                        Minimum mean quality score (default: 20).
  -n MIN_SEQUENCES, --min-sequences MIN_SEQUENCES
                        Minimum number of sequences (default: 100).
  --overrep-threshold OVERREP_THRESHOLD
                        Percentage threshold to flag a sequence as overrepresented 
                        (default: 0.1).
```

---

## `compare` Command

**Purpose**: Perform comparative analysis across multiple MetaQuest results to identify differential abundance and community patterns.

### Command Syntax
```bash
metaquest compare -i <dir1> <dir2> ... -m <metadata.tsv> -o <output_dir>
```

**Arguments:**
- `-i, --inputs` *(required)*: Space-separated list of MetaQuest output directories
- `-m, --metadata` *(required)*: Path to metadata file (TSV format)
- `-o, --output`: Comparison results directory (default: `comparison_results`)

**Examples:**
```bash
# Basic comparison
metaquest compare -i sample1/ sample2/ sample3/ -m metadata.tsv -o comparison/

# Multiple samples with custom output
metaquest compare -i healthy_1/ healthy_2/ diseased_1/ diseased_2/ -m clinical_metadata.tsv -o clinical_comparison/
```

**Help Output:**
```bash
$ metaquest compare --help
usage: metaquest compare [-h] -i INPUTS [INPUTS ...] -m METADATA [-o OUTPUT]

Perform comparative analysis across multiple samples.

options:
  -h, --help            show this help message and exit
  -i INPUTS [INPUTS ...], --inputs INPUTS [INPUTS ...]
                        Space-separated list of MetaQuest output directories to compare.
  -m METADATA, --metadata METADATA
                        Path to the metadata file (TSV) linking samples to groups.
  -o OUTPUT, --output OUTPUT
                        Directory to save comparison results (default: comparison_results).
```

### Metadata File Format

**Required columns:**
- `sample_id`: Must exactly match input directory names
- `group`: Group assignment for statistical comparison

**Example `metadata.tsv`:**
```tsv
sample_id	group	additional_info
healthy_sample1	Healthy	Control group
healthy_sample2	Healthy	Control group  
diseased_sample1	Diseased	Treatment group
diseased_sample2	Diseased	Treatment group
```

---

## `check` Command

**Purpose**: Verify all system dependencies and database installations.

### Command Syntax
```bash
metaquest check
```

**Example:**
```bash
metaquest check
```

**Sample Output:**
```
🧬 MetaQuest v3.3.0 | A Comprehensive Metagenomics Analysis Pipeline

Performing system-wide checks...
  -> Checking command-line tools...
  -> Checking Python packages...
  -> Checking ML model artifacts...
  -> Checking database files...

==================================================
          SYSTEM CHECK COMPLETE
==================================================

✅  Success! All dependencies and databases are correctly configured.
```

---

## 📁 Output Structure

### FASTQ Analysis Results
*Clinical-focused pipeline outputs*

```
results/
├── analysis_dashboard.html              # Interactive dashboard (green theme)
├── pathogen_summary.txt                 # Primary pathogen report
├── pathogen_results.txt                 # Detailed pathogen detection results
├── pathogen_risk_detection.html         # Interactive pathogen risk visualization
├── taxonomic_classification_report.txt  # Kraken2/Bracken classification results
├── taxonomic_classification_report.json # Taxonomic data in JSON format
├── taxonomy_krona.html                  # Interactive Krona taxonomy chart
├── taxonomic_abundance_chart.html       # Species abundance visualization
├── bracken_report.tsv                   # Bracken abundance estimation
├── bracken_report.txt                   # Bracken results (text format)
├── kraken_report.txt                    # Kraken2 classification report
├── kraken_classified.txt                # Classified sequence details
├── converted.fasta                      # FASTQ to FASTA conversion
├── functional_annotation_report.txt     # Functional analysis summary
├── functional_annotation_report.json    # Functional data in JSON format
├── swissprot_annotation.tsv             # SwissProt protein annotations
├── protein_length_distribution.html     # Protein length analysis
├── detection_method_coverage.html       # Detection method comparison
├── pathogen_detection_report.json       # Pathogen detection data
└── prokka_annotation/                   # Gene prediction results
    ├── sample.gff                       # Gene feature format
    ├── sample.gbk                       # GenBank format
    ├── sample.faa                       # Amino acid sequences
    ├── sample.ffn                       # Nucleotide sequences (genes)
    ├── sample.tsv                       # Annotation table
    └── [additional Prokka outputs]
```

### FASTA Analysis Results  
*Research-focused pipeline outputs*

```
results/
├── analysis_dashboard.html              # Interactive dashboard (blue theme)
├── blast_ml_pathogen_summary.txt        # Combined BLAST+ML pathogen report
├── blast_ml_pathogen_summary.json       # Pathogen summary in JSON format
├── blast_ml_integrated_pathogen_report.json # Integrated analysis results
├── ml_pathogen_predictions.csv          # Individual ML predictions
├── ml_pathogen_predictions.json         # ML predictions in JSON format
├── ml_pathogen_summary.json             # ML analysis summary
├── organism_comparison_data.csv         # Species comparison table
├── organism_comparison_data.json        # Comparison data in JSON format
├── blast_taxonomy_results.json          # BLAST taxonomy results
├── taxonomic_classification_report.txt  # BLAST taxonomy summary
├── taxonomic_classification_report.json # Taxonomic data in JSON format
├── taxonomy_krona.html                  # Interactive Krona taxonomy chart
├── taxonomic_abundance_chart.html       # Species abundance visualization
├── pathogen_risk_detection.html         # Interactive pathogen risk analysis
├── detection_method_coverage.html       # Method comparison visualization
├── input_statistics.json               # Input file statistics
├── validation.json                      # File validation results
├── functional_annotation_report.txt     # Functional analysis summary
├── functional_annotation_report.json    # Functional data in JSON format
├── swissprot_annotation.tsv             # SwissProt protein annotations
├── protein_length_distribution.html     # Protein length analysis
├── blast_cache/                         # BLAST result caching
│   └── blast_cache.json                 # Cached BLAST results
└── prokka_annotation/                   # Gene prediction results
    ├── sample.gff                       # Gene feature format
    ├── sample.gbk                       # GenBank format
    ├── sample.faa                       # Amino acid sequences
    ├── sample.ffn                       # Nucleotide sequences (genes)
    ├── sample.tsv                       # Annotation table
    └── [additional Prokka outputs]
```

### Comparative Analysis Results

```
comparison_results/
├── taxonomic_abundance_table.tsv        # Species abundance matrix across all samples
├── differential_abundance_report.tsv    # Statistically significant species differences
├── alpha_diversity_results.tsv          # Alpha diversity metrics for each sample
├── alpha_diversity_boxplot.html         # Alpha diversity group comparison
├── beta_diversity_pcoa.html             # Interactive PCoA plot
├── taxonomic_heatmap.html               # Interactive species abundance heatmap
└── volcano_plot.html                    # Differential abundance volcano plot
```

---

## 💡 Examples & Use Cases

### Clinical Research Workflow

```bash
# 1. Validate clinical samples
metaquest validate fastq --paired patient1_R1.fastq.gz patient1_R2.fastq.gz -q 25

# 2. Analyze individual samples
metaquest analyze fastq --paired patient1_R1.fastq.gz patient1_R2.fastq.gz -o patient1_results/
metaquest analyze fastq --paired patient2_R1.fastq.gz patient2_R2.fastq.gz -o patient2_results/
metaquest analyze fastq --paired control1_R1.fastq.gz control1_R2.fastq.gz -o control1_results/

# 3. Compare groups
metaquest compare -i patient1_results/ patient2_results/ control1_results/ -m clinical_metadata.tsv -o clinical_comparison/
```

### Environmental Metagenomics

```bash
# 1. Analyze environmental samples
metaquest analyze fastq --single soil_sample1.fastq.gz -o soil1_results/
metaquest analyze fastq --single water_sample1.fastq.gz -o water1_results/

# 2. Compare environments
metaquest compare -i soil1_results/ water1_results/ -m environment_metadata.tsv -o environment_comparison/
```

### Genome Analysis Workflow

```bash
# 1. Validate genome assembly
metaquest validate fasta assembled_genome.fasta

# 2. Full analysis with increased BLAST sampling
metaquest analyze fasta assembled_genome.fasta -s 200 -o genome_results/

# 3. Quick taxonomic-only analysis
metaquest analyze fasta genome.fasta --skip-annotation -o quick_results/
```

### Batch Processing Example

```bash
#!/bin/bash
# Process multiple samples

samples=("sample1" "sample2" "sample3")
for sample in "${samples[@]}"; do
    echo "Processing $sample..."
    metaquest analyze fastq --single ${sample}.fastq.gz -o ${sample}_results/
done

# Compare all results
metaquest compare -i sample1_results/ sample2_results/ sample3_results/ -m batch_metadata.tsv -o batch_comparison/
```

---

## 🔧 Troubleshooting

### Common Issues

#### File Not Found Errors
```bash
# Error: Input file not found: sample.fastq
# Solution: Check file path and permissions
ls -la sample.fastq
metaquest validate fastq --single sample.fastq.gz  # Use correct extension
```

#### Validation Failures
```bash
# Error: File validation failed
# Solution: Check validation thresholds
metaquest validate fastq --single sample.fastq.gz -q 15  # Lower quality threshold
metaquest validate fastq --single sample.fastq.gz --overrep-threshold 1.0  # Higher overrep threshold
```

#### Missing Dependencies
```bash
# Error: Command not found
# Solution: Run system check
metaquest check

# Install missing dependencies (see installation.md)
```

#### Insufficient Memory
```bash
# Error: Out of memory
# Solution: Reduce BLAST sample size for FASTA analysis
metaquest analyze fasta genome.fasta -s 25  # Reduce from default 50
```

### Sample Command Outputs

#### System Check Output
```bash
$ metaquest check
🧬 MetaQuest v3.3.0 | A Comprehensive Metagenomics Analysis Pipeline

Performing system-wide checks...
  -> Checking command-line tools...
  -> Checking Python packages...
  -> Checking ML model artifacts...
  -> Checking database files...

==================================================
          SYSTEM CHECK COMPLETE
==================================================
✅  Success! All dependencies and databases are correctly configured.
```

#### FASTQ Validation Output
```bash
$ metaquest validate fastq --interleaved gut.fastq.gz
🧬 MetaQuest v3.3.0 | A Comprehensive Metagenomics Analysis Pipeline

🧬 MetaQuest - File Validation

--- Validating file 1/1: gut.fastq.gz ---

============================================================
🔍 METAQUEST FILE VALIDATION
============================================================
File: gut.fastq.gz
Type: FASTQ
Time: 2025-08-27 14:08:30
============================================================

⏳ Analyzing file contents...
  Processing: Done!

📁 FILE INFORMATION
────────────────────────────────────────
  Size:        5.33 MB
  MD5:         212808449d188c0b...
  Compressed:  Yes (gzip)

🔬 CONTAMINATION & DUPLICATION
────────────────────────────────────────
  Adapter Content:    0.00% of reads sampled

  Flagged as Overrepresented (>0.10%) - Top 10:
  Sequence                                 Percentage
  ---------------------------------------- ----------
  ACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAA... 13.8187%
  ACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAA... 9.5217%
  ACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAA... 5.7296%
  ACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAA... 2.4509%
  ACTCCTACGGGAGGCAGCAGTAGGGAATCTTCGGCAA... 1.4007%

  RECOMMENDATIONS:
  - ℹ️ Overrepresented sequences found. This may indicate PCR duplication.
    Consider investigating or using a deduplication tool.

🧬 SEQUENCE STATISTICS
────────────────────────────────────────
  Total Sequences:  235,304
  Total Bases:      58,120,088
  Length Range:     243 - 251 bp
  Mean Length:      247.0 bp
  Median Length:    247.0 bp
  N50 Length:       251.0 bp
  GC Content:       52.4%

📊 QUALITY STATISTICS
────────────────────────────────────────
  Encoding:         Sanger/Illumina 1.8+ (Phred+33)
  Mean Quality:     30.0
  Quality Range:    30 - 30
  Q20 Bases:        100.0%
  Q30 Bases:        100.0%

🔍 QUALITY ASSESSMENT
────────────────────────────────────────
  Note: All bases have uniform quality (30)
  ✅ All quality metrics passed

============================================================
⚠️ FILE VALIDATION: PASSED WITH WARNINGS
   The file is technically valid, but please review the 
   recommendations before proceeding with analysis.
============================================================

✅ File validation successful. Ready for analysis.
```

#### FASTA Validation Output
```bash
$ metaquest validate fasta sequence.fasta
🧬 MetaQuest v3.3.0 | A Comprehensive Metagenomics Analysis Pipeline

🧬 MetaQuest - File Validation

--- Validating file 1/1: sequence.fasta ---

============================================================
🔍 METAQUEST FILE VALIDATION
============================================================
File: sequence.fasta
Type: FASTA
Time: 2025-08-27 14:09:18
============================================================

⏳ Analyzing file contents...
  Processing: . Done!

📁 FILE INFORMATION
────────────────────────────────────────
  Size:        0.01 MB
  MD5:         42e1352c3ae16582...
  Compressed:  No

🧬 SEQUENCE STATISTICS
────────────────────────────────────────
  Type:             Single sequence (likely genome/chromosome)
  Total Sequences:  1
  Total Bases:      6,477
  Length Range:     6,477 - 6,477
  Mean Length:      6477.0
  Median Length:    6477.0

🧪 COMPOSITION ANALYSIS
────────────────────────────────────────
  GC Content:       61.2%

🔍 QUALITY ASSESSMENT
────────────────────────────────────────
  ✅ All quality checks passed

============================================================
✅ FILE VALIDATION: PASSED
============================================================

✅ File validation successful. Ready for analysis.
```

#### Comparative Analysis Output
```bash
$ metaquest compare -i results/healthy* results/diseased* -m metadata.tsv -o final_comparison
🧬 MetaQuest v3.3.0 | A Comprehensive Metagenomics Analysis Pipeline

🔬 Initializing Comparative Analysis...
-> Loading metadata from: metadata.tsv
-> Found 6 samples in metadata.
-> Aggregating taxonomic profiles...
  ✓ Processed: healthy_1
  ✓ Processed: healthy_2
  ✓ Processed: healthy_3
  ✓ Processed: diseased_1
  ✓ Processed: diseased_2
  ✓ Processed: diseased_3
-> Cleaning and validating abundance table...
✅ Cleaned abundance table created with 819 species across 6 samples.
-> Calculating Alpha Diversity...
  -> Alpha diversity calculations complete.
-> Calculating Beta Diversity (Bray-Curtis)...
-> Running statistical tests for differential abundance...
  ✓ Found 40 significantly different species (p < 0.1).
    -> Full report saved to: final_comparison/differential_abundance_report.tsv
-> Performing Principal Coordinate Analysis (PCoA)...
-> Generating Visualizations...
  ✓ PCoA plot created: beta_diversity_pcoa.html
  ✓ Heatmap visualization created: taxonomic_heatmap.html
  ✓ Alpha diversity box plot created: alpha_diversity_boxplot.html
  ✓ Volcano plot created: volcano_plot.html
```

### Getting Help

```bash
# General help
metaquest --help

# Command-specific help
metaquest analyze --help
metaquest analyze fastq --help
metaquest validate fasta --help
metaquest compare --help
```

### Support Resources

- **Installation Guide**: See `installation.md`
- **GitHub Issues**: Report bugs and request features
- **Documentation**: Check README.md for overview

---


*For technical support and updates, visit the MetaQuest GitHub repository.*
