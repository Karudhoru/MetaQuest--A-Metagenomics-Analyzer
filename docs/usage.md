# MetaQuest Usage Guide

**Version 4.0.0** | Metagenomics Analysis Pipeline

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
- [Advanced Features](#advanced-features)
  - [Annotation Controls](#annotation-controls)
  - [Debug Mode](#debug-mode)
- [Output Structure](#output-structure)
- [Examples & Use Cases](#examples--use-cases)
- [Troubleshooting](#troubleshooting)

---

## 🔬 Overview

MetaQuest provides **analysis-specific workflows** optimized for different input types and use cases:

| Analysis Type | Optimized For | Key Features |
|---------------|---------------|-------------|
| **FASTQ Analysis** | Clinical samples | Rapid Kraken2/Bracken classification, pathogen screening, COG+SwissProt annotation |
| **FASTA Analysis** | Research applications | High-accuracy BLAST classification, ML pathogen prediction, dual-database annotation |
| **Comparative Analysis** | Multi-sample studies | Advanced statistical testing, ML biomarker discovery, publication-ready visualizations |
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
# FASTQ analysis with default settings
metaquest analyze fastq --single sample.fastq.gz -o results/

# FASTA analysis with enhanced annotation
metaquest analyze fasta genome.fasta -o results/

# Debug mode for troubleshooting
metaquest --debug analyze fastq --single sample.fastq.gz -o debug_results/
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
| `metaquest compare` | Multi-sample comparison with advanced statistics | MetaQuest results |
| `metaquest check` | System dependencies | None |
| `metaquest --version` | Version information | None |
| `metaquest --help` | General help | None |

### Global Options

```bash
--debug          Full diagnostic output (shows all tool commands and output)
-v, --version    Display version information
-h, --help       Show help message
```

---

## 🔧 Individual Commands

## `analyze` Command

**Purpose**: Run the complete MetaQuest analysis pipeline on a single sample.

### Command Syntax
```bash
metaquest [--debug] analyze <type> <input_files> [options]
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

**Annotation Control Options:**
- `--no-filter-contigs`: Disable contig filtering before Prokka annotation
- `--min-contig-length <int>`: Minimum contig length for filtering in bp (default: `1000`)
- `--no-kill-tbl2asn`: Disable auto-kill of long-running tbl2asn process
- `--tbl2asn-timeout <int>`: Timeout for tbl2asn process in seconds (default: `300`)
- `--annotation-threads <int>`: Number of threads for functional annotation (default: `8`)

**Examples:**
```bash
# Basic analysis with default settings
metaquest analyze fasta genome.fasta -o fasta_results/

# Custom BLAST sample size
metaquest analyze fasta genome.fasta -s 100 -o results/

# Annotate all contigs (no filtering)
metaquest analyze fasta assembly.fasta --no-filter-contigs -o all_contigs/

# Custom contig length threshold
metaquest analyze fasta assembly.fasta --min-contig-length 500 -o custom_filter/

# Extended timeout for large assemblies
metaquest analyze fasta large_assembly.fasta --tbl2asn-timeout 600 -o large_results/

# More threads for faster annotation
metaquest analyze fasta genome.fasta --annotation-threads 16 -o fast_results/

# Disable tbl2asn auto-kill (unlimited runtime)
metaquest analyze fasta genome.fasta --no-kill-tbl2asn -o unlimited/

# Debug mode for troubleshooting
metaquest --debug analyze fasta genome.fasta -o debug_results/
```

**Help Output:**
```bash
$ metaquest analyze fasta --help
usage: metaquest analyze fasta [-h] [-o DIR] [--skip-validation] [--skip-annotation]
                              [-s NUM] [--no-filter-contigs] [--min-contig-length BP]
                              [--no-kill-tbl2asn] [--tbl2asn-timeout SEC]
                              [--annotation-threads NUM] input_file

Analyze FASTA sequences

positional arguments:
  input_file            Path to the FASTA file

options:
  -h, --help            show this help message and exit
  -o DIR, --output DIR  Output directory (default: results)
  --skip-validation     Skip input validation (not recommended)
  --skip-annotation     Skip annotation for faster taxonomic-only analysis
  -s NUM, --blast-sample-size NUM
                        Number of sequences to BLAST (default: 50)

Annotation Control Options:
  Fine-tune the functional annotation process

  --no-filter-contigs   Skip contig filtering before Prokka annotation
  --min-contig-length BP
                        Minimum contig length for filtering (default: 1000)
  --no-kill-tbl2asn     Don't auto-kill long-running tbl2asn process
  --tbl2asn-timeout SEC
                        Timeout for tbl2asn process in seconds (default: 300)
  --annotation-threads NUM
                        Number of threads for functional annotation (default: 8)
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
- `--overrep-threshold <float>`: Overrepresentation threshold (default: `0.1`)

**Annotation Control Options:**
- `--no-filter-contigs`: Disable contig filtering before Prokka annotation
- `--min-contig-length <int>`: Minimum contig length for filtering in bp (default: `1000`)
- `--no-kill-tbl2asn`: Disable auto-kill of long-running tbl2asn process
- `--tbl2asn-timeout <int>`: Timeout for tbl2asn process in seconds (default: `300`)
- `--annotation-threads <int>`: Number of threads for functional annotation (default: `8`)

**Examples:**
```bash
# Single-end analysis with default settings
metaquest analyze fastq --single reads.fastq.gz -o single_results/

# Paired-end analysis with custom quality threshold
metaquest analyze fastq --paired R1.fastq.gz R2.fastq.gz -q 25 -o paired_results/

# Interleaved analysis skipping validation
metaquest analyze fastq --interleaved interleaved.fastq.gz --skip-validation -o fast_results/

# Custom contig filtering
metaquest analyze fastq --single reads.fq --min-contig-length 500 -o custom_results/

# More threads for annotation
metaquest analyze fastq --paired R1.fq R2.fq --annotation-threads 16 -o threaded_results/

# Extended tbl2asn timeout
metaquest analyze fastq --single reads.fq --tbl2asn-timeout 600 -o extended_results/

# Skip annotation for fast taxonomic-only analysis
metaquest analyze fastq --single reads.fq --skip-annotation -o taxonomy_only/

# Debug mode with full output
metaquest --debug analyze fastq --single reads.fq -o debug_results/
```

**Help Output:**
```bash
$ metaquest analyze fastq --help
usage: metaquest analyze fastq [-h] [-o DIR] [--skip-validation] [--skip-annotation]
                              [-q SCORE] [-n NUM] [--overrep-threshold PCT]
                              [--no-filter-contigs] [--min-contig-length BP]
                              [--no-kill-tbl2asn] [--tbl2asn-timeout SEC]
                              [--annotation-threads NUM]
                              (--single READS.fastq | --paired R1.fastq R2.fastq | 
                               --interleaved INTERLEAVED.fastq)

Analyze FASTQ sequences

options:
  -h, --help            show this help message and exit
  -o DIR, --output DIR  Output directory (default: results)
  --skip-validation     Skip input validation (not recommended)
  --skip-annotation     Skip annotation for faster taxonomic-only analysis
  --single READS.fastq  Single-end FASTQ file
  --paired R1.fastq R2.fastq
                        Paired-end FASTQ files
  --interleaved INTERLEAVED.fastq
                        Interleaved paired-end FASTQ file

FASTQ Validation Options:
  Configure quality control parameters for FASTQ files

  -q SCORE, --min-quality SCORE
                        Minimum mean quality score (default: 20)
  -n NUM, --min-sequences NUM
                        Minimum number of sequences (default: 100)
  --overrep-threshold PCT
                        Threshold to flag overrepresented sequences (default: 0.1)

Annotation Control Options:
  Fine-tune the functional annotation process

  --no-filter-contigs   Skip contig filtering before Prokka annotation
  --min-contig-length BP
                        Minimum contig length for filtering (default: 1000)
  --no-kill-tbl2asn     Don't auto-kill long-running tbl2asn process
  --tbl2asn-timeout SEC
                        Timeout for tbl2asn process in seconds (default: 300)
  --annotation-threads NUM
                        Number of threads for functional annotation (default: 8)
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
usage: metaquest validate fastq [-h] [-q SCORE] [-n NUM] [--overrep-threshold PCT]
                                (--single READS.fastq | --paired R1.fastq R2.fastq | 
                                 --interleaved INTERLEAVED.fastq)

Validate FASTQ file(s)

options:
  -h, --help            show this help message and exit
  --single READS.fastq  Single-end FASTQ file
  --paired R1.fastq R2.fastq
                        Paired-end FASTQ files (R1 and R2)
  --interleaved INTERLEAVED.fastq
                        Interleaved paired-end FASTQ file

FASTQ Validation Options:
  Configure quality control parameters for FASTQ files

  -q SCORE, --min-quality SCORE
                        Minimum mean quality score (default: 20)
  -n NUM, --min-sequences NUM
                        Minimum number of sequences (default: 100)
  --overrep-threshold PCT
                        Threshold to flag overrepresented sequences (default: 0.1)
```

---

## `compare` Command

**Purpose**: Perform advanced comparative analysis across multiple MetaQuest results with statistical testing, machine learning biomarker discovery, and publication-ready visualizations.

### Command Syntax
```bash
metaquest compare -i <dir1> <dir2> ... -m <metadata.tsv> -o <output_dir>
```

**Arguments:**
- `-i, --inputs` *(required)*: Space-separated list of MetaQuest output directories
- `-m, --metadata` *(required)*: Path to metadata file (TSV format)
- `-o, --output`: Comparison results directory (default: `comparison_results`)

**Statistical Analysis Features:**
- **Alpha Diversity**: Shannon, Simpson, Chao1, and Observed Species metrics with Mann-Whitney U tests
- **Beta Diversity**: Bray-Curtis dissimilarity with PERMANOVA and ANOSIM tests
- **Differential Abundance**: Multiple testing correction (FDR, Bonferroni) with volcano plot visualization
- **Machine Learning**: Random Forest biomarker discovery with cross-validation and feature importance

**Examples:**
```bash
# Basic comparison
metaquest compare -i sample1/ sample2/ sample3/ -m metadata.tsv -o comparison/

# Multiple samples with custom output
metaquest compare -i healthy_1/ healthy_2/ diseased_1/ diseased_2/ -m clinical_metadata.tsv -o clinical_comparison/

# Large cohort comparison (12+ samples)
metaquest compare -i healthy_*/  diseased_*/ -m large_study_metadata.tsv -o cohort_comparison/
```

**Help Output:**
```bash
$ metaquest compare --help
usage: metaquest compare [-h] -i DIR [DIR ...] -m FILE [-o DIR]

Perform comparative analysis across samples

options:
  -h, --help            show this help message and exit
  -i DIR [DIR ...], --inputs DIR [DIR ...]
                        MetaQuest output directories to compare
  -m FILE, --metadata FILE
                        Metadata file (TSV) linking samples to groups
  -o DIR, --output DIR  Output directory (default: comparison_results)
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

### Statistical Test Interpretations

**Alpha Diversity p-values:**
- p < 0.001: Highly significant (***)
- p < 0.01: Significant (**)
- p < 0.05: Significant (*)
- p < 0.1: Trend/marginally significant (†)

**ANOSIM R-values:**
- R > 0.75: Well separated groups
- R > 0.5: Overlapping but clearly different
- R > 0.25: Barely separable
- R ≈ 0: Indistinguishable groups

**Machine Learning Cross-validation:**
- >90% accuracy: Excellent separation
- 80-90% accuracy: Good separation
- 70-80% accuracy: Moderate separation
- <70% accuracy: Poor separation

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
```bash
🧬 MetaQuest v4.0.0 | An Integrated Metagenomics Analysis Pipeline

══════════════════════════════════════════════════════════════════════
                        SYSTEM CHECK
══════════════════════════════════════════════════════════════════════

Checking command-line tools...
  ✓ kraken2 found
  ✓ bracken found
  ✓ diamond found
  ✓ prokka found
  ✓ blastn found

Checking Python packages...
  ✓ numpy found
  ✓ pandas found
  ✓ scikit-learn found

Checking databases...
  ✓ Kraken2 database found
  ✓ COG database found
  ✓ SwissProt database found
  ✓ Pathogen database found

✓ System check completed successfully
```

---

## Advanced Features

### Annotation Controls

MetaQuest v4.0.0 introduces fine-grained control over the functional annotation process.

#### Contig Filtering

**Default Behavior:**
- Filters contigs shorter than 1000bp before annotation
- Reduces annotation time and focuses on high-quality sequences

**Custom Filtering:**
```bash
# Use 500bp threshold
metaquest analyze fastq --single reads.fq --min-contig-length 500

# Use 2000bp threshold for high-quality only
metaquest analyze fasta assembly.fasta --min-contig-length 2000
```

**Disable Filtering:**
```bash
# Annotate all contigs regardless of length
metaquest analyze fastq --single reads.fq --no-filter-contigs
```

**When to Use:**
- **Default (1000bp)**: Most cases, balanced speed and coverage
- **Higher threshold (2000bp+)**: High-quality genomes, faster annotation
- **Lower threshold (500bp)**: Fragmented assemblies, maximize coverage
- **No filtering**: Complete annotation needed, time is not a constraint

#### tbl2asn Timeout Management

**Default Behavior:**
- 300 second (5 minute) timeout
- Automatically terminates stuck tbl2asn processes

**Custom Timeout:**
```bash
# 10 minute timeout for large datasets
metaquest analyze fasta large_assembly.fasta --tbl2asn-timeout 600

# 30 minute timeout for very large assemblies
metaquest analyze fastq --paired R1.fq R2.fq --tbl2asn-timeout 1800
```

**Disable Auto-Kill:**
```bash
# Let tbl2asn run indefinitely
metaquest analyze fasta assembly.fasta --no-kill-tbl2asn
```

**When to Use:**
- **Default (300s)**: Most analyses, prevents hanging processes
- **Extended timeout (600-1800s)**: Large assemblies, complex genomes
- **Disabled auto-kill**: Critical runs where completion is essential

#### Threading Optimization

**Default Behavior:**
- Uses 8 threads for annotation

**Custom Threading:**
```bash
# Use 16 threads for faster annotation
metaquest analyze fastq --single reads.fq --annotation-threads 16

# Use 4 threads for resource-limited systems
metaquest analyze fasta genome.fasta --annotation-threads 4

# Use all available cores (example: 32 threads)
metaquest analyze fastq --paired R1.fq R2.fq --annotation-threads 32
```

**When to Use:**
- **More threads (16-32)**: Powerful servers, faster completion needed
- **Fewer threads (4-8)**: Shared systems, limited resources
- **Match system**: Set to number of available CPU cores

#### Combined Examples

```bash
# High-throughput mode: no filtering, more threads, extended timeout
metaquest analyze fastq --single reads.fq \
  --no-filter-contigs \
  --annotation-threads 32 \
  --tbl2asn-timeout 1800 \
  -o highthroughput_results/

# Conservative mode: strict filtering, default threads, short timeout
metaquest analyze fasta assembly.fasta \
  --min-contig-length 2000 \
  --annotation-threads 8 \
  --tbl2asn-timeout 300 \
  -o conservative_results/

# Quick mode: skip annotation entirely
metaquest analyze fastq --single reads.fq \
  --skip-annotation \
  -o quick_taxonomy/
```

### Debug Mode

**Purpose**: Comprehensive diagnostic output for troubleshooting and development.

#### Enabling Debug Mode

```bash
# Add --debug flag before any command
metaquest --debug analyze fastq --single reads.fq -o debug_results/
metaquest --debug validate fasta genome.fasta
metaquest --debug compare -i s1/ s2/ -m metadata.tsv -o comparison/
```

#### What Debug Mode Shows

**1. Command Invocations:**
```bash
[DEBUG] Executing: kraken2 --db databases/kraken2_db --threads 8 --output kraken_output.txt converted.fasta
```

**2. Tool Output:**
```bash
[DEBUG] Loading database... 100%
[DEBUG] Classifying sequences... 45%
[DEBUG] 18455 sequences classified (50.98%)
```

**3. Detailed Error Messages:**
```bash
[ERROR] Prokka annotation failed
[DEBUG] Command: prokka --outdir results/prokka --prefix sample converted.fasta
[DEBUG] Exit code: 1
[DEBUG] stderr: tbl2asn process hung after 300s
[DEBUG] Solution: Try increasing --tbl2asn-timeout or use --no-kill-tbl2asn
```

**4. Performance Metrics:**
```bash
[DEBUG] Kraken2 classification completed in 2m 34s
[DEBUG] Memory usage: 8.2 GB
[DEBUG] Disk usage: 1.4 GB
```

#### Standard vs Debug Output Comparison

**Standard Mode:**
```bash
══════════════════════════════════════════════════════════════════════
                       FASTQ Analysis
══════════════════════════════════════════════════════════════════════

⠋ Running Kraken2 classification...
✓ Kraken2 classification completed
⠋ Running functional annotation...
✓ Functional annotation completed

✓ Analysis completed in 5m 23s
```

**Debug Mode:**
```bash
══════════════════════════════════════════════════════════════════════
                       FASTQ Analysis
══════════════════════════════════════════════════════════════════════

[DEBUG] Starting Kraken2 classification
[DEBUG] Command: kraken2 --db databases/kraken2_db --threads 8 --output results/kraken_classified.txt converted.fasta
[DEBUG] kraken2: Loading database... 100%
[DEBUG] kraken2: Classifying sequences... 100%
[DEBUG] kraken2: 18455 sequences (50.98%) classified
✓ Kraken2 classification completed (2m 34s)

[DEBUG] Starting functional annotation
[DEBUG] Command: prokka --outdir results/prokka --prefix sample --cpus 8 converted.fasta
[DEBUG] prokka: [12:34:56] Contigs: 8579
[DEBUG] prokka: [12:35:23] Running: prodigal
[DEBUG] prokka: [12:36:45] Running: diamond
[DEBUG] prokka: [12:38:12] Found 1179 CDS, 32 tRNA, 12 rRNA
✓ Functional annotation completed (3m 18s)

✓ Analysis completed in 5m 52s
```

#### When to Use Debug Mode

- **Troubleshooting**: Understanding why an analysis failed
- **Development**: Testing new features or modifications
- **Performance tuning**: Identifying bottlenecks
- **Bug reports**: Providing detailed information to developers
- **Learning**: Understanding what MetaQuest does internally

---

## 📁 Output Structure

### FASTQ Analysis Results
*Clinical-focused pipeline outputs*

```
results/
├── metaquest.log                        # Detailed log file (all runs)
├── analysis_dashboard.html              # Interactive dashboard (green theme)
├── 01_taxonomic_report.txt              # ⭐ NEW: Enhanced taxonomic report
│   ├── Quick overview with dominant species
│   ├── Clinically significant organisms
│   ├── Clinician view (simplified summary)
│   ├── Researcher view (detailed statistics)
│   └── Diversity metrics
├── 02_functional_report.txt             # ⭐ NEW: Comprehensive functional report
│   ├── Genomic features overview
│   ├── Annotation statistics
│   ├── COG functional categories
│   ├── Mobile genetic elements analysis
│   ├── Top annotated functions
│   └── Annotation quality assessment
├── 03_pathogen_risk_report.txt          # ⭐ NEW: Three-tier pathogen assessment
│   ├── Overall pathogenicity risk score
│   ├── Tier 1: Taxonomic pathogenicity
│   ├── Tier 2: Functional pathogenicity markers
│   ├── Tier 3: ML pathogen predictions
│   ├── Integrated risk assessment
│   └── Clinical interpretation guide
├── pathogen_summary.txt                 # Legacy pathogen report
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
├── cog_annotation.tsv                   # ⭐ NEW: COG database annotations
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
├── metaquest.log                        # Detailed log file
├── analysis_dashboard.html              # Interactive dashboard (blue theme)
├── 01_taxonomic_report.txt              # ⭐ NEW: Enhanced taxonomic report
├── 02_functional_report.txt             # ⭐ NEW: Comprehensive functional report
├── 03_pathogen_risk_report.txt          # ⭐ NEW: Three-tier pathogen assessment
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
├── input_statistics.json                # Input file statistics
├── validation.json                      # File validation results
├── functional_annotation_report.txt     # Functional analysis summary
├── functional_annotation_report.json    # Functional data in JSON format
├── swissprot_annotation.tsv             # SwissProt protein annotations
├── cog_annotation.tsv                   # ⭐ NEW: COG database annotations
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
*Advanced statistical analysis outputs*

```
comparison_results/
├── taxonomic_abundance_table.tsv        # Species abundance matrix across all samples
├── differential_abundance_report.tsv    # Statistical significance testing results
├── alpha_diversity_metrics.tsv          # Alpha diversity calculations per sample
├── alpha_diversity_statistics.tsv       # Statistical test results for alpha diversity
├── permanova_results.txt                # PERMANOVA test details
├── anosim_results.txt                   # ANOSIM test details
├── random_forest_feature_importance.tsv # ML biomarker discovery results
├── beta_diversity_pcoa.html             # Interactive PCoA plot with variance explained
├── alpha_diversity_boxplot.html         # Multi-metric alpha diversity comparison
├── taxonomic_heatmap.html               # Interactive species abundance heatmap
├── volcano_plot.html                    # Differential abundance visualization
└── abundance_barplot.html               # Mean relative abundance by group
```

### Understanding the New Text Reports

#### 01_taxonomic_report.txt Structure

```
📊 MICROBIAL COMPOSITION - QUICK OVERVIEW
- Dominant species and abundance
- Total species detected
- Pathogenic species count
- Total pathogen load

⚠️ CLINICALLY SIGNIFICANT ORGANISMS DETECTED
- List of pathogenic species with:
  • Abundance percentage
  • Read counts
  • Risk level (🔴 HIGH / 🟡 MODERATE / 🟢 LOW)
  • Clinical notes

👨‍⚕️ CLINICIAN VIEW - SIMPLIFIED SUMMARY
- Known pathogens count
- Commensal organisms count
- Key findings
- Clinical recommendations

🔬 RESEARCHER VIEW - DETAILED STATISTICS
- Top 15 species by abundance
- Genus-level composition
- Rare species summary

📈 DIVERSITY METRICS
- Shannon index
- Simpson index
- Species richness
- Interpretation
```

#### 02_functional_report.txt Structure

```
🧬 GENOMIC FEATURES OVERVIEW
- Feature composition (CDS, rRNA, tRNA)
- Key metrics (bases, contigs, gene length)
- Assembly quality indicators

📊 ANNOTATION STATISTICS
- Total CDS and annotation coverage
- Hypothetical vs functionally annotated
- Annotation quality assessment
- Average sequence identity

🔬 COG FUNCTIONAL CATEGORIES
- Functional category distribution
- Bar charts showing relative abundance
- Key functional highlights

🔄 MOBILE GENETIC ELEMENTS ANALYSIS
- Total transposase hits
- IS family distribution
- Horizontal gene transfer risk assessment

🏆 TOP ANNOTATED FUNCTIONS
- Most common functional annotations
- Gene counts per function

✅ ANNOTATION QUALITY ASSESSMENT
- Overall quality score (0-100)
- Breakdown by metric
- Interpretation
```

#### 03_pathogen_risk_report.txt Structure

```
OVERALL PATHOGENICITY RISK
- Visual risk meter
- Risk level (🔴 HIGH / 🟡 MODERATE / 🟢 LOW)
- Component scores (taxonomic, functional, ML)

🦠 TIER 1: TAXONOMIC PATHOGENICITY ASSESSMENT
- Risk score
- Pathogens detected
- Pathogen risk matrix with action items

🧬 TIER 2: FUNCTIONAL PATHOGENICITY MARKERS
- Risk score
- Virulence factors detected
- AMR genes detected
- Mobile genetic elements analysis
- Pathogenicity-associated genes table

🤖 TIER 3: MACHINE LEARNING PATHOGEN PREDICTIONS
- Risk score
- Sequence analysis summary
- High-confidence predictions

🔬 INTEGRATED RISK ASSESSMENT
- Risk score breakdown
- Multipliers applied
- Intelligent interpretation

🩺 CLINICAL INTERPRETATION GUIDE
- What the risk level means
- Recommended actions
- Clinical guidance
- Additional considerations
```

---

## 💡 Examples & Use Cases

### Clinical Research Workflow

```bash
# 1. Validate clinical samples
metaquest validate fastq --paired patient1_R1.fastq.gz patient1_R2.fastq.gz -q 25

# 2. Analyze individual samples with annotation controls
metaquest analyze fastq --paired patient1_R1.fastq.gz patient1_R2.fastq.gz \
  --annotation-threads 16 \
  -o patient1_results/

metaquest analyze fastq --paired patient2_R1.fastq.gz patient2_R2.fastq.gz \
  --annotation-threads 16 \
  -o patient2_results/

metaquest analyze fastq --paired control1_R1.fastq.gz control1_R2.fastq.gz \
  --annotation-threads 16 \
  -o control1_results/

# 3. Review detailed text reports
cat patient1_results/01_taxonomic_report.txt
cat patient1_results/02_functional_report.txt
cat patient1_results/03_pathogen_risk_report.txt

# 4. Compare groups with statistical testing
metaquest compare -i patient1_results/ patient2_results/ control1_results/ \
  -m clinical_metadata.tsv -o clinical_comparison/
```

### Large Cohort Study

```bash
# Process multiple samples efficiently
for sample in healthy_{1..6} diseased_{1..6}; do
    metaquest analyze fastq --single ${sample}.fastq.gz \
      --annotation-threads 16 \
      --tbl2asn-timeout 600 \
      -o ${sample}_results/
done

# Comprehensive statistical comparison
metaquest compare -i healthy_*_results/ diseased_*_results/ \
  -m cohort_metadata.tsv -o cohort_comparison/
```

### Environmental Metagenomics

```bash
# 1. Analyze environmental samples with custom filtering
metaquest analyze fastq --single soil_sample1.fastq.gz \
  --min-contig-length 500 \
  -o soil1_results/

metaquest analyze fastq --single water_sample1.fastq.gz \
  --min-contig-length 500 \
  -o water1_results/

# 2. Compare environments with ML biomarker discovery
metaquest compare -i soil1_results/ water1_results/ \
  -m environment_metadata.tsv -o environment_comparison/
```

### Genome Analysis Workflow

```bash
# 1. Validate genome assembly
metaquest validate fasta assembled_genome.fasta

# 2. Full analysis with increased BLAST sampling and no contig filtering
metaquest analyze fasta assembled_genome.fasta \
  -s 200 \
  --no-filter-contigs \
  --annotation-threads 16 \
  -o genome_results/

# 3. Quick taxonomic-only analysis
metaquest analyze fasta genome.fasta --skip-annotation -o quick_results/

# 4. Review comprehensive reports
less genome_results/01_taxonomic_report.txt
less genome_results/02_functional_report.txt
less genome_results/03_pathogen_risk_report.txt
```

### Troubleshooting Analysis

```bash
# Run with debug mode to diagnose issues
metaquest --debug analyze fastq --single problematic_sample.fq \
  --annotation-threads 8 \
  --tbl2asn-timeout 600 \
  -o debug_results/

# Check the detailed log
less debug_results/metaquest.log

# Try with different annotation settings
metaquest analyze fastq --single problematic_sample.fq \
  --no-kill-tbl2asn \
  --min-contig-length 2000 \
  -o alternative_settings/
```

### High-Performance Computing Workflow

```bash
# Optimize for HPC environment with many cores
metaquest analyze fastq --paired R1.fq R2.fq \
  --annotation-threads 64 \
  --tbl2asn-timeout 1800 \
  --no-filter-contigs \
  -o hpc_results/
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

#### Annotation Failures

**tbl2asn Timeout:**
```bash
# Error: tbl2asn process timed out
# Solution 1: Increase timeout
metaquest analyze fastq --single reads.fq --tbl2asn-timeout 600

# Solution 2: Disable auto-kill
metaquest analyze fastq --single reads.fq --no-kill-tbl2asn

# Solution 3: Use debug mode to see details
metaquest --debug analyze fastq --single reads.fq
```

**Contig Filtering Issues:**
```bash
# Error: Too few contigs after filtering
# Solution 1: Lower threshold
metaquest analyze fasta assembly.fasta --min-contig-length 500

# Solution 2: Disable filtering
metaquest analyze fasta assembly.fasta --no-filter-contigs
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

# Or skip annotation
metaquest analyze fastq --single reads.fq --skip-annotation
```

#### Comparative Analysis Issues
```bash
# Error: Sample ID mismatch
# Solution: Ensure metadata sample_id matches directory names exactly
cat metadata.tsv  # Check sample IDs
ls *_results/     # Check directory names

# Error: Insufficient samples for statistics
# Solution: Ensure at least 3 samples per group for robust statistics
```

### Using Debug Mode for Troubleshooting

```bash
# Enable debug mode to see full details
metaquest --debug analyze fastq --single problematic.fq -o debug_out/

# Debug output will show:
# 1. Exact commands being run
# 2. Full tool output
# 3. Error messages with context
# 4. Performance metrics

# Check the log file
cat debug_out/metaquest.log
```

### Performance Optimization Tips

```bash
# 1. Use more threads for faster annotation
metaquest analyze fastq --single reads.fq --annotation-threads 32

# 2. Filter contigs more aggressively
metaquest analyze fasta assembly.fasta --min-contig-length 2000

# 3. Skip annotation if only taxonomy needed
metaquest analyze fastq --single reads.fq --skip-annotation

# 4. Reduce BLAST sample size for faster FASTA analysis
metaquest analyze fasta genome.fasta -s 25
```

### Sample Command Outputs

#### System Check Output
```bash
$ metaquest check
🧬 MetaQuest v4.0.0 | An Integrated Metagenomics Analysis Pipeline

══════════════════════════════════════════════════════════════════════
                        SYSTEM CHECK
══════════════════════════════════════════════════════════════════════

✓ Success! All dependencies and databases are correctly configured.
```

#### Standard Analysis Output
```bash
$ metaquest analyze fastq --single reads.fq -o results/
🧬 MetaQuest v4.0.0 | An Integrated Metagenomics Analysis Pipeline

══════════════════════════════════════════════════════════════════════
                      Analysis Pipeline
══════════════════════════════════════════════════════════════════════

⠋ Verifying system dependencies...
✓ System check passed

══════════════════════════════════════════════════════════════════════
                      File Validation
══════════════════════════════════════════════════════════════════════

ℹ Validating file 1/1: reads.fq
✓ Validation passed for reads.fq
✓ Input validation passed

══════════════════════════════════════════════════════════════════════
                    Annotation Settings
══════════════════════════════════════════════════════════════════════

Contig filtering     : Enabled
Min contig length    : 1000 bp
tbl2asn auto-kill    : Enabled
tbl2asn timeout      : 300s
Annotation threads   : 8

══════════════════════════════════════════════════════════════════════
                      FASTQ Analysis
══════════════════════════════════════════════════════════════════════

ℹ Input: reads.fq
ℹ Output: results/

⠋ Running taxonomic classification...
✓ Taxonomic classification completed

⠋ Running functional annotation...
✓ Functional annotation completed

⠋ Running pathogen detection...
✓ Pathogen detection completed

⠋ Generating reports and visualizations...
✓ Reports generated

✓ Analysis completed in 8m 42s
ℹ Results: results/analysis_dashboard.html
```

#### Debug Mode Output
```bash
$ metaquest --debug analyze fastq --single reads.fq -o debug_results/
🧬 MetaQuest v4.0.0 | An Integrated Metagenomics Analysis Pipeline

══════════════════════════════════════════════════════════════════════
                      Analysis Pipeline
══════════════════════════════════════════════════════════════════════

[DEBUG] Starting system check
[DEBUG] Checking for kraken2... found at /usr/local/bin/kraken2
[DEBUG] Checking for prokka... found at /usr/local/bin/prokka
[DEBUG] Checking for diamond... found at /usr/local/bin/diamond
✓ System check passed

══════════════════════════════════════════════════════════════════════
                      File Validation
══════════════════════════════════════════════════════════════════════

[DEBUG] Validating FASTQ file: reads.fq
[DEBUG] Reading first 10000 sequences for quality check
[DEBUG] Mean quality score: 32.4 (threshold: 20)
[DEBUG] Total sequences: 18455
[DEBUG] File size: 8.2 MB
✓ Validation passed for reads.fq

══════════════════════════════════════════════════════════════════════
                    Annotation Settings
══════════════════════════════════════════════════════════════════════

Contig filtering     : Enabled
Min contig length    : 1000 bp
tbl2asn auto-kill    : Enabled
tbl2asn timeout      : 300s
Annotation threads   : 8

══════════════════════════════════════════════════════════════════════
                      FASTQ Analysis
══════════════════════════════════════════════════════════════════════

[DEBUG] Converting FASTQ to FASTA
[DEBUG] Command: seqtk seq -A reads.fq > debug_results/converted.fasta
[DEBUG] seqtk: Converted 18455 sequences
[DEBUG] Starting Kraken2 classification
[DEBUG] Command: kraken2 --db databases/kraken2_db --threads 8 --output debug_results/kraken_classified.txt --report debug_results/kraken_report.txt debug_results/converted.fasta
[DEBUG] kraken2: Loading database information... done.
[DEBUG] kraken2: 18455 sequences (50.98%) classified
[DEBUG] kraken2: Classification rate: 50.98%
✓ Taxonomic classification completed (2m 34s)

[DEBUG] Starting Bracken abundance estimation
[DEBUG] Command: bracken -d databases/kraken2_db -i debug_results/kraken_report.txt -o debug_results/bracken_report.tsv -r 150
[DEBUG] bracken: Abundance estimation complete
✓ Abundance estimation completed (23s)

[DEBUG] Starting functional annotation
[DEBUG] Filtering contigs >= 1000 bp
[DEBUG] Filtered: 8579 contigs → 2341 contigs (27.3% retained)
[DEBUG] Command: prokka --outdir debug_results/prokka_annotation --prefix sample --cpus 8 --force debug_results/filtered_contigs.fasta
[DEBUG] prokka: [12:34:56] Setting up output directories
[DEBUG] prokka: [12:34:57] Running: prodigal -i filtered_contigs.fasta
[DEBUG] prokka: [12:35:23] Found 1179 CDS
[DEBUG] prokka: [12:35:24] Running: diamond blastp --query sample.faa --db databases/swissprot_db
[DEBUG] prokka: [12:36:12] SwissProt annotations: 458 hits
[DEBUG] prokka: [12:36:13] Running: diamond blastp --query sample.faa --db databases/cog_db
[DEBUG] prokka: [12:37:45] COG annotations: 573 hits
[DEBUG] prokka: [12:37:46] Running: tbl2asn (timeout: 300s)
[DEBUG] prokka: [12:38:12] tbl2asn completed successfully
[DEBUG] prokka: [12:38:13] Writing final output files
✓ Functional annotation completed (3m 17s)

[DEBUG] Starting pathogen detection
[DEBUG] Checking against pathogen database
[DEBUG] Command: diamond blastp --query debug_results/prokka_annotation/sample.faa --db databases/pathogen_db --outfmt 6
[DEBUG] diamond: Processed 1179 sequences
[DEBUG] diamond: Found 48 pathogen-associated genes
✓ Pathogen detection completed (1m 8s)

[DEBUG] Generating taxonomic report
[DEBUG] Writing: debug_results/01_taxonomic_report.txt
[DEBUG] Generating functional report
[DEBUG] Writing: debug_results/02_functional_report.txt
[DEBUG] Generating pathogen risk report
[DEBUG] Writing: debug_results/03_pathogen_risk_report.txt
[DEBUG] Generating HTML visualizations
✓ Reports generated (34s)

✓ Analysis completed in 8m 42s
ℹ Results: debug_results/analysis_dashboard.html
[DEBUG] Total memory used: 8.2 GB
[DEBUG] Total disk used: 1.4 GB
```

#### Comparative Analysis Output
```bash
$ metaquest compare -i healthy_*/ diseased_*/ -m metadata.tsv -o comparison/
🧬 MetaQuest v4.0.0 | An Integrated Metagenomics Analysis Pipeline

══════════════════════════════════════════════════════════════════════
                    Comparative Analysis
══════════════════════════════════════════════════════════════════════

Samples              : 12
Metadata             : metadata.tsv
Output               : comparison/

🔬 Starting Enhanced Comparative Microbiome Analysis
============================================================
→ Loading metadata from: metadata.tsv
→ Found 12 samples in metadata.
→ Aggregating taxonomic profiles...
  ✓ Processed: healthy_1
  ✓ Processed: healthy_2
  ✓ Processed: healthy_3
  ✓ Processed: healthy_4
  ✓ Processed: healthy_5
  ✓ Processed: healthy_6
  ✓ Processed: diseased_1
  ✓ Processed: diseased_2
  ✓ Processed: diseased_3
  ✓ Processed: diseased_4
  ✓ Processed: diseased_5
  ✓ Processed: diseased_6

→ Cleaning and validating abundance table...
✅ Cleaned abundance table: 1053 species × 12 samples

→ Calculating Alpha Diversity...
→ Running Alpha Diversity Statistical Tests...
  ✓ shannon: p = 0.0050 (significant **)
  ✓ simpson: p = 0.0050 (significant **)
  ✓ chao1: p = 0.0450 (significant *)
  ✓ sobs: p = 0.0450 (significant *)
  ✓ Alpha diversity calculations complete.

→ Calculating Beta Diversity (Bray-Curtis)...
→ Running Beta Diversity Statistical Tests...
  ✓ PERMANOVA: F = 5.418, p = 0.0040 (significant **)
  ✓ ANOSIM: R = 0.496, p = 0.0180 (overlapping but clearly different)

→ Running differential abundance analysis...
  → Comparing Healthy (n=6) vs Diseased (n=6)
  ✓ Differential abundance results:
    • FDR significant (q < 0.05): 0
    • Bonferroni significant: 0
    • Nominally significant (p < 0.05): 150
    • Trend level (p < 0.10): 241

→ Running Machine Learning for Biomarker Discovery...
  ✓ Random Forest 3-fold CV Accuracy: 0.83 ± 0.12
  ✓ Top discriminative species:
    → 'Lachnoanaerobaculum umeaense' - Importance: 0.0300
    → 'Butyrivibrio proteoclasticus' - Importance: 0.0300

→ Generating enhanced visualizations...
  ✓ PCoA plot created: beta_diversity_pcoa.html
  ✓ Heatmap visualization created: taxonomic_heatmap.html
  ✓ Alpha diversity box plot created: alpha_diversity_boxplot.html
  ✓ Volcano plot created: volcano_plot.html
  ✓ Abundance bar plot created: abundance_barplot.html

🎉 Enhanced comparative analysis completed!
📁 Results saved to: comparison/

✓ Comparison completed in 3m 15s
ℹ Results: comparison/
```

---

### Getting Help

```bash
# General help
metaquest --help

# Command-specific help
metaquest analyze --help
metaquest analyze fastq --help
metaquest analyze fasta --help
metaquest validate --help
metaquest compare --help

# Check version
metaquest --version
```

### Reading the Enhanced Reports

#### Tips for Clinicians

**Start with 01_taxonomic_report.txt:**
1. Check the "Quick Overview" for dominant species
2. Review "Clinically Significant Organisms" for pathogens
3. Read "Clinician View" for actionable summary
4. Note clinical recommendations

**Review 03_pathogen_risk_report.txt:**
1. Check overall risk score at top
2. Review Tier 1 for known pathogens
3. Check Tier 2 for resistance/virulence genes
4. Read clinical interpretation guide at bottom

#### Tips for Researchers

**Start with 01_taxonomic_report.txt:**
1. Review "Researcher View" for detailed statistics
2. Examine genus-level composition
3. Check diversity metrics
4. Identify rare species

**Review 02_functional_report.txt:**
1. Analyze COG functional categories
2. Check mobile genetic element activity
3. Review annotation quality scores
4. Examine top annotated functions

**Review 03_pathogen_risk_report.txt:**
1. Analyze all three risk tiers
2. Review functional pathogenicity markers
3. Check ML predictions
4. Examine integrated risk assessment

### Best Practices

**1. Always Validate First:**
```bash
metaquest validate fastq --single reads.fq
# Then run analysis only if validation passes
```

**2. Use Appropriate Annotation Settings:**
```bash
# For draft assemblies (fragmented)
metaquest analyze fasta draft.fasta --min-contig-length 500

# For high-quality assemblies
metaquest analyze fasta complete.fasta --min-contig-length 2000

# For time-critical analyses
metaquest analyze fastq --single reads.fq --skip-annotation
```

**3. Optimize Threading:**
```bash
# Match your CPU cores
metaquest analyze fastq --single reads.fq --annotation-threads $(nproc)
```

**4. Use Debug Mode for Issues:**
```bash
# If something fails, run with --debug
metaquest --debug analyze fastq --single reads.fq -o debug_out/
```

**5. Check Logs:**
```bash
# Always review the log file after analysis
cat results/metaquest.log
```

**6. Compare Multiple Samples:**
```bash
# Always use comparative analysis for multi-sample studies
metaquest compare -i sample_*/ -m metadata.tsv -o comparison/
```

---

### Support Resources

- **Installation Guide**: See `installation.md`
- **GitHub Issues**: Report bugs and request features
- **Documentation**: Check README.md for overview
- **Log Files**: Check `metaquest.log` in output directory for details

---

### Interpreting Report Quality Scores

#### Annotation Quality (02_functional_report.txt)

**Overall Quality Score: 0-100**
- **90-100**: Excellent - High coverage, high identity, comprehensive annotation
- **70-89**: Good - Reliable functional annotation with good coverage
- **50-69**: Moderate - Acceptable annotation, some gaps expected
- **<50**: Poor - Limited annotation, consider re-running with different settings

**Annotation Coverage:**
- **>60%**: Excellent coverage
- **40-60%**: Good coverage
- **30-40%**: Moderate coverage (⚠️)
- **<30%**: Poor coverage

**Average Identity:**
- **>80%**: High similarity to known proteins
- **70-80%**: Good similarity ✅
- **60-70%**: Moderate similarity
- **<60%**: Low similarity

#### Risk Scores (03_pathogen_risk_report.txt)

**Overall Pathogenicity Risk: 0-100**
- **75-100** 🔴 HIGH: Immediate attention required
- **50-74** 🟡 MODERATE: Clinical correlation recommended
- **25-49** 🟠 LOW-MODERATE: Monitor for symptoms
- **0-24** 🟢 LOW: Routine observation

---

*For technical support and updates, visit the MetaQuest GitHub repository.*

**Last Updated**: October 2025 - Version 4.0.0
> **Legacy document:** This guide describes the pre-0.1 experimental feature
> set. Custom pathogen, ML, HMM, ESM, island, Bayesian-risk, clinical, and
> advanced comparative claims are not part of the stable runtime. Use the root
> README for the current command behavior while this guide is rewritten.
