# MetaQuest Usage Guide

**Version 3.2.2** | Enhanced Metagenomics Analysis Pipeline

---

## Overview

MetaQuest provides **analysis-specific workflows** optimized for different input types and use cases, featuring:

- **FASTQ Analysis**: Clinical-focused pipeline with rapid Kraken2/Bracken classification
- **FASTA Analysis**: Research-focused pipeline with high-accuracy BLAST classification
- **Advanced Validation**: Comprehensive quality control with detailed statistics and customizable thresholds
- **Dual-Method Pathogen Detection**: Multi-source screening with clinical risk assessment

---

## Quick Start

### 1. Check Your Installation (Recommended)
```bash
metaquest check
```

### 2. File Validation (Recommended)

```bash
# Validate a single-end FASTQ file
metaquest validate fastq --single SRR33675829.fastq.gz

# Validate with a custom quality and overrepresentation threshold
metaquest validate fastq --paired reads_R1.fastq.gz reads_R2.fastq.gz -q 25 --overrep-threshold 0.5
```

### 3. Run an Analysis

```bash
# FASTQ analysis with default settings
metaquest analyze fastq --single SRR33675829.fastq.gz -o fastq_results/

# FASTA analysis with a larger BLAST sample size
metaquest analyze fasta sequence.fasta -o fasta_results/ -s 200
```

---

## Command Reference

### Main Commands

| Command | Purpose |
|---|---|
| `metaquest analyze` | Run the full analysis pipeline |
| `metaquest validate` | Validate files without running analysis |
| `metaquest check` | Verify all dependencies and databases |
| `metaquest -v, --version` | Display version information |

### The `analyze` Command

This is the main command for running your analysis. It requires you to specify the data type (`fasta` or `fastq`) and provide the input files.

#### **FASTA Analysis**

Optimized for assembled genomes and research applications.

```bash
metaquest analyze fasta <input_file> [options]
```

- **`<input_file>`** *(required)*: Path to your input FASTA file
- **Options:**
  - `-o, --output <dir>`: Output directory (default: `results`)
  - `-s, --blast-sample-size <int>`: Number of sequences to BLAST (default: `50`)
  - `--skip-validation`: Skip the pre-analysis file check

#### **FASTQ Analysis**

Optimized for clinical samples with rapid classification.

```bash
metaquest analyze fastq <mode> [options]
```

- **Modes** (choose one):
  - `--single <file.fastq>`: For single-end reads
  - `--paired <r1.fastq> <r2.fastq>`: For paired-end reads
  - `--interleaved <interleaved.fastq>`: For interleaved reads
- **Options:**
  - `-o, --output <dir>`: Output directory (default: `results`)
  - `-q, --min-quality <int>`: Minimum average quality for validation (default: `20`)
  - `-n, --min-sequences <int>`: Minimum sequence count for validation (default: `100`)
  - `--overrep-threshold <float>`: Percentage to flag a sequence as overrepresented (default: `0.1`)
  - `--skip-validation`: Skip the pre-analysis file check

### The `validate` Command

Checks a file without running the full analysis.

```bash
metaquest validate <type> [input_files/flags] [options]
```

- **`<type>`**: `fasta` or `fastq`
- **Options (for `fastq` type):**
  - `-q, --min-quality <int>`: Minimum average quality (default: `20`)
  - `-n, --min-sequences <int>`: Minimum sequence count (default: `100`)
  - `--overrep-threshold <float>`: Overrepresentation threshold (default: `0.1`)
- **Examples:**
  ```bash
  metaquest validate fasta my_genome.fasta
  metaquest validate fastq --paired r1.fq r2.fq -q 25
  ```

---

## Output Structure

### FASTQ Analysis Results
*Clinical-focused pipeline outputs*

```
fastq_results/
├── analysis_dashboard.html              # 🌐 Interactive dashboard (green theme)
├── pathogen_summary.txt                 # 🎯 Primary pathogen report
├── amr_hits.txt                         # 💊 Antimicrobial resistance genes
├── bracken_report.tsv                   # 📊 Species abundance estimation
├── converted.fasta                      # 🔄 FASTQ to FASTA conversion
├── functional_annotation_report.txt     # 🧪 Functional analysis
├── kraken_report.txt                    # 🔬 Taxonomic classification
├── prokka_annotation/                   # 🧬 Gene prediction results
└── [additional analysis files]
```

### FASTA Analysis Results  
*Research-focused pipeline outputs*

```
fasta_results/
├── analysis_dashboard.html              # 🌐 Interactive dashboard (blue theme)
├── blast_ml_pathogen_summary.txt        # 🎯 Primary pathogen report
├── blast_taxonomy_results.json          # 🔬 BLAST taxonomy (raw JSON)
├── functional_annotation_report.txt     # 🧪 Functional analysis
├── ml_pathogen_predictions.csv          # 🤖 Individual ML predictions
├── ml_pathogen_summary.json             # 🤖 ML analysis summary
├── prokka_annotation/                   # 🧬 Gene prediction results
└── [additional analysis files]
```

---

## Validation Reports

### FASTQ Validation Example

```
🧬 MetaQuest v3.2.2 | A Comprehensive Metagenomics Analysis Pipeline

🧬 MetaQuest - File Validation

--- Validating file 1/1: gut.fastq.gz ---

============================================================
🔍 METAQUEST FILE VALIDATION
============================================================
File: gut.fastq.gz
Type: FASTQ
Time: 2025-08-23 23:45:34
============================================================

⏳ Analyzing file contents...
  Processing: ........... Done!

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
  CGGACTACCGGGGTTTCTAATCCTGTTCGCTACCCAT... 0.6851%
  CGGACTACAGGGGTTTCTAATCCTGTTCGCTACCCAT... 0.6468%
  CGGACTACAAGGGTTTCTAATCCTGTTCGCTACCCAT... 0.6434%
  CGGACTACCGGGGTTTCTAATCCTGTTCGCTACCCAT... 0.6358%
  CGGACTACACGGGTTTCTAATCCTGTTCGCTACCCAT... 0.6090%

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

### FASTA Validation Example

```
🧬 MetaQuest v3.2.2 | A Comprehensive Metagenomics Analysis Pipeline

🧬 MetaQuest - File Validation

--- Validating file 1/1: sequence.fasta ---

============================================================
🔍 METAQUEST FILE VALIDATION
============================================================
File: sequence.fasta
Type: FASTA
Time: 2025-08-23 23:49:19
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

---

## Key Features

### Analysis-Specific Workflows
- **FASTQ**: Optimized for clinical samples with rapid classification
- **FASTA**: Optimized for research with high-accuracy methods

### Comprehensive Validation
- File integrity checks with MD5 hashing
- Quality score analysis and thresholds
- Sequence statistics and composition analysis
- Adapter contamination detection
- Overrepresented sequence identification

### Dual-Method Pathogen Detection
- Multi-source pathogen screening
- Clinical risk assessment integration
- Machine learning pathogen predictions

### Zero-Redundancy Output
- Analysis-specific file structures
- No duplicate or unnecessary files
- Clear primary report identification

---

*For technical support and updates, visit the MetaQuest documentation.*