# MetaQuest Usage Guide

**Version 3.2.1** | Enhanced Metagenomics Analysis Pipeline

---

## Overview

MetaQuest provides **analysis-specific workflows** optimized for different input types and use cases, featuring:

- **FASTQ Analysis**: Clinical-focused pipeline with rapid Kraken2/Bracken classification
- **FASTA Analysis**: Research-focused pipeline with high-accuracy BLAST classification  
- **Advanced Validation**: Comprehensive quality control with detailed statistics
- **Dual-Method Pathogen Detection**: Multi-source screening with clinical risk assessment

---

## Quick Start

### 1. Environment Setup
```bash
conda activate metaquest
```

### 2. File Validation (Recommended)
```bash
# Single-end FASTQ validation
metaquest validate fastq --single SRR33675829.fastq.gz

# FASTA validation with statistics
metaquest validate fasta sequence.fasta

# Custom quality thresholds
metaquest validate fastq --paired reads_R1.fastq.gz reads_R2.fastq.gz -q 25 -n 1000
```

### 3. Run Analysis
```bash
# Clinical FASTQ analysis
metaquest analyze fastq --single SRR33675829.fastq.gz -o fastq_results/

# Research FASTA analysis  
metaquest analyze fasta sequence.fasta -o fasta_results/
```

---

## Command Reference

### Main Commands

| Command | Purpose |
|---------|---------|
| `metaquest analyze` | Run full analysis pipeline |
| `metaquest validate` | Validate files without analysis |
| `metaquest check` | Verify dependencies and databases |
| `metaquest --version` | Display version information |

### FASTA Analysis

Optimized for assembled genomes and research applications.

```bash
metaquest analyze fasta <input_file> [options]
```

**Parameters:**
- `<input_file>` *(required)*: Path to input FASTA file

**Options:**
- `-o, --output <dir>`: Output directory (default: `results`)
- `-s, --blast-sample-size <int>`: NCBI BLAST sample size (default: 50)
- `--skip-validation`: Skip initial file validation

### FASTQ Analysis

Optimized for clinical samples with rapid classification.

```bash
metaquest analyze fastq <mode> [options]
```

**Modes:**
- `--single <file.fastq>`: Single-end reads
- `--paired <r1.fastq> <r2.fastq>`: Paired-end reads  
- `--interleaved <interleaved.fastq>`: Interleaved reads

**Options:**
- `-o, --output <dir>`: Output directory (default: `results`)
- `-q, --min-quality <int>`: Minimum average quality score (default: 20)
- `-n, --min-sequences <int>`: Minimum sequence count (default: 100)
- `--skip-validation`: Skip initial file validation

### File Validation

```bash
metaquest validate <type> [input_files] [options]
```

**Examples:**
```bash
metaquest validate fasta my_genome.fasta
metaquest validate fastq --paired r1.fq r2.fq -q 25
```

---

## Validation Reports

### FASTQ Validation Example

```
============================================================
🔍 METAQUEST FILE VALIDATION
============================================================
File: SRR33675829.fastq.gz
Type: FASTQ
Time: 2025-08-23 01:39:29
============================================================

📁 FILE INFORMATION
────────────────────────────────────────
  Size:        0.61 MB
  MD5:         a708be745aed5cd8...
  Compressed:  Yes (gzip)

🧬 SEQUENCE STATISTICS
────────────────────────────────────────
  Total Sequences:  19,643
  Total Bases:      8,479,628
  Length Range:     428 - 436 bp
  Mean Length:      431.7 bp
  Median Length:    432.0 bp
  N50 Length:       432.0 bp
  GC Content:       46.6%

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

⚙️  VALIDATION CRITERIA
────────────────────────────────────────
  Minimum sequences:  100 ✅ (Found: 19,643)
  Minimum quality:    20 ✅ (Found: 30.0)

============================================================
✅ FILE VALIDATION: PASSED
============================================================
```

### FASTA Validation Example

```
============================================================
🔍 METAQUEST FILE VALIDATION
============================================================
File: sequence.fasta
Type: FASTA
Time: 2025-08-23 01:39:29
============================================================

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

⚙️  VALIDATION CRITERIA
────────────────────────────────────────
  Minimum sequences:  1 ✅ (Found: 1)
  Unique IDs:         Required ✅ (All unique)

============================================================
✅ FILE VALIDATION: PASSED
============================================================
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

## Key Features

### Analysis-Specific Workflows
- **FASTQ**: Optimized for clinical samples with rapid classification
- **FASTA**: Optimized for research with high-accuracy methods

### Comprehensive Validation
- File integrity checks with MD5 hashing
- Quality score analysis and thresholds
- Sequence statistics and composition analysis

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