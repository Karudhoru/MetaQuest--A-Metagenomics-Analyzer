Based on our comprehensive work throughout the conversation and the major pathogen analysis completion (v3.2.0), here's the **complete rewritten usage.md** that reflects all the changes, enhanced features, and actual output files:

# MetaQuest Usage Guide

A comprehensive guide to using MetaQuest's enhanced metagenomics analysis pipeline with dual-method pathogen detection and clinical risk assessment.

## Overview

MetaQuest v3.2.0 provides **analysis-specific workflows** optimized for different input types and use cases:

- **FASTQ Analysis**: Clinical-focused with rapid Kraken2/Bracken classification and multi-source pathogen screening
- **FASTA Analysis**: Research-focused with high-accuracy BLAST classification and integrated ML pathogen predictions

## Quick Start

### Enhanced Analysis Workflows

```bash
# Activate environment
conda activate metagenomics_app

### FASTQ Analysis (Recommended for clinical samples)
# Multi-source pathogen detection with clinical risk assessment
metaquest examples/example.fastq.gz -t fastq -o fastq_results/

### FASTA Analysis (Optimized for assembled genomes)  
# BLAST taxonomy + ML pathogen predictions with separated reporting
metaquest examples/sequence.fasta -t fasta -o fasta_results/
```

### View Enhanced Results

```bash
### FASTQ Results - Clinical Workflow
# Open green-themed FASTQ dashboard
firefox fastq_results/analysis_dashboard.html

# View THE definitive pathogen report with clinical recommendations
less fastq_results/pathogen_summary.txt

# Check interactive visualizations
firefox fastq_results/pathogen_risk_detection.html
firefox fastq_results/taxonomic_abundance_chart.html

### FASTA Results - Research Workflow
# Open blue-themed FASTA dashboard
firefox fasta_results/analysis_dashboard.html

# View THE definitive integrated BLAST+ML report
less fasta_results/blast_ml_pathogen_summary.txt

# Check individual ML predictions with confidence scores
less fasta_results/ml_pathogen_predictions.csv
```

## Command Line Options

### Enhanced Global Options

```bash
metaquest [input_files] -t {fastq|fasta} -o OUTPUT_DIR [options]
```

**Required Arguments:**
- `-t, --type {fasta,fastq}`: Specify input file type (required)
- `-o, --output OUTPUT_DIR`: Output directory (required)

**Enhanced Global Options:**
- `-h, --help`: Show comprehensive help with examples
- `--threads INT`: Number of threads (default: auto-detect)
- `--version`: Show MetaQuest version and feature status
- `--check-system`: Verify installation and database status

### FASTQ Input Options (Clinical Focus)

**Single-end FASTQ:**
```bash
metaquest sample.fastq.gz -t fastq -o results/
metaquest -r sample.fastq.gz -t fastq -o results/
```

**Paired-end FASTQ:**
```bash
metaquest sample_R1.fastq.gz sample_R2.fastq.gz -t fastq -o results/
metaquest -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -t fastq -o results/
```

**Interleaved FASTQ:**
```bash
metaquest -i sample_interleaved.fastq.gz -t fastq -o results/
```

**FASTQ-Specific Options:**
- `-r, --reads`: Single-end FASTQ file
- `-1, --reads1`: First paired-end FASTQ file (R1)  
- `-2, --reads2`: Second paired-end FASTQ file (R2)
- `-i, --interleaved`: Interleaved paired-end FASTQ file
- `--clinical-mode`: Enhanced clinical risk assessment (default: enabled)

### FASTA Input Options (Research Focus)

**Single FASTA file:**
```bash
metaquest genome.fasta -t fasta -o fa_results/
metaquest sequences.fasta -t fasta -o fa_results/
```

## Output Files

MetaQuest generates analysis-specific output files optimized for each input type with **zero redundancy**:

### Directory Structure

#### FASTQ Analysis Results *(Clinical-focused with traditional pathogen screening)*

```
fastq_results/
├── **analysis_dashboard.html**                 # 🌐 Green-themed FASTQ dashboard  
├── **pathogen_summary.txt**                   # 🎯 THE definitive FASTQ pathogen report
├── amr_hits.txt                               # 💊 AMR gene detection results
├── bracken_report.tsv                         # 📊 Bracken abundance estimation (TSV)
├── bracken_report.txt                         # 📊 Bracken abundance estimation (TXT)
├── converted.fasta                            # 🔄 Converted FASTQ to FASTA
├── detection_method_coverage.html             # 📈 Pathogen detection method coverage
├── functional_annotation_report.json          # 🧪 Functional analysis (JSON)
├── functional_annotation_report.txt           # 🧪 Functional analysis summary
├── kraken_classified.txt                      # 🔬 Kraken2 classified reads
├── kraken_report.txt                          # 🔬 Kraken2 taxonomic classification
├── pathogen_detection_report.json             # 📋 Enhanced pathogen data (JSON)
├── pathogen_results.txt                       # 🦠 Traditional pathogen screening
├── pathogen_risk_detection.html               # ⚠️ Interactive pathogen risk chart
├── prokka_annotation/                         # 🧬 Complete gene prediction results
│   ├── sample.err                             # Prokka error log
│   ├── sample.faa                             # Protein sequences (FASTA)
│   ├── sample.ffn                             # Gene sequences (FASTA)
│   ├── sample.fna                             # Nucleotide sequences
│   ├── sample.fsa                             # Contig sequences
│   ├── sample.gbk                             # GenBank format annotation
│   ├── sample.gff                             # Gene feature format
│   ├── sample.log                             # Annotation process log
│   ├── sample.sqn                             # Sequin format
│   ├── sample.tbl                             # Feature table
│   ├── sample.tsv                             # Tab-separated annotations
│   └── sample.txt                             # Annotation summary
├── protein_length_distribution.html           # 📏 Protein quality analysis
├── swissprot_annotation.tsv                  # 🔬 SwissProt functional annotations
├── taxonomic_abundance_chart.html             # 📈 Interactive abundance visualization
├── taxonomic_classification_report.json       # 📋 Taxonomic data (JSON)
├── taxonomic_classification_report.txt        # 📊 Comprehensive taxonomic report
├── taxonomy_krona.html                        # 🌐 Krona taxonomic visualization
└── virulence_hits.txt                         # ⚔️ Virulence factor detection results
```

#### FASTA Analysis Results *(Research-focused with BLAST+ML integration)*

```
fa_results/
├── **analysis_dashboard.html**                # 🌐 Blue-themed FASTA dashboard
├── **blast_ml_pathogen_summary.txt**         # 🎯 THE definitive FASTA pathogen report
├── blast_cache/                              # 💾 BLAST results caching
│   └── blast_cache.json                      # Cached BLAST data
├── blast_ml_integrated_pathogen_report.json  # 📋 Integrated BLAST+ML data (JSON)
├── blast_taxonomy_results.json               # 🔬 BLAST taxonomy results (JSON)
├── detection_method_coverage.html            # 📈 Detection method coverage analysis
├── functional_annotation_report.json         # 🧪 Functional analysis (JSON)
├── functional_annotation_report.txt          # 🧪 Functional analysis summary
├── ml_pathogen_predictions.csv               # 🤖 Individual ML protein predictions
├── ml_pathogen_predictions.json              # 🤖 ML predictions (JSON)
├── ml_pathogen_summary.json                  # 🤖 ML analysis summary
├── organism_comparison_data.csv               # 📊 BLAST organism comparison (CSV)
├── organism_comparison_data.json              # 📊 BLAST organism comparison (JSON)
├── pathogen_risk_detection.html              # ⚠️ Interactive pathogen risk chart
├── prokka_annotation/                        # 🧬 Complete gene prediction results
│   ├── sample.err                            # Prokka error log
│   ├── sample.faa                            # Protein sequences (FASTA)
│   ├── sample.ffn                            # Gene sequences (FASTA)
│   ├── sample.fna                            # Nucleotide sequences
│   ├── sample.fsa                            # Contig sequences
│   ├── sample.gbk                            # GenBank format annotation
│   ├── sample.gff                            # Gene feature format
│   ├── sample.log                            # Annotation process log
│   ├── sample.sqn                            # Sequin format
│   ├── sample.tbl                            # Feature table
│   ├── sample.tsv                            # Tab-separated annotations
│   └── sample.txt                            # Annotation summary
├── protein_length_distribution.html          # 📏 Protein quality analysis
├── swissprot_annotation.tsv                 # 🔬 SwissProt functional annotations
├── taxonomic_classification_report.json      # 📋 Taxonomic data (JSON)
└── taxonomic_classification_report.txt       # 📊 BLAST-based taxonomic report
```

### 🎯 **THE Definitive Reports** *(Zero Redundancy)*

#### FASTQ Analysis - Single Source of Truth
- **pathogen_summary.txt**: THE comprehensive pathogen report with clinical risk assessment, emergency protocols, and multi-source detection evidence

#### FASTA Analysis - Single Source of Truth  
- **blast_ml_pathogen_summary.txt**: THE integrated report with separated BLAST taxonomy findings and ML pathogenicity predictions

### Key Analysis-Specific Features

#### FASTQ Analysis Advantages
- **🏥 Clinical Focus**: Emergency protocols and risk stratification (CRITICAL/HIGH/MEDIUM/LOW)
- **📊 Multi-source Detection**: Combines Bracken taxonomy, custom databases, AMR/VF screening
- **⚡ Rapid Processing**: Kraken2/Bracken for fast clinical decision support
- **🎨 Green Dashboard**: Clinical workflow optimized interface
- **📈 Abundance Estimation**: Accurate species-level quantification

#### FASTA Analysis Advantages  
- **🔬 Research Focus**: High-accuracy BLAST classification with quality metrics
- **🤖 ML Integration**: Separated BLAST taxonomy and ML pathogenicity predictions
- **📈 Confidence Scoring**: ML provides pathogenicity confidence for individual proteins
- **🎨 Blue Dashboard**: Research workflow optimized interface
- **🎯 Quality Assessment**: Detailed sequence identity and E-value metrics

## Files for Display in Documentation

### **FASTQ Analysis - Key Files to Showcase**

#### **1. pathogen_summary.txt** *(THE definitive FASTQ report)*
```
COMPREHENSIVE PATHOGEN DETECTION SUMMARY
========================================

Overall Risk Assessment: CRITICAL
Total Pathogens Detected: 15
High Risk Pathogens: 6
Medium Risk Pathogens: 5
Analysis Date: 2025-07-27 15:26:36

🚨 HIGH RISK PATHOGENS DETECTED:
• Salmonella enterica
  Risk Level: HIGH
  Detection Methods: bracken, taxonomy, sequence
  Abundance: 0.045%
  Clinical Priority: IMMEDIATE ATTENTION REQUIRED

• Escherichia coli
  Risk Level: HIGH
  Detection Methods: bracken, sequence  
  Abundance: 0.032%
  Clinical Priority: IMMEDIATE ATTENTION REQUIRED

CLINICAL RECOMMENDATIONS:
🚨 CRITICAL RISK - EMERGENCY PROTOCOLS REQUIRED
• Multiple high-risk pathogens detected across sources
• Implement immediate containment measures
• Emergency infectious disease consultation
```

#### **2. analysis_dashboard.html** *(Green-themed FASTQ interface)*
- Interactive clinical workflow dashboard
- Real-time pathogen risk visualization
- Multi-source detection coverage analysis
- Emergency protocol recommendations

#### **3. pathogen_risk_detection.html** *(Interactive risk chart)*
- Color-coded pathogen risk levels (Red=HIGH, Orange=MEDIUM)
- Clinical priority indicators
- Detection method validation

#### **4. taxonomic_abundance_chart.html** *(Species composition)*
- Top 15 abundant taxa visualization
- Interactive bar charts with abundance percentages
- Kraken2/Bracken integration

### **FASTA Analysis - Key Files to Showcase**

#### **1. blast_ml_pathogen_summary.txt** *(THE definitive FASTA report)*
```
INTEGRATED BLAST+ML PATHOGEN DETECTION SUMMARY
==============================================

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