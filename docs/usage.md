# MetaQuest Usage Guide

A comprehensive guide to using MetaQuest's enhanced metagenomics analysis pipeline with advanced file validation, dual-method pathogen detection, and clinical risk assessment.

## Overview

MetaQuest v3.2.1 provides **analysis-specific workflows** optimized for different input types and use cases, with comprehensive file validation:

- **FASTQ Analysis**: Clinical-focused with rapid Kraken2/Bracken classification and multi-source pathogen screening
- **FASTA Analysis**: Research-focused with high-accuracy BLAST classification and integrated ML pathogen predictions
- **Advanced Validation**: Comprehensive quality control for both input types with detailed statistics and quality thresholds

## Quick Start

### File Validation (Recommended First Step)

```bash
# Activate environment
conda activate metagenomics_app

# Validate FASTQ file before analysis
metaquest SRR33675829.fastq.gz -t fastq --validate-only

# Validate FASTA file with detailed statistics
metaquest sequence.fasta -t fasta  --validate-only

# Check file with custom quality thresholds
metaquest reads.fastq.gz -t fastq --min-quality 25 --min-sequences 1000 --validate-only
```

### Enhanced Analysis Workflows

```bash
### FASTQ Analysis (Recommended for clinical samples)
# Multi-source pathogen detection with clinical risk assessment
metaquest SRR33675829.fastq.gz -t fastq -o fastq_results/

### FASTA Analysis (Optimized for assembled genomes)  
# BLAST taxonomy + ML pathogen predictions with separated reporting
metaquest sequence.fasta -t fasta -o fasta_results/

### Skip validation for trusted files (advanced users)
metaquest trusted_file.fastq.gz -t fastq -o results/ --skip-validation
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

### File Validation Options (New in v3.2.1)

**Validation Commands:**
- `--validate-only`: Validate file and show statistics without running analysis
- `--skip-validation`: Skip file validation (not recommended)
- `--min-quality INT`: Minimum mean quality score for FASTQ files (default: 20)
- `--min-sequences INT`: Minimum number of sequences required (default: 100)

**Validation Examples:**
```bash
# Basic validation
metaquest sample.fastq.gz -t fastq  --validate-only

# Custom quality thresholds
metaquest sample.fastq.gz -t fastq --min-quality 25 --min-sequences 5000  --validate-only

# FASTA validation with statistics
metaquest genome.fasta -t fasta  --validate-only
```

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

## File Validation Features

### FASTQ Validation Output Example

```
============================================================
🔍 METAQUEST FILE VALIDATION
============================================================
File: SRR33675829.fastq.gz
Type: FASTQ
Time: 2025-08-03 18:47:44
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

### FASTA Validation Output Example

```
============================================================
🔍 METAQUEST FILE VALIDATION
============================================================
File: sequence.fasta
Type: FASTA
Time: 2025-08-03 18:48:16
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

### Validation Features by File Type

#### FASTQ Validation Capabilities
- **Quality Score Analysis**: Automatic detection of quality encoding (Phred+33/Phred+64)
- **Statistical Metrics**: Mean, median, N50 length calculations with quality score distributions
- **Quality Thresholds**: Q20/Q30 base percentage analysis with customizable thresholds
- **Contamination Detection**: Adapter sequence screening in first 1000 reads
- **Sequence Composition**: GC content analysis and ambiguous base detection
- **File Integrity**: MD5 checksum calculation and compression detection

#### FASTA Validation Capabilities
- **Sequence Type Detection**: Automatic classification as protein or nucleotide sequences
- **Composition Analysis**: GC content statistics with per-sequence variance analysis
- **Quality Assessment**: Duplicate ID detection, gap analysis, and N-base counting
- **Length Distribution**: Comprehensive length statistics with coefficient of variation
- **Category Classification**: Intelligent sequence categorization (genome, contigs, genes, etc.)
- **Format Validation**: Strict FASTA format compliance checking

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
- **📊 Quality Control**: Comprehensive FASTQ quality analysis and validation

#### FASTA Analysis Advantages  
- **🔬 Research Focus**: High-accuracy BLAST classification with quality metrics
- **🤖 ML Integration**: Separated BLAST taxonomy and ML pathogenicity predictions
- **📈 Confidence Scoring**: ML provides pathogenicity confidence for individual proteins
- **🎨 Blue Dashboard**: Research workflow optimized interface
- **🎯 Quality Assessment**: Detailed sequence identity and E-value metrics
- **🧪 Composition Analysis**: Advanced sequence statistics and diversity metrics

## Advanced Validation Workflows

### Quality Control Pipelines

#### High-Throughput FASTQ Validation
```bash
# Validate multiple FASTQ files
for file in *.fastq.gz; do
    metaquest --validate-only "$file" -t fastq --min-quality 25
done

# Batch validation with strict quality control
metaquest --validate-only high_quality.fastq.gz -t fastq \
    --min-quality 30 --min-sequences 10000
```

#### Research FASTA Validation
```bash
# Validate assembled genome
metaquest --validate-only genome_assembly.fasta -t fasta

# Check metagenomic contigs
metaquest --validate-only contigs.fasta -t fasta --min-sequences 100
```

### Integration with Analysis Workflows

#### Recommended Workflow Pattern
```bash
# Step 1: Validate input file
metaquest --validate-only sample.fastq.gz -t fastq

# Step 2: Review validation output
# Check for quality warnings and adjust parameters if needed

# Step 3: Run full analysis (validation runs automatically)
metaquest sample.fastq.gz -t fastq -o results/

# Alternative: Skip validation for trusted files
metaquest trusted_sample.fastq.gz -t fastq -o results/ --skip-validation
```

## Sample Output Files

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

#### **2. Bracken Abundance Report (FASTQ Input)**
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
```

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

MACHINE LEARNING PATHOGENICITY PREDICTIONS:
===========================================
High Confidence Pathogenic Proteins: 8
Medium Confidence Pathogenic Proteins: 12
Low Confidence Pathogenic Proteins: 5

INTEGRATED RISK ASSESSMENT:
==========================
BLAST Taxonomy Risk: HIGH (Salmonella enterica dominant)
ML Pathogenicity Risk: MEDIUM (8 high-confidence pathogenic proteins)
Overall Assessment: HIGH RISK
```

#### **2. ML Pathogen Predictions (FASTA Input)**
```csv
# ml_pathogen_predictions.csv format
protein_id,pathogenicity_score,confidence,prediction,organism_context
protein_001,0.85,HIGH,PATHOGENIC,Salmonella enterica
protein_002,0.72,MEDIUM,PATHOGENIC,Escherichia coli
protein_003,0.45,LOW,NON_PATHOGENIC,Cloning vector
protein_004,0.91,HIGH,PATHOGENIC,Salmonella enterica
protein_005,0.38,LOW,NON_PATHOGENIC,Unknown
```

#### **3. Organism Comparison Data (FASTA Input)**
```csv
# organism_comparison_data.csv format
organism,total_hits,sequences_with_hits,avg_hits_per_sequence,max_identity,avg_identity
Salmonella enterica,10,1,10.0,98.5,95.2
Escherichia coli,6,1,6.0,96.8,92.1
Klebsiella pneumoniae,2,1,2.0,89.3,87.6
Cloning vector,2,1,2.0,100.0,98.9
```

### **Prokka Annotation Results (Both FASTQ and FASTA)**

#### **Functional Annotation Summary**
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

#### **Annotation Statistics**
```txt
# sample.txt format
organism: Genus species strain 
contigs: 19643
bases: 8479628
CDS: 44150
rRNA: 12
tRNA: 67
tmRNA: 1
```

## Validation Best Practices

### When to Use Validation

#### **Always Validate When:**
- Working with new datasets or unknown quality files
- Preparing data for clinical applications
- Quality control is critical for downstream analysis
- Files come from external sources or collaborators

#### **Consider Skipping Validation When:**
- Files have been previously validated and haven't changed
- Running batch processing on trusted datasets
- Speed is critical and file quality is guaranteed
- Advanced users with proven data quality pipelines

### Quality Threshold Guidelines

#### **FASTQ Quality Recommendations:**
- **Minimum Quality Score**: 20 (basic), 25 (recommended), 30 (high-quality)
- **Minimum Sequences**: 100 (testing), 1,000 (basic analysis), 10,000+ (clinical)
- **Q30 Percentage**: >80% (recommended), >90% (high-quality)

#### **FASTA Quality Indicators:**
- **Sequence Count**: Depends on application (1 for genomes, 100+ for contigs)
- **Length Distribution**: Check coefficient of variation for expected uniformity
- **GC Content**: Should be within expected range for organism type
- **Duplicate IDs**: Should be zero for proper analysis

### Troubleshooting Validation Issues

#### **Common FASTQ Issues:**
```bash
# Low quality scores
metaquest low_quality.fastq.gz -t fastq --min-quality 15 --validate-only

# Insufficient sequences
metaquest small_file.fastq.gz -t fastq --min-sequences 50 --validate-only

# Quality encoding issues
# Check validation output for encoding detection
```

#### **Common FASTA Issues:**
```bash
# Duplicate sequence IDs (will be auto-renamed)
metaquest duplicates.fasta -t fasta --validate-only

# Mixed sequence types
# Check validation output for sequence type detection

# Empty or corrupted files
# Validation will detect and report format issues
```

## Interactive Features

### Dashboard Navigation

#### **FASTQ Dashboard Features (Green Theme):**
- **Clinical Risk Overview**: Color-coded pathogen risk levels
- **Abundance Charts**: Interactive species composition
- **Quality Metrics**: FASTQ-specific quality visualizations
- **Detection Coverage**: Multi-method pathogen detection analysis

#### **FASTA Dashboard Features (Blue Theme):**
- **BLAST Results**: Organism similarity and identity metrics
- **ML Predictions**: Pathogenicity confidence scoring
- **Composition Analysis**: GC content and length distributions
- **Integrated Assessment**: Combined BLAST+ML risk evaluation

### Real-time Analysis Features

#### **Progress Monitoring:**
- File validation progress with detailed statistics
- Analysis step completion indicators
- Real-time quality metric updates
- Error detection and reporting

#### **Result Interpretation:**
- Color-coded risk levels (Red=HIGH, Orange=MEDIUM, Green=LOW)
- Clinical priority indicators with emergency protocols
- Confidence scoring for ML predictions
- Multi-method validation for pathogen detection

## Advanced Usage Patterns

### Batch Processing

#### **Multiple File Validation:**
```bash
#!/bin/bash
# Validate multiple FASTQ files
for file in data/*.fastq.gz; do
    echo "Validating $file..."
    metaquest --validate-only "$file" -t fastq --min-quality 25
    if [ $? -eq 0 ]; then
        echo "✅ $file passed validation"
        # Run analysis
        metaquest "$file" -t fastq -o "results/$(basename $file .fastq.gz)/"
    else
        echo "❌ $file failed validation"
    fi
done
```

#### **Quality-Controlled Pipeline:**
```bash
#!/bin/bash
# High-quality analysis pipeline
QUALITY_THRESHOLD=30
MIN_SEQUENCES=10000

# Validate with strict criteria
metaquest --validate-only "$1" -t fastq \
    --min-quality $QUALITY_THRESHOLD \
    --min-sequences $MIN_SEQUENCES

if [ $? -eq 0 ]; then
    echo "Starting high-quality analysis..."
    metaquest "$1" -t fastq -o high_quality_results/
else
    echo "File does not meet high-quality criteria"
    exit 1
fi
```

### Integration with External Tools

#### **Pre-processing Integration:**
```bash
# Quality filtering before MetaQuest
fastp -i raw.fastq.gz -o filtered.fastq.gz -q 25 -l 100

# Validate filtered file
metaquest filtered.fastq.gz -t fastq --validate-only 

# Run MetaQuest analysis
metaquest filtered.fastq.gz -t fastq -o results/
```

## Getting Help

### Common Command Patterns

```bash
# Show version and feature status
metaquest --version

# Get comprehensive help
metaquest --help

# Check system dependencies
metaquest --check-only 

# Validate file with custom settings
metaquest [file] -t [type] --min-quality [score] --min-sequences [count] --validate-only 

# Run analysis with validation
metaquest [file] -t [type] -o [output_dir]

# Skip validation (advanced)
metaquest [file] -t [type] -o [output_dir] --skip-validation
```

### Support Resources

For issues not covered in this guide:
- Check the [Installation Guide](installation.md) for setup issues
- Review validation output for specific quality concerns
- Contact the development team for clinical interpretation guidance
- Visit the project repository for bug reports and feature requests

---

This completes the comprehensive MetaQuest usage guide with advanced file validation features. The validation system ensures data quality before analysis begins, providing detailed statistics and quality assessment for both FASTQ and FASTA inputs.