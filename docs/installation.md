# MetaQuest Installation Guide

This guide covers the complete installation process for MetaQuest, a comprehensive metagenomics analysis pipeline with enhanced pathogen detection capabilities.

## Table of Contents
- [System Requirements](#system-requirements)
- [Installation Methods](#installation-methods)
- [Database Setup](#database-setup)
- [Wrapper Script Creation](#wrapper-script-creation)
- [Installation Verification](#installation-verification)
- [Troubleshooting](#troubleshooting)

## System Requirements

### Hardware Requirements
- **RAM**: Minimum 16GB, recommended 32GB+ for comprehensive pathogen analysis
- **Storage**: At least 15GB free space for databases, ML models, and results
- **CPU**: Multi-core processor recommended (8+ cores for optimal performance)
- **GPU**: Optional, for enhanced ML pathogen prediction performance

### Software Requirements
- Linux or macOS operating system
- Conda/Miniconda package manager (version 4.9+ recommended)
- Python 3.8+ (automatically installed with conda environment)
- Internet connection for database downloads (~10GB total)
- Bash shell version 4.0+

## Installation Methods

### Method 1: Conda Environment (Recommended)

1. **Clone or download MetaQuest**
   ```bash
   # If using git
   git clone https://github.com/Karudhoru/MetaQuest--A-Metagenomics-Analyzer.git
   cd metaquest_v3
   
   # Or extract from downloaded archive
   tar -xzf metaquest_v3.tar.gz
   cd metaquest_v3
   ```

2. **Create conda environment**
   ```bash
   conda env create -f environment/environment.yml
   ```

3. **Activate the environment**
   ```bash
   conda activate metagenomics_app
   ```

4. **Install MetaQuest package**
   ```bash
   cd metaquest_v3
   pip install -e .
   ```

### Method 2: Manual Installation

If you prefer not to use conda, install dependencies manually:

```bash
# Install Python dependencies
pip install -r requirements.txt

# Install bioinformatics tools (using conda or system package manager)
conda install -c bioconda diamond kraken2 seqtk krona prokka bracken taxonkit hmmer blast prodigal
```

## Database Setup

MetaQuest requires several databases for comprehensive analysis including the new ML pathogen prediction capabilities. Run the setup script to download and prepare all necessary databases:

```bash
# Make setup script executable
chmod +x scripts/setup_databases.sh

# Run database setup (this may take 45-90 minutes depending on connection)
./scripts/setup_databases.sh
```

### What gets downloaded:
- **Kraken2 Mini Database** (~8GB): For high-accuracy taxonomic classification
- **Bracken Database**: For abundance estimation (included with Kraken2)
- **NCBI BLAST Database** (~8GB): For FASTA sequence classification
- **CARD Database** (~500MB): For antimicrobial resistance gene detection
- **VFDB Database** (~200MB): For virulence factor identification
- **Custom Pathogen Database**: Built from CARD and pathogen-specific sequences
- **ML Model Files** (~1GB): Pre-trained pathogen prediction models (when available)
- **SwissProt Database** (~2GB): For functional protein annotation

### Enhanced Database Setup Features

#### Automatic Database Validation
```bash
# The setup script now includes validation
./scripts/setup_databases.sh --validate

# Check database integrity
./scripts/check_databases.sh
```

#### Selective Database Installation
```bash
# Install only specific databases
./scripts/setup_databases.sh --kraken-only
./scripts/setup_databases.sh --pathogen-only
./scripts/setup_databases.sh --ml-models
```

### Manual Database Setup

If the automatic setup fails, you can download databases manually:

```bash
# Create database directory
mkdir -p databases

# Download Kraken2 Standard Database
cd databases
mkdir kraken2_db
cd kraken2_db
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240605.tar.gz
tar -xzf k2_standard_20240605.tar.gz

# Download CARD database
cd ../
mkdir -p pathogen_db
cd pathogen_db
wget https://card.mcmaster.ca/latest/data
tar -xjf data

# Download VFDB
wget http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas

# Download SwissProt database
cd ../
mkdir swissprot_db
cd swissprot_db
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz

# Build DIAMOND database
diamond makedb --in ../pathogen_db/nucleotide_fasta_protein_homolog_model.fasta --db pathogen_db
```

### Database Setup Flowchart

```txt
Start
│
├─ Check Kraken2 Standard DB → If missing → Download Standard Database
│
├─ Check CARD DB → If missing → Download & extract
│
├─ Check VFDB → If missing → Download
│
├─ Check SwissProt DB → If missing → Download & build
│
├─ Check BLAST DB → If missing → Download NCBI databases
│
├─ Check ML Models → If missing → Download pre-trained models
│
└─ Check DIAMOND index → If missing → Create from CARD
│
End
```

## Wrapper Script Creation

To make MetaQuest easily accessible from the command line, create a wrapper script:

### Step 1: Create the wrapper script file

1. In your project directory (where you have `MetaQuest` unpacked), create a new file named `metaquest` (no extension):

```bash
# Create the wrapper script
touch metaquest
nano metaquest # Use any one to create file
```

2. Edit the file to have exactly these contents:

```bash
#!/usr/bin/env bash
# metaquest wrapper to invoke the installed package
python -m metagenomics.cli "$@"
```

This means whenever you run `metaquest`, it will invoke your package's `cli.py` module.

### Step 2: Make the script executable

In your shell:

```bash
chmod +x metaquest
```

Now `./metaquest --help` in that directory should print your CLI help.

### Step 3: Copy it into your conda env's bin folder

When you activate your `metagenomics_app` environment, `which python` will point at something like:
```
/root/miniconda3/envs/metagenomics_app/bin/python
```

We want its `bin/` directory. Run:

```bash
ENV_BIN=$(dirname "$(which python)")
echo "Copying wrapper into $ENV_BIN"
cp metaquest "${ENV_BIN}/metaquest"
```

Now your `metaquest` script lives alongside `python`, `pip`, etc.

### Step 4: Test the wrapper script

Still in your shell (with the env activated):

```bash
# Check if metaquest is in PATH
which metaquest
# should output something like /root/miniconda3/envs/metagenomics_app/bin/metaquest

# Test the help command
metaquest --help

# Test version command
metaquest --version  # Should show v3.2.0+
```

### Alternative: Direct Installation Method

If you prefer not to create a wrapper script, you can also use the package directly:

```bash
# Install with console script entry point
pip install -e .

# This should make metagenomics-cli available directly
metagenomics-cli --help
```

## Installation Verification

After completing the installation, verify everything works correctly with the enhanced features:

### 1. Check Environment Activation
```bash
conda activate metagenomics_app
which python
# Should point to the conda environment
```

### 2. Verify Package Installation
```bash
python -c "import metagenomics; print('MetaQuest imported successfully')"
python -c "from metagenomics.reporting import PathogenReporter; print('Enhanced pathogen reporting available')"
```

### 3. Test Enhanced CLI Access
```bash
# Using wrapper script
metaquest --help
metaquest --version  # Should show v3.2.0+

# Test analysis type detection
metaquest --check-system

# Or using the package directly
python -m metagenomics.cli --help

# Or using entry point (if installed)
metagenomics-cli --help
```

### 4. Verify Enhanced Database Files
```bash
ls -la databases/
# Should show:
# - kraken2_db/ (Standard Kraken2 database)
# - pathogen_db/ (CARD, VFDB, custom pathogen databases)
# - blast_db/ (NCBI BLAST databases)
# - ml_models/ (ML pathogen prediction models, if available)
# - swissprot_db/ (SwissProt database)
```

### 5. Test Enhanced Analysis Capabilities
```bash
# Test FASTQ analysis with pathogen detection
metaquest examples/example.fastq.gz -t fastq -o test_fastq_results/

# Test FASTA analysis with BLAST+ML integration
metaquest examples/example.fasta -t fasta -o test_fasta_results/

# Verify enhanced outputs
ls test_fastq_results/pathogen_summary.txt
ls test_fasta_results/blast_ml_pathogen_summary.txt

# Check if results were generated with enhanced features
ls -la test_fastq_results/analysis_dashboard.html
ls -la test_fasta_results/ml_pathogen_predictions.csv
```

## Configuration

### Custom Database Paths

If you want to use custom database locations, edit the configuration file:

```bash
# Edit the config file
nano metaquest_v3/metagenomics/config.py
```

Update the database paths:
```python
# Custom database paths
KRAKEN2_DB = "/custom/path/to/kraken2_db"
PATHOGEN_DB = "/custom/path/to/pathogen_db"
BLAST_DB = "/custom/path/to/blast_db"
ML_MODELS_DB = "/custom/path/to/ml_models"
SWISSPROT_DB = "/custom/path/to/swissprot_db"
```

### Environment Variables

You can also set environment variables to override default paths:

```bash
export METAQUEST_KRAKEN_DB="/path/to/kraken2_db"
export METAQUEST_PATHOGEN_DB="/path/to/pathogen_db"
export METAQUEST_BLAST_DB="/path/to/blast_db"
export METAQUEST_ML_MODELS="/path/to/ml_models"
export METAQUEST_SWISSPROT_DB="/path/to/swissprot_db"
```

## Troubleshooting

### Common Installation Issues

1. **Conda environment creation fails**
   ```bash
   # Try updating conda first
   conda update conda
   conda env create -f environment/environment.yml
   ```

2. **Database download fails**
   ```bash
   # Check internet connection and try manual download
   wget --spider https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240605.tar.gz
   ```

3. **Permission denied on setup script**
   ```bash
   chmod +x scripts/setup_databases.sh
   ls -la scripts/setup_databases.sh  # Verify permissions
   ```

4. **Wrapper script not found**
   ```bash
   # Ensure you're in the correct environment
   conda activate metagenomics_app
   
   # Check if wrapper was copied correctly
   ls -la $(dirname "$(which python)")/metaquest
   ```

5. **Package import errors**
   ```bash
   # Reinstall in development mode
   pip install -e . --force-reinstall
   ```

6. **Enhanced pathogen reporting not available**
   ```bash
   # Check if all required packages are installed
   python -c "import metagenomics.reporting; print('Reporting module available')"
   
   # Verify ML dependencies
   python -c "import sklearn, joblib; print('ML dependencies available')"
   ```

### Verification Commands

Use these commands to verify your installation:

```bash
# Check conda environment
conda info --envs

# Check installed packages
conda list | grep -E "(kraken2|diamond|seqtk|blast|prokka|bracken)"

# Check Python packages
pip list | grep -E "(pandas|plotly|biopython|scikit-learn|joblib)"

# Verify database files
find databases/ -name "*.k2d" -o -name "*.dmnd" -o -name "*.fas" -o -name "*.fasta"

# Check ML models
find databases/ml_models/ -name "*.pkl" -o -name "*.joblib" 2>/dev/null || echo "ML models not found (optional)"
```

### Complete Installation Example

Here's a complete installation workflow:

```bash
# 1. Create and activate environment
conda env create -f environment/environment.yml
conda activate metagenomics_app

# 2. Install MetaQuest package
cd metaquest_v3
pip install -e .

# 3. Setup databases with validation
chmod +x scripts/setup_databases.sh
./scripts/setup_databases.sh --validate

# 4. Create wrapper script
cat > metaquest << 'EOF'
#!/usr/bin/env bash
# metaquest wrapper to invoke the installed package
+ "$(dirname "$(which python)")/python" -m metagenomics.cli "$@"
EOF

# 5. Make executable and copy to bin
chmod +x metaquest
ENV_BIN=$(dirname "$(which python)")
cp metaquest "${ENV_BIN}/metaquest"

# 6. Test enhanced installation
metaquest --help
metaquest --version
metaquest --check-system
```

### Uninstalling MetaQuest

To completely remove MetaQuest:

```bash
# Remove conda environment
conda deactivate
conda env remove -n metagenomics_app

# Remove databases (if desired)
rm -rf databases/

# Remove any custom configurations
rm -f ~/.metaquest_config
```

## Next Steps

After successful installation, proceed to the [Usage Guide](usage.md) to learn how to:
- Run enhanced FASTQ and FASTA analyses
- Understand the new clinical vs research workflows
- Interpret enhanced output files including ML predictions
- Use the new interactive dashboards
- Troubleshoot analysis issues

## Support

If you encounter issues during installation:
1. Check the troubleshooting section above
2. Verify all system requirements are met (especially increased RAM/storage)
3. Ensure you have sufficient disk space for enhanced databases (~15GB)
4. Check internet connectivity for database downloads
5. Verify ML dependencies are properly installed for enhanced pathogen prediction