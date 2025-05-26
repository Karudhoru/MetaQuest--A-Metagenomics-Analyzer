# MetaQuest Installation Guide

This guide covers the complete installation process for MetaQuest, a comprehensive metagenomics analysis pipeline.

## Table of Contents
- [System Requirements](#system-requirements)
- [Installation Methods](#installation-methods)
- [Database Setup](#database-setup)
- [Wrapper Script Creation](#wrapper-script-creation)
- [Installation Verification](#installation-verification)
- [Troubleshooting](#troubleshooting)

## System Requirements

### Hardware Requirements
- **RAM**: Minimum 8GB, recommended 32GB+ for large datasets
- **Storage**: At least 50GB free space for databases and results
- **CPU**: Multi-core processor recommended (4+ cores)

### Software Requirements
- Linux or macOS operating system
- Conda/Miniconda package manager
- Internet connection for database downloads
- Bash shell

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

MetaQuest requires several databases for comprehensive analysis. Run the setup script to download and prepare all necessary databases:

```bash
# Make setup script executable
chmod +x scripts/setup_databases.sh

# Run database setup (this may take 30-60 minutes)
./scripts/setup_databases.sh
```

### What gets downloaded:
- **MiniKraken2 Database** (~8GB): For taxonomic classification
- **CARD Database**: For antimicrobial resistance gene detection
- **VFDB Database**: For virulence factor identification
- **Custom pathogen database**: Built from CARD data

### Manual Database Setup

If the automatic setup fails, you can download databases manually:

```bash
# Create database directory
mkdir -p databases

# Download MiniKraken2
cd databases
wget https://genome-idx.s3.amazonaws.com/kraken/minikraken_8GB_202003.tgz
tar -xzf minikraken_8GB_202003.tgz

# Download CARD database
mkdir -p pathogen_db
cd pathogen_db
wget https://card.mcmaster.ca/latest/data
tar -xjf data

# Download VFDB
wget http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas

# Build DIAMOND database
diamond makedb --in nucleotide_fasta_protein_homolog_model.fasta --db pathogen_db
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

After completing the installation, verify everything works correctly:

### 1. Check Environment Activation
```bash
conda activate metagenomics_app
which python
# Should point to the conda environment
```

### 2. Verify Package Installation
```bash
python -c "import metagenomics; print('MetaQuest imported successfully')"
```

### 3. Test CLI Access
```bash
# Using wrapper script
metaquest --help

# Or using the package directly
python -m metagenomics.cli --help

# Or using entry point (if installed)
metagenomics-cli --help
```

### 4. Check Database Files
```bash
ls -la databases/
# Should show:
# - hash.k2d and taxo.k2d (Kraken2 files)
# - pathogen_db/ directory with CARD and VFDB files
```

### 5. Test with Example Data
```bash
# Run a quick test with example data
metaquest examples/example.fastq.gz -t fastq -o test_results/

# Check if results were generated
ls -la test_results/
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
```

### Environment Variables

You can also set environment variables to override default paths:

```bash
export METAQUEST_KRAKEN_DB="/path/to/kraken2_db"
export METAQUEST_PATHOGEN_DB="/path/to/pathogen_db"
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
   wget --spider https://genome-idx.s3.amazonaws.com/kraken/minikraken_8GB_202003.tgz
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

### Verification Commands

Use these commands to verify your installation:

```bash
# Check conda environment
conda info --envs

# Check installed packages
conda list | grep -E "(kraken2|diamond|seqtk)"

# Check Python packages
pip list | grep -E "(pandas|plotly|biopython)"

# Verify database files
find databases/ -name "*.k2d" -o -name "*.dmnd" -o -name "*.fas"
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

# 3. Setup databases
chmod +x scripts/setup_databases.sh
./scripts/setup_databases.sh

# 4. Create wrapper script
cat > metaquest << 'EOF'
#!/usr/bin/env bash
# metaquest wrapper to invoke the installed package
python -m metagenomics.cli "$@"
EOF

# 5. Make executable and copy to bin
chmod +x metaquest
ENV_BIN=$(dirname "$(which python)")
cp metaquest "${ENV_BIN}/metaquest"

# 6. Test installation
metaquest --help
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
- Run your first analysis
- Understand command-line options
- Interpret output files
- Troubleshoot analysis issues

## Support

If you encounter issues during installation:
1. Check the troubleshooting section above
2. Verify all system requirements are met
3. Ensure you have sufficient disk space and memory
4. Check internet connectivity for database downloads