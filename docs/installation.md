# MetaQuest Installation Guide

This guide covers the complete installation process for MetaQuest, a comprehensive metagenomics analysis pipeline with enhanced pathogen detection capabilities and machine learning integration.

## Table of Contents
- [System Requirements](#system-requirements)
- [Installation Methods](#installation-methods)
- [Database Setup](#database-setup)
- [ML Model Setup](#ml-model-setup)
- [Package Installation](#package-installation)
- [Installation Verification](#installation-verification)
- [Troubleshooting](#troubleshooting)

## System Requirements

### Hardware Requirements
- **RAM**: Minimum 8GB, recommended 16GB+ for comprehensive pathogen analysis and ML predictions
- **Storage**: At least 20GB free space for databases, ML models, and results
- **CPU**: Multi-core processor recommended (8+ cores for optimal performance)
- **GPU**: Optional, for enhanced ML pathogen prediction performance

### Software Requirements
- Linux or macOS operating system
- Conda/Miniconda package manager (version 4.9+ recommended)
- Python 3.10+ (automatically installed with conda environment)
- Internet connection for database downloads (~15GB total)
- Bash shell version 4.0+

## Installation Methods

### Method 1: Conda Environment + Package Installation (Recommended)

1. **Clone or download MetaQuest**
```bash
   # If using git
   git clone https://github.com/Karudhoru/MetaQuest.git
   cd MetaQuest
   
   # Or extract from downloaded archive
   tar -xzf MetaQuest.tar.gz
   cd MetaQuest
```

2. **Create conda environment**
   ```bash
   conda env create -f environment.yml
   ```

3. **Activate the environment**
```bash
   conda activate metaquest
```

4. **Install MetaQuest package**
```bash
   # Install in development mode (recommended for contributors)
   pip install -e .
   
   # Or install normally
   pip install .
```

### Method 2: Direct Package Installation

If you have a compatible Python environment with required dependencies:

```bash
# Install directly from source
pip install -e .

# Install required bioinformatics tools
conda install -c bioconda diamond kraken2 seqtk krona prokka bracken taxonkit hmmer blast prodigal
```

## Package Installation

The new modular structure requires proper package installation:

### Understanding the New Structure

```
MetaQuest/
├── src/metaquest/              # Main package source
│   ├── cli.py                  # Command-line interface
│   ├── config.py              # Configuration management
│   ├── core/                  # Core analysis modules
│   ├── io/                    # Input/output handling
│   ├── ml/                    # Machine learning components
│   ├── reporting/             # Report generation
│   └── visualization/         # Data visualization
├── setup.py                   # Package setup configuration
├── pyproject.toml            # Modern Python packaging
└── environment.yml           # Conda environment
```

### Installation Process

1. **Install the package**
```bash
   cd MetaQuest
   pip install -e .
```

2. **Verify installation**
```bash
   # Check package installation
   python -c "import metaquest; print('MetaQuest package installed successfully')"
   
   # Check CLI availability
   metaquest --help
```

3. **Installing the Wrapper Script**
```bash
   # Navigate to your project root
   cd /mnt/d/GIT/MetaQuest

   # Update the wrapper script content
   cat > metaquest << 'EOF'
   #!/usr/bin/env bash
   # metaquest wrapper to invoke the installed package
   python -m metaquest.cli "$@"
   EOF

   # Make it executable
   chmod +x metaquest

   # Copy to conda environment bin
   ENV_BIN=$(dirname "$(which python)")
   cp metaquest "${ENV_BIN}/metaquest"
```

3. **Alternative CLI access methods**
```bash
   # Method 1: Direct package execution
   python -m metaquest.cli --help
   
   # Method 2: Installed entry point
   metaquest --help
   
   # Method 3: Full module path
   python -c "from metaquest.cli import main; main()" --help
```

## Database Setup

MetaQuest requires several databases for comprehensive analysis including ML model artifacts. Run the setup script to download and prepare all necessary databases:

```bash
# Make setup script executable
chmod +x scripts/setup_databases.sh

# Run database setup (this may take 60-120 minutes depending on connection)
./scripts/setup_databases.sh
```

### What gets downloaded:
- **Kraken2 Standard Database** (~8GB): For high-accuracy taxonomic classification
- **Bracken Database**: For abundance estimation (included with Kraken2)
- **NCBI BLAST Database** (~8GB): For FASTA sequence classification
- **CARD Database** (~500MB): For antimicrobial resistance gene detection
- **VFDB Database** (~200MB): For virulence factor identification
- **Custom Pathogen Database**: Built from CARD and pathogen-specific sequences
- **ML Model Files** (~1GB): Pre-trained pathogen prediction models
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
./scripts/setup_databases.sh --all
```

## ML Model Setup

The machine learning components require additional setup:

### Pre-trained Model Artifacts

MetaQuest includes pre-trained models in `src/metaquest/ml/model_artifacts/`:

```bash
# Verify ML model files exist
ls -la src/metaquest/ml/model_artifacts/
# Should show:
# - best_model.pkl          (Main classification model)
# - scaler.pkl             (Feature scaling model)
# - feature_selector.pkl   (Feature selection model)
# - feature_names.pkl      (Selected feature names)
# - all_feature_names.pkl  (Complete feature vocabulary)
```

### ML Dependencies Installation

```bash
# Ensure ML dependencies are installed
pip install scikit-learn joblib pandas numpy

# Verify ML components
python -c "from metaquest.ml.pathogen_predictor import PathogenPredictor; print('ML components available')"
python -c "from metaquest.ml.feature_extractor import FeatureExtractor; print('Feature extraction available')"
```

### Custom Model Training (Optional)

To train custom models with your own data:

```bash
# Access the ML training utilities
python -c "from metaquest.ml.pathogen_predictor import PathogenPredictor; pp = PathogenPredictor(); print('Ready for custom training')"
```

## Configuration

### Default Configuration

The package uses configuration from `src/metaquest/config.py`:

```python
# Default database paths
KRAKEN2_DB = "databases/kraken2_db"
PATHOGEN_DB = "databases/pathogen_db"
BLAST_DB = "databases/blast_db"
ML_MODELS_DB = "src/metaquest/ml/model_artifacts"
SWISSPROT_DB = "databases/swissprot_db"
```

### Custom Configuration

You can override default settings:

#### Method 1: Environment Variables
```bash
export METAQUEST_KRAKEN_DB="/custom/path/to/kraken2_db"
export METAQUEST_PATHOGEN_DB="/custom/path/to/pathogen_db"
export METAQUEST_BLAST_DB="/custom/path/to/blast_db"
export METAQUEST_ML_MODELS="/custom/path/to/ml_models"
export METAQUEST_SWISSPROT_DB="/custom/path/to/swissprot_db"
```

#### Method 2: Configuration File
```bash
# Create custom config file
cat > ~/.metaquest_config.py << EOF
KRAKEN2_DB = "/custom/path/to/kraken2_db"
PATHOGEN_DB = "/custom/path/to/pathogen_db"
ML_MODELS_DB = "/custom/path/to/ml_models"
EOF
```

## Installation Verification

After completing the installation, verify everything works correctly:

### 1. Check Environment
```bash
conda activate metaquest
which python
# Should point to the conda environment
python --version
# Should show Python 3.10+
```

### 2. Verify Package Installation
```bash
# Check package import
python -c "import metaquest; print(f'MetaQuest version: {metaquest.__version__}')"

# Check all modules
python -c "
from metaquest.core import analysis
from metaquest.io import file_validator
from metaquest.ml import pathogen_predictor
from metaquest.reporting import reporting
from metaquest.visualization import visualization
print('All modules imported successfully')
"
```

### 3. Test CLI Access
```bash
# Test command-line interface
metaquest --help
metaquest --version

# Test analysis capabilities
metaquest --check-system

# Test ML integration
python -c "from metaquest.ml.pathogen_predictor import PathogenPredictor; pp = PathogenPredictor(); print('ML pipeline ready')"
```

### 4. Verify Database Files
```bash
# Check database structure
ls -la databases/
# Should show:
# - kraken2_db/           (Kraken2 database files)
# - pathogen_db/          (CARD, VFDB databases)
# - blast_db/             (NCBI BLAST databases)
# - swissprot_db/         (SwissProt database)

# Check ML models
ls -la src/metaquest/ml/model_artifacts/
# Should show .pkl model files
```

### 5. Test Complete Analysis Pipeline
```bash
# Create test data directory
mkdir -p test_data

# Test FASTQ analysis (if you have test data)
metaquest test_data/example.fastq.gz -t fastq -o test_fastq_results/

# Test FASTA analysis with ML
metaquest test_data/example.fasta -t fasta -o test_fasta_results/

# Verify outputs
ls -la test_fastq_results/pathogen_summary.txt
ls -la test_fasta_results/blast_ml_pathogen_summary.txt
ls -la test_fasta_results/ml_pathogen_predictions.csv
```

### 6. Test ML Components Separately
```bash
# Test feature extraction
python -c "
from metaquest.ml.feature_extractor import FeatureExtractor
fe = FeatureExtractor()
print('Feature extractor initialized')
"

# Test pathogen prediction
python -c "
from metaquest.ml.pathogen_predictor import PathogenPredictor
pp = PathogenPredictor()
print('Pathogen predictor ready')
print(f'Model loaded: {pp.model is not None}')
"
```

## Troubleshooting

### Common Installation Issues

1. **Package import fails**
```bash
   # Reinstall in development mode
   pip install -e . --force-reinstall
   
   # Check Python path
   python -c "import sys; print(sys.path)"
   
   # Verify package location
   pip show metaquest
```

2. **CLI command not found**
```bash
   # Check if entry point was installed
   pip show metaquest | grep Entry-points
   
   # Try alternative access methods
   python -m metaquest.cli --help
   
   # Reinstall with entry points
   pip install -e . --force-reinstall
```

3. **ML components not working**
```bash
   # Check ML dependencies
   python -c "import sklearn, joblib, pandas; print('ML dependencies OK')"
   
   # Verify model files exist
   find src/metaquest/ml/model_artifacts/ -name "*.pkl"
   
   # Test ML import
   python -c "from metaquest.ml.pathogen_predictor import PathogenPredictor; print('ML import OK')"
```

4. **Database setup fails**
```bash
   # Check available disk space
   df -h
   
   # Test internet connectivity
   wget --spider https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240605.tar.gz
   
   # Try manual database setup
   mkdir -p databases && cd databases
```

5. **Module not found errors**
```bash
   # Ensure you're in the correct environment
   conda activate metaquest
   
   # Check package installation
   pip list | grep metaquest
   
   # Verify directory structure
   ls -la src/metaquest/
```

6. **Permission issues**
```bash
   # Fix script permissions
   chmod +x scripts/*.sh
   
   # Check directory permissions
   ls -la databases/
```

### Verification Commands

Use these commands to verify your installation:

```bash
# Environment verification
conda info --envs
conda list | head -20

# Package verification
python -c "import metaquest; print('Package OK')"
python -c "from metaquest.cli import main; print('CLI OK')"
python -c "from metaquest.ml.pathogen_predictor import PathogenPredictor; print('ML OK')"

# Database verification
find databases/ -name "*.k2d" -o -name "*.dmnd" -o -name "*.fas" -o -name "*.fasta" | wc -l

# ML model verification
find src/metaquest/ml/model_artifacts/ -name "*.pkl" | wc -l
```

### Complete Installation Example

Here's a complete installation workflow:

```bash
# 1. Clone and setup
git clone https://github.com/Karudhoru/MetaQuest.git
cd MetaQuest

# 2. Create environment
conda env create -f environment.yml
conda activate metaquest

# 3. Install package
pip install -e .

# 4. Setup databases
chmod +x scripts/setup_databases.sh
./scripts/setup_databases.sh --all

# 5. Verify installation
metaquest --version
python -c "from metaquest.ml.pathogen_predictor import PathogenPredictor; print('Complete installation verified')"

# 6. Test with sample data
metaquest [input_file] [type] --check-only
```

### Development Installation

For developers working on MetaQuest:

```bash
# Install with development dependencies
pip install -e ".[dev]"

# Install pre-commit hooks
pre-commit install

# Run tests
python -m pytest tests/

# Check code style
flake8 src/metaquest/
black src/metaquest/
```

## Uninstalling MetaQuest

To completely remove MetaQuest:

```bash
# Remove package
pip uninstall metaquest

# Remove conda environment
conda deactivate
conda env remove -n metaquest

# Remove databases (if desired)
rm -rf databases/

# Remove any custom configurations
rm -f ~/.metaquest_config.py
```

## Next Steps

After successful installation, proceed to the [Usage Guide](usage.md) to learn how to:
- Use the modular MetaQuest package
- Run enhanced FASTQ and FASTA analyses with ML integration
- Understand the new clinical vs research workflows
- Interpret enhanced output files including ML predictions
- Use the new interactive dashboards
- Troubleshoot analysis issues
- Customize ML model training

## Support

If you encounter issues during installation:
1. Check the troubleshooting section above
2. Verify all system requirements are met (especially increased RAM/storage for ML)
3. Ensure you have sufficient disk space for enhanced databases (~20GB)
4. Check internet connectivity for database downloads
5. Verify ML dependencies are properly installed
6. Check that the package structure matches the expected layout
7. Ensure entry points are properly installed for CLI access

For additional help:
- **GitHub Issues**: Report bugs and installation problems
- **Documentation**: Check the usage guide for detailed examples
- **Community**: Join discussions for installation tips and troubleshooting