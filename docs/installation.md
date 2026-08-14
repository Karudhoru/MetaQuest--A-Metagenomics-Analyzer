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
├── setup.py
├── requirements.txt
├── environment/
│   └── environment.yml
├── scripts/
│   └── setup_databases.sh
├── docs/
│   ├── installation.md
│   └── annotation.md
└── src/metaquest/
    ├── cli.py                          # CLI entry point
    ├── config.py                       # Central configuration
    ├── config/
    │   └── pathogen_config.json        # Pathogen DB config
    ├── core/
    │   ├── analysis.py
    │   ├── assembly.py
    │   ├── taxonomic_analysis.py
    │   ├── pathogen_analysis.py
    │   ├── functional_analysis.py
    │   └── comparative_analysis.py
    ├── io/
    │   ├── file_validator.py
    │   ├── output_formatter.py
    │   ├── data_loaders.py
    │   ├── text_parsers.py
    │   └── utils.py
    ├── ml/
    │   ├── feature_extractor.py
    │   ├── pathogen_predictor.py
    │   └── model_artifacts/
    ├── reporting/
    │   ├── main_reporter.py
    │   ├── base_reporter.py
    │   ├── validation_engine.py
    │   ├── taxonomy_reporter.py
    │   ├── functional_reporter.py
    │   ├── pathogen_risk_reporter.py
    │   └── risk_scoring.py
    └── visualization/
        ├── main_visualizer.py
        ├── dashboard.py
        ├── taxonomic_visualizer.py
        ├── functional_visualizer.py
        ├── pathogenic_visualizer.py
        └── compare_visuals.py
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

3. **Verifying the entry point**
```bash
   # pip install -e . registers the entry point automatically
   # Verify it is on your PATH:
   which metaquest
   metaquest --help
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
chmod +x scripts/setup_databases.sh

# Download everything
./scripts/setup_databases.sh --all

# Or selectively:
./scripts/setup_databases.sh --kraken       # MiniKraken2 (~8GB)
./scripts/setup_databases.sh --swissprot   # SwissProt + COG combined DB
./scripts/setup_databases.sh --pathogen    # Build DIAMOND pathogen DB from your FASTA

# Then build the custom pathogen marker database
python scripts/custom_pathogen_db.py
```

### What gets downloaded

| Database | Size | Purpose |
|----------|------|---------|
| MiniKraken2 | ~8 GB | Taxonomic classification |
| SwissProt + COG | ~3 GB | Functional annotation |
| Pathogen markers | ~500 MB | Pathogen detection (built from CARD + VFDB) |

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
metaquest -h
metaquest -v

# Test analysis capabilities
metaquest check

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
python scripts/custom_pathogen_db.py

# 5. Verify installation
metaquest --version

# 6. Check all external tools are available
metaquest check
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
> **Legacy document:** Dependency and database instructions in this guide cover
> the earlier experimental pipeline. They are not the final minimal MetaQuest
> 0.1 installation contract. Use the root README for current scope.
