"""
MetaQuest - A Comprehensive Metagenomics Analysis Pipeline
============================================================

MetaQuest is an integrated bioinformatics pipeline for metagenomic data analysis,
providing taxonomic classification, pathogen detection, machine learning integration,
and statistical analysis capabilities.

Author: MetaQuest Development Team
License: TBD
Python: >=3.8
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read long description from README
readme_file = Path(__file__).parent / "README.md"
with open(readme_file, "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Read requirements
requirements_file = Path(__file__).parent / "requirements.txt"
with open(requirements_file, "r", encoding="utf-8") as fh:
    requirements = [
        line.strip() 
        for line in fh 
        if line.strip() and not line.startswith("#")
    ]

# Package metadata
PACKAGE_NAME = "MetaQuest"
VERSION = "4.0.0"
DESCRIPTION = "Comprehensive metagenomics analysis pipeline with taxonomic classification, pathogen detection, and ML integration"
AUTHOR = "MetaQuest Development Team"
AUTHOR_EMAIL = "devpatelcambay@gmail.com"
URL = "https://github.com/Karudhoru/MetaQuest--A-Metagenomics-Analyzer"
LICENSE = "TBD"
PYTHON_REQUIRES = ">=3.8"

# Classifiers for PyPI
CLASSIFIERS = [
    # Development Status
    "Development Status :: 5 - Production/Stable",
    
    # Intended Audience
    "Intended Audience :: Science/Research",
    "Intended Audience :: Healthcare Industry",
    "Intended Audience :: Developers",
    
    # License
    "License :: OSI Approved :: MIT License",
    
    # Operating Systems
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS",
    "Operating System :: Unix",
    
    # Programming Languages
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
    
    # Topics
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Topic :: Scientific/Engineering :: Medical Science Apps.",
    "Topic :: Software Development :: Libraries :: Python Modules",
    
    # Natural Language
    "Natural Language :: English",
    
    # Environment
    "Environment :: Console",
]

# Keywords for discoverability
KEYWORDS = [
    "metagenomics",
    "bioinformatics",
    "pathogen-detection",
    "taxonomic-classification",
    "machine-learning",
    "microbiome",
    "genomics",
    "fastq",
    "fasta",
    "clinical-diagnostics",
]

# Project URLs for additional resources
PROJECT_URLS = {
    "Bug Reports": "https://github.com/Karudhoru/MetaQuest--A-Metagenomics-Analyzer/issues",
    "Documentation": "https://github.com/Karudhoru/MetaQuest--A-Metagenomics-Analyzer/blob/main/README.md",
    "Source Code": "https://github.com/Karudhoru/MetaQuest--A-Metagenomics-Analyzer",
    "Changelog": "https://github.com/Karudhoru/MetaQuest--A-Metagenomics-Analyzer/releases",
}

# Optional dependencies for advanced features
EXTRAS_REQUIRE = {
    "dev": [
        "pytest>=7.0.0",
        "pytest-cov>=4.0.0",
        "black>=22.0.0",
        "flake8>=5.0.0",
        "mypy>=0.990",
        "pre-commit>=2.20.0",
    ],
    "docs": [
        "sphinx>=5.0.0",
        "sphinx-rtd-theme>=1.0.0",
        "sphinx-autodoc-typehints>=1.19.0",
    ],
    "ml": [
        "scikit-learn>=1.2.0",
        "pandas>=1.5.0",
        "numpy>=1.23.0",
        "matplotlib>=3.6.0",
        "seaborn>=0.12.0",
        "xgboost>=1.7.0",
        "lightgbm>=3.3.0",
        "catboost>=1.0.0"
    ],
}

# Add 'all' option that includes everything
EXTRAS_REQUIRE["all"] = list(set(sum(EXTRAS_REQUIRE.values(), [])))

setup(
    # Basic package information
    name=PACKAGE_NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type="text/markdown",
    
    # Author information
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    
    # URLs
    url=URL,
    project_urls=PROJECT_URLS,
    
    # License
    license=LICENSE,
    
    # Package discovery
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    
    # Package data
    include_package_data=True,
    package_data={
        "metaquest": [
            "ml/model_artifacts/*"
        ],
    },
    
    # Dependencies
    python_requires=PYTHON_REQUIRES,
    install_requires=requirements,
    extras_require=EXTRAS_REQUIRE,
    
    # Entry points
    entry_points={
        "console_scripts": [
            "metaquest=metaquest.cli:main",
        ],
    },
    
    # PyPI metadata
    classifiers=CLASSIFIERS,
    keywords=KEYWORDS,
    
    # Additional metadata
    platforms=["Linux", "macOS", "Unix"],
    zip_safe=False,  # Important for package data access
)