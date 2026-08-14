"""
MetaQuest - A Comprehensive Metagenomics Analysis Pipeline
============================================================

MetaQuest is a research-use pipeline for taxonomic classification, metagenomic
assembly, gene prediction, and functional annotation.

Author: MetaQuest Development Team
License: MIT
Python: >=3.10,<3.13
"""

from setuptools import setup, find_packages
from pathlib import Path
import sys

# Add src to path to import version from config
sys.path.insert(0, str(Path(__file__).parent / "src"))
from metaquest.config import __version__, __release_date__

# Read long description from README
readme_file = Path(__file__).parent / "README.md"
if readme_file.exists():
    with open(readme_file, "r", encoding="utf-8") as fh:
        long_description = fh.read()
else:
    long_description = __doc__

# Read requirements
requirements_file = Path(__file__).parent / "requirements.txt"
if requirements_file.exists():
    with open(requirements_file, "r", encoding="utf-8") as fh:
        requirements = [
            line.strip() 
            for line in fh 
            if line.strip() and not line.startswith("#")
        ]
else:
    # Fallback to minimal requirements
    requirements = [
        "biopython>=1.81,<2.0",
        "pandas>=1.5.0,<3.0",
        "numpy>=1.23.0,<2.0",
    ]

# Package metadata
PACKAGE_NAME = "metaquest"
VERSION = __version__
DESCRIPTION = "Research-use short-read metagenomics analysis pipeline"
AUTHOR = "Dev Patel"
AUTHOR_EMAIL = "devpatelcambay@gmail.com"
URL = "https://github.com/Karudhoru/MetaQuest--A-Metagenomics-Analyzer"
LICENSE = "MIT"
PYTHON_REQUIRES = ">=3.10,<3.13"

MAINTAINER = "MetaQuest Development Team"
MAINTAINER_EMAIL = "devpatelcambay@gmail.com"

# Classifiers for PyPI
CLASSIFIERS = [
    # Development Status
    "Development Status :: 3 - Alpha",
    
    # Intended Audience
    "Intended Audience :: Science/Research",
    "Intended Audience :: Developers",
    
    # License
    "License :: OSI Approved :: MIT License",
    
    # Operating Systems
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS",
    "Operating System :: Unix",
    
    # Programming Languages
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
    "Programming Language :: Python :: 3 :: Only",
    
    # Topics
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Topic :: Software Development :: Libraries :: Python Modules",
    
    # Natural Language
    "Natural Language :: English",
    
    # Environment
    "Environment :: Console",
    
    # Additional relevant classifiers
    "Framework :: Matplotlib",
    "Typing :: Typed",
]

# Keywords for discoverability
KEYWORDS = [
    "metagenomics",
    "bioinformatics",
    "taxonomic-classification",
    "microbiome-analysis",
    "genomics",
    "ngs",
    "fastq",
    "kraken2",
    "bracken",
    "comparative-metagenomics",
]

# Project URLs for additional resources
PROJECT_URLS = {
    "Homepage": URL,
    "Bug Reports": f"{URL}/issues",
    "Documentation": f"{URL}/blob/main/README.md",
    "Source Code": URL,
    "Changelog": f"{URL}/releases",
    "Download": f"{URL}/releases/latest",
    # ADD when ready:
    # "Paper": "https://doi.org/10.xxxx/xxxxx",
    # "Tutorial": "https://metaquest.readthedocs.io",
}

# Optional dependencies for advanced features
EXTRAS_REQUIRE = {
    "dev": [
        "pytest>=7.0.0",
        "pytest-cov>=4.0.0",
        "pytest-timeout>=2.1.0",
        "black>=22.0.0",
        "flake8>=5.0.0",
        "mypy>=0.990",
        "pre-commit>=2.20.0",
        "ipython>=8.0.0",
    ],
    "docs": [
        "sphinx>=5.0.0",
        "sphinx-rtd-theme>=1.0.0",
        "sphinx-autodoc-typehints>=1.19.0",
        "myst-parser>=1.0.0",
    ],
    "test": [
        "pytest>=7.0.0",
        "pytest-cov>=4.0.0",
        "pytest-timeout>=2.1.0",
        "coverage>=7.0.0",
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
    maintainer=MAINTAINER,
    maintainer_email=MAINTAINER_EMAIL,
    
    # URLs
    url=URL,
    project_urls=PROJECT_URLS,
    download_url=f"{URL}/archive/v{VERSION}.tar.gz",
    
    # License
    license=LICENSE,
    
    # Package discovery
    package_dir={"": "src"},
    packages=find_packages(where="src", exclude=["tests", "tests.*"]),
    
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
    keywords=", ".join(KEYWORDS),
    
    # Additional metadata
    platforms=["Linux", "macOS", "Unix"],
    zip_safe=False,  # Important for package data access
    
    # Test suite
    test_suite="tests",
    tests_require=EXTRAS_REQUIRE["test"],
)

if __name__ == "__main__":
    # ============================================================================
    # POST-INSTALL CHECKS (Only run during direct setup.py execution)
    # ============================================================================
    import subprocess

    def check_external_tools():
        """Check if external bioinformatics tools are installed."""
        required_tools = ["kraken2", "bracken", "megahit", "diamond"]
        missing_tools = []

        for tool in required_tools:
            try:
                subprocess.run([tool, "--version"],
                             capture_output=True,
                             check=True,
                             timeout=5)
            except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
                missing_tools.append(tool)

        if missing_tools:
            print("\n⚠️  WARNING: Missing external tools:")
            for tool in missing_tools:
                print(f"  - {tool}")
            print("\nInstall via conda:")
            print(f"  conda install -c bioconda {' '.join(missing_tools)}")
            print("\nOr use: metaquest check")

    if "install" in sys.argv or "develop" in sys.argv:
        try:
            check_external_tools()
        except Exception:
            pass  # Don't fail install if check fails
