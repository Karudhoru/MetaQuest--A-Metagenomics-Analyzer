"""
MetaQuest Utilities Module - CLEANED v2.0
==========================================

Core utility functions for assembly, validation, and system checks.
Progress parsing removed - using spinners instead.
"""

import subprocess
import os
from pathlib import Path
import gc
import re
import importlib
from .output_formatter import get_formatter

formatter = get_formatter()

def _check_command(name, version_cmd="--version"):
    """Helper to check for a single command-line tool."""
    try:
        if name == 'seqtk':
            subprocess.run([name], capture_output=True, text=True, check=False)
            return True
        
        cmd = name.split() + [version_cmd] if version_cmd else name.split()
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False

def run_system_check(
    formatter,
    *,
    config=None,
    taxonomy_only: bool = False,
    require_interleaved: bool = False,
):
    """
    Runs a comprehensive check of all dependencies using the formatter for output.
    """
    errors = []
    warnings = []
    
    formatter.info("Performing system-wide checks...")

    # Command-Line Tools
    formatter.substep("Checking command-line tools...")
    tools = {'kraken2': '--version', 'bracken': '--version'}
    if not taxonomy_only:
        tools.update({
            'megahit': '--version',
            'diamond': 'version',
        })
    if require_interleaved:
        tools['reformat.sh'] = None
    for tool, cmd in tools.items():
        if not _check_command(tool, cmd):
            errors.append(f"Tool not found: '{tool}'. Please install it, e.g., via 'conda install -c bioconda {tool}'.")

    # Python Packages
    formatter.substep("Checking Python packages...")
    packages = {
        'Bio': 'biopython',
        'pandas': 'pandas',
        'numpy': 'numpy',
        'plotly': 'plotly',
        'matplotlib': 'matplotlib',
        'pyrodigal': 'pyrodigal',
    }
    for pkg, install_name in packages.items():
        try:
            importlib.import_module(pkg)
        except ImportError:
            errors.append(f"Python package not found: '{install_name}'. Please install it, e.g., 'pip install {install_name}'.")

    # Databases
    formatter.substep("Checking database files...")
    if config is None:
        from metaquest.settings import get_config
        config = get_config()
    required_dbs = {"Kraken2 DB": config.databases.kraken_db / "hash.k2d"}
    if not taxonomy_only:
        required_dbs["Functional annotation DB"] = config.databases.swissprot_cog
    for name, path in required_dbs.items():
        if not path.exists():
            errors.append(f"Required database not found: {name} (expected at {path}).")

    # Final Report
    if warnings:
        for warning in warnings:
            formatter.warning(warning)
    
    if errors:
        formatter.error(
            "System check failed due to missing critical dependencies.",
            solutions=errors + ["Please resolve the issues above before running an analysis."]
        )
        raise SystemExit(1)

def get_tool_version(tool_name: str) -> str:
    """
    Dynamically gets the version of a command-line tool.
    
    Args:
        tool_name: The name of the tool (e.g., 'kraken2', 'diamond').
        
    Returns:
        The version string or 'N/A' if not found.
    """
    version_commands = {
        'kraken2': ['kraken2', '--version'],
        'diamond': ['diamond', '--version'],
        'spades': ['spades.py', '--version'],
        'megahit': ['megahit', '--version']
    }
    
    cmd = version_commands.get(tool_name.lower())
    if not cmd:
        return "N/A"
        
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=10)
        output = result.stdout + result.stderr
        
        match = re.search(r'(\d+\.\d+(\.\d+)?)', output)
        if match:
            return match.group(1)
        return "Unknown"
    except (FileNotFoundError, subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return "Not Found"


def split_interleaved(interleaved_fastq: str, output_dir: Path, formatter) -> list:
    """
    Splits an interleaved FASTQ into two files using reformat.sh.
    """
    r1 = output_dir / "split_R1.fastq.gz"
    r2 = output_dir / "split_R2.fastq.gz"

    cmd = [
        "reformat.sh", f"in={interleaved_fastq}",
        f"out1={r1}", f"out2={r2}", "overwrite=true",
    ]
    
    returncode, _, stderr = formatter.run_subprocess(cmd, "Splitting interleaved file", show_command=True)

    if returncode != 0:
        formatter.error(
            "reformat.sh failed to split the file.",
            solutions=["Ensure BBTools is installed and the input file is a valid FASTQ."]
        )
        raise RuntimeError(f"reformat.sh failed. Stderr: {stderr}")
    
    return [str(r1), str(r2)]

def explicit_cleanup(*objects):
    """Explicitly delete objects and force garbage collection."""
    for obj in objects:
        if obj is not None:
            del obj
    gc.collect()
