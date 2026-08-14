"""
MetaQuest Configuration Module
===============================

Centralized configuration management for the MetaQuest metagenomic analysis pipeline.
All database paths, analysis parameters, thresholds, and pathogen definitions are 
defined here for consistency across modules.

Version: 2.0.0-alpha.1
Release Date: 2026-08-14
"""

import os
from pathlib import Path
from typing import Dict, Any, Optional
import sys
import hashlib
import json
from datetime import datetime

# ============================================================================
# VERSION AND METADATA
# ============================================================================

__version__ = "2.0.0-alpha.1"
__release_date__ = "2026-08-14"
__pipeline_modules__ = [
    "core.analysis",
    "core.taxonomic_analysis", 
    "core.pathogen_analysis",
    "core.functional_analysis",
    "core.comparative_analysis"
]

# ============================================================================
# PATHOGEN CONFIGURATION - Central Data Source
# ============================================================================

def _load_pathogen_config() -> Dict[str, Any]:
    """
    Load pathogen configuration from JSON file.
    """
    # Use relative path from this file location
    config_path = Path(__file__).parent / "pathogen_config.json"
    
    if not config_path.exists():
        # Fallback for development environments or alternative layouts
        config_path = Path(__file__).parent / "config" / "pathogen_config.json"

    if not config_path.exists():
        raise FileNotFoundError(
            f"Critical configuration file missing: {config_path}\n"
            f"Please ensure pathogen_config.json exists in the package config directory."
        )
    
    try:
        with open(config_path, 'r') as f:
            return json.load(f)
    except json.JSONDecodeError as e:
        raise RuntimeError(f"Failed to parse pathogen_config.json: {e}")
    except Exception as e:
        raise RuntimeError(f"Failed to load pathogen configuration: {e}")


# Load pathogen configuration at module initialization
PATHOGEN_CONFIG = _load_pathogen_config()

# ============================================================================
# PATHOGEN CONFIGURATION - Convenience Accessors
# ============================================================================
# (These remain unchanged as they are just dictionary lookups)
CRITICAL_MOTIFS = PATHOGEN_CONFIG.get('critical_motifs', {})
COMMENSAL_SPECIES = PATHOGEN_CONFIG.get('commensal_species', [])
PATHOGENIC_GENERA = PATHOGEN_CONFIG.get('pathogenic_genera', {})
GENE_RISK_SCORES = PATHOGEN_CONFIG.get('gene_risk_scores', {})
CLINICAL_NOTES = PATHOGEN_CONFIG.get('clinical_notes', {})
RISK_WEIGHTS = PATHOGEN_CONFIG.get('risk_weights', {})
RISK_MULTIPLIERS = PATHOGEN_CONFIG.get('risk_multipliers', {})
CONFIDENCE_THRESHOLDS = PATHOGEN_CONFIG.get('confidence_thresholds', {})
FILTERING_THRESHOLDS = PATHOGEN_CONFIG.get('filtering_thresholds', {})
RISK_LEVEL_THRESHOLDS = PATHOGEN_CONFIG.get('risk_level_thresholds', {})
DIVERSITY_THRESHOLDS = PATHOGEN_CONFIG.get('diversity_interpretation', {})

# ============================================================================
# DATABASE PATHS - FIXED PORTABILITY
# ============================================================================

# Base database directory (can be overridden by environment variable)
DB_DIR = Path("./databases")

# Kraken2 taxonomic classification database
KRAKEN_DB = DB_DIR 

# Custom pathogen marker database (Diamond format)
PATHOGEN_DB_CUSTOM = DB_DIR / "metaquest_pathogen_markers.dmnd"

# SwissProt + COG functional annotation database (Diamond format)
SWISSPROT_DB = DB_DIR / "SwissProt_COG_db.dmnd"

# Machine learning model artifacts
MODEL_ARTIFACTS_DIR = Path(__file__).parent / "ml" / "model_artifacts"


# ============================================================================
# MACHINE LEARNING CONFIGURATION
# ============================================================================

ML_CONFIG = {
    'enable_ml_prediction': False,
    'confidence_threshold': 0.7,
    'batch_size': 100,
    'feature_cache': True,
    'random_seed': 42,
    'model_type': 'random_forest',
    'n_estimators': 100,
}

# ============================================================================
# ASSEMBLY CONFIGURATION
# ============================================================================

ASSEMBLY_CONFIG = {
    'default_assembler': 'megahit',
    'alternative_assembler': 'spades',
    'min_contig_length': 500,
    'memory_limit_gb': 8,
    'threads': 8,
    'kmer_min': 21,
    'kmer_max': 141,
    'kmer_step': 12,
}

# ============================================================================
# ANNOTATION CONFIGURATION
# ============================================================================

ANNOTATION_CONFIG = {
    'tool': 'prokka',
    'threads': 8,
    'evalue': 1e-6,
    'kill_tbl2asn': True,
    'tbl2asn_timeout': 30,
    'rfam': False,
    'note': 'Prokka default minimum contig length (200bp) is used'
}

# ============================================================================
# PATHOGEN DETECTION CONFIGURATION
# ============================================================================

PATHOGEN_CONFIG_PARAMS = {
    'min_fragment_length': 80,
    'min_identity': 80.0,
    'min_query_coverage': 0.7,
    'min_subject_coverage': 0.7,
    'confidence_threshold_high': 70,
    'confidence_threshold_medium': 45,
    'evalue_threshold': 1e-5,
}

# ============================================================================
# COMPARATIVE ANALYSIS CONFIGURATION
# ============================================================================

COMPARATIVE_CONFIG = {
    'default_normalization': 'tss',
    'normalization_options': ['tss', 'clr', 'none'],
    'min_prevalence': 0.10,
    'random_seed': 42,
    'permutations': 999,
    'alpha_level': 0.05,
    'fdr_method': 'fdr_bh',
}

# ============================================================================
# VISUALIZATION CONFIGURATION
# ============================================================================

TAXONOMIC_COLORS = ['#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3', 
                   '#FDB462', '#B3DE69', '#FCCDE5', '#D9D9D9', '#BC80BD']
FUNCTIONAL_COLORS = ['#FBB4AE', '#B3CDE3', '#CCEBC5', '#DECBE4', '#FED9A6',
                     '#FFE5CC', '#E5D8BD', '#FDDAEC', '#F2F2F2']

PATHOGEN_COLORS = {
    'Pathogenic': '#FF6B6B',
    'Non-pathogenic': '#4ECDC4',
    'Unknown': '#95A5A6',
    'AMR': '#E74C3C',
    'Virulence': '#F39C12',
    'High': '#E74C3C',
    'Moderate': '#F39C12',
    'Low': '#2ECC71',
}

RISK_LEVEL_COLORS = {
    'High': '#DC3545',
    'Moderate': '#FFC107',
    'Low': '#28A745',
    'Unknown': '#6C757D',
}

# ============================================================================
# BIOSAFETY LEVEL (BSL) DECLARATIONS
# Source: CDC/NIH Biosafety in Microbiological and Biomedical Laboratories (BMBL)
# ============================================================================

BSL_LEVELS = {
    # BSL-3 Pathogens
    'Yersinia pestis': '3',
    'Bacillus anthracis': '3',
    'Francisella tularensis': '3',
    'Brucella abortus': '3',
    'Brucella melitensis': '3',
    'Brucella suis': '3',
    'Burkholderia mallei': '3',
    'Burkholderia pseudomallei': '3',
    'Mycobacterium tuberculosis': '3',
    'Salmonella enterica subsp. enterica serovar Typhi': '3',
    'Rickettsia rickettsii': '3',
    'Coxiella burnetii': '3',
    'Clostridium botulinum': '3',

    # BSL-2 Pathogens
    'Escherichia coli': '2',
    'Klebsiella pneumoniae': '2',
    'Klebsiella variicola': '2',
    'Klebsiella aerogenes': '2',
    'Klebsiella quasipneumoniae': '2',
    'Klebsiella oxytoca': '2',
    'Salmonella enterica': '2',
    'Salmonella bongori': '2',
    'Shigella dysenteriae': '2',
    'Shigella flexneri': '2',
    'Shigella sonnei': '2',
    'Shigella boydii': '2',
    'Staphylococcus aureus': '2',
    'Staphylococcus haemolyticus': '2',
    'Pseudomonas aeruginosa': '2',
    'Acinetobacter baumannii': '2',
    'Enterococcus faecium': '2',
    'Vibrio cholerae': '2',
    'Listeria monocytogenes': '2',
    'Neisseria meningitidis': '2',
    'Streptococcus pneumoniae': '2',
    'Streptococcus pyogenes': '2',
    'Clostridioides difficile': '2',
    'Campylobacter jejuni': '2',
    'Helicobacter pylori': '2',
    'Burkholderia dolosa': '2',
    'Aeromonas hydrophila': '2',
    'Aeromonas caviae': '2',
    'Proteus mirabilis': '2',
    'Providencia rettgeri': '2',
    'Morganella morganii': '2',
    'Serratia marcescens': '2',
    'Enterobacter cloacae': '2',
    'Enterobacter hormaechei': '2',
    'Citrobacter freundii': '2',
    'Cronobacter sakazakii': '2',
    'Yersinia enterocolitica': '2',
    'Edwardsiella ictaluri': '2',
    'Bacteroides fragilis': '2',
    'Prevotella copri': '2',

    # BSL-1 Organisms (Environmental / Low Risk)
    'Cupriavidus taiwanensis': '1',
    'Bacillus subtilis': '1',
    'Pseudomonas fluorescens': '1',
    'Lactobacillus acidophilus': '1',
    'Bifidobacterium longum': '1',
    'Faecalibacterium prausnitzii': '1',
    'Phocaeicola dorei': '1',
    'Phocaeicola vulgatus': '1',
    'Bacteroides thetaiotaomicron': '1',
    'Bacteroides ovatus': '1',
    'Alistipes onderdonkii': '1',
    'Parabacteroides distasonis': '1',
    'Lachnoclostridium sp. YL32': '1',
    'Enterocloster clostridioformis': '1',
    'Blautia obeum': '1',
    'Roseburia intestinalis': '1',
    'Veillonella parvula': '1',
}

# ============================================================================
# DATABASE AND CONTACT CONFIGURATION
# Update KRAKEN_DB_VERSION to match your actual database build.
# ============================================================================

KRAKEN_DB_VERSION = "Standard-16 (Built: 2024-01-10)"
CONTACT_EMAIL = "metaquest-support@example.org"

# ============================================================================
# MEMORY OPTIMIZATION CONFIGURATION
# ============================================================================

MEMORY_CONFIG = {
    'target_system': '8-16GB',
    'max_sequences_in_memory': 10000,
    'chunk_size': 1000,
    'explicit_cleanup': True,
    'streaming_enabled': True,
}

# ============================================================================
# REPORTING CONFIGURATION
# ============================================================================

REPORT_CONFIG = {
    'output_formats': ['txt', 'json', 'csv'],
    'include_plots': True,
    'include_tables': True,
    'confidence_warnings': True,
    'max_table_rows': 100,
    'decimal_precision': 2,
}

# ============================================================================
# QUALITY CONTROL THRESHOLDS
# ============================================================================

QC_THRESHOLDS = {
    'min_read_quality': 20,
    'min_read_length': 50,
    'min_classification_rate': 0.5,
    'min_species_level_rate': 0.3,
    'max_unclassified_rate': 0.5,
    'min_contig_coverage': 2.0,
}

# ============================================================================
# LOGGING CONFIGURATION
# ============================================================================

LOGGING_CONFIG = {
    'level': 'INFO',
    'format': '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    'date_format': '%Y-%m-%d %H:%M:%S',
    'log_to_file': True,
    'log_to_console': True,
}

# ============================================================================
# PIPELINE REPRODUCIBILITY FUNCTIONS
# ============================================================================

def get_pipeline_hash() -> str:
    """
    Calculate hash of critical pipeline modules for reproducibility.
    """
    try:
        # Assuming this file is inside src/metaquest/
        base_path = Path(__file__).parent
        
        hash_obj = hashlib.sha256()
        
        for module_name in __pipeline_modules__:
            # Convert python module dot notation to path
            # e.g., core.analysis -> core/analysis.py
            rel_path = module_name.replace('.', '/') + '.py'
            module_path = base_path / rel_path
            
            if module_path.exists():
                with open(module_path, 'rb') as f:
                    hash_obj.update(f.read())
        
        return hash_obj.hexdigest()[:16]
    except Exception:
        return "unknown"


def get_analysis_metadata() -> Dict[str, Any]:
    """
    Capture complete analysis environment for reproducibility.
    """
    import numpy as np
    import pandas as pd
    import scipy
    
    try:
        import sklearn
        sklearn_version = sklearn.__version__
    except ImportError:
        sklearn_version = "not installed"
    
    try:
        from Bio import __version__ as bio_version
    except ImportError:
        bio_version = "not installed"
    
    metadata = {
        'metaquest_version': __version__,
        'pipeline_hash': get_pipeline_hash(),
        'release_date': __release_date__,
        'timestamp': datetime.now().isoformat(),
        'environment': {
            'python_version': sys.version.split()[0],
            'python_full': sys.version,
            'platform': sys.platform,
            'numpy_version': np.__version__,
            'pandas_version': pd.__version__,
            'scipy_version': scipy.__version__,
            'sklearn_version': sklearn_version,
            'biopython_version': bio_version,
        },
        'pipeline_modules': __pipeline_modules__,
        'configuration': {
            'ml_enabled': ML_CONFIG['enable_ml_prediction'],
            'assembler': ASSEMBLY_CONFIG['default_assembler'],
            'min_contig_length': ASSEMBLY_CONFIG['min_contig_length'],
            'threads': ASSEMBLY_CONFIG['threads'],
        }
    }
    
    return metadata


def save_analysis_metadata(output_dir: Path, parameters: Optional[Dict] = None) -> Path:
    """
    Save complete analysis metadata to output directory.
    """
    metadata = get_analysis_metadata()
    
    if parameters:
        metadata['parameters'] = parameters
    
    metadata_file = output_dir / "analysis_metadata.json"
    
    with open(metadata_file, 'w') as f:
        json.dump(metadata, f, indent=2)
    
    return metadata_file

# ============================================================================
# CONFIGURATION VALIDATION
# ============================================================================

def validate_config() -> tuple[bool, list]:
    """
    Validate configuration parameters on import.
    """
    errors = []
    
    # Validate numeric parameters
    if MEMORY_CONFIG['max_sequences_in_memory'] < 1000:
        errors.append("max_sequences_in_memory must be >= 1000")
    
    if ASSEMBLY_CONFIG['memory_limit_gb'] < 4:
        errors.append("memory_limit_gb should be >= 4 for reliable assembly")
    
    if ML_CONFIG['confidence_threshold'] < 0 or ML_CONFIG['confidence_threshold'] > 1:
        errors.append("ML confidence_threshold must be between 0 and 1")
    
    # Validate database paths exist 
    # NOTE: We only warn now, because DB setup might happen later
    if not DB_DIR.exists():
        pass # Created dynamically or waiting for setup-db
    
    # Validate pathogen config completeness
    required_pathogen_keys = [
        'critical_motifs', 'commensal_species', 'pathogenic_genera',
        'risk_weights', 'gene_risk_scores'
    ]
    
    for key in required_pathogen_keys:
        if key not in PATHOGEN_CONFIG:
            errors.append(f"Missing required pathogen config section: {key}")
    
    # Validate risk weights sum to 1.0
    if 'risk_weights' in PATHOGEN_CONFIG:
        weight_sum = sum(PATHOGEN_CONFIG['risk_weights'].values())
        if abs(weight_sum - 1.0) > 0.01:
            errors.append(f"Risk weights must sum to 1.0 (currently {weight_sum:.3f})")
    
    return len(errors) == 0, errors


def initialize_config(verbose: bool = False) -> bool:
    """
    Initialize and validate configuration.
    """
    is_valid, errors = validate_config()
    
    if verbose:
        print("=" * 70)
        print(f"MetaQuest v{__version__} Configuration")
        print("=" * 70)
        print(f"Release Date:   {__release_date__}")
        print(f"Pipeline Hash:  {get_pipeline_hash()}")
        print(f"Config Status:  {'VALID' if is_valid else 'INVALID'}")
        print("=" * 70)
        print()
    
    if not is_valid:
        print("[!]  Configuration Validation Errors:")
        for error in errors:
            print(f"  ✗ {error}")
        print()
        
    return is_valid


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def get_config_summary() -> Dict[str, Any]:
    """
    Get summary of current configuration for logging/reporting.
    """
    return {
        'version': __version__,
        'databases': {
            'kraken': str(KRAKEN_DB),
            'pathogen': str(PATHOGEN_DB_CUSTOM),
            'swissprot': str(SWISSPROT_DB),
        },
        'thresholds': {
            'min_identity': PATHOGEN_CONFIG_PARAMS['min_identity'],
            'min_coverage': PATHOGEN_CONFIG_PARAMS['min_query_coverage'],
            'confidence_high': CONFIDENCE_THRESHOLDS['high_confidence'],
        },
        'risk_weights': RISK_WEIGHTS,
        'ml_enabled': ML_CONFIG['enable_ml_prediction'],
    }


def export_config(output_file: Path) -> None:
    """
    Export current configuration to JSON file for reproducibility.
    """
    config_export = {
        'version': __version__,
        'timestamp': datetime.now().isoformat(),
        'pipeline_hash': get_pipeline_hash(),
        'databases': {
            'db_dir': str(DB_DIR),
            'kraken_db': str(KRAKEN_DB),
            'pathogen_db': str(PATHOGEN_DB_CUSTOM),
            'swissprot_db': str(SWISSPROT_DB),
        },
        'ml_config': ML_CONFIG,
        'assembly_config': ASSEMBLY_CONFIG,
        'annotation_config': ANNOTATION_CONFIG,
        'pathogen_params': PATHOGEN_CONFIG_PARAMS,
        'risk_weights': RISK_WEIGHTS,
        'risk_multipliers': RISK_MULTIPLIERS,
        'thresholds': CONFIDENCE_THRESHOLDS,
    }
    
    with open(output_file, 'w') as f:
        json.dump(config_export, f, indent=2)

# ============================================================================
# AUTO-VALIDATION ON IMPORT
# ============================================================================

# Run validation on import but don't fail (just warn)
_config_valid, _config_errors = validate_config()

if not _config_valid and __name__ != "__main__":
    import warnings
    warnings.warn(
        f"MetaQuest configuration has {len(_config_errors)} validation error(s). "
        f"Run metaquest.config.initialize_config(verbose=True) for details.",
        UserWarning
    )

# ============================================================================
# MODULE EXPORTS
# ============================================================================

__all__ = [
    # Version info
    '__version__',
    '__release_date__',
    
    # Pathogen configuration
    'PATHOGEN_CONFIG',
    'CRITICAL_MOTIFS',
    'COMMENSAL_SPECIES',
    'PATHOGENIC_GENERA',
    'GENE_RISK_SCORES',
    'CLINICAL_NOTES',
    'RISK_WEIGHTS',
    'RISK_MULTIPLIERS',
    'CONFIDENCE_THRESHOLDS',
    'FILTERING_THRESHOLDS',
    'RISK_LEVEL_THRESHOLDS',
    'DIVERSITY_THRESHOLDS',
    
    # Database paths
    'DB_DIR',
    'KRAKEN_DB',
    'PATHOGEN_DB_CUSTOM',
    'SWISSPROT_DB',
    'MODEL_ARTIFACTS_DIR',
    
    # Configuration dictionaries
    'ML_CONFIG',
    'ASSEMBLY_CONFIG',
    'ANNOTATION_CONFIG',
    'PATHOGEN_CONFIG_PARAMS',
    'COMPARATIVE_CONFIG',
    'MEMORY_CONFIG',
    'REPORT_CONFIG',
    'QC_THRESHOLDS',
    
    # Colors
    'TAXONOMIC_COLORS',
    'FUNCTIONAL_COLORS',
    'PATHOGEN_COLORS',
    'RISK_LEVEL_COLORS',
    
    # Functions
    'get_pipeline_hash',
    'get_analysis_metadata',
    'save_analysis_metadata',
    'validate_config',
    'initialize_config',
    'get_config_summary',
    'export_config',

    # Biosafety and reporting
    'BSL_LEVELS',
    'KRAKEN_DB_VERSION',
    'CONTACT_EMAIL',
]
