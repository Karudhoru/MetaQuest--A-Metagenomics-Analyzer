"""
MetaQuest CLI - Command Line Interface (FASTQ-only version)
============================================================

A professional command-line interface for the MetaQuest metagenomics pipeline.
Now exclusively supports FASTQ input for streamlined workflow.

Author: MetaQuest Development Team
"""

import argparse
import sys
import time
from pathlib import Path
from textwrap import dedent
from datetime import datetime

from .exceptions import MetaQuestError
from .io.output_formatter import Colors, OutputFormatter, set_formatter
from .config import __version__

# ============================================================================
# APP METADATA
# ============================================================================

__app_name__ = "MetaQuest"
__tagline__ = "An Integrated Metagenomics Analysis Pipeline"
__author__ = "MetaQuest Development Team"


# ============================================================================
# ARGUMENT GROUP SETUP FUNCTIONS
# ============================================================================

def setup_fastq_validation_args(parser):
    """
    Add FASTQ-specific validation arguments to a parser.
    
    Configures quality control thresholds for FASTQ validation including:
    - Minimum quality scores
    - Sequence count requirements
    - Overrepresentation detection
    """
    val_group = parser.add_argument_group(
        'FASTQ Validation Options',
        'Configure quality control parameters for FASTQ input files'
    )
    
    val_group.add_argument(
        '-q', '--min-quality',
        type=int,
        default=20,
        metavar='SCORE',
        help="Minimum mean Phred quality score threshold (default: 20)"
    )
    
    val_group.add_argument(
        '-n', '--min-sequences',
        type=int,
        default=100,
        metavar='NUM',
        help="Minimum number of sequences required for analysis (default: 100)"
    )
    
    val_group.add_argument(
        '--overrep-threshold',
        type=float,
        default=0.1,
        metavar='PCT',
        help="Fraction threshold for flagging overrepresented sequences (default: 0.1)"
    )

def setup_annotation_args(parser):
    """
    Add functional annotation control arguments to a parser.
    
    Configures the thread allocation for provisional functional annotation.
    
    Pyrodigal runs in metagenomic mode on contigs of at least 200 bp.
    """
    annot_group = parser.add_argument_group(
        'Annotation Control Options',
        'Control Pyrodigal gene prediction and functional annotation'
    )
    
    annot_group.add_argument(
        '--annotation-threads',
        type=int,
        default=8,
        metavar='NUM',
        help="Number of CPU threads for functional annotation (default: 8)"
    )

# ============================================================================
# VALIDATION HANDLER
# ============================================================================

def handle_validation(file_paths, args, formatter):
    """
    Run comprehensive FASTQ file validation and return validation status.
    
    Args:
        file_paths: List of FASTQ file paths to validate
        args: Command line arguments containing validation parameters
        formatter: OutputFormatter instance for consistent logging
        
    Returns:
        bool: True if all files are valid, False if any file fails validation
    """
    formatter.section_header("FILE VALIDATION")
    
    from .io.file_validator import FileValidator

    validator = FileValidator()
    
    # Configure FASTQ validation parameters
    validator.quality_threshold = args.min_quality
    validator.min_sequences = args.min_sequences
    validator.overrep_threshold = args.overrep_threshold
    
    formatter.info("Validation Configuration:")
    formatter.result({
        'Quality threshold': f"Q{args.min_quality}",
        'Minimum sequences': f"{args.min_sequences:,}",
        'Overrepresentation cutoff': f"{args.overrep_threshold:.1%}"
    }, indent=2)

    all_valid = True
    total_files = len(file_paths)
    
    # Validate each file individually
    for i, file_path in enumerate(file_paths, 1):
        if total_files > 1:
            formatter.info(f"Validating file {i}/{total_files}: {Path(file_path).name}")
        else:
            formatter.info(f"Validating file: {Path(file_path).name}")
        
        is_valid, report = validator.validate_and_analyze(file_path, 'fastq')
        
        if not is_valid:
            all_valid = False
            formatter.error(f"Validation failed: {file_path}")
        else:
            formatter.success(f"Validation passed: {file_path}")
    
    print()  # Add spacing
    
    if total_files > 1:
        formatter.info(f"Validation Summary: {total_files} file(s) processed")
    
    if all_valid:
        formatter.success("All files validated successfully")
        formatter.info("Files meet quality standards and are ready for analysis")
    else:
        formatter.error("Validation failed for one or more files")
        formatter.info("Please address quality issues before proceeding")
    
    return all_valid


# ============================================================================
# CUSTOM HELP FORMATTER
# ============================================================================

class CustomHelpFormatter(argparse.RawDescriptionHelpFormatter):
    """
    Enhanced help text formatter with improved readability.
    
    Features:
    - Extended help position for better alignment
    - Consistent column widths
    - Cleaner option string formatting
    """
    
    def __init__(self, prog):
        super().__init__(prog, max_help_position=35, width=100)
    
    def _format_action_invocation(self, action):
        """Format command-line option strings with consistent styling."""
        if not action.option_strings:
            metavar, = self._metavar_formatter(action, action.dest)(1)
            return metavar
        else:
            parts = []
            if action.nargs == 0:
                parts.extend(action.option_strings)
            else:
                default = action.dest.upper()
                args_string = self._format_args(action, default)
                for option_string in action.option_strings:
                    parts.append(f'{option_string}')
                parts[-1] += f' {args_string}'
            return ', '.join(parts)


# ============================================================================
# PARSER CREATION
# ============================================================================

def create_parser():
    """
    Create and configure the main argument parser with all subcommands.
    
    Returns:
        argparse.ArgumentParser: Fully configured parser with all commands and options
    """
    
    description = dedent(f"""
    ═══════════════════════════════════════════════════════════════════════
                           MetaQuest v{__version__}
              {__tagline__}
    ═══════════════════════════════════════════════════════════════════════
    
    OVERVIEW:
        MetaQuest is a research-use pipeline for short-read metagenomic FASTQ
        analysis. The stable surface currently provides input validation,
        taxonomic profiling, provisional functional annotation, and basic
        comparative analysis.
    
    KEY CAPABILITIES:
        • Taxonomic Profiling       - Kraken2/Bracken classification
        • Functional Analysis        - Pyrodigal + provisional DIAMOND workflow
        • Comparative Studies        - Basic group comparison across samples
        • Descriptive Reports        - Text, JSON, tables, and visualizations
    
    QUICK START EXAMPLES:
        # Verify system configuration
        metaquest check
        
        # Validate input files before analysis
        metaquest validate --single reads.fastq.gz
        
        # Run complete analysis pipeline
        metaquest analyze --single reads.fastq.gz -o results/
        
        # Compare multiple samples statistically
        metaquest compare -m metadata.tsv -i sample1/ sample2/ sample3/ -o comparison_results/
    
    COMMON WORKFLOWS:
        Single-end FASTQ:
            metaquest analyze --single reads.fq -o output/
        
    Paired-end FASTQ:
        metaquest analyze --paired R1.fq R2.fq -o output/
        
        Fast taxonomic profiling only:
            metaquest analyze --single reads.fq --skip-annotation
    
    DOCUMENTATION:
        GitHub: https://github.com/Karudhoru/MetaQuest--A-Metagenomics-Analyzer
    """)
    
    parser = argparse.ArgumentParser(
        description=description,
        prog="metaquest",
        formatter_class=CustomHelpFormatter,
        epilog=f"Developed by {__author__}"
    )
    
    # Global options
    parser.add_argument(
        '-v', '--version',
        action='version',
        version=f'MetaQuest version {__version__}'
    )
    
    parser.add_argument(
        '--debug',
        action='store_true',
        help='Enable debug mode with full diagnostic output (includes all commands and stderr)'
    )
    output_mode = parser.add_mutually_exclusive_group()
    output_mode.add_argument(
        '--quiet',
        action='store_true',
        help='Show errors only'
    )
    parser.add_argument(
        '--no-color',
        action='store_true',
        help='Disable ANSI colors and terminal animation'
    )
    output_mode.add_argument(
        '--verbose',
        action='store_true',
        help='Show additional progress and diagnostic context'
    )

    parser.add_argument(
        '--config',
        metavar='FILE',
        help='Path to YAML configuration file (overrides defaults). Generate one with: metaquest init-config'
    )
    
    # Create subcommand parsers
    subparsers = parser.add_subparsers(
        dest='command',
        required=True,
        title='Available Commands',
        description='Choose a command to execute',
        metavar='COMMAND'
    )

    # ========================================================================
    # CHECK COMMAND
    # ========================================================================
    parser_check = subparsers.add_parser(
        'check',
        help='Verify system dependencies and database status',
        description=dedent("""
            System Dependency Verification
            ==============================
            
            Performs comprehensive validation of:
                • Stable command-line tools (Kraken2, Bracken, MEGAHIT,
                  Pyrodigal, and DIAMOND)
                • Stable Python package dependencies
                • Reference database files and indices
            
            Run this before your first analysis to ensure proper configuration.
            Identifies missing dependencies with installation suggestions.
        """),
        formatter_class=CustomHelpFormatter
    )
    parser_check.add_argument(
        '--db-dir',
        type=Path,
        default=None,
        metavar='DIR',
        help='Database root override'
    )

    # ========================================================================
    # VALIDATE COMMAND
    # ========================================================================
    parser_validate = subparsers.add_parser(
        'validate',
        help='Validate FASTQ input files without running analysis',
        description=dedent("""
            FASTQ File Validation
            =====================
            
            Quality control checks for FASTQ files:
                • Format validation (proper headers, sequences)
                • Quality score analysis
                • Sequence count verification
                • Overrepresentation detection
                • File integrity checks
            
            Recommended to run before analysis to catch issues early.
        """),
        formatter_class=CustomHelpFormatter
    )
    
    validate_mode = parser_validate.add_mutually_exclusive_group(required=True)
    validate_mode.add_argument(
        '--single',
        metavar='READS.fastq',
        help="Single-end FASTQ file"
    )
    validate_mode.add_argument(
        '--paired',
        nargs=2,
        metavar=('R1.fastq', 'R2.fastq'),
        help="Paired-end FASTQ files"
    )
    validate_mode.add_argument(
        '--interleaved',
        metavar='INTERLEAVED.fastq',
        help="Interleaved paired-end FASTQ file"
    )
    setup_fastq_validation_args(parser_validate)

    # ========================================================================
    # ANALYZE COMMAND
    # ========================================================================
    parser_analyze = subparsers.add_parser(
        'analyze',
        help='Run complete FASTQ analysis pipeline',
        description=dedent("""
            Complete Metagenomics Analysis Pipeline
            ========================================
            
            Executes the full MetaQuest pipeline on FASTQ reads:
            
            Pipeline Steps:
                1. Taxonomic Classification (Kraken2 + Bracken)
                2. Metagenomic Assembly (MEGAHIT)
                3. Gene Prediction (Pyrodigal metagenomic mode)
                4. Provisional Functional Annotation (DIAMOND database search)
                5. Descriptive Report Generation
            
            Use --skip-annotation for rapid taxonomic classification only.
        """),
        formatter_class=CustomHelpFormatter
    )
    
    # Add all standard arguments
    parser_analyze.add_argument(
        '-o', '--output',
        default='results',
        metavar='DIR',
        help="Output directory for results (default: results)"
    )
    parser_analyze.add_argument(
        '--db-dir',
        type=Path,
        default=None,
        metavar='DIR',
        help='Database root override'
    )
    parser_analyze.add_argument(
        '--skip-validation',
        action='store_true',
        help="Skip input file validation (not recommended for production use)"
    )
    parser_analyze.add_argument(
        '--skip-annotation',
        action='store_true',
        help="Skip annotation and functional analysis for faster taxonomic-only results"
    )
    
    # Add argument groups
    setup_annotation_args(parser_analyze)
    setup_fastq_validation_args(parser_analyze)
    
    fastq_mode = parser_analyze.add_mutually_exclusive_group(required=True)
    fastq_mode.add_argument(
        '--single',
        metavar='READS.fastq',
        help="Single-end FASTQ file"
    )
    fastq_mode.add_argument(
        '--paired',
        nargs=2,
        metavar=('R1.fastq', 'R2.fastq'),
        help="Paired-end FASTQ files (forward and reverse reads)"
    )
    fastq_mode.add_argument(
        '--interleaved',
        metavar='INTERLEAVED.fastq',
        help="Interleaved paired-end FASTQ file"
    )

    # ========================================================================
    # COMPARE COMMAND
    # ========================================================================
    parser_compare = subparsers.add_parser(
        'compare',
        help='Perform comparative analysis across samples',
        description=dedent("""
            Comparative Metagenomic Analysis
            =================================
            
            Statistical comparison of multiple MetaQuest results:
            
            Diversity Analysis:
                • Alpha diversity metrics (Shannon, Simpson, Chao1)
                • Rarefaction curves
                • Species richness estimation
            
            Beta Diversity:
                • Bray-Curtis dissimilarity
                • Principal Coordinate Analysis (PCoA)
                • PERMANOVA testing
            
            Differential Abundance:
                • Statistical testing between groups
                • Effect size calculation
                • Multiple testing correction
            
            Requirements:
                • Multiple MetaQuest output directories
                • Metadata file (TSV) with sample-to-group mappings
            
            Metadata Format:
                sample_id    group    [optional_columns]
                sample1      control  ...
                sample2      control  ...
                sample3      treatment ...
        """),
        formatter_class=CustomHelpFormatter
    )
    parser_compare.add_argument(
        '-i', '--inputs',
        nargs='+',
        required=True,
        metavar='DIR',
        help="MetaQuest output directories to compare (space-separated list)"
    )
    parser_compare.add_argument(
        '-m', '--metadata',
        required=True,
        metavar='FILE',
        help="Metadata file (TSV format) linking sample IDs to experimental groups"
    )
    parser_compare.add_argument(
        '-o', '--output',
        default='comparison_results',
        metavar='DIR',
        help="Output directory for comparison results (default: comparison_results)"
    )

    # ========================================================================
    # SETUP-DB COMMAND
    # ========================================================================
    parser_setup_db = subparsers.add_parser(
        'setup-db',
        help='Download and configure required reference databases',
        description=dedent("""
            Database Setup Utility
            ======================
            
            Installs versioned reference data outside the source repository.
            Use --list to inspect availability and size before downloading.

            Example:
                metaquest setup-db --db-dir /data/metaquest --database taxonomy
        """),
        formatter_class=CustomHelpFormatter
    )

    # ========================================================================
    # INIT-CONFIG COMMAND
    # ========================================================================
    parser_init_config = subparsers.add_parser(
        'init-config',
        help='Generate a default configuration YAML file',
        description=dedent("""
            Generate Default Configuration
            ===============================

            Creates a metaquest.yaml configuration file with all default values.
            Edit this file to customize pipeline behavior without CLI flags.

            Usage:
                metaquest init-config
                metaquest init-config -o my_config.yaml
        """),
        formatter_class=CustomHelpFormatter
    )
    parser_init_config.add_argument(
        '-o', '--output',
        default='metaquest.yaml',
        metavar='FILE',
        help='Output path for the config file (default: metaquest.yaml)'
    )

    db_options = parser_setup_db.add_argument_group('Database Selection')
    db_options.add_argument(
        '--db-dir',
        type=Path,
        default=None,
        metavar='DIR',
        help='Database root (default: METAQUEST_DB_DIR, config, or ./databases)'
    )
    db_options.add_argument(
        '--database',
        choices=('taxonomy', 'functional', 'amr', 'virulence'),
        metavar='NAME',
        help='Install one database'
    )
    db_options.add_argument(
        '--list',
        action='store_true',
        help='List database releases, sizes, and installation status'
    )
    db_options.add_argument(
        '--force',
        action='store_true',
        help='Replace an existing invalid or older installation'
    )

    return parser

# ============================================================================
# MAIN ENTRY POINT
# ============================================================================

def main():
    """
    Main entry point for the MetaQuest CLI.
    
    Orchestrates:
        • Argument parsing and validation
        • Output formatter initialization
        • Command execution routing
        • Error handling and user feedback
        • Runtime tracking and reporting
    """
    
    parser = create_parser()
    args = parser.parse_args()
    
    # Determine verbosity level from command-line flags
    verbosity = (
        'debug' if args.debug else
        'verbose' if args.verbose else
        'silent' if args.quiet else
        'standard'
    )
    
    # Setup log file path if output directory specified (not for init-config)
    log_file = None
    if hasattr(args, 'output') and args.output and args.command not in ('init-config', 'check'):
        log_file = Path(args.output) / 'metaquest.log'
    
    # Initialize global output formatter
    formatter = OutputFormatter(verbosity=verbosity, log_file=log_file)
    if args.no_color:
        formatter.colors_enabled = False
        Colors.disable()
    set_formatter(formatter)
    
    # Show banner (unless version or help requested)
    if not any(arg in sys.argv for arg in ['-v', '--version', '-h', '--help']):
        formatter.banner(__app_name__, __version__, __tagline__)

    start_time = time.time()

    try:
        # ====================================================================
        # COMMAND: CHECK
        # ====================================================================
        if args.command == 'check':
            from .io.utils import run_system_check

            formatter.section_header("SYSTEM DEPENDENCY CHECK")
            formatter.info("Verifying installation of required tools and databases...")
            
            from .settings import load_config
            run_system_check(
                formatter,
                config=load_config(
                    Path(args.config) if args.config else None,
                    db_dir=args.db_dir,
                ),
            )
            
            formatter.success("System check completed successfully")
            formatter.info("All required dependencies are properly configured")
            sys.exit(0)

        # ====================================================================
        # FILE PATH VALIDATION (for analyze and validate commands)
        # ====================================================================
        if args.command in ['analyze', 'validate']:
            file_paths = []
            
            if args.single:
                file_paths = [args.single]
            elif args.paired:
                file_paths = args.paired
            elif args.interleaved:
                file_paths = [args.interleaved]

            # Check file existence before proceeding
            for f in file_paths:
                if not Path(f).exists():
                    formatter.error(
                        f"Input file not found: {f}",
                        solutions=[
                            "Check the file path for typos",
                            "Verify the file exists in the current directory",
                            "Use absolute path instead of relative path",
                            "Check file permissions"
                        ]
                    )
                    sys.exit(1)

        # ====================================================================
        # COMMAND: VALIDATE
        # ====================================================================
        if args.command == 'validate':
            is_valid = handle_validation(file_paths, args, formatter)
            
            if is_valid:
                formatter.success("All files validated successfully")
                formatter.info("Files meet quality standards and are ready for analysis")
            else:
                formatter.error("Validation failed for one or more files")
                formatter.info("Please address the quality issues before proceeding with analysis")
            
            sys.exit(0 if is_valid else 1)
        
        # ====================================================================
        # COMMAND: ANALYZE
        # ====================================================================
        elif args.command == 'analyze':
            from .core.analysis import run_analysis
            from .io.utils import run_system_check

            formatter.section_header("ANALYSIS PIPELINE INITIALIZATION")

            # System check
            formatter.info("Performing system dependency check...")
            with formatter.spinner("Verifying required tools and databases"):
                from .settings import load_config
                check_config = load_config(
                    Path(args.config) if args.config else None,
                    db_dir=args.db_dir,
                )
                run_system_check(
                    formatter,
                    config=check_config,
                    taxonomy_only=args.skip_annotation,
                    require_interleaved=bool(args.interleaved),
                )
            formatter.success("System check passed - all dependencies available")

            # File validation (unless skipped)
            if not args.skip_validation:
                is_valid = handle_validation(file_paths, args, formatter)
                if not is_valid:
                    formatter.error(
                        "Analysis aborted due to validation failure",
                        solutions=[
                            "Fix validation errors reported above",
                            "Improve input file quality if needed",
                            "Use --skip-validation flag to bypass (not recommended)"
                        ]
                    )
                    sys.exit(1)
                formatter.success("Input validation completed - files meet quality standards")
            else:
                formatter.warning("Skipping file validation (--skip-validation flag detected)")
                formatter.info("Proceeding with unvalidated data")
            
            formatter.section_header("ASSEMBLY CONFIGURATION")
            formatter.result({
                'Read type': 'Single-end' if args.single else 'Paired-end' if args.paired else 'Interleaved',
                'Assembler': 'MEGAHIT'
            })

            # Display annotation settings if not skipping annotation
            if not args.skip_annotation:
                formatter.section_header("ANNOTATION SETTINGS")
                formatter.result({
                    'Gene prediction': 'Pyrodigal metagenomic mode',
                    'Min contig length': '200 bp',
                    'Annotation threads': str(args.annotation_threads)
                })
                
            # Run analysis pipeline
            formatter.section_header("FASTQ ANALYSIS")
            formatter.info(f"Input: {', '.join(file_paths)}")
            formatter.info(f"Output: {args.output}")
            formatter.info("Initiating analysis pipeline...\n")
            
            run_analysis(file_paths, args.output, args)

            # Show elapsed time and results location
            elapsed = time.time() - start_time
            formatter.success(f"Analysis completed in {formatter._format_time(elapsed)}")
            formatter.info(f"Results: {Path(args.output) / 'analysis_summary.json'}")

        # ====================================================================
        # COMMAND: COMPARE
        # ====================================================================
        elif args.command == 'compare':
            from .core.comparative_analysis import run_comparison

            formatter.section_header("COMPARATIVE ANALYSIS")
            
            formatter.result({
                'Samples to compare': str(len(args.inputs)),
                'Metadata file': args.metadata,
                'Output directory': args.output
            })
            
            formatter.info("Initiating comparative statistical analysis...")
            run_comparison(args.inputs, args.metadata, args.output)
            
            elapsed = time.time() - start_time
            formatter.success(f"Comparison completed in {formatter._format_time(elapsed)}")
            formatter.info(f"Results: {args.output}")

        # ====================================================================
        # COMMAND: SETUP-DB
        # ====================================================================
        elif args.command == 'setup-db':
            from .database_manager import DATABASES, database_rows, install_database
            from .settings import load_config

            formatter.section_header("DATABASE SETUP")

            configured = load_config(Path(args.config)) if args.config else load_config()
            db_root = args.db_dir or configured.databases.base_dir
            formatter.result({
                'Database root': str(Path(db_root).expanduser().resolve()),
            })

            if args.list:
                for row in database_rows(db_root):
                    formatter.info(
                        f"{row['Database']}: {row['Status']} | "
                        f"{row['Release']} | {row['Size']}"
                    )
                sys.exit(0)

            if not args.database:
                formatter.error(
                    "No database selected",
                    solutions=[
                        "Use --list to inspect available databases",
                        "Use --database taxonomy to install Kraken2 Standard-8",
                    ]
                )
                sys.exit(1)

            selected = DATABASES[args.database]
            if selected.download_gb is not None:
                formatter.warning(
                    f"The {args.database} download is approximately "
                    f"{selected.download_gb:g} GB"
                )
            target = install_database(
                args.database,
                db_root,
                force=args.force,
                notify=formatter.info,
            )
            elapsed = time.time() - start_time
            formatter.success(
                f"Database setup completed in {formatter._format_time(elapsed)}"
            )
            formatter.info(f"Installed at: {target}")

        # ====================================================================
        # COMMAND: INIT-CONFIG
        # ====================================================================
        elif args.command == 'init-config':
            import shutil
            from .settings import _DEFAULT_CONFIG_PATH

            output_path = Path(args.output)
            if output_path.exists():
                formatter.warning(f"File already exists: {output_path}")
                formatter.info("Overwriting with default configuration")

            shutil.copy2(_DEFAULT_CONFIG_PATH, output_path)
            formatter.success(f"Configuration file created: {output_path}")
            formatter.info("Edit this file to customize pipeline behavior.")
            formatter.info("Usage: metaquest analyze --config metaquest.yaml --single reads.fq -o results/")

    except KeyboardInterrupt:
        formatter.warning("Operation interrupted by user (Ctrl+C)")
        sys.exit(130)

    except MetaQuestError as e:
        formatter.error(str(e))
        if e.suggestions:
            formatter.info("Suggestions:")
            for s in e.suggestions:
                formatter.info(f"  - {s}", indent=2)
        sys.exit(1)

    except Exception as e:
        formatter.error(
            f"Unexpected error: {e}",
            solutions=[
                "Check the log file for details",
                "Verify input file integrity",
                "Ensure sufficient disk space and memory",
                "Report issue at GitHub if problem persists"
            ],
            doc_link="https://github.com/Karudhoru/MetaQuest--A-Metagenomics-Analyzer/issues"
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
