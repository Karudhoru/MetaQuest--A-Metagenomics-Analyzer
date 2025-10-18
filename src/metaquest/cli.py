"""
MetaQuest CLI - Command Line Interface
=======================================

A professional command-line interface for the MetaQuest metagenomics pipeline.
Provides intuitive commands for validation, analysis, and comparative studies.

Author: MetaQuest Development Team
"""

import argparse
import sys
import time
from pathlib import Path
from textwrap import dedent
from datetime import datetime

from .core.analysis import run_analysis
from .core.comparative_analysis import run_comparison 
from .io.file_validator import FileValidator
from .io.utils import run_system_check
from .io.output_formatter import OutputFormatter, set_formatter

# Version and metadata
__version__ = "4.0.0"
__app_name__ = "MetaQuest"
__tagline__ = "An Integrated Metagenomics Analysis Pipeline"
__author__ = "MetaQuest Development Team"


def setup_fastq_validation_args(parser):
    """Add FASTQ-specific validation arguments to a parser."""
    val_group = parser.add_argument_group(
        'FASTQ Validation Options',
        'Configure quality control parameters for FASTQ files'
    )
    val_group.add_argument(
        '-q', '--min-quality',
        type=int,
        default=20,
        metavar='SCORE',
        help="Minimum mean quality score (default: 20)"
    )
    val_group.add_argument(
        '-n', '--min-sequences',
        type=int,
        default=100,
        metavar='NUM',
        help="Minimum number of sequences (default: 100)"
    )
    val_group.add_argument(
        '--overrep-threshold',
        type=float,
        default=0.1,
        metavar='PCT',
        help="Threshold to flag overrepresented sequences (default: 0.1)"
    )


def setup_annotation_args(parser):
    """Add annotation control arguments to a parser."""
    annot_group = parser.add_argument_group(
        'Annotation Control Options',
        'Fine-tune the functional annotation process'
    )
    annot_group.add_argument(
        '--no-filter-contigs',
        dest='filter_contigs',
        action='store_false',
        default=True,
        help="Skip contig filtering before Prokka annotation"
    )
    annot_group.add_argument(
        '--min-contig-length',
        type=int,
        default=1000,
        metavar='BP',
        help="Minimum contig length for filtering (default: 1000)"
    )
    annot_group.add_argument(
        '--no-kill-tbl2asn',
        dest='kill_tbl2asn',
        action='store_false',
        default=True,
        help="Don't auto-kill long-running tbl2asn process"
    )
    annot_group.add_argument(
        '--tbl2asn-timeout',
        type=int,
        default=300,
        metavar='SEC',
        help="Timeout for tbl2asn process in seconds (default: 300)"
    )
    annot_group.add_argument(
        '--annotation-threads',
        type=int,
        default=8,
        metavar='NUM',
        help="Number of threads for functional annotation (default: 8)"
    )


def handle_validation(file_paths, file_type, args, formatter):
    """
    Run file validation and return validation status.
    
    Args:
        file_paths: List of file paths to validate
        file_type: Type of files ('fasta' or 'fastq')
        args: Command line arguments
        formatter: OutputFormatter instance
        
    Returns:
        bool: True if all files are valid, False otherwise
    """
    formatter.section_header("File Validation")
    
    validator = FileValidator()
    
    if file_type == 'fastq':
        validator.quality_threshold = args.min_quality
        validator.min_sequences = args.min_sequences
        validator.overrep_threshold = args.overrep_threshold

    all_valid = True
    total_files = len(file_paths)
    
    for i, file_path in enumerate(file_paths, 1):
        formatter.info(f"Validating file {i}/{total_files}: {file_path}")
        is_valid, _ = validator.validate_and_analyze(file_path, file_type)
        
        if not is_valid:
            all_valid = False
            formatter.error(f"Validation failed for {file_path}")
        else:
            formatter.success(f"Validation passed for {file_path}")
    
    return all_valid


class CustomHelpFormatter(argparse.RawDescriptionHelpFormatter):
    """Custom formatter for better help text formatting."""
    
    def __init__(self, prog):
        super().__init__(prog, max_help_position=35, width=100)
    
    def _format_action_invocation(self, action):
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


def create_parser():
    """Create and configure the argument parser."""
    
    description = dedent(f"""
    MetaQuest v{__version__} - {__tagline__}
    
    A comprehensive pipeline for metagenomic analysis featuring:
      * Taxonomic classification (FASTQ/FASTA)
      * Pathogen detection with risk assessment
      * Machine learning-based predictions
      * Statistical comparative analysis
      * Interactive HTML reports
    
    Quick Start:
      metaquest check                              # Verify system setup
      metaquest validate fastq --single reads.fq   # Validate input files
      metaquest analyze fastq --single reads.fq    # Run full analysis
      metaquest compare -i sample1/ sample2/ -m metadata.tsv
    
    Documentation:
      GitHub: https://github.com/Karudhoru/MetaQuest--A-Metagenomics-Analyzer
    """)
    
    parser = argparse.ArgumentParser(
        description=description,
        prog="metaquest",
        formatter_class=CustomHelpFormatter,
        epilog=f"Developed by {__author__}"
    )
    
    parser.add_argument(
        '-v', '--version',
        action='version',
        version=f'MetaQuest version {__version__}'
    )
    
    # Global verbosity options
    parser.add_argument(
        '--debug',
        action='store_true',
        help='Full diagnostic output (includes all tool commands and output)'
    )
    
    subparsers = parser.add_subparsers(
        dest='command',
        required=True,
        title='Available Commands',
        description='Choose a command to run',
        metavar='COMMAND'
    )

    # ========== CHECK COMMAND ==========
    parser_check = subparsers.add_parser(
        'check',
        help='Verify system dependencies and database status',
        description='Check if all required dependencies and databases are properly configured.',
        formatter_class=CustomHelpFormatter
    )

    # ========== VALIDATE COMMAND ==========
    parser_validate = subparsers.add_parser(
        'validate',
        help='Validate input files without running analysis',
        description='Perform quality control checks on input files.',
        formatter_class=CustomHelpFormatter
    )
    
    validate_subparsers = parser_validate.add_subparsers(
        dest='type',
        required=True,
        title='File Types',
        metavar='TYPE'
    )

    # Validate FASTA
    validate_fasta = validate_subparsers.add_parser(
        'fasta',
        help='Validate FASTA file(s)',
        formatter_class=CustomHelpFormatter
    )
    validate_fasta.add_argument(
        'input_file',
        help="Path to the FASTA file"
    )

    # Validate FASTQ
    validate_fastq = validate_subparsers.add_parser(
        'fastq',
        help='Validate FASTQ file(s)',
        formatter_class=CustomHelpFormatter
    )
    validate_mode = validate_fastq.add_mutually_exclusive_group(required=True)
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
    setup_fastq_validation_args(validate_fastq)

    # ========== ANALYZE COMMAND ==========
    parser_analyze = subparsers.add_parser(
        'analyze',
        help='Run complete analysis pipeline',
        description='Execute the full MetaQuest analysis pipeline on your data.',
        formatter_class=CustomHelpFormatter
    )
    
    # Shared arguments for analysis
    analysis_parent = argparse.ArgumentParser(add_help=False)
    analysis_parent.add_argument(
        '-o', '--output',
        default='results',
        metavar='DIR',
        help="Output directory (default: results)"
    )
    analysis_parent.add_argument(
        '--skip-validation',
        action='store_true',
        help="Skip input validation (not recommended)"
    )
    analysis_parent.add_argument(
        '--skip-annotation',
        action='store_true',
        help="Skip annotation for faster taxonomic-only analysis"
    )
    
    # Add annotation control arguments to shared parent
    setup_annotation_args(analysis_parent)
    
    analysis_subparsers = parser_analyze.add_subparsers(
        dest='type',
        required=True,
        title='Analysis Types',
        metavar='TYPE'
    )
    
    # Analyze FASTA
    analyze_fasta = analysis_subparsers.add_parser(
        'fasta',
        help='Analyze FASTA sequences',
        parents=[analysis_parent],
        formatter_class=CustomHelpFormatter
    )
    analyze_fasta.add_argument(
        'input_file',
        help="Path to the FASTA file"
    )
    analyze_fasta.add_argument(
        '-s', '--blast-sample-size',
        type=int,
        default=50,
        metavar='NUM',
        help="Number of sequences to BLAST (default: 50)"
    )
    
    # Analyze FASTQ
    analyze_fastq = analysis_subparsers.add_parser(
        'fastq',
        help='Analyze FASTQ sequences',
        parents=[analysis_parent],
        formatter_class=CustomHelpFormatter
    )
    setup_fastq_validation_args(analyze_fastq)
    
    fastq_mode = analyze_fastq.add_mutually_exclusive_group(required=True)
    fastq_mode.add_argument(
        '--single',
        metavar='READS.fastq',
        help="Single-end FASTQ file"
    )
    fastq_mode.add_argument(
        '--paired',
        nargs=2,
        metavar=('R1.fastq', 'R2.fastq'),
        help="Paired-end FASTQ files"
    )
    fastq_mode.add_argument(
        '--interleaved',
        metavar='INTERLEAVED.fastq',
        help="Interleaved paired-end FASTQ file"
    )

    # ========== COMPARE COMMAND ==========
    parser_compare = subparsers.add_parser(
        'compare',
        help='Perform comparative analysis across samples',
        description='Statistical comparison of multiple MetaQuest analysis results.',
        formatter_class=CustomHelpFormatter
    )
    parser_compare.add_argument(
        '-i', '--inputs',
        nargs='+',
        required=True,
        metavar='DIR',
        help="MetaQuest output directories to compare"
    )
    parser_compare.add_argument(
        '-m', '--metadata',
        required=True,
        metavar='FILE',
        help="Metadata file (TSV) linking samples to groups"
    )
    parser_compare.add_argument(
        '-o', '--output',
        default='comparison_results',
        metavar='DIR',
        help="Output directory (default: comparison_results)"
    )

    return parser


def main():
    """Main entry point for the MetaQuest CLI."""
    
    parser = create_parser()
    args = parser.parse_args()
    
    # Determine verbosity level
    verbosity = 'debug' if args.debug else 'standard'
    
    # Setup log file path if output directory is specified
    log_file = None
    if hasattr(args, 'output') and args.output:
        log_file = Path(args.output) / 'metaquest.log'
    
    # Initialize global formatter
    formatter = OutputFormatter(verbosity=verbosity, log_file=log_file)
    set_formatter(formatter)
    
    # Show banner (unless --version or --help)
    if not any(arg in sys.argv for arg in ['-v', '--version', '-h', '--help']):
        formatter.banner(__app_name__, __version__, __tagline__)

    start_time = time.time()

    try:
        # ========== COMMAND: CHECK ==========
        if args.command == 'check':
            formatter.section_header("System Check")
            run_system_check(formatter) # Pass formatter here
            formatter.success("System check completed successfully")
            sys.exit(0)

        # ========== VALIDATE FILE PATHS ==========
        if args.command in ['analyze', 'validate']:
            file_paths = []
            
            if args.type == 'fasta':
                file_paths = [args.input_file]
            elif args.type == 'fastq':
                if args.single:
                    file_paths = [args.single]
                elif args.paired:
                    file_paths = args.paired
                elif args.interleaved:
                    file_paths = [args.interleaved]

            # Check file existence
            for f in file_paths:
                if not Path(f).exists():
                    formatter.error(
                        f"Input file not found: {f}",
                        solutions=[
                            "Check the file path for typos",
                            "Verify the file exists in the current directory",
                            "Use absolute path instead of relative path"
                        ]
                    )
                    sys.exit(1)

        # ========== COMMAND: VALIDATE ==========
        if args.command == 'validate':
            is_valid = handle_validation(file_paths, args.type, args, formatter)
            
            if is_valid:
                formatter.success("All files validated successfully")
                formatter.info("Files are ready for analysis")
            else:
                formatter.error("Validation failed for one or more files")
                formatter.info("Please fix the issues before proceeding with analysis")
            
            sys.exit(0 if is_valid else 1)
        
        # ========== COMMAND: ANALYZE ==========
        elif args.command == 'analyze':
            formatter.section_header("Analysis Pipeline")

            # System check
            with formatter.spinner("Verifying system dependencies"):
                run_system_check(formatter) # Pass formatter here
            formatter.success("System check passed")

            # File validation (unless skipped)
            if not args.skip_validation:
                is_valid = handle_validation(file_paths, args.type, args, formatter)
                if not is_valid:
                    formatter.error(
                        "Analysis aborted due to validation failure",
                        solutions=[
                            "Fix validation errors and try again",
                            "Use --skip-validation flag to bypass (not recommended)"
                        ]
                    )
                    sys.exit(1)
                formatter.success("Input validation passed")
            else:
                formatter.warning("Skipping file validation (--skip-validation flag)")

            # Display annotation settings if not skipping annotation
            if not args.skip_annotation:
                formatter.section_header("Annotation Settings")
                formatter.result({
                    'Contig filtering': 'Enabled' if args.filter_contigs else 'Disabled',
                    'Min contig length': f"{args.min_contig_length} bp" if args.filter_contigs else 'N/A',
                    'tbl2asn auto-kill': 'Enabled' if args.kill_tbl2asn else 'Disabled',
                    'tbl2asn timeout': f"{args.tbl2asn_timeout}s" if args.kill_tbl2asn else 'N/A',
                    'Annotation threads': str(args.annotation_threads)
                })

            # Run analysis
            formatter.section_header(f"{args.type.upper()} Analysis")
            formatter.info(f"Input: {', '.join(file_paths)}")
            formatter.info(f"Output: {args.output}")
            
            run_analysis(file_paths, args.type, args.output, args)

            # Show elapsed time
            elapsed = time.time() - start_time
            formatter.success(f"Analysis completed in {formatter._format_time(elapsed)}")
            formatter.info(f"Results: {Path(args.output) / 'analysis_dashboard.html'}")

        # ========== COMMAND: COMPARE ==========
        elif args.command == 'compare':
            formatter.section_header("Comparative Analysis")
            
            formatter.result({
                'Samples': str(len(args.inputs)),
                'Metadata': args.metadata,
                'Output': args.output
            })
            
            run_comparison(args.inputs, args.metadata, args.output)
            
            elapsed = time.time() - start_time
            formatter.success(f"Comparison completed in {formatter._format_time(elapsed)}")
            formatter.info(f"Results: {args.output}")

    except KeyboardInterrupt:
        formatter.warning("Operation interrupted by user")
        sys.exit(130)
    
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