import argparse
import sys
import os
from pathlib import Path
from .core.analysis import run_analysis
from .core.comparative_analysis import run_comparison 
from .io.file_validator import FileValidator
from .io.utils import run_system_check

# Version information
__version__ = "3.3.0" 
__app_name__ = "MetaQuest"

def _display_header():
    """Displays the MetaQuest tagline header."""
    tagline = "A Comprehensive Metagenomics Analysis Pipeline"
    print(f"\n🧬 {__app_name__} v{__version__} | {tagline}\n")

def setup_fastq_validation_args(parser):
    """Helper to add FASTQ-specific validation arguments to a parser."""
    val_group = parser.add_argument_group('FASTQ Validation Options')
    val_group.add_argument('-q', '--min-quality', type=int, default=20, help="Minimum mean quality score (default: 20).")
    val_group.add_argument('-n', '--min-sequences', type=int, default=100, help="Minimum number of sequences (default: 100).")
    val_group.add_argument('--overrep-threshold', type=float, default=0.1, help="Percentage threshold to flag a sequence as overrepresented (default: 0.1).")

def handle_validation(file_paths, file_type, args):
    """Helper function to run the file validator and return validation status."""
    print(f"🧬 {__app_name__} - File Validation")
    validator = FileValidator()
    
    if file_type == 'fastq':
        validator.quality_threshold = args.min_quality
        validator.min_sequences = args.min_sequences
        validator.overrep_threshold = args.overrep_threshold

    all_valid = True
    for i, file_path in enumerate(file_paths):
        print(f"\n--- Validating file {i+1}/{len(file_paths)}: {file_path} ---")
        is_valid, _ = validator.validate_and_analyze(file_path, file_type)
        if not is_valid:
            all_valid = False
    
    return all_valid

def main():
    """Main entry point for the MetaQuest CLI."""
    if not any(arg in sys.argv for arg in ['-v', '--version', '-h', '--help']):
        _display_header()

    parser = argparse.ArgumentParser(
        description=f"{__app_name__} v{__version__} - A Metagenomics Analysis Pipeline",
        prog="metaquest"
    )
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s {__version__}')
    
    subparsers = parser.add_subparsers(dest='command', required=True, help="Available commands")

    # --- Command: analyze ---
    parser_analyze = subparsers.add_parser('analyze', help="Run the full analysis pipeline on a single sample.")
    analysis_parent_parser = argparse.ArgumentParser(add_help=False)
    analysis_parent_parser.add_argument('-o', '--output', default='results', help="Output directory name (default: results).")
    analysis_parent_parser.add_argument('--skip-validation', action='store_true', help="Skip input file validation.")
    analysis_parent_parser.add_argument('--skip-annotation', action='store_true', help="Skip functional and pathogen annotation steps for a faster taxonomic-only analysis.")
    
    analysis_subparsers = parser_analyze.add_subparsers(dest='type', required=True, help="Input data type")
    analyze_fasta = analysis_subparsers.add_parser('fasta', help='Analyze a single FASTA file.', parents=[analysis_parent_parser])
    analyze_fasta.add_argument('input_file', help="Path to the input FASTA file.")
    analyze_fasta.add_argument('-s', '--blast-sample-size', type=int, default=50, help="Number of sequences to BLAST (default: 50).")
    
    analyze_fastq = analysis_subparsers.add_parser('fastq', help='Analyze FASTQ files.', parents=[analysis_parent_parser])
    setup_fastq_validation_args(analyze_fastq)
    fastq_mode = analyze_fastq.add_mutually_exclusive_group(required=True)
    fastq_mode.add_argument('--single', metavar='READS.fastq', help="Single-end FASTQ file.")
    fastq_mode.add_argument('--paired', nargs=2, metavar=('R1.fastq', 'R2.fastq'), help="Paired-end FASTQ files.")
    fastq_mode.add_argument('--interleaved', metavar='INTERLEAVED.fastq', help="Interleaved paired-end FASTQ file.")

    # --- Command: validate ---
    parser_validate = subparsers.add_parser('validate', help="Validate input file(s) without running analysis.")
    validate_subparsers = parser_validate.add_subparsers(dest='type', required=True, help="Input data type")

    validate_fasta = validate_subparsers.add_parser('fasta', help='Validate a single FASTA file.')
    validate_fasta.add_argument('input_file', help="Path to the input FASTA file.")

    validate_fastq = validate_subparsers.add_parser('fastq', help='Validate FASTQ files.')
    validate_mode = validate_fastq.add_mutually_exclusive_group(required=True)
    validate_mode.add_argument('--single', metavar='READS.fastq', help="Single-end FASTQ file.")
    validate_mode.add_argument('--paired', nargs=2, metavar=('R1.fastq', 'R2.fastq'), help="Paired-end FASTQ files (R1 and R2).")
    validate_mode.add_argument('--interleaved', metavar='INTERLEAVED.fastq', help="Interleaved paired-end FASTQ file.")
    setup_fastq_validation_args(validate_fastq)

    # --- Command: check ---
    parser_check = subparsers.add_parser('check', help="Check all dependencies and database status.")
    
    # --- NEW: Command: compare ---
    parser_compare = subparsers.add_parser('compare', help="Perform comparative analysis across multiple samples.")
    parser_compare.add_argument('-i', '--inputs', nargs='+', required=True, help="Space-separated list of MetaQuest output directories to compare.")
    parser_compare.add_argument('-m', '--metadata', required=True, help="Path to the metadata file (TSV) linking samples to groups.")
    parser_compare.add_argument('-o', '--output', default='comparison_results', help="Directory to save comparison results (default: comparison_results).")


    args = parser.parse_args()

    # --- Command Execution Logic ---
    try:
        if args.command == 'check':
            run_system_check()
            sys.exit(0)

        # Logic for 'analyze' and 'validate' commands (now handles file paths)
        if args.command in ['analyze', 'validate']:
            file_paths = []
            if args.type == 'fasta':
                file_paths = [args.input_file]
            elif args.type == 'fastq':
                if args.single: file_paths = [args.single]
                elif args.paired: file_paths = args.paired
                elif args.interleaved: file_paths = [args.interleaved]

            for f in file_paths:
                if not Path(f).exists():
                    parser.error(f"Input file not found: {f}")

        if args.command == 'validate':
            is_valid = handle_validation(file_paths, args.type, args)
            if is_valid:
                print("\n✅ File validation successful. Ready for analysis.")
            else:
                print("\n❌ File validation failed.")
            sys.exit(0 if is_valid else 1)
        
        elif args.command == 'analyze':
            print("🔍 Initializing analysis environment...")
            run_system_check()

            if not args.skip_validation:
                is_valid = handle_validation(file_paths, args.type, args)
                if not is_valid:
                    print("\n❌ Analysis aborted due to validation failure.")
                    sys.exit(1)
                print("\n✅ File validation passed! Proceeding with analysis...")
            else:
                print("\n⚠️ Warning: Skipping file validation (--skip-validation flag used).")

            print(f"\n🚀 Starting {args.type.upper()} analysis on: {', '.join(file_paths)}")
            run_analysis(file_paths, args.type, args.output)

            print(f"\n🎉 Analysis complete! Results saved to '{args.output}'")

        elif args.command == 'compare':
            print("🔬 Initializing Comparative Analysis...")
            run_comparison(args.inputs, args.metadata, args.output)

    except SystemExit as e:
        if e.code != 0:
            print(f"\nHalting due to system check failure.")
        sys.exit(e.code)
    except Exception as e:
        print(f"\n❌ An unexpected error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()