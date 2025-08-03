import argparse
import subprocess
import sys
import os
from pathlib import Path
from .file_validator import FileValidator

# Version information
__version__ = "3.2.1"
__app_name__ = "MetaQuest"

def main():
    # Suppress startup messages for version/help commands
    if '--version' in sys.argv or '--help' in sys.argv or '-h' in sys.argv or '--validate-only' in sys.argv:
        os.environ['METAQUEST_QUIET'] = '1'
    
    parser = argparse.ArgumentParser(
        description="Metagenomics Analysis Pipeline",
        prog="metaquest"
    )
    
    # Add version argument
    parser.add_argument('--version', action='version', 
                       version=f'{__app_name__} v{__version__}')
    
    # Quick check for version flag to avoid unnecessary imports
    if '--version' in sys.argv:
        parser.parse_args()  # This will print version and exit
        return
    
    # Import modules only when needed (after version check)
    from .analysis import run_analysis
    from .utils import check_dependencies, check_database_status
    
    # Add the input file argument for FASTA
    parser.add_argument('input', nargs='?', help="Input FASTA/FASTQ file")
    
    # For FASTA we still take a single input; for FASTQ allow single- or paired-end
    parser.add_argument('-t', '--type', required=True, choices=['fasta', 'fastq'],
                        help="Input file type")
    fq_group = parser.add_mutually_exclusive_group()
    fq_group.add_argument('-r','--reads', nargs=1, metavar=('R1.fastq',),
                          help="Single-end FASTQ file")
    fq_group.add_argument('-1','--reads1', nargs=1, metavar=('R1.fastq',),
                          help="Paired-end FASTQ: R1 file") 
    parser.add_argument('-2','--reads2', nargs=1, metavar=('R2.fastq',),
                            help="Paired-end FASTQ: R2 file (with --reads1)")
    fq_group.add_argument('-i','--interleaved', nargs=1, metavar=('reads.interleaved.fastq',),
                          help="Interleaved paired-end FASTQ in one file")
    parser.add_argument('-o', '--output', default='results', help="Output directory")
    parser.add_argument('--check-only', action='store_true', 
                       help="Only check dependencies and databases")
    
    # Add validation-specific arguments
    parser.add_argument('--validate-only', action='store_true',
                       help="Only validate input file and show statistics (no analysis)")
    parser.add_argument('--min-quality', type=int, default=20,
                       help="Minimum mean quality score for FASTQ files (default: 20)")
    parser.add_argument('--min-sequences', type=int, default=100,
                       help="Minimum number of sequences required (default: 100)")
    parser.add_argument('--skip-validation', action='store_true',
                       help="Skip file validation (not recommended)")
    
    args = parser.parse_args()
    
    # Handle validation-only mode
    if args.validate_only:
        print(f"🧬 {__app_name__} v{__version__} - File Validation Mode")
        
        # Create validator with custom thresholds if provided
        validator = FileValidator()
        if args.min_quality:
            validator.quality_threshold = args.min_quality
        if args.min_sequences:
            validator.min_sequences = args.min_sequences
        
        # Determine input file and type
        input_file = None
        if args.type == 'fastq':
            if args.reads:
                input_file = args.reads[0]
            elif args.interleaved:
                input_file = args.interleaved[0]
            else:
                input_file = args.reads1[0] if args.reads1 else None
        else:  # fasta
            input_file = args.input
        
        if not input_file:
            print("❌ Error: No input file specified")
            exit(1)
        
        # Run validation
        is_valid, stats = validator.validate_and_analyze(input_file, args.type)
        
        if is_valid:
            print("\n✅ File is ready for analysis!")
            print(f"💡 Run without --validate-only to start analysis")
        else:
            print("\n❌ File failed validation")
            print(f"💡 Please check the quality warnings above")
        
        exit(0 if is_valid else 1)
    
    # Normal analysis mode
    try:
        print(f"🧬 {__app_name__} v{__version__} - Metagenomics Analysis Pipeline")
        print("🔍 Initializing analysis environment...")
        check_dependencies()
        check_database_status()
        
        if args.check_only:
            print("\n✅ All systems ready!")
            exit(0)
        
        # Determine input files based on type
        if args.type == 'fastq':
            # verify single or paired
            if args.reads:
                r1 = args.reads[0]
                if not Path(r1).exists(): 
                    raise FileNotFoundError(f"Input file not found: {r1}")
                reads = [r1]
                input_file_for_validation = r1
            elif args.interleaved:
                # mark for de-interleaving
                reads = args.interleaved
                input_file_for_validation = args.interleaved[0]
            else:
                r1 = args.reads1 and args.reads1[0]
                r2 = args.reads2 and args.reads2[0]
                if not (r1 and r2):
                    raise ValueError("Must provide both --reads1 and --reads2 for paired-end")
                for f in (r1,r2):
                    if not Path(f).exists(): 
                        raise FileNotFoundError(f"Input file not found: {f}")
                reads = [r1, r2]
                input_file_for_validation = r1  # Validate R1 for paired-end
        else:  # FASTA type
            if not args.input:
                raise ValueError("Input FASTA file is required for FASTA analysis")
            if not Path(args.input).exists():
                raise FileNotFoundError(f"Input file not found: {args.input}")
            reads = [args.input]
            input_file_for_validation = args.input
        
        # Validate input file unless skipped
        if not args.skip_validation:
            print("\n🔍 Validating input file...")
            validator = FileValidator()
            if args.min_quality:
                validator.quality_threshold = args.min_quality
            if args.min_sequences:
                validator.min_sequences = args.min_sequences
            
            is_valid, file_stats = validator.validate_and_analyze(input_file_for_validation, args.type)
            
            if not is_valid:
                print("\n❌ Analysis aborted due to validation failure.")
                print("💡 Tips:")
                print("   - For FASTQ: Check quality scores and sequence count")
                print("   - For FASTA: Ensure proper format and unique IDs")
                print("   - Use --validate-only flag to check files without analysis")
                print("   - Use --skip-validation to bypass validation (not recommended)")
                exit(1)
            
            print("\n✅ File validation passed! Proceeding with analysis...")
        else:
            print("\n⚠️  Warning: Skipping file validation (--skip-validation flag used)")
            print("   This may lead to unexpected errors during analysis")
        
        # Start analysis
        if args.type == 'fastq':
            print(f"\n🚀 Starting FASTQ analysis of {reads}")
            run_analysis(reads, 'fastq', args.output)
        else:  # FASTA type
            print(f"\n🚀 Starting FASTA analysis of {args.input}")
            run_analysis(reads, 'fasta', args.output)
        
        print(f"\n🎉 Analysis complete! Results saved to {args.output}")
        
        # Print summary of outputs
        output_path = Path(args.output)
        print("\n📊 Generated files:")
        total_size = 0
        file_count = 0
        for file in sorted(output_path.rglob("*")):
            if file.is_file():
                size_mb = file.stat().st_size / (1024*1024)
                total_size += size_mb
                file_count += 1
                # Show only important files
                if size_mb > 0.1 or file.suffix in ['.html', '.txt', '.json']:
                    print(f"  📁 {file.relative_to(output_path)} ({size_mb:.1f} MB)")
        
        print(f"\n📈 Summary: {file_count} files, {total_size:.1f} MB total")
        print(f"🌐 Open {output_path}/analysis_dashboard.html to explore results")
                
    except subprocess.CalledProcessError as e:
        print(f"\n❌ Error running command: {e.cmd}")
        print("Check that all dependencies are installed and databases are available")
        exit(1)
    except FileNotFoundError as e:
        print(f"\n❌ File error: {e}")
        exit(1)
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        exit(1)

if __name__ == "__main__":
    main()