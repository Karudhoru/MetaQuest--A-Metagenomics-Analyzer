import argparse
import subprocess
import sys
import os
from pathlib import Path

# Version information
__version__ = "3.2.0"
__app_name__ = "MetaQuest"

def main():
    # Suppress startup messages for version/help commands
    if '--version' in sys.argv or '--help' in sys.argv or '-h' in sys.argv:
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
    parser.add_argument('input', nargs='?', help="Input FASTA file")
    
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
    
    args = parser.parse_args()
    
    try:
        print(f"🧬 {__app_name__} v{__version__} - Metagenomics Analysis Pipeline")
        print("🔍 Initializing analysis environment...")
        check_dependencies()
        check_database_status()
        
        if args.check_only:
            print("\n✅ All systems ready!")
            exit(0)
        
        if args.type == 'fastq':
            # verify single or paired
            if args.reads:
                r1 = args.reads[0]
                if not Path(r1).exists(): raise FileNotFoundError(r1)
                reads = [r1]
            elif args.interleaved:
                # mark for de-interleaving
                reads = args.interleaved
            else:
                r1 = args.reads1 and args.reads1[0]
                r2 = args.reads2 and args.reads2[0]
                if not (r1 and r2):
                    raise ValueError("Must provide both --reads1 and --reads2 for paired-end")
                for f in (r1,r2):
                    if not Path(f).exists(): raise FileNotFoundError(f)
                reads = [r1, r2]
            print(f"\n🚀 Starting FASTQ analysis of {reads}")
            run_analysis(reads, 'fastq', args.output)
        else:  # FASTA type
            if not args.input:
                raise ValueError("Input FASTA file is required for FASTA analysis")
            if not Path(args.input).exists():
                raise FileNotFoundError(f"Input file not found: {args.input}")
            print(f"\n🚀 Starting FASTA analysis of {args.input}")
            run_analysis([args.input], 'fasta', args.output)
        
        print(f"\n🎉 Analysis complete! Results saved to {args.output}")
        
        # Print summary of outputs
        output_path = Path(args.output)
        print("\nGenerated files:")
        for file in sorted(output_path.rglob("*")):
            if file.is_file():
                size_mb = file.stat().st_size / (1024*1024)
                print(f"  📁 {file.relative_to(output_path)} ({size_mb:.1f} MB)")
                
    except subprocess.CalledProcessError as e:
        print(f"\n❌ Error running command: {e.cmd}")
        print("Check that all dependencies are installed and databases are available")
    except FileNotFoundError as e:
        print(f"\n❌ File error: {e}")
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()