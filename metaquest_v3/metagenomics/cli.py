import argparse
import subprocess
from pathlib import Path
from .analysis import run_analysis
from .utils import check_dependencies, check_database_status

def main():
    parser = argparse.ArgumentParser(description="Metagenomics Analysis Pipeline")
    parser.add_argument('input', help="Input file path (FASTA or FASTQ)")
    parser.add_argument('-t', '--type', required=True, choices=['fasta', 'fastq'],
                       help="Input file type")
    parser.add_argument('-o', '--output', default='results', help="Output directory")
    parser.add_argument('--check-only', action='store_true', 
                       help="Only check dependencies and databases")
    
    args = parser.parse_args()
    
    try:
        print("=== Metagenomics Analysis Pipeline ===")
        print("Checking dependencies and databases...")
        check_dependencies()
        check_database_status()
        
        if args.check_only:
            print("\n✓ All checks passed!")
            exit(0)
        
        if not Path(args.input).exists():
            raise FileNotFoundError(f"Input file not found: {args.input}")
        
        print(f"\nStarting analysis of {args.input}")
        run_analysis(args.input, args.type, args.output)
        
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
