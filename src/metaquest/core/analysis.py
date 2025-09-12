import subprocess
import datetime
import os
import pandas as pd
from pathlib import Path
from Bio import SeqIO
import numpy as np
import json
from .taxonomic_analysis import run_kraken, run_bracken, run_fasta_blast_taxonomy
from .pathogen_analysis import run_pathogen_scan
from .functional_analysis import run_prokka, run_swissprot_annotation
from ..config import *
from ..io.file_validator import FileValidator
from ..io.utils import convert_fastq_to_fasta, split_interleaved
from ..visualization.dashboard import create_dashboard
from ..reporting.main_reporter import (
    PathogenReporter,
    generate_taxonomic_report, 
    generate_functional_report,
    generate_fastq_pathogen_report,
    generate_fasta_ml_pathogen_report,
    generate_sequence_hits_report
)
from ..visualization.visualization import (
    create_taxonomic_visualizations,
    create_pathogen_visualizations, 
    create_functional_visualizations
)

# ML pathogen prediction integration
ML_PATHOGEN_AVAILABLE = False
try:
    from ..ml.pathogen_predictor import run_ml_pathogen_prediction
    ML_PATHOGEN_AVAILABLE = True
except ImportError as e:
    ML_PATHOGEN_AVAILABLE = False
    print(f"⚠️ PathogenPredictor not available: {e}")
    print(" Analysis will continue with traditional pathogen detection methods")

def should_use_ml_prediction(prokka_dir):
    """Determine if ML prediction is appropriate based on sequence lengths"""
    try:
        prokka_path = Path(prokka_dir)
        protein_files = list(prokka_path.glob("*.faa"))
        
        if not protein_files:
            return False
        
        from Bio import SeqIO
        sequences = list(SeqIO.parse(protein_files[0], "fasta"))
        
        if not sequences:
            return False
        
        lengths = [len(seq.seq) for seq in sequences]
        avg_length = np.mean(lengths)
        min_length_threshold = 200  
        
        print(f"📏 Protein length analysis:")
        print(f"   • Average protein length: {avg_length:.1f} amino acids")
        print(f"   • Total proteins: {len(sequences)}")
        
        if avg_length >= min_length_threshold:
            print(f"   ✅ Suitable for ML prediction (avg ≥ {min_length_threshold} aa)")
            return True
        else:
            print(f"   ⚠️ Too short for ML prediction (avg < {min_length_threshold} aa)")
            print(f"   📝 ML model trained on 500-1000 aa proteins")
            return False
            
    except Exception as e:
        print(f"⚠️ Could not analyze protein lengths: {e}")
        return False

def run_analysis(input_file, file_type, output_dir, cli_args=None):
    """Main analysis controller that accepts CLI arguments and calls the correct pipeline."""
    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)

    if file_type == 'fastq':
        analyze_fastq(input_file, output_dir_path, cli_args)
    else:
        # For FASTA, ensure we pass the actual file path, not a list
        fasta_path = input_file[0] if isinstance(input_file, list) else input_file
        analyze_fasta(Path(fasta_path), output_dir_path, cli_args)
    
    # --- This is the correct, final call for the dashboard ---
    print("\nGenerating final analysis dashboard...")
    create_dashboard(analysis_type=file_type, output_dir=output_dir_path)
    
    print(f"\n🎉 Analysis complete! Open {output_dir_path / 'analysis_dashboard.html'} to explore results.")

def analyze_fastq(reads, output_dir: Path, args):
    """Process FASTQ files with streamlined pathogen reporting and optional annotation skip"""
    print("\n=== FASTQ Analysis Pipeline ===")
    
    fasta_path = None
    kraken_report = None
    bracken_report = None
    
    try:
        # Handle interleaved files
        if hasattr(args, 'interleaved') and args.interleaved:
            print("Splitting interleaved FASTQ...")
            reads = split_interleaved(reads[0], output_dir)
        elif hasattr(args, 'single') and args.single:
            print("Processing single-end FASTQ...")
            # Ensure reads is a list for consistency
            reads = [reads] if isinstance(reads, str) else reads
        
        # --- 1. Taxonomic classification (Always runs) ---
        print("\n--- Running Taxonomic Analysis ---")
        print("1. Running Kraken2 classification...")
        kraken_report_path = run_kraken(reads, output_dir)
        
        print("2. Running Bracken abundance estimation...")
        bracken_report = run_bracken(kraken_report_path, output_dir)
        
        # Generate taxonomic reports and visualizations
        print("3. Generating taxonomic reports and visualizations...")
        generate_taxonomic_report(output_dir, bracken_data=bracken_report)
        create_taxonomic_visualizations(output_dir, bracken_report)
        print("✓ Taxonomic analysis completed")
        
    except Exception as e:
        print(f"⚠️ Taxonomic analysis failed: {str(e)}")
    
    # --- 2. Convert to FASTA (needed for annotation steps) ---
    if args and not args.skip_annotation:
        try:
            print("4. Converting to FASTA...")
            fasta_path = convert_fastq_to_fasta(reads if len(reads) > 1 else reads[0], output_dir)
        except Exception as e:
            print(f"⚠️ FASTA conversion failed: {str(e)}")
            print(" Skipping analyses that require FASTA format.")
            fasta_path = None
    else:
        print("4. ⏭️ Skipping FASTA conversion (--skip-annotation flag detected)")
    
    # --- 3. Conditionally run annotation steps ---
    if args and not args.skip_annotation:
        if fasta_path:
            # Enhanced pathogen screening
            pathogen_results = None
            try:
                print("\n--- Running Pathogen & Functional Analysis (this may take a while) ---")
                print("5. Running comprehensive pathogen screening...")
                pathogen_results = run_pathogen_scan(
                    fasta_path, output_dir, 
                    bracken_results=bracken_report, 
                    taxonomy_results=None
                )
                
                print("6. Generating Comprehensive Pathogen Hits Report...")
                generate_sequence_hits_report(output_dir, 'fastq', pathogen_results, "Pathogen Screening Hits")
                print("  -> Pathogen screening complete.")
     
                
            except Exception as e:
                print(f"⚠️ Pathogen analysis failed: {str(e)}")
            
            # Functional annotation with ML pathogen prediction
            prokka_dir = None
            ml_results = None
            ml_summary = None
            
            try:
                print("7. Running gene prediction...")
                prokka_dir = run_prokka(fasta_path, output_dir)
                protein_file = prokka_dir / "sample.faa"
                gff_file = prokka_dir / "sample.gff"
                
                protein_files = list(Path(prokka_dir).glob("*.faa"))
                if not protein_files:
                    print("⚠️ Warning: No protein FASTA files found from Prokka. Skipping functional annotation.")
                elif all(os.path.getsize(pf) == 0 for pf in protein_files):
                    print("⚠️ Warning: All protein files are empty. Skipping functional annotation.")
                else:
                    print("8. Running functional annotation...")
                    swissprot_results = run_swissprot_annotation(prokka_dir, output_dir)
                    
                    print("9. Generating functional analysis reports...")
                    generate_functional_report(output_dir, prokka_results=prokka_dir, swissprot_results=swissprot_results)
                    create_functional_visualizations(output_dir, prokka_results=prokka_dir, swissprot_results=swissprot_results)
                    print("  -> Functional annotation complete.")
                    
                    # ML-based pathogen prediction
                    if ML_PATHOGEN_AVAILABLE:
                        if should_use_ml_prediction(prokka_dir):
                            try:
                                print("10. Running ML-based pathogen prediction...")
                                ml_results, ml_summary = run_ml_pathogen_prediction(prokka_dir, output_dir)
                                
                                if ml_results and ml_summary:
                                    print(f"✓ ML pathogen prediction completed:")
                                    print(f" 📊 {ml_summary['pathogenic_predictions']}/{ml_summary['total_sequences_analyzed']} proteins predicted as pathogenic")
                                    print(f" ⭐ {ml_summary['high_confidence_predictions']} high-confidence predictions")
                                    print(f" 🎯 Average confidence: {ml_summary['average_confidence']:.3f}")
                                else:
                                    print("⚠️ ML pathogen prediction produced no results")
                                    
                            except Exception as e:
                                print(f"⚠️ ML pathogen prediction failed: {str(e)}")
                        else:
                            print("10. ⏭️ Skipping ML prediction - protein sequences too short")
                    else:
                        print("10. ℹ️ ML pathogen predictor not available - skipping ML analysis")
                        
            except Exception as e:
                print(f"⚠️ Functional analysis failed: {str(e)}")
            
            # STREAMLINED pathogen reporting for FASTQ
            try:
                print("11. Generating comprehensive pathogen analysis for FASTQ...")
                
                # Extract pathogens directly from Bracken data
                bracken_pathogens = []
                if bracken_report and Path(bracken_report).exists():
                    reporter = PathogenReporter(Path(output_dir), analysis_mode='fastq')
                    bracken_pathogens = reporter._extract_pathogens_from_bracken(bracken_report)
                
                # Generate the comprehensive FASTQ pathogen report
                generate_fastq_pathogen_report(
                    output_dir, 
                    bracken_pathogens=bracken_pathogens,
                    taxonomy_pathogens=[], 
                    sequence_pathogens=[]
                )

                # Generate pathogen visualizations
                report_file = output_dir / "pathogen_detection_report.json"
                if report_file.exists():
                    create_pathogen_visualizations(output_dir, traditional_data=str(report_file))
                
                print("✓ Comprehensive FASTQ pathogen analysis completed")
                        
            except Exception as e:
                print(f"⚠️ Pathogen report generation failed: {str(e)}")
                
    else:
        print("\n--skip-annotation flag detected. Skipping pathogen and functional analysis.")
        print("  -> Only taxonomic classification will be performed for faster analysis.")
    
    # --- 4. Generate visualizations (Always runs) ---
    print("\n--- Generating Visualizations ---")
    if bracken_report:
        create_taxonomic_visualizations(output_dir, bracken_report)
        print("✓ Taxonomic visualizations completed")
    
    print(f"\n📁 All results saved to: {output_dir}")
    print("🔍 Key files to review:")
    print(" • taxonomic_classification_report.txt - Species identification and abundance")
    if args and not args.skip_annotation:
        print(" • pathogen_summary.txt - Comprehensive pathogen analysis with clinical recommendations")
        print(" • functional_annotation_report.txt - Gene prediction and annotation")
    print(" • analysis_dashboard.html - Interactive dashboard")

def analyze_fasta(fasta_path: Path, output_dir: Path, args):
    """Process FASTA files with BLAST taxonomy + ML and optional annotation skip"""
    print("\n=== FASTA Analysis Pipeline ===")
    
    blast_results_data = None
    
    try:
        # --- 1. Taxonomic classification (Always runs) ---
        print("\n--- Running Taxonomic Analysis ---")
        print("1. Running BLAST taxonomic classification...")
        
        # Use blast_sample_size if provided in args, otherwise default to 50
        max_sequences = getattr(args, 'blast_sample_size', 50) if args else 50
        
        blast_results_file = run_fasta_blast_taxonomy(fasta_path, output_dir, database="nt", max_sequences=max_sequences)
        
        if blast_results_file and Path(blast_results_file).exists():
            with open(blast_results_file, 'r') as f:
                blast_results_data = json.load(f)
            print(f"✓ Loaded {len(blast_results_data)} BLAST results")
            
            print("2. Generating comprehensive taxonomic reports...")
            generate_taxonomic_report(output_dir, blast_data=blast_results_data)
            
            print("3. Creating BLAST-based taxonomic visualizations...")
            create_taxonomic_visualizations(output_dir, blast_results_data)
            print("✓ Comprehensive BLAST taxonomic analysis completed")
        else:
            print("⚠️ BLAST analysis produced no results file.")

    except Exception as e:
        print(f"⚠️ BLAST taxonomic analysis failed: {str(e)}")
    
    # --- 2. Conditionally run annotation steps ---
    if args and not args.skip_annotation:
        print("\n--- Running Functional & ML Analysis (this may take a while) ---")
        print("4. ⏭️ Skipping traditional pathogen screening for FASTA analysis (using BLAST+ML instead).")
        
        # Functional annotation with ML pathogen prediction
        ml_results = None
        ml_summary = None
        
        try:
            print("5. Running gene prediction...")
            prokka_dir = run_prokka(fasta_path, output_dir)
            
            protein_files = list(Path(prokka_dir).glob("*.faa"))
            if not protein_files or all(os.path.getsize(pf) == 0 for pf in protein_files):
                print("⚠️ No valid proteins predicted. Skipping functional and ML analysis.")
            else:
                print("6. Running functional annotation...")
                swissprot_results = run_swissprot_annotation(prokka_dir, output_dir)
                
                print("7. Generating functional analysis reports...")
                generate_functional_report(output_dir, prokka_results=prokka_dir, swissprot_results=swissprot_results)
                create_functional_visualizations(output_dir, prokka_results=prokka_dir, swissprot_results=swissprot_results)
                print("  -> Functional annotation complete.")
                
                # ML pathogen prediction
                if ML_PATHOGEN_AVAILABLE and should_use_ml_prediction(prokka_dir):
                    print("8. Running ML-based pathogen prediction...")
                    ml_results, ml_summary = run_ml_pathogen_prediction(prokka_dir, output_dir)
                    
                    if ml_results and ml_summary:
                        print(f"✓ ML pathogen prediction completed:")
                        print(f" 📊 {ml_summary['pathogenic_predictions']}/{ml_summary['total_sequences_analyzed']} proteins predicted as pathogenic")
                        print(f" ⭐ {ml_summary['high_confidence_predictions']} high-confidence predictions")
                        print(f" 🎯 Average confidence: {ml_summary['average_confidence']:.3f}")
                else:
                    print("8. ⏭️ Skipping ML prediction.")
                    
        except Exception as e:
            print(f"⚠️ Functional/ML analysis failed: {str(e)}")
        
        # INTEGRATED BLAST+ML pathogen reporting
        try:
            print("9. Generating integrated BLAST+ML pathogen analysis...")
            
            # Extract pathogenic species from BLAST results
            blast_taxonomy_pathogens = []
            if blast_results_data:
                reporter = PathogenReporter(Path(output_dir), analysis_mode='fasta')
                blast_taxonomy_pathogens = reporter._extract_pathogens_from_blast_taxonomy(blast_results_data)

            # Generate the integrated FASTA+ML pathogen report if we have any evidence
            if blast_taxonomy_pathogens or (ml_results and ml_summary):
                generate_fasta_ml_pathogen_report(
                    output_dir, 
                    blast_taxonomy_pathogens=blast_taxonomy_pathogens,
                    ml_results=ml_results,
                    ml_summary=ml_summary
                )
                
                # Generate visualizations for the integrated report
                report_file = output_dir / "blast_ml_integrated_pathogen_report.json"
                if report_file.exists():
                    create_pathogen_visualizations(output_dir, traditional_data=str(report_file))
                    
                print("✓ Integrated BLAST+ML pathogen analysis completed")
            else:
                print("⚠️ No pathogenic evidence found from either BLAST or ML analysis. No report generated.")
                    
        except Exception as e:
            print(f"⚠️ Integrated pathogen report generation failed: {str(e)}")
            
    else:
        print("\n--skip-annotation flag detected. Skipping functional and ML pathogen analysis.")
        print("  -> Only BLAST taxonomic classification will be performed for faster analysis.")
    
    # --- 3. Generate visualizations ---
    print("\n--- Generating Visualizations ---")
    if blast_results_data:
        create_taxonomic_visualizations(output_dir, blast_results_data)
        print("✓ BLAST taxonomic visualizations completed")
    
    print(f"\n📁 All results saved to: {output_dir}")
    print("🔍 Key files to review:")
    print(" • taxonomic_classification_report.txt - Comprehensive BLAST-based species identification")
    if args and not args.skip_annotation:
        print(" • blast_ml_pathogen_summary.txt - Integrated BLAST+ML pathogen analysis")
        print(" • functional_annotation_report.txt - Gene prediction and annotation quality")
    print(" • analysis_dashboard.html - Interactive dashboard")