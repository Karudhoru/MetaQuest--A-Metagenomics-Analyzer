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
from .functional_analysis import run_prokka, run_functional_annotation
from ..config import *
from ..io.file_validator import FileValidator
from ..io.utils import assemble_reads_to_fasta, split_interleaved, should_use_ml_prediction, extract_pathogens_from_bracken, extract_pathogens_from_blast_taxonomy
from ..visualization.dashboard import create_dashboard
from ..reporting.main_reporter import (
    generate_taxonomic_report, 
    generate_functional_report,
    generate_fastq_pathogen_report,
    generate_fasta_ml_pathogen_report,
    generate_sequence_hits_report
)
from ..visualization.main_visualizer import (
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
    print("   Analysis will continue with traditional pathogen detection methods")

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
    
    # Generate final dashboard
    print("\n" + "="*80)
    print("GENERATING FINAL ANALYSIS DASHBOARD")
    print("="*80)
    create_dashboard(analysis_type=file_type, output_dir=output_dir_path)
    
    print(f"\n🎉 Analysis complete! Open {output_dir_path / 'analysis_dashboard.html'} to explore results.")


def analyze_fastq(reads, output_dir: Path, args):
    """Process FASTQ files with comprehensive v4.0.0 reporting and visualization"""
    print("\n" + "="*80)
    print("FASTQ ANALYSIS PIPELINE v4.0.0")
    print("="*80)
    
    fasta_path = None
    kraken_report = None
    bracken_report = None
    
    # ========================================================================
    # STEP 1: TAXONOMIC CLASSIFICATION (Always runs)
    # ========================================================================
    try:
        print("\n" + "="*80)
        print("STEP 1: TAXONOMIC CLASSIFICATION")
        print("="*80)
        
        # Handle interleaved files
        if hasattr(args, 'interleaved') and args.interleaved:
            print("→ Splitting interleaved FASTQ...")
            reads = split_interleaved(reads[0], output_dir)
        elif hasattr(args, 'single') and args.single:
            print("→ Processing single-end FASTQ...")
            reads = [reads] if isinstance(reads, str) else reads
        
        print("\n1. Running Kraken2 classification...")
        kraken_report_path = run_kraken(reads, output_dir)
        
        print("2. Running Bracken abundance estimation...")
        bracken_report = run_bracken(kraken_report_path, output_dir)
        
        print("3. Generating taxonomic reports...")
        generate_taxonomic_report(str(output_dir), bracken_data=bracken_report)
        
        print("4. Creating taxonomic visualizations...")
        create_taxonomic_visualizations(output_dir, bracken_report)
        
        print("\n✓ Taxonomic classification completed successfully")
        
    except Exception as e:
        print(f"\n✗ Taxonomic analysis failed: {str(e)}")
        import traceback
        traceback.print_exc()
    
    # ========================================================================
    # STEP 2: FASTA CONVERSION (Conditional)
    # ========================================================================
    if args and not args.skip_annotation:
        try:
            print("\n" + "="*80)
            print("STEP 2: METAGENOMIC ASSEMBLY")
            print("="*80)
            print("→ Converting reads to FASTA format...")
            fasta_path = assemble_reads_to_fasta(reads, output_dir)
            print(f"✓ FASTA conversion completed: {fasta_path}")
        except Exception as e:
            print(f"✗ FASTA conversion failed: {str(e)}")
            fasta_path = None
    else:
        print("\n⏭️ Skipping FASTA conversion (--skip-annotation flag)")
    
    # ========================================================================
    # STEP 3: PATHOGEN & FUNCTIONAL ANALYSIS (Conditional)
    # ========================================================================
    if args and not args.skip_annotation and fasta_path:
        
        # --- Pathogen Screening ---
        pathogen_results = None
        try:
            print("\n" + "="*80)
            print("STEP 3: PATHOGEN SCREENING")
            print("="*80)
            print("→ Running comprehensive pathogen database search...")
            pathogen_results = run_pathogen_scan(
                fasta_path, output_dir, 
                bracken_results=bracken_report, 
                taxonomy_results=None
            )
            
            if pathogen_results and Path(pathogen_results).exists():
                print("→ Generating pathogen hits summary report...")
                generate_sequence_hits_report(
                    str(output_dir), 
                    'Pathogen Database Search', 
                    pathogen_results, 
                    "Pathogen Screening Results"
                )
                print("✓ Pathogen screening completed")
            
        except Exception as e:
            print(f"✗ Pathogen screening failed: {str(e)}")
        
        # --- Functional Annotation ---
        prokka_dir = None
        ml_results = None
        ml_summary = None
        
        try:
            print("\n" + "="*80)
            print("STEP 4: FUNCTIONAL ANNOTATION")
            print("="*80)
            
            print("1. Running gene prediction (Prokka)...")
            prokka_dir = run_prokka(fasta_path, output_dir)
            
            protein_files = list(Path(prokka_dir).glob("*.faa"))
            if not protein_files:
                print("⚠️ No protein files found - skipping functional annotation")
            elif all(os.path.getsize(pf) == 0 for pf in protein_files):
                print("⚠️ All protein files are empty - skipping functional annotation")
            else:
                print("2. Running SwissProt annotation...")
                swissprot_results = run_functional_annotation(prokka_dir, output_dir)
                
                print("3. Generating functional reports...")
                generate_functional_report(
                    str(output_dir), 
                    prokka_results=str(prokka_dir), 
                    swissprot_results=swissprot_results
                )
                
                print("4. Creating functional visualizations...")
                create_functional_visualizations(
                    output_dir, 
                    prokka_results=prokka_dir, 
                    swissprot_results=swissprot_results
                )
                
                print("✓ Functional annotation completed")
                
                # --- ML Pathogen Prediction ---
                if ML_PATHOGEN_AVAILABLE:
                    if should_use_ml_prediction(prokka_dir):
                        try:
                            print("\n5. Running ML-based pathogen prediction...")
                            ml_results, ml_summary = run_ml_pathogen_prediction(prokka_dir, output_dir)
                            
                            if ml_results and ml_summary:
                                print(f"✓ ML pathogen prediction completed:")
                                print(f"   📊 {ml_summary['pathogenic_predictions']}/{ml_summary['total_sequences_analyzed']} proteins predicted as pathogenic")
                                print(f"   ⭐ {ml_summary['high_confidence_predictions']} high-confidence predictions")
                                print(f"   🎯 Average confidence: {ml_summary['average_confidence']:.3f}")
                            else:
                                print("⚠️ ML prediction produced no results")
                                
                        except Exception as e:
                            print(f"⚠️ ML prediction failed: {str(e)}")
                    else:
                        print("\n5. ⏭️ Skipping ML prediction - proteins too short")
                else:
                    print("\n5. ℹ️ ML predictor not available - skipping")
                    
        except Exception as e:
            print(f"\n✗ Functional analysis failed: {str(e)}")
            import traceback
            traceback.print_exc()
        
        # --- Comprehensive Pathogen Report ---
        try:
            print("\n" + "="*80)
            print("STEP 5: COMPREHENSIVE PATHOGEN ANALYSIS")
            print("="*80)
            
            print("→ Extracting pathogenic organisms from Bracken...")
            bracken_pathogens = []
            if bracken_report and Path(bracken_report).exists():
                bracken_pathogens = extract_pathogens_from_bracken(bracken_report)
                print(f"   ✓ Extracted {len(bracken_pathogens)} pathogenic organisms")
            
            print("→ Generating comprehensive FASTQ pathogen report...")
            pathogen_report_result = generate_fastq_pathogen_report(
                str(output_dir), 
                bracken_pathogens=bracken_pathogens,
                taxonomy_pathogens=[], 
                sequence_pathogens=[]
            )

            print("→ Creating pathogen visualizations...")
            report_file = output_dir / "pathogen_detection_report.json"
            if report_file.exists():
                create_pathogen_visualizations(
                    output_dir, 
                    traditional_data=str(report_file)
                )
            
            print("\n✓ Comprehensive pathogen analysis completed")
            print(f"   Risk Level: {pathogen_report_result.get('risk_level', 'UNKNOWN')}")
            print(f"   Pathogens Detected: {pathogen_report_result.get('pathogen_count', 0)}")
            print(f"   Critical Pathogens: {pathogen_report_result.get('critical_pathogens', 0)}")
                    
        except Exception as e:
            print(f"\n✗ Pathogen report generation failed: {str(e)}")
            import traceback
            traceback.print_exc()
            
    else:
        if args and args.skip_annotation:
            print("\n" + "="*80)
            print("⏭️ SKIPPING ANNOTATION STEPS")
            print("="*80)
            print("--skip-annotation flag detected")
            print("Only taxonomic classification will be performed")
    
    # ========================================================================
    # FINAL SUMMARY
    # ========================================================================
    print("\n" + "="*80)
    print("FASTQ ANALYSIS SUMMARY")
    print("="*80)
    print(f"📁 Output Directory: {output_dir}")
    print("\n📋 Key Files Generated:")
    print("   • taxonomic_classification_report.txt - Species identification")
    print("   • taxonomic_abundance_chart.html - Interactive abundance visualization")
    
    if args and not args.skip_annotation:
        print("   • pathogen_detection_report.txt - Comprehensive pathogen analysis")
        print("   • pathogen_risk_assessment.html - Risk visualization")
        print("   • functional_annotation_report.txt - Gene prediction results")
        print("   • annotation_quality_dashboard.html - Quality metrics")
    
    print("   • analysis_dashboard.html - Complete interactive dashboard")
    print("\n✓ FASTQ pipeline completed successfully")


def analyze_fasta(fasta_path: Path, output_dir: Path, args):
    """Process FASTA files with comprehensive v4.0.0 reporting and visualization"""
    print("\n" + "="*80)
    print("FASTA ANALYSIS PIPELINE v4.0.0")
    print("="*80)
    
    blast_results_data = None
    
    # ========================================================================
    # STEP 1: BLAST TAXONOMIC CLASSIFICATION (Always runs)
    # ========================================================================
    try:
        print("\n" + "="*80)
        print("STEP 1: BLAST TAXONOMIC CLASSIFICATION")
        print("="*80)
        
        max_sequences = getattr(args, 'blast_sample_size', 50) if args else 50
        print(f"→ Running BLAST taxonomy (analyzing {max_sequences} sequences)...")
        
        blast_results_file = run_fasta_blast_taxonomy(
            fasta_path, output_dir, 
            database="nt", 
            max_sequences=max_sequences
        )
        
        if blast_results_file and Path(blast_results_file).exists():
            with open(blast_results_file, 'r') as f:
                blast_results_data = json.load(f)
            print(f"✓ Loaded {len(blast_results_data)} BLAST results")
            
            print("→ Generating comprehensive taxonomic reports...")
            generate_taxonomic_report(str(output_dir), blast_data=blast_results_data)
            
            print("→ Creating BLAST taxonomic visualizations...")
            create_taxonomic_visualizations(output_dir, blast_results_data)
            
            print("\n✓ BLAST taxonomic analysis completed successfully")
        else:
            print("⚠️ BLAST analysis produced no results")

    except Exception as e:
        print(f"\n✗ BLAST taxonomic analysis failed: {str(e)}")
        import traceback
        traceback.print_exc()
    
    # ========================================================================
    # STEP 2: FUNCTIONAL & ML ANALYSIS (Conditional)
    # ========================================================================
    if args and not args.skip_annotation:
        
        print("\n" + "="*80)
        print("STEP 2: FUNCTIONAL & ML ANALYSIS")
        print("="*80)
        print("ℹ️ Skipping traditional pathogen screening (using BLAST+ML approach)")
        
        ml_results = None
        ml_summary = None
        
        try:
            print("\n1. Running gene prediction (Prokka)...")
            prokka_dir = run_prokka(fasta_path, output_dir)
            
            protein_files = list(Path(prokka_dir).glob("*.faa"))
            if not protein_files or all(os.path.getsize(pf) == 0 for pf in protein_files):
                print("⚠️ No valid proteins predicted - skipping functional/ML analysis")
            else:
                print("2. Running SwissProt annotation...")
                swissprot_results = run_functional_annotation(prokka_dir, output_dir)
                
                print("3. Generating functional reports...")
                generate_functional_report(
                    str(output_dir), 
                    prokka_results=str(prokka_dir), 
                    swissprot_results=swissprot_results
                )
                
                print("4. Creating functional visualizations...")
                create_functional_visualizations(
                    output_dir, 
                    prokka_results=prokka_dir, 
                    swissprot_results=swissprot_results
                )
                
                print("✓ Functional annotation completed")
                
                # ML pathogen prediction
                if ML_PATHOGEN_AVAILABLE and should_use_ml_prediction(prokka_dir):
                    print("\n5. Running ML-based pathogen prediction...")
                    ml_results, ml_summary = run_ml_pathogen_prediction(prokka_dir, output_dir)
                    
                    if ml_results and ml_summary:
                        print(f"✓ ML pathogen prediction completed:")
                        print(f"   📊 {ml_summary['pathogenic_predictions']}/{ml_summary['total_sequences_analyzed']} proteins predicted as pathogenic")
                        print(f"   ⭐ {ml_summary['high_confidence_predictions']} high-confidence predictions")
                        print(f"   🎯 Average confidence: {ml_summary['average_confidence']:.3f}")
                else:
                    print("\n5. ⏭️ Skipping ML prediction")
                    
        except Exception as e:
            print(f"\n✗ Functional/ML analysis failed: {str(e)}")
            import traceback
            traceback.print_exc()
        
        # --- Integrated BLAST+ML Pathogen Report ---
        try:
            print("\n" + "="*80)
            print("STEP 3: INTEGRATED BLAST+ML PATHOGEN ANALYSIS")
            print("="*80)
            
            print("→ Extracting pathogenic organisms from BLAST results...")
            blast_taxonomy_pathogens = []
            if blast_results_data:
                blast_taxonomy_pathogens = extract_pathogens_from_blast_taxonomy(blast_results_data)
                print(f"   ✓ Extracted {len(blast_taxonomy_pathogens)} pathogenic organisms")

            # Generate report if we have any evidence
            if blast_taxonomy_pathogens or (ml_results and ml_summary):
                print("→ Generating integrated BLAST+ML pathogen report...")
                pathogen_report_result = generate_fasta_ml_pathogen_report(
                    str(output_dir), 
                    blast_taxonomy_pathogens=blast_taxonomy_pathogens,
                    ml_results=ml_results,
                    ml_summary=ml_summary
                )
                
                print("→ Creating pathogen visualizations...")
                report_file = output_dir / "blast_ml_pathogen_report.json"
                if report_file.exists():
                    create_pathogen_visualizations(
                        output_dir, 
                        traditional_data=str(report_file)
                    )
                
                print("\n✓ Integrated BLAST+ML pathogen analysis completed")
                print(f"   Risk Level: {pathogen_report_result.get('risk_level', 'UNKNOWN')}")
                print(f"   BLAST Pathogens: {pathogen_report_result.get('blast_pathogens', 0)}")
                print(f"   ML Pathogenic Proteins: {pathogen_report_result.get('ml_pathogenic_proteins', 0)}")
            else:
                print("⚠️ No pathogenic evidence found - no report generated")
                    
        except Exception as e:
            print(f"\n✗ Integrated pathogen report generation failed: {str(e)}")
            import traceback
            traceback.print_exc()
            
    else:
        if args and args.skip_annotation:
            print("\n" + "="*80)
            print("⏭️ SKIPPING ANNOTATION STEPS")
            print("="*80)
            print("--skip-annotation flag detected")
            print("Only BLAST taxonomic classification will be performed")
    
    # ========================================================================
    # FINAL SUMMARY
    # ========================================================================
    print("\n" + "="*80)
    print("FASTA ANALYSIS SUMMARY")
    print("="*80)
    print(f"📁 Output Directory: {output_dir}")
    print("\n📋 Key Files Generated:")
    print("   • taxonomic_classification_report.txt - BLAST-based taxonomy")
    print("   • taxonomic_abundance_chart.html - BLAST hit distribution")
    
    if args and not args.skip_annotation:
        print("   • blast_ml_pathogen_report.txt - Integrated pathogen analysis")
        print("   • pathogen_risk_assessment.html - Risk visualization")
        print("   • functional_annotation_report.txt - Gene prediction results")
        print("   • annotation_quality_dashboard.html - Quality metrics")
    
    print("   • analysis_dashboard.html - Complete interactive dashboard")
    print("\n✓ FASTA pipeline completed successfully")