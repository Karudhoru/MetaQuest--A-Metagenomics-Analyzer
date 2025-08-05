import subprocess
import datetime
import os
if not os.getenv('METAQUEST_QUIET'):
    print("🧬 MetaQuest ML pathogen predictor ready")
import pandas as pd
from pathlib import Path
from Bio import SeqIO
import numpy as np
import json
from .taxonomic_analysis import run_kraken, run_bracken, run_fasta_blast_taxonomy
from .pathogen_analysis import run_pathogen_scan, run_antimicrobial_resistance_scan, run_virulence_factor_scan
from .functional_analysis import run_prokka, run_swissprot_annotation
from ..config import *
from ..io.file_validator import FileValidator
from ..io.utils import check_dependencies, convert_fastq_to_fasta, split_interleaved
from .taxonomic_analysis import run_kraken, run_bracken, run_fasta_blast_taxonomy
from .pathogen_analysis import run_pathogen_scan, run_antimicrobial_resistance_scan, run_virulence_factor_scan
from .functional_analysis import run_prokka, run_swissprot_annotation
from ..reporting.reporting import (
    PathogenReporter,
    generate_taxonomic_report, 
    generate_pathogen_report, 
    generate_functional_report,
    generate_pathogen_summary_report
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
    print("✅ PathogenPredictor with integrated ML functionality loaded successfully")
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

def run_analysis(input_file, file_type, output_dir):
    """Main analysis controller with validation"""
    try:
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True)
        
        # Handle validation for different input types
        validator = FileValidator()
        
        if file_type == 'fastq':
            # For FASTQ, input_file might be a list (paired-end) or single file
            if isinstance(input_file, list):
                # Validate the first file (R1 for paired-end, or single file)
                validation_file = input_file[0]
                print(f"🔍 Validating primary FASTQ file: {validation_file}")
                is_valid, file_stats = validator.validate_and_analyze(validation_file, file_type)
                
                # For paired-end, also validate R2 if present
                if len(input_file) > 1:
                    print(f"🔍 Validating secondary FASTQ file: {input_file[1]}")
                    is_valid_r2, file_stats_r2 = validator.validate_and_analyze(input_file[1], file_type)
                    is_valid = is_valid and is_valid_r2
                    # Combine stats (use R1 stats as primary, note R2 in metadata)
                    file_stats['paired_end'] = True
                    file_stats['r2_file'] = input_file[1]
                    file_stats['r2_sequences'] = file_stats_r2.get('total_sequences', 0)
            else:
                # Single FASTQ file
                validation_file = input_file
                is_valid, file_stats = validator.validate_and_analyze(validation_file, file_type)
                file_stats['paired_end'] = False
        else:
            # For FASTA, input_file should be a single file (already handled in CLI)
            validation_file = input_file[0] if isinstance(input_file, list) else input_file
            is_valid, file_stats = validator.validate_and_analyze(validation_file, file_type)
        
        if not is_valid:
            print("\n❌ Analysis aborted due to validation failure.")
            print("💡 Tips:")
            print("   - For FASTQ: Check quality scores and sequence count")
            print("   - For FASTA: Ensure proper format and unique IDs")
            print("   - Use --validate-only flag to check files without analysis")
            return False
        
        # Save statistics for pipeline use
        import json
        stats_file = output_dir / "input_statistics.json"
        with open(stats_file, 'w') as f:
            json.dump(file_stats, f, indent=2)
        
        print(f"\n🚀 Starting {file_type.upper()} analysis pipeline...")
        print(f"📁 Output directory: {output_dir}\n")
        
        # Continue with existing analysis using the original input_file format
        if file_type == 'fastq':
            analyze_fastq(input_file, output_dir)
        else:
            # For FASTA, ensure we pass the actual file path, not a list
            fasta_path = input_file[0] if isinstance(input_file, list) else input_file
            analyze_fasta(fasta_path, output_dir)
        
        print("\nGenerating analysis dashboard...")
        create_analysis_dashboard(output_dir)
        
        print(f"\n🎉 Analysis complete! Open {output_dir}/analysis_dashboard.html to explore results")
        
    except Exception as e:
        print(f"❌ Analysis failed: {str(e)}")
        raise

def analyze_fastq(reads, output_dir):
    """Process FASTQ files with streamlined pathogen reporting - NO REDUNDANCY"""
    print("\n=== FASTQ Analysis Pipeline ===")
    
    fasta_path = None
    kraken_report = None
    bracken_report = None
    pathogen_results = None
    
    try:
        # Handle interleaved files
        if len(reads) == 1 and reads[0].endswith(('.fq', '.fastq')):
            print("Splitting interleaved FASTQ...")
            reads = split_interleaved(reads[0], output_dir)
        
        # Taxonomic classification
        print("1. Running Kraken2 classification...")
        kraken_report = run_kraken(reads, output_dir)
        
        print("2. Running Bracken abundance estimation...")
        bracken_report = run_bracken(kraken_report, output_dir)
        
        # Generate taxonomic reports and visualizations
        print("3. Generating taxonomic reports and visualizations...")
        generate_taxonomic_report(output_dir, bracken_data=bracken_report)
        create_taxonomic_visualizations(output_dir, bracken_report)
        print("✓ Taxonomic analysis completed")
        
    except Exception as e:
        print(f"⚠️ Taxonomic analysis failed: {str(e)}")
    
    try:
        print("4. Converting to FASTA...")
        fasta_path = convert_fastq_to_fasta(reads if len(reads) > 1 else reads[0], output_dir)
    except Exception as e:
        print(f"⚠️ FASTA conversion failed: {str(e)}")
        print(" Skipping analyses that require FASTA format.")
        return
    
    # Enhanced pathogen screening (only if FASTA conversion succeeded)
    if fasta_path:
        try:
            print("5. Running comprehensive pathogen screening...")
            pathogen_results = run_pathogen_scan(
                fasta_path, output_dir, 
                bracken_results=bracken_report, 
                taxonomy_results=None
            )
            
            print("6. Scanning for AMR genes...")
            run_antimicrobial_resistance_scan(fasta_path, output_dir)
            
            print("7. Scanning for virulence factors...")
            run_virulence_factor_scan(fasta_path, output_dir)
            
        except Exception as e:
            print(f"⚠️ Pathogen analysis failed: {str(e)}")
    
    # Functional annotation with ML pathogen prediction
    prokka_dir = None
    ml_results = None
    ml_summary = None
    
    try:
        print("8. Running gene prediction...")
        prokka_dir = run_prokka(fasta_path, output_dir)
        
        protein_files = list(Path(prokka_dir).glob("*.faa"))
        if not protein_files:
            print("⚠️ Warning: No protein FASTA files found from Prokka. Skipping functional annotation.")
        elif all(os.path.getsize(pf) == 0 for pf in protein_files):
            print("⚠️ Warning: All protein files are empty. Skipping functional annotation.")
        else:
            print("9. Running functional annotation...")
            swissprot_results = run_swissprot_annotation(prokka_dir, output_dir)
            
            print("10. Generating functional analysis reports...")
            generate_functional_report(output_dir, prokka_results=prokka_dir, swissprot_results=swissprot_results)
            create_functional_visualizations(output_dir, prokka_results=prokka_dir, swissprot_results=swissprot_results)
            print("✓ Functional annotation completed")
            
            # ML-based pathogen prediction
            if ML_PATHOGEN_AVAILABLE:
                if should_use_ml_prediction(prokka_dir):
                    try:
                        print("11. Running ML-based pathogen prediction...")
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
                    print("11. ⏭️ Skipping ML prediction - protein sequences too short")
            else:
                print("11. ℹ️ ML pathogen predictor not available - skipping ML analysis")
                
    except Exception as e:
        print(f"⚠️ Functional analysis failed: {str(e)}")
    
    # STREAMLINED pathogen reporting - ONLY ONE COMPREHENSIVE REPORT
    try:
        print("12. Generating comprehensive pathogen analysis...")
        
        # Extract pathogens directly from Bracken data
        bracken_pathogens = []
        taxonomy_pathogens = []
        sequence_pathogens = []
        
        if bracken_report and Path(bracken_report).exists():
            try:
                temp_reporter = PathogenReporter(Path(output_dir))
                bracken_pathogens = temp_reporter._extract_pathogens_from_bracken(bracken_report)
                print(f"   ✓ Extracted {len(bracken_pathogens)} pathogens from Bracken classification")
            except Exception as e:
                print(f"   ⚠️ Could not extract pathogens from Bracken data: {e}")
        
        # Generate THE ONLY comprehensive pathogen report (pathogen_summary.txt + JSON)
        comprehensive_report = generate_pathogen_summary_report(
            output_dir, 
            bracken_pathogens=bracken_pathogens,
            taxonomy_pathogens=taxonomy_pathogens, 
            sequence_pathogens=sequence_pathogens
        )

        # Generate pathogen visualizations
        try:
            comprehensive_report_file = output_dir / "pathogen_detection_report.json"
            if comprehensive_report_file.exists():
                create_pathogen_visualizations(output_dir, traditional_data=str(comprehensive_report_file))
        except Exception as e:
            print(f"⚠️ Pathogen visualization failed: {e}")
        
        # ONLY generate additional JSON for ML integration when available
        if ml_results and ml_summary:
            generate_pathogen_report(
                output_dir, 
                traditional_results=comprehensive_report,
                ml_results=ml_results, 
                ml_summary=ml_summary
            )
            print("✓ Enhanced pathogen analysis with ML integration completed")
        else:
            print("✓ Comprehensive pathogen analysis completed")
                
    except Exception as e:
        print(f"⚠️ Pathogen report generation failed: {str(e)}")
    
    print(f"\n📁 All results saved to: {output_dir}")
    print("🔍 Key files to review:")
    print(" • taxonomic_classification_report.txt - Species identification and abundance")
    print(" • pathogen_summary.txt - Comprehensive pathogen analysis with clinical recommendations")  # THE ONLY ONE
    print(" • functional_annotation_report.txt - Gene prediction and annotation")
    if ml_results:
        print(" • pathogen_detection_report.json - Enhanced ML+traditional pathogen analysis")
    print(" • analysis_dashboard.html - Interactive dashboard")

def analyze_fasta(fasta_path, output_dir):
    """Process FASTA files with BLAST taxonomy + ML only - NO traditional pathogen screening"""
    print("\n=== FASTA Analysis Pipeline ===")
    
    blast_results_data = None
    ml_results = None
    ml_summary = None
    
    try:
        print("1. Running BLAST taxonomic classification...")
        blast_results_file = run_fasta_blast_taxonomy(fasta_path, output_dir, database="nt", max_sequences=50)
        
        if blast_results_file and Path(blast_results_file).exists():
            try:
                with open(blast_results_file, 'r') as f:
                    blast_results_data = json.load(f)
                print(f"✓ Loaded {len(blast_results_data)} BLAST results")
                
                print("2. Generating comprehensive taxonomic reports and visualizations...")
                generated_reports = generate_taxonomic_report(output_dir, blast_data=blast_results_data)
                print("✓ Comprehensive BLAST taxonomic analysis completed")
                
                # Create taxonomic visualizations directly from BLAST data (NO kraken-style report)
                print("3. Creating BLAST-based taxonomic visualizations...")
                create_taxonomic_visualizations(output_dir, blast_results_data)
                print("✓ BLAST taxonomic visualizations created")
                    
            except Exception as e:
                print(f"⚠️ Could not load BLAST results: {e}")
                
    except Exception as e:
        print(f"⚠️ BLAST taxonomic analysis failed: {str(e)}")
    
    # SKIP TRADITIONAL PATHOGEN SCREENING FOR FASTA
    print("4. ⏭️ Skipping traditional pathogen screening for FASTA analysis")
    print("   🎯 Using BLAST taxonomy + ML predictions only")
    
    # Functional annotation with ML pathogen prediction
    try:
        print("5. Running gene prediction...")
        prokka_dir = run_prokka(fasta_path, output_dir)
        
        protein_files = list(Path(prokka_dir).glob("*.faa"))
        if not protein_files:
            print("⚠️ Warning: No protein FASTA files found from Prokka. Skipping functional annotation.")
        elif all(os.path.getsize(pf) == 0 for pf in protein_files):
            print("⚠️ Warning: All protein files are empty. Skipping functional annotation.")
        else:
            print("6. Running functional annotation...")
            swissprot_results = run_swissprot_annotation(prokka_dir, output_dir)
            
            print("7. Generating functional analysis reports...")
            generate_functional_report(output_dir, prokka_results=prokka_dir, swissprot_results=swissprot_results)
            create_functional_visualizations(output_dir, prokka_results=prokka_dir, swissprot_results=swissprot_results)
            
            # ML pathogen prediction
            if ML_PATHOGEN_AVAILABLE:
                if should_use_ml_prediction(prokka_dir):
                    try:
                        print("8. Running ML-based pathogen prediction...")
                        ml_results, ml_summary = run_ml_pathogen_prediction(prokka_dir, output_dir)
                        
                        if ml_results and ml_summary:
                            print(f"✓ ML pathogen prediction completed:")
                            print(f" 📊 {ml_summary['pathogenic_predictions']}/{ml_summary['total_sequences_analyzed']} proteins predicted as pathogenic")
                            print(f" ⭐ {ml_summary['high_confidence_predictions']} high-confidence predictions")
                        else:
                            print("⚠️ ML pathogen prediction produced no results")
                            
                    except Exception as e:
                        print(f"⚠️ ML pathogen prediction failed: {str(e)}")
                else:
                    print("8. ⏭️ Skipping ML prediction - protein sequences too short")
            else:
                print("8. ℹ️ ML pathogen predictor not available - skipping ML analysis")
                
    except Exception as e:
        print(f"⚠️ Functional analysis failed: {str(e)}")
    
    # INTEGRATED BLAST+ML pathogen reporting (REQUIRES MISSING METHODS)
    try:
        print("9. Generating integrated BLAST+ML pathogen analysis...")
        
        # Extract pathogens from BLAST taxonomy - REQUIRES _extract_pathogens_from_blast_taxonomy()
        blast_taxonomy_pathogens = []
        if blast_results_data:
            try:
                temp_reporter = PathogenReporter(Path(output_dir))
                blast_taxonomy_pathogens = temp_reporter._extract_pathogens_from_blast_taxonomy(blast_results_data)
                print(f"   ✓ Extracted {len(blast_taxonomy_pathogens)} pathogens from BLAST taxonomy")
            except Exception as e:
                print(f"   ⚠️ Could not extract pathogens from BLAST taxonomy: {e}")
        
        # Generate INTEGRATED BLAST+ML pathogen report - REQUIRES generate_blast_ml_integrated_report()
        if blast_taxonomy_pathogens or (ml_results and ml_summary):
            integrated_report = temp_reporter.generate_blast_ml_integrated_report(
                output_dir, 
                blast_taxonomy_pathogens=blast_taxonomy_pathogens,
                ml_results=ml_results,
                ml_summary=ml_summary
            )
            
            # Generate visualizations for integrated report
            try:
                integrated_report_file = output_dir / "blast_ml_integrated_pathogen_report.json"
                if integrated_report_file.exists():
                    create_pathogen_visualizations(output_dir, traditional_data=str(integrated_report_file))
            except Exception as e:
                print(f"⚠️ Pathogen visualization failed: {e}")
                
            print("✓ Integrated BLAST+ML pathogen analysis completed")
        else:
            print("⚠️ No pathogenic evidence found in BLAST or ML analysis")
                
    except Exception as e:
        print(f"⚠️ Integrated pathogen analysis failed: {str(e)}")
    
    print(f"\n📁 All results saved to: {output_dir}")
    print("🔍 Key files to review:")
    print(" • taxonomic_classification_report.txt - Comprehensive BLAST-based species identification")
    print(" • blast_ml_pathogen_summary.txt - Integrated BLAST+ML pathogen analysis")  # INTEGRATED REPORT
    print(" • functional_annotation_report.txt - Gene prediction and annotation quality")
    if ml_results:
        print(" • ml_pathogen_predictions.csv - ML pathogen predictions data")
    print(" • analysis_dashboard.html - Interactive dashboard")

def create_analysis_dashboard(output_dir):
    """Create analysis-specific dashboards - fixed for actual file names"""
    try:
        # Fixed detection logic based on your actual FASTQ file names
        is_fastq_analysis = (
            (output_dir / "bracken_report.tsv").exists() or  # ACTUAL Bracken file
            (output_dir / "bracken_report.txt").exists() or  # Alternative Bracken file
            (output_dir / "kraken_report.txt").exists() or   # ACTUAL Kraken file
            (output_dir / "pathogen_summary.txt").exists()   # FASTQ-specific pathogen report
        )
        
        is_fasta_analysis = (
            (output_dir / "blast_taxonomy_results.json").exists() or
            (output_dir / "comprehensive_blast_taxonomy_summary.txt").exists() or
            (output_dir / "blast_ml_pathogen_summary.txt").exists()
        )
        
        print(f"🔍 Analysis type detection:")
        print(f"   • bracken_report.tsv: {(output_dir / 'bracken_report.tsv').exists()}")
        print(f"   • kraken_report.txt: {(output_dir / 'kraken_report.txt').exists()}")
        print(f"   • pathogen_summary.txt: {(output_dir / 'pathogen_summary.txt').exists()}")
        print(f"   • FASTQ analysis: {is_fastq_analysis}")
        print(f"   • FASTA analysis: {is_fasta_analysis}")
        
        if is_fastq_analysis:
            print("📊 Detected FASTQ analysis - creating FASTQ dashboard")
            create_fastq_dashboard(output_dir)
        elif is_fasta_analysis:
            print("📊 Detected FASTA analysis - creating FASTA dashboard")
            create_fasta_dashboard(output_dir)
        else:
            print("⚠️ Could not determine analysis type - creating generic dashboard")
            create_generic_dashboard(output_dir)
            
    except Exception as e:
        print(f"⚠️ Dashboard creation failed: {str(e)}")


def create_fastq_dashboard(output_dir):
    """Create FASTQ-specific dashboard with Bracken + traditional pathogen screening"""
    print("Creating FASTQ analysis dashboard...")
    
    # Check FASTQ-specific files
    files_exist = {
        'taxonomic_report': (output_dir / "taxonomic_classification_report.txt").exists(),
        'pathogen_summary': (output_dir / "pathogen_summary.txt").exists(),  # ONLY THIS ONE for FASTQ
        'functional_report': (output_dir / "functional_annotation_report.txt").exists(),
        'ml_pathogen_json': (output_dir / "pathogen_detection_report.json").exists(),
        'taxonomy_abundance': (output_dir / "taxonomic_abundance_chart.html").exists(),
        'taxonomy_krona': (output_dir / "taxonomy_krona.html").exists(),
        'pathogen_risk': (output_dir / "pathogen_risk_detection.html").exists(),
        'detection_coverage': (output_dir / "detection_method_coverage.html").exists(),
        'annotation_quality': (output_dir / "annotation_quality_analysis.html").exists(),
        'protein_length': (output_dir / "protein_length_distribution.html").exists(),
    }
    
    has_ml = files_exist['ml_pathogen_json']
    ml_status = "✅ Available" if ML_PATHOGEN_AVAILABLE else "⚠️ Not available"
    
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>MetaQuest FASTQ Analysis Dashboard</title>
    <style>
        body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 0; padding: 20px; background: linear-gradient(135deg, #4CAF50 0%, #45a049 100%); color: #333; }}
        .container {{ max-width: 1200px; margin: 0 auto; background: white; border-radius: 15px; box-shadow: 0 8px 32px rgba(0,0,0,0.1); padding: 30px; }}
        .header {{ text-align: center; margin-bottom: 40px; }}
        .header h1 {{ color: #2c3e50; margin: 0; font-size: 2.5em; font-weight: 300; }}
        .header .analysis-type {{ background: #4CAF50; color: white; padding: 5px 15px; border-radius: 20px; font-size: 0.9em; display: inline-block; margin-top: 10px; }}
        .header p {{ color: #7f8c8d; font-size: 1.1em; margin: 10px 0; }}
        .section {{ margin: 30px 0; padding: 25px; border-radius: 10px; background: #f8f9fa; border-left: 5px solid #4CAF50; }}
        .section h2 {{ color: #2c3e50; margin-top: 0; font-size: 1.5em; }}
        .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; margin: 20px 0; }}
        .card {{ background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); border-top: 3px solid #4CAF50; }}
        .card h3 {{ margin-top: 0; color: #2c3e50; }}
        .card a {{ color: #4CAF50; text-decoration: none; font-weight: 500; }}
        .card a:hover {{ text-decoration: underline; }}
        .metadata {{ background: #e8f5e8; padding: 15px; border-radius: 8px; margin: 20px 0; border-left: 4px solid #4CAF50; }}
        .highlight-fastq {{ background: linear-gradient(120deg, #a8e6a3 0%, #88d498 100%); padding: 20px; border-radius: 8px; margin: 20px 0; }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 MetaQuest Analysis Dashboard</h1>
            <div class="analysis-type">FASTQ Analysis</div>
            <p>Comprehensive metagenomics with Kraken2/Bracken + traditional pathogen screening</p>
        </div>

        <div class="metadata">
            <strong>Analysis Type:</strong> FASTQ (Read-based)<br>
            <strong>Taxonomic Method:</strong> Kraken2/Bracken<br>
            <strong>Pathogen Detection:</strong> Multi-source traditional screening<br>
            <strong>Generated:</strong> {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
        </div>

        <div class="section">
            <h2>📊 Taxonomic Classification (Kraken2/Bracken)</h2>
            <p>Fast, accurate taxonomic profiling of metagenomic reads with abundance estimation.</p>
            <div class="grid">"""
    
    if files_exist['taxonomic_report']:
        html_content += f"""
                <div class="card">
                    <h3>📋 Classification Report</h3>
                    <p>Comprehensive taxonomic analysis with diversity metrics and top taxa.</p>
                    <a href="taxonomic_classification_report.txt">View Report</a>
                </div>"""
    
    if files_exist['taxonomy_abundance']:
        html_content += f"""
                <div class="card">
                    <h3>📈 Abundance Chart</h3>
                    <p>Interactive visualization of top 15 most abundant taxa.</p>
                    <a href="taxonomic_abundance_chart.html">View Chart</a>
                </div>"""
    
    if files_exist['taxonomy_krona']:
        html_content += f"""
                <div class="card">
                    <h3>🌐 Krona Plot</h3>
                    <p>Hierarchical taxonomic composition visualization.</p>
                    <a href="taxonomy_krona.html">View Krona</a>
                </div>"""
    
    html_content += """
            </div>
        </div>"""
    
    # FASTQ-specific pathogen section
    html_content += f"""
        <div class="section">
            <h2>🦠 Pathogen Detection (Multi-Source Traditional)</h2>
            <p>Comprehensive pathogen screening using Bracken taxonomy + custom databases + AMR/VF analysis.</p>
            <div class="grid">"""
    
    if files_exist['pathogen_summary']:
        html_content += f"""
                <div class="card">
                    <h3>📝 Comprehensive Pathogen Analysis</h3>
                    <p>THE definitive pathogen report with clinical risk assessment from multi-source detection.</p>
                    <a href="pathogen_summary.txt">View Analysis</a>
                </div>"""
    
    if files_exist['pathogen_risk']:
        html_content += f"""
                <div class="card">
                    <h3>⚠️ Risk Detection Chart</h3>
                    <p>Interactive visualization of detected pathogens by risk level.</p>
                    <a href="pathogen_risk_detection.html">View Chart</a>
                </div>"""
    
    if files_exist['detection_coverage']:
        html_content += f"""
                <div class="card">
                    <h3>🎯 Detection Method Coverage</h3>
                    <p>Multi-method pathogen detection coverage analysis.</p>
                    <a href="detection_method_coverage.html">View Coverage</a>
                </div>"""
    
    html_content += """
            </div>
        </div>"""
    
    # Functional analysis section
    html_content += create_functional_section(files_exist)
    
    # FASTQ report index
    html_content += f"""
        <div class="section">
            <h2>📚 FASTQ Analysis Reports</h2>
            <ul>"""
    
    fastq_reports = [
        ("taxonomic_classification_report.txt", "📊 Taxonomic Classification (Kraken2/Bracken)"),
        ("pathogen_summary.txt", "📝 THE Definitive Pathogen Analysis (Multi-source)"),
        ("functional_annotation_report.txt", "🧪 Functional Annotation Report"),
        ("taxonomic_classification_report.json", "📋 Taxonomic Data Export"),
        ("functional_annotation_report.json", "🔬 Functional Data Export")
    ]
    
    for filename, description in fastq_reports:
        if (output_dir / filename).exists():
            html_content += f'<li><a href="{filename}">{description}</a></li>'
    
    html_content += """
            </ul>
        </div>
    </div>
</body>
</html>"""
    
    # Save FASTQ dashboard
    dashboard_path = output_dir / "analysis_dashboard.html"
    with open(dashboard_path, 'w') as f:
        f.write(html_content)
    
    print(f"✓ FASTQ analysis dashboard created: {dashboard_path}")

def create_fasta_dashboard(output_dir):
    """Create FASTA-specific dashboard with BLAST taxonomy + ML predictions"""
    print("Creating FASTA analysis dashboard...")
    
    # Check FASTA-specific files
    files_exist = {
        'taxonomic_report': (output_dir / "taxonomic_classification_report.txt").exists(),
        'blast_ml_pathogen_summary': (output_dir / "blast_ml_pathogen_summary.txt").exists(),  # FASTA-specific
        'functional_report': (output_dir / "functional_annotation_report.txt").exists(),
        'ml_predictions_csv': (output_dir / "ml_pathogen_predictions.csv").exists(),
        'taxonomy_abundance': (output_dir / "taxonomic_abundance_chart.html").exists(),
        'pathogen_risk': (output_dir / "pathogen_risk_detection.html").exists(),
        'annotation_quality': (output_dir / "annotation_quality_analysis.html").exists(),
        'protein_length': (output_dir / "protein_length_distribution.html").exists(),
        'organism_comparison': (output_dir / "organism_comparison_data.csv").exists(),
    }
    
    has_ml = files_exist['ml_predictions_csv']
    ml_status = "✅ Available" if ML_PATHOGEN_AVAILABLE else "⚠️ Not available"
    
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>MetaQuest FASTA Analysis Dashboard</title>
    <style>
        body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 0; padding: 20px; background: linear-gradient(135deg, #2196F3 0%, #1976D2 100%); color: #333; }}
        .container {{ max-width: 1200px; margin: 0 auto; background: white; border-radius: 15px; box-shadow: 0 8px 32px rgba(0,0,0,0.1); padding: 30px; }}
        .header {{ text-align: center; margin-bottom: 40px; }}
        .header h1 {{ color: #2c3e50; margin: 0; font-size: 2.5em; font-weight: 300; }}
        .header .analysis-type {{ background: #2196F3; color: white; padding: 5px 15px; border-radius: 20px; font-size: 0.9em; display: inline-block; margin-top: 10px; }}
        .header p {{ color: #7f8c8d; font-size: 1.1em; margin: 10px 0; }}
        .section {{ margin: 30px 0; padding: 25px; border-radius: 10px; background: #f8f9fa; border-left: 5px solid #2196F3; }}
        .section h2 {{ color: #2c3e50; margin-top: 0; font-size: 1.5em; }}
        .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; margin: 20px 0; }}
        .card {{ background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); border-top: 3px solid #2196F3; }}
        .card h3 {{ margin-top: 0; color: #2c3e50; }}
        .card a {{ color: #2196F3; text-decoration: none; font-weight: 500; }}
        .card a:hover {{ text-decoration: underline; }}
        .metadata {{ background: #e3f2fd; padding: 15px; border-radius: 8px; margin: 20px 0; border-left: 4px solid #2196F3; }}
        .two-section {{ background: linear-gradient(120deg, #e1f5fe 0%, #b3e5fc 100%); padding: 20px; border-radius: 8px; margin: 20px 0; }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 MetaQuest Analysis Dashboard</h1>
            <div class="analysis-type">FASTA Analysis</div>
            <p>BLAST-based taxonomic classification with ML pathogen predictions</p>
        </div>

        <div class="metadata">
            <strong>Analysis Type:</strong> FASTA (Sequence-based)<br>
            <strong>Taxonomic Method:</strong> BLAST against NCBI database<br>
            <strong>Pathogen Detection:</strong> BLAST taxonomy + ML predictions (NO traditional screening)<br>
            <strong>ML Status:</strong> {ml_status}<br>
            <strong>Generated:</strong> {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
        </div>

        <div class="section">
            <h2>📊 Taxonomic Classification (BLAST)</h2>
            <p>High-quality sequence-based taxonomic identification against NCBI database.</p>
            <div class="grid">"""
    
    if files_exist['taxonomic_report']:
        html_content += f"""
                <div class="card">
                    <h3>📋 BLAST Classification Report</h3>
                    <p>Comprehensive BLAST-based taxonomic analysis with quality metrics.</p>
                    <a href="taxonomic_classification_report.txt">View Report</a>
                </div>"""
    
    if files_exist['organism_comparison']:
        html_content += f"""
                <div class="card">
                    <h3>📊 Organism Comparison</h3>
                    <p>Detailed organism hit statistics and comparison data.</p>
                    <a href="organism_comparison_data.csv">View Data</a>
                </div>"""
    
    if files_exist['taxonomy_abundance']:
        html_content += f"""
                <div class="card">
                    <h3>📈 Abundance Visualization</h3>
                    <p>Interactive BLAST-based taxonomic abundance chart.</p>
                    <a href="taxonomic_abundance_chart.html">View Chart</a>
                </div>"""
    
    html_content += """
            </div>
        </div>"""
    
    # FASTA-specific integrated pathogen section
    html_content += f"""
        <div class="two-section">
            <h2>🦠 Integrated BLAST+ML Pathogen Analysis</h2>
            <p>Dual-method pathogen detection: BLAST taxonomy identification + ML pathogenicity predictions.</p>
            <div class="grid">"""
    
    if files_exist['blast_ml_pathogen_summary']:
        html_content += f"""
                <div class="card">
                    <h3>📝 Integrated BLAST+ML Report</h3>
                    <p>THE definitive report with separated BLAST taxonomy and ML prediction sections.</p>
                    <a href="blast_ml_pathogen_summary.txt">View Integrated Analysis</a>
                </div>"""
    
    if files_exist['ml_predictions_csv']:
        html_content += f"""
                <div class="card">
                    <h3>🤖 ML Pathogen Predictions</h3>
                    <p>Individual protein pathogenicity predictions with confidence scores (CSV format).</p>
                    <a href="ml_pathogen_predictions.csv">View Predictions</a>
                </div>"""
    
    if files_exist['pathogen_risk']:
        html_content += f"""
                <div class="card">
                    <h3>⚠️ Integrated Risk Chart</h3>
                    <p>Combined BLAST+ML pathogen risk visualization.</p>
                    <a href="pathogen_risk_detection.html">View Chart</a>
                </div>"""
    
    html_content += """
            </div>
        </div>"""
    
    # Functional analysis section
    html_content += create_functional_section(files_exist)
    
    # FASTA report index
    html_content += f"""
        <div class="section">
            <h2>📚 FASTA Analysis Reports</h2>
            <ul>"""
    
    fasta_reports = [
        ("taxonomic_classification_report.txt", "📊 BLAST Taxonomic Classification"),
        ("functional_annotation_report.txt", "🧪 Functional Annotation Report"),
        ("ml_pathogen_predictions.csv", "🤖 ML Pathogen Predictions (CSV)"),
        ("organism_comparison_data.csv", "📊 BLAST Organism Comparison Data")
    ]
    
    for filename, description in fasta_reports:
        if (output_dir / filename).exists():
            html_content += f'<li><a href="{filename}">{description}</a></li>'
    
    html_content += """
            </ul>
        </div>
    </div>
</body>
</html>"""
    
    # Save FASTA dashboard
    dashboard_path = output_dir / "analysis_dashboard.html"
    with open(dashboard_path, 'w') as f:
        f.write(html_content)
    
    print(f"✓ FASTA analysis dashboard created: {dashboard_path}")

def create_functional_section(files_exist):
    """Common functional analysis section for both dashboards"""
    html_content = f"""
        <div class="section">
            <h2>🧪 Functional Analysis</h2>
            <p>Gene prediction, protein annotation, and functional characterization.</p>
            <div class="grid">"""
    
    if files_exist['functional_report']:
        html_content += f"""
                <div class="card">
                    <h3>🔬 Annotation Report</h3>
                    <p>Comprehensive gene prediction and functional annotation analysis.</p>
                    <a href="functional_annotation_report.txt">View Report</a>
                </div>"""
    
    if files_exist.get('annotation_quality'):
        html_content += f"""
                <div class="card">
                    <h3>📊 Quality Analysis</h3>
                    <p>SwissProt annotation quality metrics and distribution analysis.</p>
                    <a href="annotation_quality_analysis.html">View Analysis</a>
                </div>"""
    
    if files_exist.get('protein_length'):
        html_content += f"""
                <div class="card">
                    <h3>📏 Protein Length Distribution</h3>
                    <p>Analysis of predicted protein lengths with quality indicators.</p>
                    <a href="protein_length_distribution.html">View Distribution</a>
                </div>"""
    
    html_content += """
            </div>
        </div>"""
    
    return html_content

def create_generic_dashboard(output_dir):
    """Fallback dashboard when analysis type cannot be determined"""
    dashboard_path = output_dir / "analysis_dashboard.html"
    with open(dashboard_path, 'w') as f:
        f.write("""<!DOCTYPE html>
<html><head><title>MetaQuest Analysis</title></head>
<body><h1>Analysis Complete</h1><p>Please check individual report files.</p></body>
</html>""")
    print(f"✓ Generic dashboard created: {dashboard_path}")
