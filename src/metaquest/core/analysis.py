"""
MetaQuest Core Analysis Module v4.0.0
======================================

Main analysis orchestrator with incremental reporting and visualization.
Each analysis step immediately generates its reports and visualizations.
"""

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
from ..ml.pathogen_predictor import run_ml_pathogen_prediction
from ..io.utils import (
    assemble_reads_to_fasta,
    split_interleaved,
    should_use_ml_prediction,
    extract_pathogens_from_bracken,
    extract_pathogens_from_blast_taxonomy
)
from ..io.output_formatter import get_formatter
from ..visualization.dashboard import create_dashboard

# Import NEW v4.0.0 reporting system
from ..reporting.main_reporter import MainReporter
from ..visualization.main_visualizer import (
    create_taxonomic_visualizations,
    create_pathogen_visualizations,
    create_functional_visualizations
)


def run_analysis(input_file, file_type, output_dir, cli_args=None):
    """Main analysis controller that accepts CLI arguments and calls the correct pipeline."""
    fmt = get_formatter()

    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)

    if file_type == 'fastq':
        analyze_fastq(input_file, output_dir_path, cli_args)
    else:
        fasta_path = input_file[0] if isinstance(input_file, list) else input_file
        analyze_fasta(Path(fasta_path), output_dir_path, cli_args)

    # Generate final dashboard
    fmt.section_header("GENERATING FINAL DASHBOARD")
    with fmt.spinner("Creating interactive analysis dashboard"):
        create_dashboard(analysis_type=file_type, output_dir=output_dir_path)
    fmt.success("Dashboard generated successfully")

    fmt.success(f"Analysis complete - view results at: {output_dir_path / 'analysis_dashboard.html'}")


def analyze_fastq(reads, output_dir: Path, args):
    """
    Process FASTQ files with incremental reporting.
    Each analysis step immediately generates reports and visualizations.
    """

    fmt = get_formatter()
    
    # Initialize main reporter
    reporter = MainReporter(output_dir, view_mode='both')

    # Display configuration
    fmt.info(f"Output Directory: {output_dir}")
    input_type = 'Interleaved' if (hasattr(args, 'interleaved') and args.interleaved) else \
                 'Paired-end' if (hasattr(args, 'paired') and args.paired) else 'Single-end'
    fmt.info(f"Input Type: {input_type}")

    fasta_path = None
    total_steps = 1 if (args and args.skip_annotation) else 4

    # ========================================================================
    # STEP 1: TAXONOMIC CLASSIFICATION → IMMEDIATE REPORTING
    # ========================================================================
    bracken_file = None
    try:
        fmt.step_header(1, total_steps, "TAXONOMIC CLASSIFICATION")

        # Handle input format
        if hasattr(args, 'interleaved') and args.interleaved:
            fmt.operation("Splitting interleaved FASTQ file")
            reads = split_interleaved(reads[0], output_dir, fmt)
            fmt.success("Interleaved file split successfully")
        elif hasattr(args, 'single') and args.single:
            reads = [reads] if isinstance(reads, str) else reads

        # Run Kraken/Bracken
        kraken_report_path = run_kraken(reads, output_dir)
        bracken_file = run_bracken(kraken_report_path, output_dir)

        # ►►► IMMEDIATE REPORTING & VISUALIZATION ◄◄◄
        fmt.section_header("► Generating Taxonomic Reports")
        
        # Calculate preliminary taxonomic risk
        from ..reporting.risk_scoring import RiskScorer
        pathogen_organism_db = Path("databases/pathogen_organisms.txt")
        if not pathogen_organism_db.exists():
            pathogen_organism_db = _create_default_pathogen_db(output_dir)
        risk_scorer = RiskScorer(pathogen_organism_file=pathogen_organism_db)
        bracken_df = pd.read_csv(bracken_file, sep='\t')
        taxonomic_risk = risk_scorer.calculate_taxonomic_risk(bracken_df)
        
        # Build preliminary risk data
        preliminary_risk_data = {
            'taxonomic': taxonomic_risk,
            'functional': {'score': 0, 'details': {}},
            'ml': {'score': 0, 'details': []},
            'integrated': {
                'final_score': taxonomic_risk['score'],
                'risk_level': _get_risk_level(taxonomic_risk['score']),
                'tier_scores': {
                    'taxonomic': taxonomic_risk['score'],
                    'functional': 0,
                    'ml': 0
                }
            }
        }
        
        taxonomic_report = reporter.generate_taxonomy_report(
            bracken_file=Path(bracken_file),
            risk_data=preliminary_risk_data
        )
        
        # Save standalone report
        with open(output_dir / "01_taxonomic_report.txt", 'w') as f:
            f.write(taxonomic_report)
        fmt.success("✓ Taxonomic report saved")

        # Generate visualizations
        with fmt.spinner("Creating taxonomic visualizations"):
            viz_files = create_taxonomic_visualizations(output_dir, bracken_file)
            for viz_name, viz_path in viz_files.items():
                fmt.info(f"  • {viz_name}: {Path(viz_path).name}", indent=2)
        
        fmt.stage_complete("Taxonomic Classification → Reports Generated")

    except Exception as e:
        fmt.error(f"Taxonomic analysis failed: {str(e)}", 
                 solutions=["Verify Kraken2/Bracken installation and database paths.", 
                           "Check input FASTQ file integrity."])
        raise

    # Skip annotation steps if flag set
    if args and args.skip_annotation:
        fmt.section_header("SKIPPING ANNOTATION & PATHOGEN STEPS")
        fmt.info("--skip-annotation flag detected")
        fmt.info("Only taxonomic classification was performed")
        _display_completion_summary(output_dir, fmt, "partial")
        return

    # ========================================================================
    # STEP 2: METAGENOMIC ASSEMBLY
    # ========================================================================
    try:
        fmt.step_header(2, total_steps, "METAGENOMIC ASSEMBLY")
        with fmt.spinner("Assembling reads into contigs (FASTA)"):
            fasta_path = assemble_reads_to_fasta(reads, output_dir, fmt)
        fmt.success(f"Assembly complete: {fasta_path.name}")
        
        # Save assembly stats for reporting
        _save_assembly_stats(fasta_path, output_dir)
        fmt.stage_complete("Assembly → Stats Recorded")
        
    except Exception as e:
        fmt.error(f"Assembly failed: {str(e)}", 
                 solutions=["Ensure sufficient memory for assembly.", 
                           "Verify input FASTQ files are not corrupted."])
        fasta_path = None
        return

    # ========================================================================
    # STEP 3: FUNCTIONAL ANNOTATION → IMMEDIATE REPORTING
    # ========================================================================
    prokka_dir = None
    swissprot_results = None
    
    try:
        fmt.step_header(3, total_steps, "FUNCTIONAL ANNOTATION")
        
        # Run Prokka
        filter_contigs = getattr(args, 'filter_contigs', True)
        min_contig_length = getattr(args, 'min_contig_length', 500)
        kill_tbl2asn = getattr(args, 'kill_tbl2asn', True)
        tbl2asn_timeout = getattr(args, 'tbl2asn_timeout', 30)

        prokka_dir = run_prokka(fasta_path, output_dir, 
                               filter_contigs_flag=filter_contigs, 
                               min_contig_length=min_contig_length, 
                               kill_tbl2asn=kill_tbl2asn, 
                               tbl2asn_timeout=tbl2asn_timeout)

        # Check if proteins were predicted
        protein_files = list(Path(prokka_dir).glob("*.faa"))
        if not protein_files or all(os.path.getsize(pf) == 0 for pf in protein_files):
            fmt.warning("No proteins predicted - skipping functional annotation")
        else:
            # Run SwissProt annotation
            annotation_threads = getattr(args, 'annotation_threads', 8)
            swissprot_results = run_functional_annotation(prokka_dir, output_dir, 
                                                          threads=annotation_threads)

            # ►►► IMMEDIATE REPORTING & VISUALIZATION ◄◄◄
            fmt.section_header("► Generating Functional Reports")
            
            # Calculate functional risk
            functional_df = pd.read_csv(swissprot_results, sep='\t', header=None)

            # *** FIX: Add column names ***
            ann_cols = ['query_id', 'subject_id', 'identity', 'length', 'mismatches', 
                        'gaps', 'q_start', 'q_end', 's_start', 's_end', 'evalue', 
                        'bitscore', 'description']
            functional_df.columns = ann_cols[:len(functional_df.columns)]

            pathogen_hits_file = output_dir / "pathogen_results.txt"

            if pathogen_hits_file.exists():
                pathogen_df = pd.read_csv(pathogen_hits_file, sep='\t', header=None)
                pathogen_df.columns = ann_cols[:len(pathogen_df.columns)]
                
                # *** FIX: Calculate total_cds and pass it ***
                total_cds = len(functional_df['query_id'].unique())
                functional_risk = risk_scorer.calculate_functional_risk(
                    functional_df, 
                    pathogen_df,
                    total_cds=total_cds  # ← ADD THIS!
                )
            else:
                functional_risk = {'score': 0, 'details': {}}
            
            # Update risk data with functional component
            preliminary_risk_data['functional'] = functional_risk
            preliminary_risk_data['integrated']['tier_scores']['functional'] = functional_risk['score']
            preliminary_risk_data['integrated']['final_score'] = (
                preliminary_risk_data['taxonomic']['score'] * 0.4 +
                functional_risk['score'] * 0.3 +
                preliminary_risk_data['ml']['score'] * 0.3
            )
            preliminary_risk_data['integrated']['risk_level'] = _get_risk_level(
                preliminary_risk_data['integrated']['final_score']
            )
            
            sample_info_file = output_dir / "sample.txt"
            functional_report = reporter.generate_functional_report(
                sample_info_file=sample_info_file,
                functional_annotations_file=Path(swissprot_results),
                risk_data=preliminary_risk_data
            )
            
            # Save standalone report
            with open(output_dir / "02_functional_report.txt", 'w') as f:
                f.write(functional_report)
            fmt.success("✓ Functional report saved")

            # Generate visualizations
            with fmt.spinner("Creating functional visualizations"):
                viz_files = create_functional_visualizations(
                    output_dir, 
                    prokka_results=prokka_dir, 
                    swissprot_results=swissprot_results
                )
                for viz_name, viz_path in viz_files.items():
                    fmt.info(f"  • {viz_name}: {Path(viz_path).name}", indent=2)
            
        fmt.stage_complete("Functional Annotation → Reports Generated")
        
    except Exception as e:
        fmt.error(f"Functional analysis failed: {str(e)}", 
                 solutions=["Verify Prokka installation.", 
                           "Check DIAMOND database configuration."])
        import traceback
        traceback.print_exc()

    # CRITICAL FIX: In analysis.py around line 240-260
    # Replace the pathogen detection section with this:

    # ========================================================================
    # STEP 4: PATHOGEN DETECTION → IMMEDIATE REPORTING
    # ========================================================================
    ml_results = None
    ml_summary = None
    pathogen_hits_file = None
    protein_file = prokka_dir / "sample.faa" 

    try:
        fmt.step_header(4, total_steps, "PATHOGEN DETECTION & RISK ASSESSMENT")

        # 4A: Sequence Database Scan
        fmt.section_header("Sequence-Based Pathogen Screening")
        pathogen_hits_file = run_pathogen_scan(
            protein_file, 
            output_dir,
            bracken_results=Path(bracken_file),
            taxonomy_results=None
        )
        
        if pathogen_hits_file and pathogen_hits_file.exists():
            fmt.success(f"Database scan complete: {pathogen_hits_file.name}")
        else:
            fmt.info("No hits in pathogen sequence database")

        # 4B: Machine Learning Prediction
        if prokka_dir and should_use_ml_prediction(prokka_dir, fmt):
            fmt.section_header("ML-Based Pathogen Prediction v2.1")
            try:
                with fmt.spinner("Running ML pathogen prediction"):
                    ml_results, ml_summary = run_ml_pathogen_prediction(prokka_dir, output_dir)
                if ml_results and ml_summary:
                    fmt.success("ML prediction complete")
                else:
                    fmt.warning("ML prediction produced no results")
            except Exception as e:
                fmt.warning(f"ML prediction failed: {str(e)}")
        else:
            fmt.info("ML prediction skipped (insufficient protein data)")

        # ►►► IMMEDIATE PATHOGEN REPORTING & VISUALIZATION ◄◄◄
        fmt.section_header("► Generating Pathogen Risk Reports")
        
        # Calculate comprehensive risk assessment with all three tiers
        from ..reporting.risk_scoring import calculate_all_risks
        
        pathogen_organism_db = Path("databases/pathogen_organisms.txt")
        if not pathogen_organism_db.exists():
            pathogen_organism_db = _create_default_pathogen_db(output_dir)
        
        ml_predictions_file = output_dir / "ml_pathogen_predictions.json"
        
        risk_data = calculate_all_risks(
            bracken_file=Path(bracken_file),
            functional_file=Path(swissprot_results) if swissprot_results else None,
            pathogen_hits_file=pathogen_hits_file,
            ml_predictions_file=ml_predictions_file if ml_predictions_file.exists() else None,
            pathogen_organism_file=pathogen_organism_db
        )
        
        # Generate pathogen risk report using MainReporter
        pathogen_report = reporter.generate_pathogen_report(
            risk_data=risk_data,
            pathogen_hits_file=pathogen_hits_file if pathogen_hits_file else output_dir / "pathogen_results.txt",
            ml_predictions_file=ml_predictions_file
        )
                
        # Save standalone report
        with open(output_dir / "03_pathogen_risk_report.txt", 'w') as f:
            f.write(pathogen_report)
        fmt.success("✓ Pathogen risk report saved")

        pathogen_detections = _build_detailed_pathogen_detections(
                    risk_data, 
                    bracken_file, 
                    pathogen_hits_file
                )
                
        # Save detailed pathogen detections for visualizer
        pathogen_detections_file = output_dir / "pathogen_detections.json"
        with open(pathogen_detections_file, 'w') as f:
            json.dump(pathogen_detections, f, indent=2)                
            fmt.success(f"✓ Detailed pathogen detections saved: {pathogen_detections_file.name}")

        # Generate visualizations using DETAILED detections file
        with fmt.spinner("Creating pathogen visualizations"):
            viz_files = create_pathogen_visualizations(
                output_dir,
                traditional_data=str(pathogen_detections_file)  # ← Use detailed file!
            )
            for viz_name, viz_path in viz_files.items():
                fmt.info(f"  • {viz_name}: {Path(viz_path).name}", indent=2)
        
        # Display risk summary
        _display_risk_summary(risk_data, fmt)
        
        fmt.stage_complete("Pathogen Detection → Reports Generated")

    except Exception as e:
        fmt.error(f"Pathogen analysis failed: {str(e)}", 
                solutions=["Check pathogen database files.", 
                        "Verify ML model files are present."])
        import traceback
        traceback.print_exc()

    # ========================================================================
    # COMPREHENSIVE INTEGRATED REPORT (ALL COMPONENTS)
    # ========================================================================
    try:
        fmt.section_header("═══ COMPREHENSIVE INTEGRATED REPORT ═══")
        fmt.operation("Compiling all analyses into comprehensive report")
        
        sample_info_file = output_dir / "sample.txt"
        functional_file = Path(swissprot_results) if swissprot_results else output_dir / "functional_annotations.tsv"
        
        # Use MainReporter to generate complete integrated report
        comprehensive_report = reporter.generate_report(
            bracken_file=Path(bracken_file),
            sample_info_file=sample_info_file,
            functional_annotations_file=functional_file,
            pathogen_hits_file=pathogen_hits_file if pathogen_hits_file else output_dir / "pathogen_results.txt",
            ml_predictions_file=output_dir / "ml_pathogen_predictions.json",
            pathogen_organism_file=pathogen_organism_db
        )
        
        fmt.success("✓ Comprehensive integrated report generated")
        fmt.info(f"Main report: {output_dir / 'comprehensive_report.txt'}")
        
        # Export data tables using MainReporter
        reporter.export_tables(
            bracken_file=Path(bracken_file),
            annotation_file=functional_file,
            risk_data=risk_data
        )
        fmt.success("✓ Data tables exported")
        
    except Exception as e:
        fmt.warning(f"Comprehensive report generation failed: {str(e)}")
        fmt.info("Individual module reports are still available")
        import traceback
        traceback.print_exc()

    # ========================================================================
    # FINAL SUMMARY
    # ========================================================================
    _display_completion_summary(output_dir, fmt, "complete")


# ═══════════════════════════════════════════════════════════════
# NEW HELPER FUNCTION: Build detailed pathogen detections
# ═══════════════════════════════════════════════════════════════
def _build_detailed_pathogen_detections(risk_data: dict, bracken_file: Path, pathogen_hits_file: Path) -> dict:
    """
    Build detailed pathogen detection data structure for visualization.
    This creates the format expected by pathogenic_visualizer.py
    """
    import pandas as pd
    
    # Read Bracken data
    bracken_df = pd.read_csv(bracken_file, sep='\t')
    
    # Build detailed detection structure
    detections = {
        'data': {
            'summary': {
                'analysis_type': 'FASTQ Metagenomics',
                'total_species': len(bracken_df),
                'pathogen_count': len(risk_data['taxonomic']['pathogens_detected']),
                'overall_risk_assessment': risk_data['integrated']['risk_level']
            },
            'pathogen_detections': {}
        }
    }
    
    # Add each pathogen with full details
    for pathogen_info in risk_data['taxonomic']['pathogens_detected']:
        organism_name = pathogen_info['species']
        
        # Find this organism in Bracken results
        organism_row = bracken_df[bracken_df['name'] == organism_name]
        
        if not organism_row.empty:
            abundance = float(organism_row['fraction_total_reads'].iloc[0]) * 100
            reads = int(organism_row['new_est_reads'].iloc[0])
        else:
            abundance = 0.0
            reads = 0
        
        # Determine risk level from pathogen info
        risk_level = pathogen_info.get('risk_level', 'MEDIUM').upper()
        
        # Determine WHO priority (if available in pathogen database)
        who_priority = _map_risk_to_who_priority(risk_level)
        
        # Build detection methods list
        detection_methods = []
        if reads > 0:
            detection_methods.append('Kraken2/Bracken')
        if pathogen_hits_file and pathogen_hits_file.exists():
            detection_methods.append('Sequence Database')
        
        # Calculate confidence score based on reads and detection methods
        confidence_score = min(1.0, (reads / 100) * len(detection_methods) * 0.1)
        
        detections['data']['pathogen_detections'][organism_name] = {
            'risk_level': risk_level,
            'abundance_percentage': round(abundance, 4),
            'sequence_hits': reads,
            'sequence_identity': 95.0,  # Default for Bracken
            'confidence_score': round(confidence_score, 3),
            'detection_methods': detection_methods,
            'who_priority': who_priority,
            'clinical_notes': _get_clinical_notes(organism_name, risk_level)
        }
    
    return detections

def _map_risk_to_who_priority(risk_level: str) -> str:
    """Map risk level to WHO priority classification"""
    mapping = {
        'HIGH': 'Critical Priority',
        'CRITICAL': 'Critical Priority',
        'MEDIUM': 'High Priority',
        'MODERATE': 'High Priority',
        'LOW': 'Not Listed'
    }
    return mapping.get(risk_level, 'Not Listed')


def _get_clinical_notes(organism_name: str, risk_level: str) -> str:
    """Get clinical notes for common pathogens"""
    notes = {
        'Clostridioides difficile': 'Antibiotic-associated diarrhea, pseudomembranous colitis',
        'Clostridium difficile': 'Antibiotic-associated diarrhea, pseudomembranous colitis',
        'Staphylococcus aureus': 'Skin infections, pneumonia, septicemia. Check for MRSA',
        'Escherichia coli': 'Gastrointestinal infections, UTIs. May carry virulence factors',
        'Salmonella enterica': 'Gastroenteritis, typhoid fever',
        'Streptococcus pneumoniae': 'Pneumonia, meningitis, otitis media',
        'Pseudomonas aeruginosa': 'Nosocomial infections, immunocompromised patients',
        'Klebsiella pneumoniae': 'Pneumonia, UTIs. Check for carbapenem resistance'
    }
    
    default_notes = {
        'HIGH': 'High-risk pathogen - requires immediate clinical attention',
        'CRITICAL': 'Critical pathogen - immediate intervention required',
        'MEDIUM': 'Opportunistic pathogen - monitor clinically',
        'MODERATE': 'Opportunistic pathogen - monitor clinically',
        'LOW': 'Opportunistic pathogen - monitor clinically'
    }
    
    return notes.get(organism_name, default_notes.get(risk_level, 'Monitor for symptoms'))


def analyze_fasta(fasta_path: Path, output_dir: Path, args):
    """
    Process FASTA files with incremental reporting.
    Each analysis step immediately generates reports and visualizations.
    """
    
    fmt = get_formatter()
    
    # Initialize main reporter
    reporter = MainReporter(output_dir, view_mode='both')
    
    fmt.info(f"Output Directory: {output_dir}")
    fmt.info(f"Input File: {fasta_path.name}")
    
    total_steps = 1 if (args and args.skip_annotation) else 3
    
    # ========================================================================
    # STEP 1: BLAST TAXONOMIC CLASSIFICATION → IMMEDIATE REPORTING
    # ========================================================================
    blast_results_data = None
    try:
        fmt.step_header(1, total_steps, "BLAST TAXONOMIC CLASSIFICATION")
        
        max_sequences = getattr(args, 'blast_sample_size', 50) if args else 50
        fmt.operation(f"Running BLAST taxonomy (sampling {max_sequences} sequences)")
        
        blast_results_file = run_fasta_blast_taxonomy(
            fasta_path, output_dir, 
            database="nt", 
            max_sequences=max_sequences
        )
        
        if blast_results_file and Path(blast_results_file).exists():
            with open(blast_results_file, 'r') as f:
                blast_results_data = json.load(f)
            fmt.success(f"BLAST complete: {len(blast_results_data)} results")
            
            # ►►► IMMEDIATE REPORTING & VISUALIZATION ◄◄◄
            fmt.section_header("► Generating BLAST Taxonomic Reports")
            
            # Note: For FASTA with BLAST, we don't have Bracken-style taxonomy
            # So we generate a simplified taxonomy report from BLAST results
            # This could be enhanced in future to convert BLAST to Bracken-like format
            
            fmt.info("BLAST taxonomy reporting (simplified mode)")
            fmt.info("For full taxonomic reports, use FASTQ input with Kraken2/Bracken")
            
            # Generate visualizations
            with fmt.spinner("Creating BLAST visualizations"):
                viz_files = create_taxonomic_visualizations(output_dir, blast_results_data)
                for viz_name, viz_path in viz_files.items():
                    fmt.info(f"  • {viz_name}: {Path(viz_path).name}", indent=2)
            
            fmt.stage_complete("BLAST Taxonomy → Visualizations Generated")
        else:
            fmt.warning("BLAST analysis produced no results")

    except Exception as e:
        fmt.error(f"BLAST analysis failed: {str(e)}",
                 solutions=["Check internet connection", 
                           "Verify NCBI BLAST service status"])
        import traceback
        traceback.print_exc()
    
    # Skip annotation if flag set
    if args and args.skip_annotation:
        fmt.section_header("SKIPPING ANNOTATION STEPS")
        _display_completion_summary(output_dir, fmt, "partial")
        return
    
    # ========================================================================
    # STEP 2: FUNCTIONAL ANNOTATION → IMMEDIATE REPORTING
    # ========================================================================
    prokka_dir = None
    swissprot_results = None
    
    try:
        fmt.step_header(2, total_steps, "FUNCTIONAL ANNOTATION")
        
        filter_contigs = getattr(args, 'filter_contigs', True)
        min_contig_length = getattr(args, 'min_contig_length', 1000)
        kill_tbl2asn = getattr(args, 'kill_tbl2asn', True)
        tbl2asn_timeout = getattr(args, 'tbl2asn_timeout', 300)
        
        prokka_dir = run_prokka(
            fasta_path, output_dir,
            filter_contigs_flag=filter_contigs,
            min_contig_length=min_contig_length,
            kill_tbl2asn=kill_tbl2asn,
            tbl2asn_timeout=tbl2asn_timeout
        )
        
        protein_files = list(Path(prokka_dir).glob("*.faa"))
        if protein_files and any(os.path.getsize(pf) > 0 for pf in protein_files):
            annotation_threads = getattr(args, 'annotation_threads', 8)
            swissprot_results = run_functional_annotation(
                prokka_dir, output_dir, threads=annotation_threads
            )
            
            # ►►► IMMEDIATE REPORTING & VISUALIZATION ◄◄◄
            fmt.section_header("► Generating Functional Reports")
            
            sample_info_file = output_dir / "sample.txt"
            
            # Build preliminary risk data for functional report
            preliminary_risk_data = {
                'taxonomic': {'score': 0},
                'functional': {'score': 0, 'details': {}},
                'ml': {'score': 0},
                'integrated': {'final_score': 0, 'risk_level': 'Pending'}
            }
            
            functional_report = reporter.generate_functional_report(
                sample_info_file=sample_info_file,
                functional_annotations_file=Path(swissprot_results),
                risk_data=preliminary_risk_data
            )
            
            with open(output_dir / "02_functional_report.txt", 'w') as f:
                f.write(functional_report)
            fmt.success("✓ Functional report saved")
            
            with fmt.spinner("Creating functional visualizations"):
                viz_files = create_functional_visualizations(
                    output_dir, prokka_dir, swissprot_results
                )
                for viz_name, viz_path in viz_files.items():
                    fmt.info(f"  • {viz_name}: {Path(viz_path).name}", indent=2)
            
            fmt.stage_complete("Functional Annotation → Reports Generated")
        else:
            fmt.warning("No proteins predicted")
            
    except Exception as e:
        fmt.error(f"Functional analysis failed: {str(e)}")
        import traceback
        traceback.print_exc()

    # ========================================================================
    # STEP 3: ML PATHOGEN PREDICTION → IMMEDIATE REPORTING
    # ========================================================================
    try:
        fmt.step_header(3, total_steps, "ML PATHOGEN PREDICTION")
        
        if prokka_dir and should_use_ml_prediction(prokka_dir, fmt):
            with fmt.spinner("Running ML pathogen prediction"):
                ml_results, ml_summary = run_ml_pathogen_prediction(prokka_dir, output_dir)
            
            if ml_results and ml_summary:
                # ►►► IMMEDIATE REPORTING ◄◄◄
                fmt.section_header("► Generating ML Prediction Report")
                
                # For FASTA, we use BLAST pathogens + ML results
                blast_pathogens = extract_pathogens_from_blast_taxonomy(blast_results_data) if blast_results_data else []
                
                # Build risk data for ML report
                ml_risk_data = {
                    'taxonomic': {'score': 0, 'pathogens_detected': blast_pathogens},
                    'functional': {'score': 0},
                    'ml': {
                        'score': ml_summary.get('pathogenic_percentage', 0),
                        'high_confidence_pathogenic': ml_summary.get('high_confidence_pathogenic', 0),
                        'predictions': ml_results
                    },
                    'integrated': {
                        'final_score': ml_summary.get('pathogenic_percentage', 0),
                        'risk_level': _get_risk_level(ml_summary.get('pathogenic_percentage', 0))
                    }
                }
                
                # Generate report using MainReporter
                ml_report = reporter.generate_pathogen_report(
                    risk_data=ml_risk_data,
                    pathogen_hits_file=output_dir / "pathogen_results.txt",
                    ml_predictions_file=output_dir / "ml_pathogen_predictions.json"
                )
                
                with open(output_dir / "03_ml_pathogen_report.txt", 'w') as f:
                    f.write(ml_report)
                
                fmt.success("✓ ML prediction report saved")
                fmt.stage_complete("ML Prediction → Reports Generated")
        else:
            fmt.info("ML prediction skipped")
            
    except Exception as e:
        fmt.warning(f"ML prediction failed: {str(e)}")

    # ========================================================================
    # FINAL SUMMARY
    # ========================================================================
    _display_completion_summary(output_dir, fmt, "complete")


# ============================================================================
# HELPER FUNCTIONS (Minimal - only what MainReporter doesn't provide)
# ============================================================================

def _save_assembly_stats(fasta_path: Path, output_dir: Path):
    """Save assembly statistics to sample.txt for MainReporter"""
    from Bio import SeqIO
    
    sequences = list(SeqIO.parse(fasta_path, "fasta"))
    total_bases = sum(len(seq.seq) for seq in sequences)
    
    with open(output_dir / "sample.txt", 'w') as f:
        f.write(f"Organism: {fasta_path.stem}\n")
        f.write(f"Contigs: {len(sequences)}\n")
        f.write(f"Total bases: {total_bases}\n")
        f.write(f"Genes: TBD\n")



def _get_risk_level(score: float) -> str:
    """Convert risk score to risk level"""
    if score >= 70:
        return 'High'
    elif score >= 40:
        return 'Moderate'
    else:
        return 'Low'


def _display_risk_summary(risk_data: dict, fmt):
    """Display key risk assessment findings"""
    integrated = risk_data['integrated']
    
    risk_emoji = {
        'High': '🔴',
        'Moderate': '🟡',
        'Low': '🟢'
    }.get(integrated['risk_level'], '⚪')
    
    fmt.section_header("RISK ASSESSMENT SUMMARY")
    fmt.info(f"Overall Risk: {risk_emoji} {integrated['risk_level']} ({integrated['final_score']:.0f}/100)")
    
    tier_scores = integrated.get('tier_scores', {})
    if tier_scores:
        fmt.info("\nTier Breakdown:")
        for tier, score in tier_scores.items():
            fmt.info(f"  • {tier.title()}: {score:.0f}/100", indent=2)


def _create_default_pathogen_db(output_dir: Path) -> Path:
    """Create default pathogen organism database"""
    db_file = output_dir / "pathogen_organisms_default.txt"
    
    default_pathogens = """# Pathogenic organism names
Brucella melitensis	high
Brucella abortus	high
Brucella suis	high
Salmonella enterica	high
Escherichia coli	high
Clostridioides difficile	high
Clostridium difficile	high
Staphylococcus aureus	high
Streptococcus pneumoniae	high
Streptococcus thermophilus	medium
Streptococcus salivarius	medium
Streptococcus pyogenes	high
Mycobacterium tuberculosis	high
Pseudomonas aeruginosa	high
Vibrio cholerae	high
Yersinia pestis	high
Francisella tularensis	high
Bacillus anthracis	high
Listeria monocytogenes	high
Klebsiella pneumoniae	medium
Acinetobacter baumannii	medium
Enterococcus faecium	medium
Campylobacter jejuni	medium
Campylobacter gracilis	medium
Helicobacter pylori	medium
Shigella flexneri	medium
Neisseria meningitidis	medium
Haemophilus influenzae	medium
Burkholderia pseudomallei	medium
Legionella pneumophila	medium
Chlamydia trachomatis	medium
Rickettsia rickettsii	medium
Coxiella burnetii	medium
Erysipelatoclostridium ramosum	low
Clostridium scindens	low
Clostridium innocuum	low
Lachnoclostridium	low
Borrelia burgdorferi	low
Treponema pallidum	low
Mycoplasma pneumoniae	low
Bartonella henselae	low
"""
    
    with open(db_file, 'w') as f:
        f.write(default_pathogens)
    
    return db_file


def _display_completion_summary(output_dir: Path, fmt, analysis_type: str):
    """Display final completion summary"""
    fmt.section_header("ANALYSIS COMPLETE")
    fmt.info(f"Output Directory: {output_dir}")
    
    if analysis_type == "partial":
        fmt.info("Status: Partial analysis (annotation skipped)")
        report_files = ["01_taxonomic_report.txt"]
    else:
        fmt.info("Status: Complete analysis")
        report_files = [
            "01_taxonomic_report.txt",
            "02_functional_report.txt", 
            "03_pathogen_risk_report.txt",
            "comprehensive_report.txt"
        ]
    
    # Check which files actually exist
    existing_files = []
    for report_file in report_files:
        if (output_dir / report_file).exists():
            existing_files.append(report_file)
    
    if existing_files:
        fmt.file_list("Generated Reports", existing_files)
    
    fmt.success("Pipeline completed successfully")




def _map_risk_to_who_priority(risk_level: str) -> str:
    """Map risk level to WHO priority classification"""
    mapping = {
        'HIGH': 'Critical Priority',
        'CRITICAL': 'Critical Priority',
        'MEDIUM': 'High Priority',
        'MODERATE': 'High Priority',
        'LOW': 'Not Listed'
    }
    return mapping.get(risk_level, 'Not Listed')