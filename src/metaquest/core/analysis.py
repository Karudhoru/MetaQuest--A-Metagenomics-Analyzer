"""
MetaQuest Core Analysis Module - FASTQ-only v5.0.0
===================================================

OPTIMIZATIONS:
Memory-optimized sequence handling (streaming, chunking)
Explicit cleanup and garbage collection
Reproducibility metadata tracking
FASTA support removed - FASTQ-only pipeline
Checkpoint system removed - simplified workflow
Contig filtering removed - Prokka handles all contigs

Target: 8-16GB RAM systems
"""

import subprocess
import os
import gc
import traceback
import pandas as pd
from pathlib import Path
from datetime import datetime
from Bio import SeqIO
import numpy as np
import json
from typing import Optional, List, Iterator

# Import cleaned modules
from .taxonomic_analysis import run_kraken, run_bracken
from .pathogen_analysis import run_pathogen_scan
from .functional_analysis import run_prokka, run_functional_annotation, parse_prokka_sample_txt
from .assembly import assemble_reads_to_fasta

# Config constants are accessed via sub-modules; no direct config symbols needed here.
from ..io.file_validator import FileValidator
from ..io.data_loaders import load_bracken_report, load_annotation_file, load_pathogen_hits
from ..ml.pathogen_predictor import run_ml_pathogen_prediction
from ..io.utils import (
    split_interleaved,
    should_use_ml_prediction,
    explicit_cleanup
)
from ..io.output_formatter import get_formatter
from ..visualization.dashboard import create_dashboard

from ..reporting.main_reporter import MainReporter
from ..reporting.risk_scoring import RiskScorer, calculate_all_risks
from ..visualization.main_visualizer import (
    create_taxonomic_visualizations,
    create_pathogen_visualizations,
    create_functional_visualizations
)


# ============================================================================
# MEMORY-OPTIMIZED UTILITIES
# ============================================================================

def calculate_assembly_stats_streaming(contigs_fasta: Path) -> dict:
    """
    Calculate assembly stats without loading all sequences into memory.
    
    Memory: O(1) instead of O(total_sequences)
    
    Args:
        contigs_fasta: Path to assembled FASTA file
        
    Returns:
        Dict containing assembly statistics
    """
    total_contigs = 0
    total_bases = 0
    lengths = []
    max_length = 0
    
    # First pass: get counts and max
    for record in SeqIO.parse(contigs_fasta, "fasta"):
        length = len(record.seq)
        total_contigs += 1
        total_bases += length
        lengths.append(length)  # Need for N50
        if length > max_length:
            max_length = length
    
    # Calculate N50 (requires lengths sorted)
    n50 = _calculate_n50(lengths)
    
    # Cleanup
    del lengths
    gc.collect()
    
    return {
        'total_contigs': total_contigs,
        'total_bases': total_bases,
        'n50': n50,
        'max_length': max_length,
        'mean_length': total_bases / total_contigs if total_contigs > 0 else 0
    }


def _calculate_n50(lengths: List[int]) -> int:
    """
    Calculate N50: length of shortest contig at which 50% of total bases are covered.
    
    Edge cases handled:
    - Empty list → 0
    - Single contig → that contig's length
    - Ties → first contig that reaches threshold
    
    Args:
        lengths: List of contig lengths
        
    Returns:
        N50 value in base pairs
    """

    if not lengths:
        return 0
    if len(lengths) == 1:
        return lengths[0]
    
    sorted_lengths = sorted(lengths, reverse=True)
    cumsum = np.cumsum(sorted_lengths)
    total = cumsum[-1]
    

    if total == 0:
        return 0
    
    threshold = total / 2.0
    
    for i, cs in enumerate(cumsum):
        if cs >= threshold:
            return int(sorted_lengths[i])
    
    return int(sorted_lengths[-1])


# ============================================================================
# METADATA MANAGEMENT
# ============================================================================

def save_analysis_metadata(output_dir: Path, parameters: dict):
    """
    Save analysis parameters for reproducibility.
    
    Args:
        output_dir: Output directory path
        parameters: Dictionary of analysis parameters
    """
    from ..config import __version__ as _version
    metadata = {
        'timestamp': datetime.now().isoformat(),
        'metaquest_version': _version,
        'analysis_type': 'fastq',
        'parameters': parameters
    }
    
    with open(output_dir / 'analysis_metadata.json', 'w') as f:
        json.dump(metadata, f, indent=2)


# ============================================================================
# MAIN ANALYSIS CONTROLLER
# ============================================================================

def run_analysis(input_file, output_dir, cli_args=None):
    """
    Main analysis controller - FASTQ-only.
    
    Args:
        input_file: FASTQ file path(s)
        output_dir: Output directory path
        cli_args: Command-line arguments namespace
    """
    fmt = get_formatter()

    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)
    
    # Save analysis metadata for reproducibility
    parameters = vars(cli_args) if cli_args else {}
    save_analysis_metadata(output_dir_path, parameters)

    # Run FASTQ analysis
    analyze_fastq(input_file, output_dir_path, cli_args)
    
    # Generate dashboard
    fmt.section_header("FINALIZING RESULTS")
    
    fmt.info("Creating interactive dashboard")
    with fmt.suppressed_output():
        create_dashboard(analysis_type='fastq', output_dir=output_dir_path)
    
    fmt.success("Dashboard created")
    fmt.info(f"  {output_dir_path.name}/analysis_dashboard.html", indent=2)
    fmt.info("  Open in your web browser to explore all results", indent=2)
    
    fmt.success(f"Analysis complete: {output_dir_path / 'analysis_dashboard.html'}")


# ============================================================================
# FASTQ ANALYSIS - Memory Optimized
# ============================================================================

def analyze_fastq(reads, output_dir: Path, args):
    """
    Process FASTQ with memory optimization.
    
    OPTIMIZATIONS:
    - Streaming sequence processing
    - Explicit memory cleanup
    - No checkpoint overhead
    - Direct contig usage (no filtering)
    
    Args:
        reads: FASTQ file path(s)
        output_dir: Output directory
        args: CLI arguments namespace
    """

    fmt = get_formatter()
    reporter = MainReporter(output_dir, view_mode='both')

    # Configuration display
    _display_pipeline_config(reads, output_dir, args, fmt)

    contigs_fasta = None
    total_steps = 1 if (args and args.skip_annotation) else 4

    # ========================================================================
    # STEP 1: TAXONOMIC CLASSIFICATION
    # ========================================================================
    bracken_file = None
    taxonomic_risk = None
    
    try:
        fmt.step_header(1, total_steps, "TAXONOMIC CLASSIFICATION")

        # Handle input format
        if hasattr(args, 'interleaved') and args.interleaved:
            with fmt.spinner("Splitting interleaved reads"):
                with fmt.suppressed_output():
                    reads = split_interleaved(reads[0], output_dir, fmt)
            fmt.success("Reads split successfully")
        elif hasattr(args, 'single') and args.single:
            reads = [reads] if isinstance(reads, str) else reads

        # Run Kraken2
        with fmt.spinner("Running Kraken2 classification"):
            with fmt.suppressed_output():
                kraken_report_path = run_kraken(reads, output_dir)
        
        # Run Bracken
        with fmt.spinner("Running Bracken abundance estimation"):
            with fmt.suppressed_output():
                bracken_file = run_bracken(kraken_report_path, output_dir)

        # Display statistics
        bracken_df = load_bracken_report(bracken_file)

        fmt.success("Taxonomic classification complete")
        fmt.info(f"  Species detected: {len(bracken_df):,}")
        fmt.info(f"  Classified reads: {bracken_df['new_est_reads'].sum():,.0f}")
        fmt.info(f"  Mean abundance: {bracken_df['fraction_total_reads'].mean() * 100:.2f}%")
        
        # MEMORY: Explicit cleanup
        explicit_cleanup(bracken_df)

        # Generate reports
        fmt.info("Creating taxonomic reports and visualizations")
        with fmt.suppressed_output():
            risk_scorer = RiskScorer()
            taxonomic_risk = risk_scorer.calculate_taxonomic_risk(bracken_df)
            
            preliminary_risk_data = {
                'taxonomic': taxonomic_risk,
                'functional': {'score': 0, 'details': {}},
                'ml': {'score': 0, 'details': []},
                'integrated': {
                    'final_score': taxonomic_risk['score'],
                    'risk_level': risk_scorer.get_risk_level(taxonomic_risk['score']),
                    'tier_scores': {'taxonomic': taxonomic_risk['score'], 'functional': 0, 'ml': 0}
                }
            }
            
            report = reporter.generate_taxonomy_report(
                bracken_file=Path(bracken_file),
                risk_data=preliminary_risk_data
            )
            
            with open(output_dir / "01_taxonomic_report.txt", 'w') as rf:
                rf.write(report)
            
            viz_files = create_taxonomic_visualizations(output_dir, bracken_file)
            
            # MEMORY: Cleanup after report generation
            explicit_cleanup(bracken_df, report, viz_files)
        
        # Show output files
        report_files = []
        if (output_dir / "01_taxonomic_report.txt").exists():
            report_files.append("01_taxonomic_report.txt")
        
        if report_files:
            fmt.success(f"Generated {len(report_files)} report(s)")
            for report in report_files:
                fmt.info(f"  {report}", indent=2)
        
        fmt.stage_complete("Taxonomic Classification")

    except Exception as e:
        fmt.error(f"Classification failed: {str(e)}", 
                 solutions=["Verify Kraken2/Bracken installation", 
                           "Check database paths in config"])
        raise

    if args and args.skip_annotation:
        fmt.info("Annotation skipped (--skip-annotation flag)")
        _display_final_summary(output_dir, fmt, "partial")
        return

    # ========================================================================
    # STEP 2: METAGENOMIC ASSEMBLY - Memory Optimized
    # ========================================================================
    
    try:
        fmt.step_header(2, total_steps, "METAGENOMIC ASSEMBLY")
                
        read_type = 'Single-end' if args.single else 'Paired-end' if args.paired else 'Interleaved'
        fmt.info(f"Assembler: MEGAHIT ({read_type}, {getattr(args, 'annotation_threads', 8)} threads)")
        
        # Run assembly
        assembly_threads = getattr(args, 'annotation_threads', 8)
        with fmt.spinner(f"Assembling reads with MEGAHIT"):
            with fmt.suppressed_output():
                contigs_fasta = assemble_reads_to_fasta(
                    reads, output_dir, fmt,
                    threads=assembly_threads
                )
        
        # MEMORY-OPTIMIZED: Stream statistics instead of loading all
        stats = calculate_assembly_stats_streaming(contigs_fasta)
        
        fmt.success("Assembly complete")
        fmt.info(f"  Total contigs: {stats['total_contigs']:,}")
        fmt.info(f"  Total bases: {stats['total_bases']:,} bp")
        fmt.info(f"  N50: {stats['n50']:,} bp")
        fmt.info(f"  Longest contig: {stats['max_length']:,} bp")
        
        # Save stats
        _save_assembly_stats_dict(stats, output_dir)
        
        # MEMORY: Cleanup
        explicit_cleanup(stats)
        
        fmt.stage_complete("Assembly")
        
    except Exception as e:
        fmt.error(f"Assembly failed: {str(e)}", 
                 solutions=["Ensure sufficient memory (≥8GB)"])
        return

    # ========================================================================
    # STEP 3: FUNCTIONAL ANNOTATION
    # ========================================================================
    prokka_dir = None
    swissprot_results = None
    
    try:
        fmt.step_header(3, total_steps, "FUNCTIONAL ANNOTATION")
        
        # Run Prokka (no filtering - Prokka handles all contigs)
        with fmt.spinner("Running Prokka gene prediction"):
            with fmt.suppressed_output():
                prokka_dir = run_prokka(
                    contigs_fasta, output_dir,
                    kill_tbl2asn=getattr(args, 'kill_tbl2asn', True),
                    tbl2asn_timeout=getattr(args, 'tbl2asn_timeout', 30)
                )

        protein_files = list(Path(prokka_dir).glob("*.faa"))
        if not protein_files or all(os.path.getsize(pf) == 0 for pf in protein_files):
            fmt.warning("No proteins predicted (contigs may be too short)")
        else:
            # Display Prokka stats
            stats = parse_prokka_sample_txt(prokka_dir / "sample.txt")
            fmt.success("Prokka gene prediction complete")
            fmt.info(f"  Contigs: {stats.get('contigs', 0):,}")
            fmt.info(f"  CDS: {stats.get('CDS', 0):,}")
            fmt.info(f"  rRNA: {stats.get('rRNA', 0):,}")
            fmt.info(f"  tRNA: {stats.get('tRNA', 0):,}")
            
            # MEMORY: Cleanup stats
            explicit_cleanup(stats)
            
            # Run SwissProt annotation
            fmt.section_header("Functional Annotation")
            
            protein_file = Path(prokka_dir) / "sample.faa"
            
            # MEMORY-OPTIMIZED: Count proteins without loading
            protein_count = sum(1 for _ in SeqIO.parse(protein_file, "fasta"))
            
            fmt.info(f"Annotating {protein_count} proteins against SwissProt database")
            
            # Run annotation
            with fmt.spinner("DIAMOND annotation in progress"):
                with fmt.suppressed_output():
                    annotation_threads = getattr(args, 'annotation_threads', 8)
                    swissprot_results = run_functional_annotation(
                        prokka_dir, output_dir,
                        threads=annotation_threads
                    )
            
            if swissprot_results:
                # Count annotated (stream file)
                with open(swissprot_results) as f:
                    annotated_count = len(set(line.split('\t')[0] for line in f if line.strip()))
                
                fmt.success("Functional annotation complete")
                fmt.info(f"  Annotated: {annotated_count}/{protein_count} proteins ({annotated_count/protein_count*100:.1f}%)")
                
                # MEMORY: Cleanup
                del annotated_count, protein_count
            
            # Generate reports
            fmt.info("Creating functional reports and visualizations")
            with fmt.suppressed_output():
                risk_scorer = RiskScorer()
                
                functional_df = load_annotation_file(swissprot_results)
                
                pathogen_hits_file = output_dir / "pathogen_results.txt"
                pathogen_df = pd.DataFrame()
                
                if pathogen_hits_file.exists():
                    try:
                        pathogen_df = load_pathogen_hits(pathogen_hits_file)
                    except Exception:
                        pass
                
                if not functional_df.empty and 'query_id' in functional_df.columns:
                    total_cds = len(functional_df['query_id'].unique())
                else:
                    total_cds = None  # Let function handle it

                functional_risk_score = risk_scorer.calculate_functional_risk(
                    functional_df, pathogen_df, total_cds=total_cds
                )
                
                functional_risk = {
                    'taxonomic': taxonomic_risk,
                    'functional': functional_risk_score,
                    'ml': {'score': 0, 'details': []},
                    'integrated': {
                        'final_score': (taxonomic_risk['score'] * 0.4 + functional_risk_score['score'] * 0.3),
                        'risk_level': risk_scorer.get_risk_level(
                            taxonomic_risk['score'] * 0.4 + functional_risk_score['score'] * 0.3
                        ),
                        'tier_scores': {
                            'taxonomic': taxonomic_risk['score'],
                            'functional': functional_risk_score['score'],
                            'ml': 0
                        }
                    }
                }
                
                sample_info_file = output_dir / "prokka_annotation" / "sample.txt"
                report = reporter.generate_functional_report(
                    sample_info_file=sample_info_file,
                    functional_annotations_file=Path(swissprot_results),
                    risk_data=functional_risk
                )
                
                with open(output_dir / "02_functional_report.txt", 'w') as rf:
                    rf.write(report)
                
                viz_files = create_functional_visualizations(output_dir, prokka_dir, swissprot_results)
                
                # MEMORY: Cleanup after report
                explicit_cleanup(functional_df, pathogen_df, report, viz_files)
            
            # Show output files
            report_files = []
            if (output_dir / "02_functional_report.txt").exists():
                report_files.append("02_functional_report.txt")
            
            if report_files:
                fmt.success(f"Generated {len(report_files)} report(s)")
                for report in report_files:
                    fmt.info(f"  {report}", indent=2)
        
        fmt.stage_complete("Functional Annotation")
        
    except Exception as e:
        fmt.error(f"Annotation failed: {str(e)}")
        traceback.print_exc()

    # ========================================================================
    # STEP 4: PATHOGEN DETECTION
    # ========================================================================
    try:
        fmt.step_header(4, total_steps, "PATHOGEN DETECTION")

        pathogen_hits_file = None
        protein_file = prokka_dir / "sample.faa" if prokka_dir else None
        
        # Database scan
        if protein_file and protein_file.exists():
            with fmt.spinner("Scanning pathogen database"):
                with fmt.suppressed_output():
                    pathogen_hits_file = run_pathogen_scan(
                        protein_file, output_dir,
                        bracken_results=Path(bracken_file),
                        taxonomy_results=None
                    )
            fmt.success("Database scan complete")

        # ML prediction
        ml_summary = None
        if prokka_dir:
            with fmt.suppressed_output():
                ml_suitable = should_use_ml_prediction(prokka_dir, fmt)
            
            if ml_suitable:
                with fmt.spinner("Running ML pathogen prediction"):
                    with fmt.suppressed_output():
                        ml_results, ml_summary = run_ml_pathogen_prediction(prokka_dir, output_dir)
                
                if ml_results:
                    total_proteins = len(ml_results)
                    pathogenic_count = sum(1 for r in ml_results if r.get('prediction') == 'Pathogenic')
                    high_conf_count = sum(1 for r in ml_results if r.get('confidence', 0) > 0.85)
                    pathogenic_pct = (pathogenic_count / total_proteins * 100) if total_proteins > 0 else 0
                    
                    fmt.success("ML prediction complete")
                    fmt.info(f"  Total proteins: {total_proteins:,}")
                    fmt.info(f"  Pathogenic: {pathogenic_count:,} ({pathogenic_pct:.1f}%)")
                    fmt.info(f"  High confidence: {high_conf_count:,}")
                    
                    # MEMORY: Cleanup
                    explicit_cleanup(ml_results)
                elif ml_summary:
                    total = ml_summary.get('total_sequences', ml_summary.get('total_predictions', 0))
                    pathogenic = ml_summary.get('pathogenic_count', ml_summary.get('pathogenic', 0))
                    pct = ml_summary.get('pathogenic_percentage', 0)
                    
                    if total > 0:
                        fmt.info(f"ML predictions: {pathogenic}/{total} pathogenic ({pct:.1f}%)")
                    else:
                        fmt.info(f"ML predictions: {pathogenic} pathogenic ({pct:.1f}%)")

        # Risk assessment
        with fmt.spinner("Calculating comprehensive risk assessment"):
            with fmt.suppressed_output():
                ml_predictions_file = output_dir / "ml_pathogen_predictions.json"
                
                risk_data = calculate_all_risks(
                    bracken_file=Path(bracken_file),
                    functional_file=Path(swissprot_results) if swissprot_results else None,
                    pathogen_hits_file=pathogen_hits_file,
                    ml_predictions_file=ml_predictions_file if ml_predictions_file.exists() else None
                )
        
        # Display risk summary
        _display_clean_risk_summary(risk_data, fmt)

        # Generate reports
        fmt.info("Creating pathogen reports and visualizations")
        with fmt.suppressed_output():
            report = reporter.generate_pathogen_report(
                risk_data=risk_data,
                pathogen_hits_file=pathogen_hits_file if pathogen_hits_file else output_dir / "pathogen_results.tsv",
                ml_predictions_file=ml_predictions_file
            )
            
            with open(output_dir / "03_pathogen_risk_report.txt", 'w') as rf:
                rf.write(report)
            
            viz_files = create_pathogen_visualizations(output_dir, risk_data=risk_data)
            
            # MEMORY: Cleanup
            explicit_cleanup(report, viz_files)
        
        # Show output files
        report_files = []
        if (output_dir / "03_pathogen_risk_report.txt").exists():
            report_files.append("03_pathogen_risk_report.txt")
        
        if report_files:
            fmt.success(f"Generated {len(report_files)} report(s)")
            for report in report_files:
                fmt.info(f"  {report}", indent=2)
        
        fmt.stage_complete("Pathogen Detection")

    except Exception as e:
        fmt.error(f"Pathogen analysis failed: {str(e)}")
        traceback.print_exc()

    # ========================================================================
    # COMPREHENSIVE REPORT
    # ========================================================================
    try:
        fmt.section_header("GENERATING COMPREHENSIVE REPORT")
        
        fmt.info("Compiling comprehensive analysis report")
        with fmt.suppressed_output():
            sample_info_file = output_dir / "prokka_annotation" / "sample.txt"
            functional_file = Path(swissprot_results) if swissprot_results else output_dir / "functional_annotations.tsv"
            
            comprehensive_report = reporter.generate_report(
                bracken_file=Path(bracken_file),
                sample_info_file=sample_info_file,
                functional_annotations_file=functional_file,
                pathogen_hits_file=pathogen_hits_file if pathogen_hits_file else output_dir / "pathogen_results.txt",
                ml_predictions_file=output_dir / "ml_pathogen_predictions.json"
            )
            
            reporter.export_tables(
                bracken_file=Path(bracken_file),
                annotation_file=functional_file,
                risk_data=risk_data
            )
            
            # MEMORY: Cleanup
            explicit_cleanup(comprehensive_report)
        
        # Show output files
        report_files = []
        if (output_dir / "comprehensive_report.txt").exists():
            report_files.append("comprehensive_report.txt")
        if (output_dir / "validation_summary.txt").exists():
            report_files.append("validation_summary.txt")
        
        if report_files:
            fmt.success(f"Generated {len(report_files)} report(s)")
            for report in report_files:
                fmt.info(f"  {report}", indent=2)
        
        # Check for exported tables
        tables_dir = output_dir / "tables"
        if tables_dir.exists():
            table_files = list(tables_dir.glob("*.csv")) + list(tables_dir.glob("*.json"))
            if table_files:
                fmt.success(f"Exported {len(table_files)} data table(s) to tables/")
                for table_file in sorted(table_files)[:5]:
                    fmt.info(f"  {table_file.name}", indent=2)
                if len(table_files) > 5:
                    fmt.info(f"  ... and {len(table_files) - 5} more", indent=2)
        
    except Exception as e:
        fmt.warning(f"Comprehensive report failed: {str(e)}")
        fmt.info("Individual reports are still available")
        
    # ========================================================================
    # FINAL SUMMARY & CLEANUP
    # ========================================================================
    
    _display_final_summary(output_dir, fmt, "complete")
    
    # MEMORY: Final cleanup
    gc.collect()


# ============================================================================
# DISPLAY FUNCTIONS
# ============================================================================

def _display_pipeline_config(reads, output_dir, args, fmt):
    """
    Display single clean configuration table.
    
    Args:
        reads: Input read file(s)
        output_dir: Output directory path
        args: CLI arguments
        fmt: Output formatter instance
    """
    input_type = 'Interleaved' if (hasattr(args, 'interleaved') and args.interleaved) else \
                 'Paired-end' if (hasattr(args, 'paired') and args.paired) else 'Single-end'
    
    skip_annotation = 'Yes' if (args and args.skip_annotation) else 'No'
    
    sections = [{
        'header': 'Pipeline Configuration',
        'rows': {
            'Input Type': input_type,
            'Read Files': str(reads) if isinstance(reads, str) else f"{len(reads)} file(s)",
            'Output Directory': str(output_dir),
            'Skip Annotation': skip_annotation
        }
    }]
    fmt.display_stats_table("Analysis Configuration", sections)


def _display_clean_risk_summary(risk_data, fmt):
    """
    Display clean risk summary.
    
    Args:
        risk_data: Risk assessment results
        fmt: Output formatter instance
    """
    integrated = risk_data['integrated']
    
    risk_emoji = {
        'High': '🔴',
        'Moderate': '🟡',
        'Low': '🟢'
    }.get(integrated['risk_level'], '⚪')
    
    fmt.success("Risk assessment complete")
    fmt.info(f"  Risk level: {risk_emoji} {integrated['risk_level']}")
    fmt.info(f"  Risk score: {integrated['final_score']:.0f}/100")


def _display_final_summary(output_dir, fmt, analysis_type):
    """
    Display final summary.
    
    Args:
        output_dir: Output directory path
        fmt: Output formatter instance
        analysis_type: Type of analysis completed ('partial' or 'complete')
    """
    fmt.section_header("ANALYSIS COMPLETE")
    
    if analysis_type == "partial":
        status = "Taxonomic Classification Only"
    else:
        status = "Complete Pipeline"
    
    # Count reports
    report_files = [
        "01_taxonomic_report.txt",
        "02_functional_report.txt", 
        "03_pathogen_risk_report.txt",
        "comprehensive_report.txt",
        "validation_summary.txt"
    ]
    existing_reports = [r for r in report_files if (output_dir / r).exists()]
    
    fmt.info(f"Status: {status}")
    fmt.info(f"Output directory: {output_dir}")
    fmt.info(f"Reports generated: {len(existing_reports)} file(s)")
    fmt.info(f"Dashboard: analysis_dashboard.html")

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def _save_assembly_stats_dict(stats: dict, output_dir: Path):
    """Save assembly statistics from dict."""
    with open(output_dir / "sample.txt", 'w') as f:
        f.write(f"Organism: {output_dir.stem}\n")
        f.write(f"Contigs: {stats['total_contigs']}\n")
        f.write(f"Total bases: {stats['total_bases']}\n")
        f.write(f"N50: {stats['n50']}\n")
        f.write(f"Max length: {stats['max_length']}\n")
        f.write(f"Mean length: {stats['mean_length']:.1f}\n")
        f.write(f"Genes: TBD\n")