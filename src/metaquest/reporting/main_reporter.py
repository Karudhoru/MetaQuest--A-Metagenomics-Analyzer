#!/usr/bin/env python3
"""
MetaQuest Professional Main Reporting Module
Author: MetaQuest Metagenomics Team

Main entry points for all MetaQuest reporting functionality.
Provides unified interface for generating taxonomic, functional, and pathogen detection reports.
"""

from pathlib import Path
from typing import Dict, Optional, List, Any
from .taxonomic_reporter import TaxonomicReporter
from .pathogen_reporter import PathogenReporter, SequenceHitsReporter
from .functional_reporter import FunctionalReporter


def generate_taxonomic_report(output_dir: str, bracken_data: Optional[Any] = None,
                              blast_data: Optional[Any] = None) -> Dict[str, str]:
    """
    Generate comprehensive taxonomic classification report.
    
    Supports both FASTQ (Kraken2/Bracken) and FASTA (BLAST) pipelines,
    or integrated dual-method analysis.
    
    Args:
        output_dir: Output directory for reports
        bracken_data: Bracken/Kraken2 output file path or DataFrame
        blast_data: BLAST results JSON file path or data structure
        
    Returns:
        Dictionary containing paths to generated reports and summary statistics
        
    Example:
        >>> results = generate_taxonomic_report(
        ...     output_dir="./reports",
        ...     bracken_data="bracken_output.txt",
        ...     blast_data="blast_results.json"
        ... )
        >>> print(f"Report saved to: {results['text_report']}")
    """
    reporter = TaxonomicReporter(Path(output_dir))
    return reporter.generate_report(bracken_data=bracken_data, blast_data=blast_data)


def generate_fastq_pathogen_report(output_dir: str, 
                                   bracken_pathogens: Optional[List] = None,
                                   taxonomy_pathogens: Optional[List] = None,
                                   sequence_pathogens: Optional[List] = None) -> Dict[str, str]:
    """
    Generate pathogen detection report for FASTQ pipeline.
    
    Integrates pathogen detections from multiple sources:
    - Kraken2/Bracken taxonomic profiling
    - BLAST taxonomy assignments
    - Pathogen-specific sequence databases
    
    Provides WHO priority pathogen classification and clinical risk assessment.
    
    Args:
        output_dir: Output directory for reports
        bracken_pathogens: List of pathogenic organisms from Bracken analysis
        taxonomy_pathogens: List of pathogenic organisms from BLAST taxonomy
        sequence_pathogens: List of pathogenic organisms from sequence database
        
    Returns:
        Dictionary containing:
        - text_report: Path to text report
        - json_report: Path to JSON report
        - risk_level: Overall clinical risk assessment (CRITICAL/HIGH/MEDIUM/LOW)
        - pathogen_count: Number of detected pathogens
        - critical_pathogens: Count of WHO critical priority pathogens
        
    Example:
        >>> bracken_pathogens = [
        ...     {'organism': 'Escherichia coli', 'abundance': 0.15, 'reads': 45000},
        ...     {'organism': 'Staphylococcus aureus', 'abundance': 0.08, 'reads': 24000}
        ... ]
        >>> results = generate_fastq_pathogen_report(
        ...     output_dir="./reports",
        ...     bracken_pathogens=bracken_pathogens
        ... )
        >>> print(f"Risk Level: {results['risk_level']}")
    """
    reporter = PathogenReporter(Path(output_dir), analysis_mode='fastq')
    return reporter.generate_fastq_report(
        bracken_pathogens=bracken_pathogens,
        taxonomy_pathogens=taxonomy_pathogens,
        sequence_pathogens=sequence_pathogens
    )


def generate_fasta_ml_pathogen_report(output_dir: str,
                                      blast_taxonomy_pathogens: Optional[List] = None,
                                      ml_results: Optional[List] = None,
                                      ml_summary: Optional[Dict] = None) -> Dict[str, str]:
    """
    Generate integrated BLAST+ML pathogen detection report for FASTA pipeline.
    
    Combines two complementary approaches:
    1. BLAST-based taxonomic pathogen identification
    2. Machine learning pathogenicity prediction
    
    Provides integrated risk assessment and clinical recommendations.
    
    Args:
        output_dir: Output directory for reports
        blast_taxonomy_pathogens: List of pathogenic organisms from BLAST taxonomy
        ml_results: Individual ML predictions for proteins (optional, for detailed analysis)
        ml_summary: Summary statistics from ML analysis containing:
            - total_sequences_analyzed: Total number of proteins analyzed
            - pathogenic_predictions: Number of pathogenic predictions
            - high_confidence_predictions: Number of high-confidence predictions
            - average_confidence: Average prediction confidence
            
    Returns:
        Dictionary containing:
        - text_report: Path to integrated report
        - json_report: Path to JSON report
        - risk_level: Overall risk assessment
        - blast_pathogens: Number of BLAST-detected pathogens
        - ml_pathogenic_proteins: Number of ML-predicted pathogenic proteins
        
    Example:
        >>> ml_summary = {
        ...     'total_sequences_analyzed': 3500,
        ...     'pathogenic_predictions': 420,
        ...     'high_confidence_predictions': 315,
        ...     'average_confidence': 0.87
        ... }
        >>> results = generate_fasta_ml_pathogen_report(
        ...     output_dir="./reports",
        ...     ml_summary=ml_summary
        ... )
        >>> print(f"ML Pathogenic Ratio: {ml_summary['pathogenic_predictions']/ml_summary['total_sequences_analyzed']:.1%}")
    """
    reporter = PathogenReporter(Path(output_dir), analysis_mode='fasta')
    return reporter.generate_fasta_ml_report(
        blast_taxonomy_pathogens=blast_taxonomy_pathogens,
        ml_results=ml_results,
        ml_summary=ml_summary
    )


def generate_functional_report(output_dir: str,
                               prokka_results: Optional[str] = None,
                               swissprot_results: Optional[str] = None) -> Dict[str, str]:
    """
    Generate comprehensive functional annotation report.
    
    Analyzes:
    - Gene prediction results (Prokka)
    - Protein functional annotations (SwissProt/UniProt)
    - Functional category distribution
    - Annotation quality metrics
    
    Args:
        output_dir: Output directory for reports
        prokka_results: Path to Prokka output directory
        swissprot_results: Path to SwissProt BLAST annotation TSV file
        
    Returns:
        Dictionary containing:
        - text_report: Path to text report
        - json_report: Path to JSON report
        - total_genes: Number of predicted genes
        - annotation_rate: Percentage of genes with functional annotations
        - functional_categories: Number of functional categories detected
        
    Example:
        >>> results = generate_functional_report(
        ...     output_dir="./reports",
        ...     prokka_results="./prokka_output",
        ...     swissprot_results="./swissprot_annotations.tsv"
        ... )
        >>> print(f"Annotation Rate: {results['annotation_rate']:.1f}%")
    """
    reporter = FunctionalReporter(Path(output_dir))
    return reporter.generate_report(
        prokka_results=prokka_results,
        swissprot_results=swissprot_results
    )


def generate_sequence_hits_report(output_dir: str, 
                                  analysis_type: str,
                                  hits_file_path: str,
                                  title: str = "Sequence Hits Summary") -> Dict[str, str]:
    """
    Generate summary report from DIAMOND/BLAST TSV output.
    
    Creates a concise, human-readable summary of sequence similarity hits,
    useful for quick review of pathogen-specific database searches.
    
    Args:
        output_dir: Output directory for reports
        analysis_type: Type of analysis (e.g., "Pathogen Database Search")
        hits_file_path: Path to DIAMOND/BLAST TSV file
        title: Custom report title
        
    Returns:
        Dictionary containing text_report path
        
    Example:
        >>> results = generate_sequence_hits_report(
        ...     output_dir="./reports",
        ...     analysis_type="VFDB Search",
        ...     hits_file_path="./vfdb_hits.tsv",
        ...     title="Virulence Factor Detection"
        ... )
    """
    reporter = SequenceHitsReporter(Path(output_dir), analysis_type)
    return reporter.generate_report(hits_file_path, title)


# Convenience function for batch report generation
def generate_all_reports(output_dir: str,
                        bracken_data: Optional[Any] = None,
                        blast_data: Optional[Any] = None,
                        prokka_results: Optional[str] = None,
                        swissprot_results: Optional[str] = None,
                        bracken_pathogens: Optional[List] = None,
                        taxonomy_pathogens: Optional[List] = None,
                        sequence_pathogens: Optional[List] = None,
                        ml_summary: Optional[Dict] = None) -> Dict[str, Dict]:
    """
    Generate all available reports in a single call.
    
    Convenience function that generates taxonomic, functional, and pathogen
    detection reports based on available input data.
    
    Args:
        output_dir: Output directory for all reports
        bracken_data: Bracken/Kraken2 data for taxonomic report
        blast_data: BLAST data for taxonomic report
        prokka_results: Prokka directory for functional report
        swissprot_results: SwissProt annotations for functional report
        bracken_pathogens: Bracken pathogens for pathogen report
        taxonomy_pathogens: BLAST taxonomy pathogens for pathogen report
        sequence_pathogens: Sequence database pathogens for pathogen report
        ml_summary: ML predictions for FASTA pathogen report
        
    Returns:
        Dictionary containing results from all generated reports:
        - taxonomic_report: Taxonomic classification results
        - functional_report: Functional annotation results
        - pathogen_report: Pathogen detection results (if applicable)
        
    Example:
        >>> results = generate_all_reports(
        ...     output_dir="./reports",
        ...     bracken_data="bracken_output.txt",
        ...     prokka_results="./prokka_output",
        ...     swissprot_results="./swissprot.tsv"
        ... )
        >>> print(f"Generated {len(results)} report types")
    """
    all_results = {}
    
    # Generate taxonomic report
    if bracken_data or blast_data:
        print("\n" + "="*80)
        print("GENERATING TAXONOMIC CLASSIFICATION REPORT")
        print("="*80)
        all_results['taxonomic_report'] = generate_taxonomic_report(
            output_dir, bracken_data, blast_data
        )
    
    # Generate functional report
    if prokka_results or swissprot_results:
        print("\n" + "="*80)
        print("GENERATING FUNCTIONAL ANNOTATION REPORT")
        print("="*80)
        all_results['functional_report'] = generate_functional_report(
            output_dir, prokka_results, swissprot_results
        )
    
    # Generate pathogen report (FASTQ pipeline)
    if bracken_pathogens or taxonomy_pathogens or sequence_pathogens:
        print("\n" + "="*80)
        print("GENERATING PATHOGEN DETECTION REPORT (FASTQ)")
        print("="*80)
        all_results['pathogen_report_fastq'] = generate_fastq_pathogen_report(
            output_dir, bracken_pathogens, taxonomy_pathogens, sequence_pathogens
        )
    
    # Generate pathogen report (FASTA pipeline with ML)
    if ml_summary:
        print("\n" + "="*80)
        print("GENERATING INTEGRATED PATHOGEN REPORT (FASTA+ML)")
        print("="*80)
        all_results['pathogen_report_fasta'] = generate_fasta_ml_pathogen_report(
            output_dir, taxonomy_pathogens, ml_summary=ml_summary
        )
    
    print("\n" + "="*80)
    print("REPORT GENERATION COMPLETE")
    print("="*80)
    print(f"Total report types generated: {len(all_results)}")
    print(f"Output directory: {output_dir}")
    
    return all_results


# Module-level documentation and version info
__version__ = "4.0.0"
__author__ = "MetaQuest Metagenomics Team"
__all__ = [
    'generate_taxonomic_report',
    'generate_fastq_pathogen_report',
    'generate_fasta_ml_pathogen_report',
    'generate_functional_report',
    'generate_sequence_hits_report',
    'generate_all_reports'
]
