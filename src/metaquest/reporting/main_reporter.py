#!/usr/bin/env python3
"""
MetaQuest Professional Main Reporting Module.
"""
from pathlib import Path
from .taxonomic_reporter import TaxonomicReporter
from .pathogen_reporter import PathogenReporter, SequenceHitsReporter
from .functional_reporter import FunctionalReporter


# Main Function to generate taxonomic report
def generate_taxonomic_report(output_dir, bracken_data=None, blast_data=None):
    """Generate professional taxonomic classification report"""
    reporter = TaxonomicReporter(Path(output_dir))
    return reporter.generate_report(bracken_data=bracken_data, blast_data=blast_data)

# Main functions for pathogen reports
def generate_fastq_pathogen_report(output_dir, bracken_pathogens=None, taxonomy_pathogens=None, sequence_pathogens=None):
    """
    Main entry point for generating the FASTQ pathogen report.
    Generates: pathogen_summary.txt + pathogen_detection_report.json
    
    Args:
        output_dir (str/Path): Output directory path
        bracken_pathogens (list): Pathogenic organisms from Bracken analysis
        taxonomy_pathogens (list): Pathogenic organisms from BLAST taxonomy
        sequence_pathogens (list): Pathogenic organisms from sequence database
    
    Returns:
        dict: Generated file paths and summary statistics
    """
    reporter = PathogenReporter(Path(output_dir), analysis_mode='fastq')
    return reporter.generate_fastq_report(
        bracken_pathogens=bracken_pathogens,
        taxonomy_pathogens=taxonomy_pathogens,
        sequence_pathogens=sequence_pathogens
    )

def generate_fasta_ml_pathogen_report(output_dir, blast_taxonomy_pathogens=None, ml_results=None, ml_summary=None):
    """
    Main entry point for generating the integrated FASTA+ML pathogen report.
    Generates: blast_ml_pathogen_summary.txt + blast_ml_integrated_pathogen_report.json
    
    Args:
        output_dir (str/Path): Output directory path
        blast_taxonomy_pathogens (list): Pathogenic organisms from BLAST taxonomy
        ml_results (list): Individual ML predictions for proteins
        ml_summary (dict): Summary statistics from ML analysis
    
    Returns:
        dict: Generated file paths and summary statistics
    """
    reporter = PathogenReporter(Path(output_dir), analysis_mode='fasta')
    return reporter.generate_fasta_ml_report(
        blast_taxonomy_pathogens=blast_taxonomy_pathogens,
        ml_results=ml_results,
        ml_summary=ml_summary
    )

# Main function for Functional Report
def generate_functional_report(output_dir, prokka_results=None, swissprot_results=None):
    """
    Main entry point for generating comprehensive functional annotation report.
    
    Args:
        output_dir (str/Path): Output directory path
        prokka_results (str/Path): Path to Prokka results directory
        swissprot_results (str/Path): Path to SwissProt annotation TSV file
    
    Returns:
        dict: Generated file paths and analysis summary
    """
    reporter = FunctionalReporter(Path(output_dir))
    return reporter.generate_report(
        prokka_results=prokka_results,
        swissprot_results=swissprot_results
    )

def generate_sequence_hits_report(output_dir: Path, analysis_type: str, hits_file_path: str, title: str):
    """Wrapper to generate a summary report from a DIAMOND/BLAST TSV file."""
    if hits_file_path:
        reporter = SequenceHitsReporter(output_dir, analysis_type)
        reporter.generate_report(hits_file_path, title)