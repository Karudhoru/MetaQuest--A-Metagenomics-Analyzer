#!/usr/bin/env python3
"""
MetaQuest Visualization Module v5.0.0
Main Visualizer Functions
"""
from pathlib import Path
import pandas as pd
from .taxonomic_visualizer import TaxonomicVisualizer  
from .pathogenic_visualizer import PathogenVisualizer
from .functional_visualizer import FunctionalVisualizer
from ..io.output_formatter import get_formatter

# ============================================================================
# Main visualization functions for pipeline integration
# ============================================================================

def create_taxonomic_visualizations(output_dir, data, **kwargs):
    """
    Create all taxonomic visualizations for v5.0.0 reports.
    
    Args:
        output_dir: Output directory for visualizations
        data: Bracken/Kraken file path or BLAST results list
    
    Returns:
        Dictionary of generated file paths
    """
    # --- ADD FORMATTER ---
    formatter = get_formatter()
    visualizer = TaxonomicVisualizer(Path(output_dir))
    generated_files = {}
    
    # --- USE FORMATTER ---
    formatter.section_header("Creating Taxonomic Visualizations")
    
    # Abundance chart
    abundance_chart = visualizer.create_abundance_chart(data)
    if abundance_chart:
        generated_files['abundance_chart'] = abundance_chart
    
    # Diversity metrics (only for Bracken/Kraken)
    diversity_chart = visualizer.create_diversity_summary(data)
    if diversity_chart:
        generated_files['diversity_metrics'] = diversity_chart
    
    # Krona plot
    krona_plot = visualizer.create_krona_plot(data)
    if krona_plot:
        generated_files['krona_plot'] = krona_plot
    
    formatter.success(f"Generated {len(generated_files)} taxonomic visualizations")
    return generated_files


def create_pathogen_visualizations(output_dir, traditional_data=None, **kwargs):
    """
    Create all pathogen visualizations for v5.0.0 reports.
    
    Args:
        output_dir: Output directory for visualizations
        traditional_data: Path to pathogen detection JSON report
    
    Returns:
        Dictionary of generated file paths
    """
    formatter = get_formatter()
    visualizer = PathogenVisualizer(Path(output_dir))
    generated_files = {}
    
    # --- USE FORMATTER ---
    formatter.section_header("Creating Pathogen Visualizations")
    
    if traditional_data:
        # Risk assessment chart
        risk_chart = visualizer.create_risk_assessment_chart(traditional_data)
        if risk_chart:
            generated_files['risk_assessment'] = risk_chart
        
        # WHO priority distribution
        who_chart = visualizer.create_who_priority_distribution(traditional_data)
        if who_chart:
            generated_files['who_priority'] = who_chart
        
        # Detection confidence
        confidence_chart = visualizer.create_detection_confidence_chart(traditional_data)
        if confidence_chart:
            generated_files['detection_confidence'] = confidence_chart
        
    formatter.success(f"Generated {len(generated_files)} pathogen visualizations")    
    return generated_files


def create_functional_visualizations(output_dir, prokka_results=None, 
                                    swissprot_results=None, **kwargs):
    """
    Create all functional visualizations for v5.0.0 reports.
    
    Args:
        output_dir: Output directory for visualizations
        prokka_results: Path to Prokka output directory
        swissprot_results: Path to SwissProt annotation file
    
    Returns:
        Dictionary of generated file paths
    """
    # --- ADD FORMATTER ---
    formatter = get_formatter()
    visualizer = FunctionalVisualizer(Path(output_dir))
    generated_files = {}
    
    # --- USE FORMATTER ---
    formatter.section_header("Creating Functional Visualizations")
    
    # Annotation quality dashboard
    if swissprot_results:
        quality_dashboard = visualizer.create_annotation_quality_dashboard(swissprot_results)
        if quality_dashboard:
            generated_files['annotation_quality'] = quality_dashboard
    
    # Protein length analysis
    if prokka_results:
        length_analysis = visualizer.create_protein_length_analysis(prokka_results)
        if length_analysis:
            generated_files['protein_length'] = length_analysis
        
        # Functional categories
        category_chart = visualizer.create_functional_category_chart(prokka_results)
        if category_chart:
            generated_files['functional_categories'] = category_chart
    
    formatter.success(f"Generated {len(generated_files)} functional visualizations")
    return generated_files

# ============================================================================
# Utility functions
# ============================================================================

def validate_visualization_data(data, data_type='taxonomic'):
    """
    Validate input data before visualization.
    
    Args:
        data: Input data (file path, list, or DataFrame)
        data_type: Type of data ('taxonomic', 'pathogen', 'functional')
    
    Returns:
        Tuple of (is_valid: bool, message: str)
    """
    if data is None:
        return False, "No data provided"
    
    if isinstance(data, (str, Path)):
        path = Path(data)
        if not path.exists():
            return False, f"File not found: {path}"
        if path.stat().st_size == 0:
            return False, f"File is empty: {path}"
        return True, "Valid file path"
    
    elif isinstance(data, list):
        if len(data) == 0:
            return False, "Empty data list"
        return True, "Valid data list"
    
    elif isinstance(data, pd.DataFrame):
        if data.empty:
            return False, "Empty DataFrame"
        return True, "Valid DataFrame"
    
    else:
        return False, f"Unsupported data type: {type(data)}"


def generate_visualization_summary(output_dir):
    """
    Generate a summary report of all visualizations in output directory.
    
    Args:
        output_dir: Directory containing visualizations
    
    Returns:
        Dictionary with visualization summary
    """
    output_path = Path(output_dir)
    html_files = list(output_path.glob("*.html"))
    
    summary = {
        'total_visualizations': len(html_files),
        'taxonomic': [],
        'pathogen': [],
        'functional': [],
        'other': []
    }
    
    # Categorize visualizations
    for html_file in html_files:
        name = html_file.name
        if 'taxonomic' in name or 'abundance' in name or 'diversity' in name or 'krona' in name:
            summary['taxonomic'].append(name)
        elif 'pathogen' in name or 'risk' in name or 'who' in name or 'confidence' in name:
            summary['pathogen'].append(name)
        elif 'functional' in name or 'annotation' in name or 'protein' in name or 'category' in name:
            summary['functional'].append(name)
        else:
            summary['other'].append(name)
    
    return summary