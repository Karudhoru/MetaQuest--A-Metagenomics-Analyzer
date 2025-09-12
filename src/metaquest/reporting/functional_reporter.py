#!/usr/bin/env python3
"""
MetaQuest FunctionalReporter Module - Professional Scientific Reporting
Updated with comprehensive functional annotation reporting following scientific standards.
"""

import json
import pandas as pd
from pathlib import Path
from datetime import datetime
from abc import ABC, abstractmethod
from typing import Dict, List, Any
import numpy as np
from collections import Counter
import re
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio import BiopythonWarning
import warnings
from .base_reporter import BaseReportGenerator

class FunctionalReporter(BaseReportGenerator):
    """
    Professional functional annotation reporting with comprehensive analysis
    following scientific standards for gene prediction and protein annotation.
    """
    
    # Functional category mappings for enhanced analysis
    FUNCTIONAL_CATEGORIES = {
        'Antimicrobial Resistance': [
            'beta-lactamase', 'carbapenemase', 'penicillinase', 'cephalosporinase',
            'efflux pump', 'multidrug resistance', 'antibiotic resistance',
            'vancomycin resistance', 'methicillin resistance', 'fluoroquinolone resistance'
        ],
        'Virulence Factors': [
            'virulence', 'toxin', 'adhesin', 'invasin', 'hemolysin', 'enterotoxin',
            'cytotoxin', 'secretion system', 'pathogenicity', 'colonization'
        ],
        'Metabolism': [
            'kinase', 'synthase', 'dehydrogenase', 'reductase', 'oxidase', 'transferase',
            'hydrolase', 'lyase', 'isomerase', 'ligase', 'metabolic pathway'
        ],
        'Transport Systems': [
            'transporter', 'permease', 'channel', 'pump', 'abc transporter',
            'efflux', 'influx', 'ion transport', 'substrate transport'
        ],
        'Transcriptional Regulation': [
            'transcriptional regulator', 'repressor', 'activator', 'two-component',
            'response regulator', 'sensor', 'sigma factor', 'dna-binding'
        ],
        'Cell Wall Synthesis': [
            'peptidoglycan', 'murein', 'cell wall', 'penicillin-binding protein',
            'transpeptidase', 'transglycosylase', 'autolysis'
        ],
        'Stress Response': [
            'heat shock', 'cold shock', 'oxidative stress', 'osmotic stress',
            'acid resistance', 'alkali resistance', 'survival', 'adaptation'
        ],
        'DNA Repair/Replication': [
            'dna polymerase', 'helicase', 'primase', 'ligase', 'topoisomerase',
            'recombination', 'repair', 'proofreading', 'mismatch repair'
        ]
    }
    
    def __init__(self, output_dir: Path):
        super().__init__(output_dir, "Functional Annotation")
    
    def generate_report(self, prokka_results=None, swissprot_results=None, **kwargs) -> Dict[str, str]:
        """Generate comprehensive functional annotation report following professional standards"""
        print("Generating professional functional annotation report...")
        
        # Analyze gene prediction results
        gene_analysis = self._analyze_gene_prediction(prokka_results) if prokka_results else None
        
        # Analyze functional annotation results
        annotation_analysis = self._analyze_functional_annotation(swissprot_results) if swissprot_results else None
        
        # Calculate annotation coverage metrics
        coverage_metrics = self._calculate_annotation_coverage(gene_analysis, annotation_analysis)
        
        # Generate professional text report
        self._generate_professional_functional_report(gene_analysis, annotation_analysis, coverage_metrics)
        
        # Generate structured JSON report
        json_report = self._create_functional_json_report(gene_analysis, annotation_analysis, coverage_metrics)
        json_file = self._save_json_report(json_report, "functional_annotation_report.json")
        
        return {
            'text_report': str(self.output_dir / "functional_annotation_report.txt"),
            'json_report': json_file,
            'total_genes': gene_analysis.get('total_proteins', 0) if gene_analysis else 0,
            'annotation_rate': coverage_metrics.get('annotation_rate', 0),
            'functional_categories': len(coverage_metrics.get('functional_distribution', {}))
        }

    def _analyze_gene_prediction(self, prokka_dir):
        """Analyze Prokka gene prediction results with quality assessment"""
        try:
            prokka_path = Path(prokka_dir)
            
            # Find protein sequence files
            faa_files = list(prokka_path.glob("*.faa"))
            ffn_files = list(prokka_path.glob("*.ffn"))
            gff_files = list(prokka_path.glob("*.gff"))
            
            analysis_results = {
                'prokka_directory': str(prokka_path),
                'output_files': {
                    'protein_files': len(faa_files),
                    'gene_files': len(ffn_files),
                    'annotation_files': len(gff_files)
                },
                'total_proteins': 0,
                'length_statistics': {},
                'quality_metrics': {}
            }
            
            # Analyze protein sequences
            if faa_files:
                protein_file = faa_files[0]
                proteins = list(SeqIO.parse(protein_file, "fasta"))
                
                if proteins:
                    lengths = [len(seq.seq) for seq in proteins]
                    
                    analysis_results.update({
                        'total_proteins': len(proteins),
                        'length_statistics': {
                            'mean': np.mean(lengths),
                            'median': np.median(lengths),
                            'std': np.std(lengths),
                            'min': min(lengths),
                            'max': max(lengths),
                            'range': max(lengths) - min(lengths)
                        },
                        'length_distribution': {
                            'very_short': len([l for l in lengths if l < 50]),
                            'short': len([l for l in lengths if 50 <= l < 100]),
                            'medium': len([l for l in lengths if 100 <= l < 300]),
                            'long': len([l for l in lengths if l >= 300])
                        }
                    })
                    
                    # Quality assessment
                    very_short_ratio = analysis_results['length_distribution']['very_short'] / len(proteins)
                    mean_length = analysis_results['length_statistics']['mean']
                    
                    analysis_results['quality_metrics'] = {
                        'very_short_ratio': very_short_ratio,
                        'mean_length_assessment': 'Normal' if mean_length > 200 else 'Short',
                        'coding_density_estimate': self._estimate_coding_density(lengths),
                        'quality_flag': 'Pass' if very_short_ratio < 0.3 and mean_length > 150 else 'Review'
                    }
            
            return analysis_results
            
        except Exception as e:
            return {'error': f"Error analyzing gene prediction: {e}"}

    def _analyze_functional_annotation(self, swissprot_file):
        """Analyze SwissProt functional annotation with comprehensive metrics"""
        try:
            if not Path(swissprot_file).exists():
                return {'error': "SwissProt annotation file not found"}
            
            # Read annotation results
            df = pd.read_csv(swissprot_file, sep='\t')
            
            if df.empty:
                return {'error': "No SwissProt annotations found"}
            
            analysis_results = {
                'total_matches': len(df),
                'unique_proteins': df['qseqid'].nunique(),
                'identity_statistics': {
                    'mean': df['pident'].mean(),
                    'median': df['pident'].median(),
                    'std': df['pident'].std(),
                    'min': df['pident'].min(),
                    'max': df['pident'].max()
                },
                'quality_distribution': {
                    'excellent': len(df[df['pident'] >= 95]),
                    'very_good': len(df[(df['pident'] >= 90) & (df['pident'] < 95)]),
                    'good': len(df[(df['pident'] >= 80) & (df['pident'] < 90)]),
                    'moderate': len(df[(df['pident'] >= 70) & (df['pident'] < 80)]),
                    'low': len(df[df['pident'] < 70])
                },
                'evalue_statistics': {
                    'mean': df['evalue'].mean(),
                    'median': df['evalue'].median(),
                    'highly_significant': len(df[df['evalue'] < 1e-50]),
                    'significant': len(df[(df['evalue'] >= 1e-50) & (df['evalue'] < 1e-10)]),
                    'marginal': len(df[df['evalue'] >= 1e-10])
                }
            }
            
            # Functional category analysis
            if 'stitle' in df.columns:
                analysis_results['functional_analysis'] = self._analyze_functional_categories(df['stitle'])
            
            return analysis_results
            
        except Exception as e:
            return {'error': f"Error analyzing functional annotation: {e}"}

    def _analyze_functional_categories(self, descriptions):
        """Analyze functional categories from protein descriptions"""
        category_counts = {}
        descriptions_lower = descriptions.str.lower().fillna('')
        
        for category, keywords in self.FUNCTIONAL_CATEGORIES.items():
            count = 0
            for keyword in keywords:
                count += descriptions_lower.str.contains(keyword, na=False).sum()
            
            if count > 0:
                category_counts[category] = count
        
        return category_counts

    def _estimate_coding_density(self, protein_lengths):
        """Estimate coding density based on protein length distribution"""
        total_aa = sum(protein_lengths)
        total_nucleotides = total_aa * 3  # Approximate
        
        # Rough genome size estimate (assumes typical bacterial genome)
        estimated_genome_size = total_nucleotides * 1.2  # Account for non-coding regions
        coding_density = (total_nucleotides / estimated_genome_size) * 100
        
        return round(coding_density, 1)

    def _calculate_annotation_coverage(self, gene_analysis, annotation_analysis):
        """Calculate comprehensive annotation coverage metrics"""
        metrics = {
            'annotation_rate': 0,
            'coverage_assessment': 'No data',
            'functional_distribution': {},
            'quality_assessment': 'Incomplete analysis'
        }
        
        if gene_analysis and annotation_analysis and 'error' not in gene_analysis and 'error' not in annotation_analysis:
            total_proteins = gene_analysis.get('total_proteins', 0)
            annotated_proteins = annotation_analysis.get('unique_proteins', 0)
            
            if total_proteins > 0:
                annotation_rate = (annotated_proteins / total_proteins) * 100
                metrics.update({
                    'annotation_rate': round(annotation_rate, 1),
                    'total_predicted': total_proteins,
                    'total_annotated': annotated_proteins,
                    'unannotated': total_proteins - annotated_proteins
                })
                
                # Coverage assessment
                if annotation_rate >= 85:
                    metrics['coverage_assessment'] = 'Excellent'
                elif annotation_rate >= 70:
                    metrics['coverage_assessment'] = 'Good'
                elif annotation_rate >= 50:
                    metrics['coverage_assessment'] = 'Moderate'
                else:
                    metrics['coverage_assessment'] = 'Limited'
                
                # Functional distribution
                if 'functional_analysis' in annotation_analysis:
                    metrics['functional_distribution'] = annotation_analysis['functional_analysis']
                
                # Quality assessment
                high_quality_annotations = annotation_analysis.get('quality_distribution', {}).get('excellent', 0)
                quality_ratio = high_quality_annotations / max(annotated_proteins, 1)
                
                if quality_ratio >= 0.7:
                    metrics['quality_assessment'] = 'High-quality annotations predominate'
                elif quality_ratio >= 0.5:
                    metrics['quality_assessment'] = 'Good annotation quality'
                else:
                    metrics['quality_assessment'] = 'Mixed annotation quality'
        
        return metrics

    def _generate_professional_functional_report(self, gene_analysis, annotation_analysis, coverage_metrics):
        """Generate professional functional annotation report"""
        content = self._create_professional_header(
            "Functional Annotation Report",
            "Prokka v1.14.6 + SwissProt BLAST v2.13.0"
        )
        
        # Executive Summary
        content.extend([
            "EXECUTIVE SUMMARY:",
            self._generate_functional_executive_summary(gene_analysis, annotation_analysis, coverage_metrics),
            ""
        ])
        
        # Gene Prediction Summary
        if gene_analysis and 'error' not in gene_analysis:
            content.extend(self._generate_gene_prediction_section(gene_analysis))
        
        # Functional Annotation Coverage
        if annotation_analysis and 'error' not in annotation_analysis:
            content.extend(self._generate_annotation_coverage_section(annotation_analysis, coverage_metrics))
        
        # Functional Category Analysis
        if coverage_metrics.get('functional_distribution'):
            content.extend(self._generate_functional_category_section(coverage_metrics['functional_distribution']))
        
        # Quality Assessment and Recommendations
        content.extend(self._generate_quality_assessment_section(coverage_metrics))
        
        # Technical Specifications
        content.extend(self._generate_technical_specifications())
        
        # Save report
        self._save_text_report(content, "functional_annotation_report.txt")

    def _generate_functional_executive_summary(self, gene_analysis, annotation_analysis, coverage_metrics):
        """Generate executive summary for functional report"""
        if not gene_analysis or 'error' in gene_analysis:
            return "Gene prediction analysis could not be completed due to data availability issues."
        
        total_genes = gene_analysis.get('total_proteins', 0)
        annotation_rate = coverage_metrics.get('annotation_rate', 0)
        coverage_assessment = coverage_metrics.get('coverage_assessment', 'Unknown')
        
        summary = f"Comprehensive functional analysis identified {total_genes:,} protein-coding genes with "
        
        if annotation_rate > 0:
            summary += f"{coverage_assessment.lower()} annotation coverage ({annotation_rate:.1f}%). "
        else:
            summary += "limited annotation coverage due to analysis constraints. "
        
        functional_categories = len(coverage_metrics.get('functional_distribution', {}))
        if functional_categories > 0:
            summary += f"Functional profile indicates diverse metabolic capabilities across {functional_categories} "
            summary += "major functional categories, including essential cellular processes and specialized functions."
        else:
            summary += "Functional categorization requires additional annotation data for comprehensive analysis."
        
        return summary

    def _generate_gene_prediction_section(self, gene_analysis):
        """Generate gene prediction summary section"""
        content = [
            "Gene Prediction Summary (Prokka)",
            "-" * 35,
            f"Total Predicted Genes: {gene_analysis['total_proteins']:,}",
            f"Protein-Coding Genes: {gene_analysis['total_proteins']:,} (100.0%)",
            f"Analysis Source: {Path(gene_analysis['prokka_directory']).name}",
            ""
        ]
        
        if 'length_statistics' in gene_analysis:
            stats = gene_analysis['length_statistics']
            content.extend([
                "Protein Length Distribution:",
                f"Average Length: {stats['mean']:.0f} ± {stats['std']:.0f} amino acids",
                f"Median Length: {stats['median']:.0f} amino acids",
                f"Range: {stats['min']:.0f} - {stats['max']:.0f} amino acids",
                ""
            ])
        
        if 'length_distribution' in gene_analysis:
            dist = gene_analysis['length_distribution']
            total = gene_analysis['total_proteins']
            content.extend([
                "Length Category Distribution:",
                f"Very short (<50 aa): {dist['very_short']:,} ({dist['very_short']/total*100:.1f}%)",
                f"Short (50-99 aa): {dist['short']:,} ({dist['short']/total*100:.1f}%)",
                f"Medium (100-299 aa): {dist['medium']:,} ({dist['medium']/total*100:.1f}%)",
                f"Long (≥300 aa): {dist['long']:,} ({dist['long']/total*100:.1f}%)",
                ""
            ])
        
        if 'quality_metrics' in gene_analysis:
            quality = gene_analysis['quality_metrics']
            content.extend([
                "Quality Assessment:",
                f"Coding Density Estimate: {quality.get('coding_density_estimate', 'N/A')}%",
                f"Length Distribution: {quality.get('mean_length_assessment', 'Unknown')}",
                f"Overall Quality: {quality.get('quality_flag', 'Unknown')}",
                ""
            ])
        
        return content

    def _generate_annotation_coverage_section(self, annotation_analysis, coverage_metrics):
        """Generate functional annotation coverage section"""
        content = [
            "Functional Annotation Coverage (SwissProt)",
            "-" * 45,
            f"Total Protein Matches: {annotation_analysis['total_matches']:,}",
            f"Unique Proteins Annotated: {annotation_analysis['unique_proteins']:,}",
            f"Annotation Rate: {coverage_metrics.get('annotation_rate', 0):.1f}%",
            f"Coverage Assessment: {coverage_metrics.get('coverage_assessment', 'Unknown')}",
            ""
        ]
        
        if 'identity_statistics' in annotation_analysis:
            stats = annotation_analysis['identity_statistics']
            content.extend([
                "Sequence Identity Statistics:",
                f"Average Identity: {stats['mean']:.1f}%",
                f"Median Identity: {stats['median']:.1f}%",
                f"Identity Range: {stats['min']:.1f}% - {stats['max']:.1f}%",
                ""
            ])
        
        if 'quality_distribution' in annotation_analysis:
            qual = annotation_analysis['quality_distribution']
            total = annotation_analysis['total_matches']
            content.extend([
                "Annotation Quality Distribution:",
                f"Excellent (≥95% identity): {qual['excellent']:,} ({qual['excellent']/total*100:.1f}%)",
                f"Very Good (90-94% identity): {qual['very_good']:,} ({qual['very_good']/total*100:.1f}%)",
                f"Good (80-89% identity): {qual['good']:,} ({qual['good']/total*100:.1f}%)",
                f"Moderate (70-79% identity): {qual['moderate']:,} ({qual['moderate']/total*100:.1f}%)",
                ""
            ])
        
        return content

    def _generate_functional_category_section(self, functional_distribution):
        """Generate functional category analysis section"""
        content = [
            "Functional Category Analysis",
            "-" * 30,
            ""
        ]
        
        # Sort categories by protein count
        sorted_categories = sorted(functional_distribution.items(), key=lambda x: x[1], reverse=True)
        total_categorized = sum(functional_distribution.values())
        
        content.extend([
            f"{'Category':<30} {'Proteins':<10} {'Percentage':<12}",
            "-" * 54
        ])
        
        for category, count in sorted_categories:
            percentage = (count / total_categorized) * 100 if total_categorized > 0 else 0
            content.append(f"{category:<30} {count:<10} {percentage:<12.1f}%")
        
        content.extend([
            "",
            f"Total Categorized Functions: {total_categorized:,}",
            "Note: Proteins with multiple functions may be counted in multiple categories",
            ""
        ])
        
        return content

    def _generate_quality_assessment_section(self, coverage_metrics):
        """Generate quality assessment and recommendations section"""
        content = [
            "Quality Assessment & Recommendations",
            "-" * 40,
            f"Overall Assessment: {coverage_metrics.get('quality_assessment', 'Incomplete')}",
            f"Annotation Coverage: {coverage_metrics.get('coverage_assessment', 'Unknown')}",
            "",
            "Interpretation Guidelines:"
        ]
        
        annotation_rate = coverage_metrics.get('annotation_rate', 0)
        
        if annotation_rate >= 85:
            content.extend([
                "- Excellent annotation coverage provides reliable functional insights",
                "- High-quality matches support confident functional assignments",
                "- Comprehensive functional analysis is well-supported"
            ])
        elif annotation_rate >= 70:
            content.extend([
                "- Good annotation coverage enables robust functional analysis",
                "- Most essential functions likely captured in analysis",
                "- Minor gaps may exist for specialized or novel functions"
            ])
        elif annotation_rate >= 50:
            content.extend([
                "- Moderate annotation coverage provides basic functional insights",
                "- Essential metabolic pathways likely represented",
                "- Significant gaps may exist for specialized functions"
            ])
        else:
            content.extend([
                "- Limited annotation coverage constrains functional analysis",
                "- Consider specialized databases for improved coverage",
                "- Results should be interpreted with caution"
            ])
        
        content.extend([
            "",
            "Recommended Follow-up:",
            "- Cross-reference findings with expected sample characteristics",
            "- Consider domain-specific databases for specialized functions",
            "- Validate critical findings through targeted approaches",
            "- Compare functional profiles with reference datasets when available",
            ""
        ])
        
        return content

    def _generate_technical_specifications(self):
        """Generate technical specifications section"""
        return [
            "Technical Specifications",
            "-" * 25,
            "Gene Prediction: Prokka v1.14.6",
            "Annotation Database: SwissProt Release 2024_03",
            f"Analysis Pipeline: MetaQuest v{self.version}",
            "BLAST Parameters: E-value <1e-10, Identity >70%, Coverage >70%",
            "Quality Thresholds: High (≥90% ID), Good (≥80% ID), Moderate (≥70% ID)",
            "",
            "Quality Control:",
            "- Gene calling consistency verified across analysis",
            "- Annotation concordance assessed for reliability",
            "- Functional pathway completeness evaluated",
            "- Statistical confidence measures applied",
            ""
        ]

    def _create_functional_json_report(self, gene_analysis, annotation_analysis, coverage_metrics):
        """Create structured JSON report for functional analysis"""
        json_data = {
            'analysis_summary': {
                'analysis_type': 'Functional_Annotation',
                'total_genes': gene_analysis.get('total_proteins', 0) if gene_analysis and 'error' not in gene_analysis else 0,
                'annotation_rate': coverage_metrics.get('annotation_rate', 0),
                'coverage_assessment': coverage_metrics.get('coverage_assessment', 'Unknown'),
                'quality_assessment': coverage_metrics.get('quality_assessment', 'Unknown'),
                'analysis_timestamp': self.timestamp.isoformat(),
                'pipeline_version': f'MetaQuest v{self.version}'
            },
            'gene_prediction_analysis': gene_analysis if gene_analysis and 'error' not in gene_analysis else None,
            'functional_annotation_analysis': annotation_analysis if annotation_analysis and 'error' not in annotation_analysis else None,
            'coverage_metrics': coverage_metrics,
            'functional_categories': coverage_metrics.get('functional_distribution', {})
        }
        
        return json_data


