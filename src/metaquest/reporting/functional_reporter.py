#!/usr/bin/env python3
"""
MetaQuest Professional Functional Annotation Reporting Module
Author: MetaQuest Metagenomics Team

Provides comprehensive functional annotation reporting with gene prediction analysis,
protein annotation coverage, and functional category profiling.
"""

import json
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Any
from collections import Counter
from Bio import SeqIO
import warnings
warnings.filterwarnings('ignore')

from .base_reporter import BaseReportGenerator


class FunctionalReporter(BaseReportGenerator):
    """
    Professional functional annotation reporting following scientific standards.
    
    Analyzes:
    - Gene prediction results (Prokka)
    - Protein functional annotation (SwissProt/UniProt)
    - Functional category distribution
    - Annotation quality and coverage metrics
    """
    
    # Functional category keywords for automated classification
    FUNCTIONAL_CATEGORIES = {
        'Antimicrobial Resistance': {
            'keywords': [
                'beta-lactamase', 'carbapenemase', 'penicillinase', 'cephalosporinase',
                'efflux pump', 'multidrug resistance', 'antibiotic resistance',
                'vancomycin resistance', 'methicillin resistance', 'fluoroquinolone resistance',
                'aminoglycoside resistance', 'tetracycline resistance'
            ],
            'priority': 'critical'
        },
        'Virulence Factors': {
            'keywords': [
                'virulence', 'toxin', 'adhesin', 'invasin', 'hemolysin', 'enterotoxin',
                'cytotoxin', 'secretion system', 'pathogenicity', 'colonization factor',
                'fimbrial', 'pilus', 'capsule biosynthesis'
            ],
            'priority': 'high'
        },
        'Central Metabolism': {
            'keywords': [
                'glycolysis', 'tca cycle', 'citric acid cycle', 'pentose phosphate',
                'gluconeogenesis', 'oxidative phosphorylation', 'atp synthase',
                'electron transport', 'respiration'
            ],
            'priority': 'medium'
        },
        'Amino Acid Metabolism': {
            'keywords': [
                'amino acid biosynthesis', 'amino acid degradation', 'aminotransferase',
                'decarboxylase', 'amino acid transporter'
            ],
            'priority': 'medium'
        },
        'Nucleotide Metabolism': {
            'keywords': [
                'nucleotide biosynthesis', 'purine', 'pyrimidine', 'ribonucleotide reductase',
                'thymidylate synthase'
            ],
            'priority': 'medium'
        },
        'Lipid Metabolism': {
            'keywords': [
                'fatty acid biosynthesis', 'fatty acid degradation', 'lipid biosynthesis',
                'phospholipid', 'lipopolysaccharide'
            ],
            'priority': 'medium'
        },
        'Transport Systems': {
            'keywords': [
                'abc transporter', 'permease', 'channel protein', 'membrane transporter',
                'ion transport', 'substrate transport', 'drug efflux'
            ],
            'priority': 'medium'
        },
        'Transcriptional Regulation': {
            'keywords': [
                'transcriptional regulator', 'repressor', 'activator', 'two-component system',
                'response regulator', 'sensor kinase', 'sigma factor', 'dna-binding protein'
            ],
            'priority': 'medium'
        },
        'Cell Wall Synthesis': {
            'keywords': [
                'peptidoglycan', 'murein', 'cell wall biosynthesis', 'penicillin-binding protein',
                'transpeptidase', 'transglycosylase', 'murein hydrolase'
            ],
            'priority': 'medium'
        },
        'Stress Response': {
            'keywords': [
                'heat shock protein', 'cold shock protein', 'oxidative stress', 'chaperone',
                'acid resistance', 'osmotic stress', 'universal stress protein'
            ],
            'priority': 'medium'
        },
        'DNA Replication & Repair': {
            'keywords': [
                'dna polymerase', 'helicase', 'primase', 'ligase', 'topoisomerase',
                'recombination', 'dna repair', 'mismatch repair', 'nucleotide excision'
            ],
            'priority': 'medium'
        },
        'Protein Synthesis': {
            'keywords': [
                'ribosomal protein', 'translation factor', 'aminoacyl-trna synthetase',
                'peptidyl transferase', 'elongation factor'
            ],
            'priority': 'medium'
        }
    }
    
    def __init__(self, output_dir: Path):
        """Initialize functional reporter."""
        super().__init__(output_dir, "Functional Annotation")
        self._gene_stats = None
        self._annotation_stats = None
    
    def generate_report(self, prokka_results: Optional[str] = None, 
                       swissprot_results: Optional[str] = None, **kwargs) -> Dict[str, str]:
        """
        Generate comprehensive functional annotation report.
        
        Args:
            prokka_results: Path to Prokka output directory
            swissprot_results: Path to SwissProt annotation TSV file
            
        Returns:
            Dictionary with paths to generated reports and summary statistics
        """
        print(f"[{self.timestamp.strftime('%H:%M:%S')}] Initiating functional annotation report generation...")
        
        # Parse input data
        gene_data = self._parse_gene_prediction(prokka_results) if prokka_results else None
        annotation_data = self._parse_functional_annotation(swissprot_results) if swissprot_results else None
        
        # Calculate coverage metrics
        coverage_metrics = self._calculate_coverage_metrics(gene_data, annotation_data)
        
        # Generate text report
        text_content = self._generate_text_report(gene_data, annotation_data, coverage_metrics)
        text_file = self._save_text_report(text_content, "functional_annotation_report.txt")
        
        # Generate JSON report
        json_data = self._generate_json_report(gene_data, annotation_data, coverage_metrics)
        json_file = self._save_json_report(json_data, "functional_annotation_report.json")
        
        print(f"[{self.timestamp.strftime('%H:%M:%S')}] Report generation complete.")
        print(f"  ✓ Text report: {text_file}")
        print(f"  ✓ JSON report: {json_file}")
        
        return {
            'text_report': text_file,
            'json_report': json_file,
            'total_genes': gene_data.get('total_proteins', 0) if gene_data else 0,
            'annotation_rate': coverage_metrics.get('annotation_rate', 0),
            'functional_categories': len(coverage_metrics.get('functional_distribution', {}))
        }
    
    def _parse_gene_prediction(self, prokka_dir: str) -> Optional[Dict]:
        """
        Parse Prokka gene prediction results.
        
        Args:
            prokka_dir: Path to Prokka output directory
            
        Returns:
            Dictionary containing gene prediction statistics
        """
        try:
            prokka_path = Path(prokka_dir)
            
            # Find output files
            faa_files = list(prokka_path.glob("*.faa"))
            ffn_files = list(prokka_path.glob("*.ffn"))
            gff_files = list(prokka_path.glob("*.gff"))
            
            if not faa_files:
                print("  ✗ No protein sequence files (.faa) found in Prokka directory")
                return None
            
            # Parse protein sequences
            protein_file = faa_files[0]
            proteins = list(SeqIO.parse(protein_file, "fasta"))
            
            if not proteins:
                print("  ✗ No proteins found in file")
                return None
            
            # Calculate statistics
            lengths = [len(seq.seq) for seq in proteins]
            
            # Length distribution
            length_dist = {
                'very_short': sum(1 for l in lengths if l < 50),
                'short': sum(1 for l in lengths if 50 <= l < 100),
                'medium': sum(1 for l in lengths if 100 <= l < 300),
                'long': sum(1 for l in lengths if l >= 300)
            }
            
            # Quality metrics
            very_short_ratio = length_dist['very_short'] / len(proteins)
            mean_length = np.mean(lengths)
            
            quality_flag = 'Pass' if very_short_ratio < 0.3 and mean_length > 150 else 'Review'
            
            gene_data = {
                'prokka_directory': str(prokka_path),
                'total_proteins': len(proteins),
                'output_files': {
                    'protein_files': len(faa_files),
                    'gene_files': len(ffn_files),
                    'annotation_files': len(gff_files)
                },
                'length_statistics': {
                    'mean': float(np.mean(lengths)),
                    'median': float(np.median(lengths)),
                    'std': float(np.std(lengths)),
                    'min': int(min(lengths)),
                    'max': int(max(lengths)),
                    'q25': float(np.percentile(lengths, 25)),
                    'q75': float(np.percentile(lengths, 75))
                },
                'length_distribution': length_dist,
                'quality_metrics': {
                    'very_short_ratio': float(very_short_ratio),
                    'mean_length': float(mean_length),
                    'quality_flag': quality_flag,
                    'coding_density_estimate': self._estimate_coding_density(lengths)
                }
            }
            
            self._gene_stats = gene_data
            
            print(f"  ✓ Parsed gene prediction data: {len(proteins):,} proteins identified")
            
            return gene_data
            
        except Exception as e:
            print(f"  ✗ Error parsing gene prediction data: {e}")
            return None
    
    def _parse_functional_annotation(self, swissprot_file: str) -> Optional[Dict]:
        """
        Parse SwissProt functional annotation results.
        
        Args:
            swissprot_file: Path to SwissProt BLAST TSV file
            
        Returns:
            Dictionary containing annotation statistics
        """
        try:
            swissprot_path = Path(swissprot_file)
            
            if not swissprot_path.exists():
                print(f"  ✗ SwissProt annotation file not found: {swissprot_file}")
                return None
            
            # Read annotation file
            df = pd.read_csv(swissprot_file, sep='\t')
            
            if df.empty:
                print("  ✗ No annotations found in SwissProt file")
                return None
            
            # Calculate statistics
            identity_stats = {
                'mean': float(df['pident'].mean()),
                'median': float(df['pident'].median()),
                'std': float(df['pident'].std()),
                'min': float(df['pident'].min()),
                'max': float(df['pident'].max())
            }
            
            # Quality distribution
            quality_dist = {
                'excellent': len(df[df['pident'] >= 95]),
                'very_good': len(df[(df['pident'] >= 90) & (df['pident'] < 95)]),
                'good': len(df[(df['pident'] >= 80) & (df['pident'] < 90)]),
                'moderate': len(df[(df['pident'] >= 70) & (df['pident'] < 80)]),
                'low': len(df[df['pident'] < 70])
            }
            
            # E-value statistics
            evalue_stats = {
                'mean': float(df['evalue'].mean()),
                'median': float(df['evalue'].median()),
                'highly_significant': len(df[df['evalue'] < 1e-50]),
                'significant': len(df[(df['evalue'] >= 1e-50) & (df['evalue'] < 1e-10)]),
                'marginal': len(df[df['evalue'] >= 1e-10])
            }
            
            # Functional category analysis
            functional_analysis = {}
            if 'stitle' in df.columns:
                functional_analysis = self._analyze_functional_categories(df['stitle'])
            
            annotation_data = {
                'total_matches': len(df),
                'unique_proteins': df['qseqid'].nunique(),
                'identity_statistics': identity_stats,
                'quality_distribution': quality_dist,
                'evalue_statistics': evalue_stats,
                'functional_analysis': functional_analysis
            }
            
            self._annotation_stats = annotation_data
            
            print(f"  ✓ Parsed annotation data: {len(df):,} matches for {df['qseqid'].nunique():,} unique proteins")
            
            return annotation_data
            
        except Exception as e:
            print(f"  ✗ Error parsing functional annotation data: {e}")
            return None
    
    def _analyze_functional_categories(self, descriptions: pd.Series) -> Dict[str, int]:
        """
        Analyze functional categories from protein descriptions.
        
        Args:
            descriptions: Series of protein descriptions
            
        Returns:
            Dictionary mapping categories to protein counts
        """
        category_counts = {}
        descriptions_lower = descriptions.str.lower().fillna('')
        
        for category, config in self.FUNCTIONAL_CATEGORIES.items():
            count = 0
            for keyword in config['keywords']:
                count += descriptions_lower.str.contains(keyword, regex=False, na=False).sum()
            
            if count > 0:
                category_counts[category] = count
        
        return category_counts
    
    def _estimate_coding_density(self, protein_lengths: List[int]) -> float:
        """
        Estimate coding density from protein length distribution.
        
        Args:
            protein_lengths: List of protein lengths in amino acids
            
        Returns:
            Estimated coding density percentage
        """
        total_aa = sum(protein_lengths)
        total_nucleotides = total_aa * 3
        
        # Rough genome size estimate (assumes typical bacterial genome)
        estimated_genome_size = total_nucleotides * 1.15  # Account for non-coding
        coding_density = (total_nucleotides / estimated_genome_size) * 100
        
        return round(coding_density, 1)
    
    def _calculate_coverage_metrics(self, gene_data: Optional[Dict], 
                                    annotation_data: Optional[Dict]) -> Dict[str, Any]:
        """
        Calculate annotation coverage metrics.
        
        Args:
            gene_data: Gene prediction statistics
            annotation_data: Annotation statistics
            
        Returns:
            Dictionary of coverage metrics
        """
        metrics = {
            'annotation_rate': 0.0,
            'coverage_assessment': 'No data',
            'total_predicted': 0,
            'total_annotated': 0,
            'unannotated': 0,
            'functional_distribution': {},
            'quality_assessment': 'Incomplete analysis'
        }
        
        if not gene_data or not annotation_data:
            return metrics
        
        total_proteins = gene_data.get('total_proteins', 0)
        annotated_proteins = annotation_data.get('unique_proteins', 0)
        
        if total_proteins > 0:
            annotation_rate = (annotated_proteins / total_proteins) * 100
            
            metrics.update({
                'annotation_rate': round(annotation_rate, 2),
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
            metrics['functional_distribution'] = annotation_data.get('functional_analysis', {})
            
            # Quality assessment
            quality_dist = annotation_data.get('quality_distribution', {})
            high_quality = quality_dist.get('excellent', 0) + quality_dist.get('very_good', 0)
            quality_ratio = high_quality / max(annotated_proteins, 1)
            
            if quality_ratio >= 0.7:
                metrics['quality_assessment'] = 'High-quality annotations predominate'
            elif quality_ratio >= 0.5:
                metrics['quality_assessment'] = 'Good annotation quality'
            elif quality_ratio >= 0.3:
                metrics['quality_assessment'] = 'Mixed annotation quality'
            else:
                metrics['quality_assessment'] = 'Lower confidence annotations'
        
        return metrics
    
    def _generate_text_report(self, gene_data: Optional[Dict], 
                             annotation_data: Optional[Dict],
                             coverage_metrics: Dict) -> List[str]:
        """Generate comprehensive text report."""
        content = []
        
        # Header
        content.extend(self._create_header(
            "FUNCTIONAL ANNOTATION REPORT",
            "Prokka v1.14.6 + SwissProt BLAST v2.14.0",
            "SwissProt Release 2024_05"
        ))
        
        # Executive Summary
        content.extend(self._create_section_header("EXECUTIVE SUMMARY", level=1))
        content.extend(self._generate_executive_summary(gene_data, annotation_data, coverage_metrics))
        
        # Gene Prediction Section
        if gene_data:
            content.extend(self._create_section_header("GENE PREDICTION ANALYSIS", level=1))
            content.extend(self._generate_gene_prediction_section(gene_data))
        
        # Functional Annotation Section
        if annotation_data:
            content.extend(self._create_section_header("FUNCTIONAL ANNOTATION ANALYSIS", level=1))
            content.extend(self._generate_annotation_section(annotation_data, coverage_metrics))
        
        # Functional Category Section
        if coverage_metrics.get('functional_distribution'):
            content.extend(self._create_section_header("FUNCTIONAL CATEGORY PROFILING", level=1))
            content.extend(self._generate_functional_category_section(coverage_metrics['functional_distribution']))
        
        # Quality Assessment
        content.extend(self._create_section_header("QUALITY ASSESSMENT", level=1))
        content.extend(self._generate_quality_assessment_section(gene_data, annotation_data, coverage_metrics))
        
        # Recommendations
        content.extend(self._create_section_header("RECOMMENDATIONS", level=1))
        content.extend(self._generate_recommendations_section(coverage_metrics))
        
        # Footer
        content.extend(self._create_footer())
        
        return content
    
    def _generate_executive_summary(self, gene_data: Optional[Dict], 
                                    annotation_data: Optional[Dict],
                                    coverage_metrics: Dict) -> List[str]:
        """Generate executive summary."""
        summary_items = []
        
        if gene_data:
            total_genes = gene_data['total_proteins']
            mean_length = gene_data['length_statistics']['mean']
            quality_flag = gene_data['quality_metrics']['quality_flag']
            
            summary_items.append(
                f"Gene prediction analysis identified {total_genes:,} protein-coding genes with an average "
                f"length of {mean_length:.0f} amino acids. Quality assessment: {quality_flag}."
            )
        
        if annotation_data and coverage_metrics.get('annotation_rate', 0) > 0:
            annotation_rate = coverage_metrics['annotation_rate']
            coverage_assessment = coverage_metrics['coverage_assessment']
            total_annotated = coverage_metrics['total_annotated']
            
            summary_items.append(
                f"Functional annotation achieved {coverage_assessment.lower()} coverage with {annotation_rate:.1f}% "
                f"of proteins annotated ({total_annotated:,} proteins assigned to known functions)."
            )
            
            # Functional categories
            func_dist = coverage_metrics.get('functional_distribution', {})
            if func_dist:
                top_category = max(func_dist.items(), key=lambda x: x[1])
                summary_items.append(
                    f"Functional profiling identified {len(func_dist)} major functional categories, with "
                    f"{top_category[0]} being the most represented ({top_category[1]} annotations)."
                )
        
        if not summary_items:
            summary_items.append("Functional annotation analysis requires input data for comprehensive assessment.")
        
        return self._create_summary_box("Key Findings", summary_items)
    
    def _generate_gene_prediction_section(self, gene_data: Dict) -> List[str]:
        """Generate gene prediction analysis section."""
        lines = []
        
        # Overview
        lines.extend(self._create_section_header("Prediction Overview", level=2))
        lines.extend([
            f"Gene Prediction Tool: Prokka v1.14.6",
            f"Total Genes Identified: {gene_data['total_proteins']:,}",
            f"Protein-Coding Genes: {gene_data['total_proteins']:,} (100%)",
            f"Analysis Directory: {Path(gene_data['prokka_directory']).name}",
            ""
        ])
        
        # Length statistics
        lines.extend(self._create_section_header("Protein Length Statistics", level=2))
        stats = gene_data['length_statistics']
        
        lines.extend([
            f"Mean Length: {stats['mean']:.1f} ± {stats['std']:.1f} amino acids",
            f"Median Length: {stats['median']:.0f} amino acids",
            f"Length Range: {stats['min']} - {stats['max']} amino acids",
            f"25th Percentile: {stats['q25']:.0f} amino acids",
            f"75th Percentile: {stats['q75']:.0f} amino acids",
            ""
        ])
        
        # Length distribution table
        lines.extend(self._create_section_header("Length Category Distribution", level=2))
        
        dist = gene_data['length_distribution']
        total = gene_data['total_proteins']
        
        headers = ["Category", "Range (aa)", "Count", "Percentage"]
        rows = [
            ["Very Short", "<50", f"{dist['very_short']:,}", f"{dist['very_short']/total*100:.2f}%"],
            ["Short", "50-99", f"{dist['short']:,}", f"{dist['short']/total*100:.2f}%"],
            ["Medium", "100-299", f"{dist['medium']:,}", f"{dist['medium']/total*100:.2f}%"],
            ["Long", "≥300", f"{dist['long']:,}", f"{dist['long']/total*100:.2f}%"]
        ]
        
        lines.extend(self._format_table(headers, rows, alignments=['left', 'left', 'right', 'right']))
        lines.append("")
        
        # Quality metrics
        lines.extend(self._create_section_header("Quality Metrics", level=2))
        quality = gene_data['quality_metrics']
        
        lines.extend([
            f"Very Short Protein Ratio: {quality['very_short_ratio']:.3f}",
            f"Mean Length Assessment: {'Normal' if quality['mean_length'] > 200 else 'Short'}",
            f"Estimated Coding Density: {quality['coding_density_estimate']}%",
            f"Overall Quality Flag: {quality['quality_flag']}",
            ""
        ])
        
        return lines
    
    def _generate_annotation_section(self, annotation_data: Dict, 
                                     coverage_metrics: Dict) -> List[str]:
        """Generate functional annotation analysis section."""
        lines = []
        
        # Coverage overview
        lines.extend(self._create_section_header("Annotation Coverage", level=2))
        lines.extend([
            f"Annotation Database: SwissProt Release 2024_05",
            f"Total Protein Matches: {annotation_data['total_matches']:,}",
            f"Unique Proteins Annotated: {annotation_data['unique_proteins']:,}",
            f"Annotation Rate: {coverage_metrics['annotation_rate']:.2f}%",
            f"Coverage Assessment: {coverage_metrics['coverage_assessment']}",
            ""
        ])
        
        # Identity statistics
        lines.extend(self._create_section_header("Sequence Identity Statistics", level=2))
        stats = annotation_data['identity_statistics']
        
        lines.extend([
            f"Average Identity: {stats['mean']:.2f}%",
            f"Median Identity: {stats['median']:.2f}%",
            f"Standard Deviation: {stats['std']:.2f}%",
            f"Identity Range: {stats['min']:.1f}% - {stats['max']:.1f}%",
            ""
        ])
        
        # Quality distribution
        lines.extend(self._create_section_header("Annotation Quality Distribution", level=2))
        
        qual = annotation_data['quality_distribution']
        total = annotation_data['total_matches']
        
        headers = ["Quality Level", "Identity Range", "Count", "Percentage"]
        rows = [
            ["Excellent", "≥95%", f"{qual['excellent']:,}", f"{qual['excellent']/total*100:.2f}%"],
            ["Very Good", "90-94%", f"{qual['very_good']:,}", f"{qual['very_good']/total*100:.2f}%"],
            ["Good", "80-89%", f"{qual['good']:,}", f"{qual['good']/total*100:.2f}%"],
            ["Moderate", "70-79%", f"{qual['moderate']:,}", f"{qual['moderate']/total*100:.2f}%"],
            ["Low", "<70%", f"{qual['low']:,}", f"{qual['low']/total*100:.2f}%"]
        ]
        
        lines.extend(self._format_table(headers, rows, alignments=['left', 'left', 'right', 'right']))
        lines.append("")
        
        # E-value statistics
        lines.extend(self._create_section_header("Statistical Significance (E-values)", level=2))
        evalue = annotation_data['evalue_statistics']
        
        lines.extend([
            f"Mean E-value: {evalue['mean']:.2e}",
            f"Median E-value: {evalue['median']:.2e}",
            f"Highly Significant (E < 1e-50): {evalue['highly_significant']:,} ({evalue['highly_significant']/total*100:.1f}%)",
            f"Significant (1e-50 ≤ E < 1e-10): {evalue['significant']:,} ({evalue['significant']/total*100:.1f}%)",
            f"Marginal (E ≥ 1e-10): {evalue['marginal']:,} ({evalue['marginal']/total*100:.1f}%)",
            ""
        ])
        
        return lines
    
    def _generate_functional_category_section(self, functional_dist: Dict[str, int]) -> List[str]:
        """Generate functional category profiling section."""
        lines = []
        
        # Sort categories by count
        sorted_categories = sorted(functional_dist.items(), key=lambda x: x[1], reverse=True)
        total_categorized = sum(functional_dist.values())
        
        lines.extend(self._create_section_header("Functional Category Overview", level=2))
        lines.extend([
            f"Total Functional Categories Detected: {len(functional_dist)}",
            f"Total Categorized Annotations: {total_categorized:,}",
            f"Note: Proteins with multiple functions may be counted in multiple categories",
            ""
        ])
        
        # Category table
        lines.extend(self._create_section_header("Category Distribution", level=2))
        
        headers = ["Rank", "Functional Category", "Protein Count", "Percentage"]
        rows = []
        
        for i, (category, count) in enumerate(sorted_categories, 1):
            percentage = (count / total_categorized * 100) if total_categorized > 0 else 0
            rows.append([
                str(i),
                category,
                f"{count:,}",
                f"{percentage:.2f}%"
            ])
        
        lines.extend(self._format_table(headers, rows, alignments=['left', 'left', 'right', 'right']))
        lines.append("")
        
        # Priority annotations
        priority_categories = {
            k: v for k, v in functional_dist.items() 
            if self.FUNCTIONAL_CATEGORIES.get(k, {}).get('priority') in ['critical', 'high']
        }
        
        if priority_categories:
            lines.extend(self._create_section_header("High-Priority Functional Categories", level=2))
            lines.extend([
                "The following categories warrant special attention:",
                ""
            ])
            
            for category, count in sorted(priority_categories.items(), key=lambda x: x[1], reverse=True):
                priority = self.FUNCTIONAL_CATEGORIES[category]['priority']
                lines.append(f"• {category} ({priority.upper()} priority): {count:,} proteins")
            
            lines.append("")
        
        return lines
    
    def _generate_quality_assessment_section(self, gene_data: Optional[Dict], 
                                             annotation_data: Optional[Dict],
                                             coverage_metrics: Dict) -> List[str]:
        """Generate quality assessment section."""
        lines = []
        
        # Overall assessment
        lines.extend(self._create_section_header("Overall Quality Assessment", level=2))
        lines.extend([
            f"Annotation Coverage: {coverage_metrics['coverage_assessment']}",
            f"Annotation Quality: {coverage_metrics['quality_assessment']}",
            ""
        ])
        
        # Interpretation
        annotation_rate = coverage_metrics.get('annotation_rate', 0)
        
        if annotation_rate >= 85:
            interpretation = ("Excellent annotation coverage provides high confidence in functional assignments. "
                            "The majority of predicted proteins have reliable functional annotations, enabling "
                            "comprehensive pathway analysis and functional interpretation.")
        elif annotation_rate >= 70:
            interpretation = ("Good annotation coverage supports robust functional analysis. Most essential "
                            "metabolic functions are likely well-represented, though some specialized or novel "
                            "functions may lack annotations.")
        elif annotation_rate >= 50:
            interpretation = ("Moderate annotation coverage provides basic functional insights. Essential "
                            "housekeeping genes and core metabolic pathways are likely captured, but significant "
                            "gaps may exist for specialized functions or novel proteins.")
        else:
            interpretation = ("Limited annotation coverage constrains comprehensive functional analysis. "
                            "Consider using additional databases or de novo functional prediction methods "
                            "for improved coverage.")
        
        lines.extend([interpretation, ""])
        
        # Validation metrics
        if gene_data and annotation_data:
            lines.extend(self._create_section_header("Validation Metrics", level=2))
            lines.extend([
                f"Gene Calling Quality: {gene_data['quality_metrics']['quality_flag']}",
                f"Mean Protein Length: {gene_data['length_statistics']['mean']:.0f} amino acids",
                f"Annotation Identity: {annotation_data['identity_statistics']['mean']:.1f}% average",
                f"High-Quality Annotations (≥90% ID): {sum([annotation_data['quality_distribution']['excellent'], annotation_data['quality_distribution']['very_good']]):,}",
                ""
            ])
        
        return lines
    
    def _generate_recommendations_section(self, coverage_metrics: Dict) -> List[str]:
        """Generate recommendations section."""
        lines = []
        recommendations = []
        
        annotation_rate = coverage_metrics.get('annotation_rate', 0)
        coverage_assessment = coverage_metrics.get('coverage_assessment', 'Unknown')
        
        # Data-driven recommendations
        if annotation_rate < 70:
            recommendations.extend([
                "Consider supplementing SwissProt annotations with additional databases (KEGG, COG, Pfam) "
                "to improve functional coverage.",
                "Evaluate unannotated proteins for novel or highly divergent sequences that may represent "
                "unique adaptations or horizontal gene transfer events."
            ])
        
        if coverage_assessment == 'Excellent':
            recommendations.append(
                "High annotation coverage enables comprehensive pathway reconstruction and metabolic modeling. "
                "Consider performing pathway enrichment analysis and comparative genomics."
            )
        
        # General best practices
        recommendations.extend([
            "Cross-reference functional annotations with domain predictions (InterProScan, Pfam) for enhanced confidence.",
            "Perform pathway completeness analysis to identify metabolic capabilities and potential auxotrophies.",
            "Validate critical functional assignments (e.g., virulence factors, AMR genes) through targeted analysis.",
            "Compare functional profiles with phylogenetically related organisms to identify unique features.",
            "Consider manual curation of key functional categories for publication-quality annotations."
        ])
        
        # Format recommendations
        for i, rec in enumerate(recommendations, 1):
            lines.append(f"{i}. {rec}")
            lines.append("")
        
        return lines
    
    def _generate_json_report(self, gene_data: Optional[Dict], 
                             annotation_data: Optional[Dict],
                             coverage_metrics: Dict) -> Dict:
        """Generate structured JSON report."""
        json_data = {
            'summary': {
                'total_genes_predicted': gene_data.get('total_proteins', 0) if gene_data else 0,
                'annotation_rate': coverage_metrics.get('annotation_rate', 0),
                'coverage_assessment': coverage_metrics.get('coverage_assessment', 'No data'),
                'quality_assessment': coverage_metrics.get('quality_assessment', 'No data'),
                'functional_categories_detected': len(coverage_metrics.get('functional_distribution', {}))
            },
            'gene_prediction': gene_data,
            'functional_annotation': annotation_data,
            'coverage_metrics': coverage_metrics
        }
        
        return json_data