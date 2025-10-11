#!/usr/bin/env python3
"""
MetaQuest Professional Taxonomic Reporting Module
Author: MetaQuest Metagenomics Team

Provides comprehensive taxonomic classification reporting with statistical analysis,
diversity metrics, and quality assessment for both FASTQ and FASTA pipelines.
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from collections import Counter
from skbio.diversity import alpha_diversity
from .base_reporter import BaseReportGenerator


class TaxonomicReporter(BaseReportGenerator):
    """
    Professional taxonomic classification reporting with comprehensive statistical analysis.
    
    Supports:
    - Kraken2/Bracken taxonomic profiling (FASTQ)
    - BLAST-based taxonomic assignment (FASTA)
    - Dual-method integrated analysis
    - Alpha diversity metrics (Shannon, Simpson, Chao1)
    - Quality assessment and recommendations
    """
    
    # Quality thresholds
    MIN_ABUNDANCE_THRESHOLD = 0.0001  # 0.01%
    HIGH_QUALITY_IDENTITY = 95.0
    MEDIUM_QUALITY_IDENTITY = 80.0
    
    def __init__(self, output_dir: Path):
        """Initialize taxonomic reporter."""
        super().__init__(output_dir, "Taxonomic Classification")
        self._diversity_data = None
        self._quality_data = None
    
    def generate_report(self, bracken_data: Optional[Any] = None, 
                       blast_data: Optional[Any] = None, **kwargs) -> Dict[str, str]:
        """
        Generate comprehensive taxonomic classification report.
        
        Args:
            bracken_data: Bracken/Kraken output file or DataFrame
            blast_data: BLAST results JSON file or data structure
            **kwargs: Additional parameters
            
        Returns:
            Dictionary with paths to generated reports and summary statistics
        """
        print(f"[{self.timestamp.strftime('%H:%M:%S')}] Initiating taxonomic classification report generation...")
        
        # Determine pipeline configuration
        pipeline_config = self._determine_pipeline(bracken_data, blast_data)
        
        # Parse input data
        bracken_parsed = self._parse_bracken_data(bracken_data) if bracken_data else None
        blast_parsed = self._parse_blast_data(blast_data) if blast_data else None
        
        # Generate text report
        text_content = self._generate_text_report(bracken_parsed, blast_parsed, pipeline_config)
        text_file = self._save_text_report(text_content, "taxonomic_classification_report.txt")
        
        # Generate JSON report
        json_data = self._generate_json_report(bracken_parsed, blast_parsed, pipeline_config)
        json_file = self._save_json_report(json_data, "taxonomic_classification_report.json")
        
        print(f"[{self.timestamp.strftime('%H:%M:%S')}] Report generation complete.")
        print(f"  ✓ Text report: {text_file}")
        print(f"  ✓ JSON report: {json_file}")
        
        return {
            'text_report': text_file,
            'json_report': json_file,
            'pipeline_mode': pipeline_config['mode'],
            'total_taxa': json_data.get('summary', {}).get('total_taxa_detected', 0)
        }
    
    def _determine_pipeline(self, bracken_data: Any, blast_data: Any) -> Dict[str, str]:
        """Determine pipeline configuration based on input data."""
        if bracken_data and blast_data:
            return {
                'mode': 'Dual-Method',
                'tool_info': 'Kraken2 v2.1.3 / Bracken v2.8 + BLAST v2.14.0',
                'database_info': 'Standard-8 (16GB) + NCBI NT (2024)'
            }
        elif bracken_data:
            return {
                'mode': 'FASTQ',
                'tool_info': 'Kraken2 v2.1.3 / Bracken v2.8',
                'database_info': 'Standard-8 (16GB Kraken2 Database)'
            }
        elif blast_data:
            return {
                'mode': 'FASTA',
                'tool_info': 'BLAST v2.14.0 + Taxonomic Assignment',
                'database_info': 'NCBI Nucleotide (NT) Database 2024'
            }
        else:
            return {
                'mode': 'Unknown',
                'tool_info': 'Not specified',
                'database_info': 'Not specified'
            }
    
    def _parse_bracken_data(self, data: Any) -> Optional[Dict]:
        """
        Parse Bracken/Kraken data into standardized format.
        
        Args:
            data: Bracken file path, DataFrame, or dict
            
        Returns:
            Parsed data dictionary or None if parsing fails
        """
        try:
            # Load data
            if isinstance(data, (str, Path)):
                df = pd.read_csv(data, sep='\t')
            elif isinstance(data, pd.DataFrame):
                df = data.copy()
            else:
                return None
            
            # Detect format (Bracken vs Kraken2)
            if 'fraction_total_reads' in df.columns:
                # Bracken format
                data_type = 'Bracken'
                df['abundance'] = df['fraction_total_reads']
                df['reads'] = df.get('new_est_reads', 0)
            else:
                # Kraken2 format
                df.columns = ['percentage', 'clade_reads', 'taxon_reads', 'rank', 'taxid', 'name']
                data_type = 'Kraken2'
                df['abundance'] = df['percentage'] / 100
                df['reads'] = df['clade_reads']
            
            # Filter significant taxa
            df_filtered = df[
                (df['abundance'] >= self.MIN_ABUNDANCE_THRESHOLD) &
                (~df['name'].str.contains('unclassified', case=False, na=False))
            ].copy()
            
            # Sort by abundance
            df_filtered = df_filtered.sort_values('abundance', ascending=False)
            
            # Calculate statistics
            total_classified = df_filtered['reads'].sum()
            species_count = len(df_filtered)
            
            parsed_data = {
                'type': data_type,
                'dataframe': df_filtered,
                'total_classified_reads': int(total_classified),
                'species_count': species_count,
                'top_taxa': df_filtered.head(20).to_dict('records'),
                'diversity_metrics': self._calculate_diversity_metrics(df_filtered['reads'].values)
            }
            
            self._diversity_data = parsed_data['diversity_metrics']
            
            print(f"  ✓ Parsed {data_type} data: {species_count} taxa, {total_classified:,} classified reads")
            
            return parsed_data
            
        except Exception as e:
            print(f"  ✗ Error parsing Bracken/Kraken data: {e}")
            return None
    
    def _parse_blast_data(self, data: Any) -> Optional[Dict]:
        """
        Parse BLAST taxonomy data into standardized format.
        
        Args:
            data: BLAST results file path, list, or dict
            
        Returns:
            Parsed data dictionary or None if parsing fails
        """
        try:
            # Load data
            if isinstance(data, (str, Path)):
                with open(data, 'r') as f:
                    blast_results = json.load(f)
            elif isinstance(data, (list, dict)):
                blast_results = data if isinstance(data, list) else [data]
            else:
                return None
            
            # Analyze results
            organism_counts = Counter()
            organism_sequences = {}
            identity_scores = []
            evalue_scores = []
            total_hits = 0
            sequences_with_hits = 0
            
            quality_dist = {
                'high': 0,    # >95% identity
                'medium': 0,  # 80-95% identity
                'low': 0      # <80% identity
            }
            
            for result in blast_results:
                if 'error' in result or not result.get('hits'):
                    continue
                
                sequences_with_hits += 1
                query_id = result.get('query_id', 'Unknown')
                
                for hit in result['hits']:
                    organism = hit.get('organism', 'Unknown')
                    identity = hit.get('identity', 0)
                    evalue = hit.get('e_value', 1.0)
                    
                    organism_counts[organism] += 1
                    total_hits += 1
                    identity_scores.append(identity)
                    evalue_scores.append(evalue)
                    
                    # Track sequences per organism
                    if organism not in organism_sequences:
                        organism_sequences[organism] = set()
                    organism_sequences[organism].add(query_id)
                    
                    # Quality categorization
                    if identity >= self.HIGH_QUALITY_IDENTITY:
                        quality_dist['high'] += 1
                    elif identity >= self.MEDIUM_QUALITY_IDENTITY:
                        quality_dist['medium'] += 1
                    else:
                        quality_dist['low'] += 1
            
            # Calculate statistics
            total_sequences = len(blast_results)
            classification_rate = (sequences_with_hits / total_sequences * 100) if total_sequences > 0 else 0
            avg_identity = np.mean(identity_scores) if identity_scores else 0
            avg_evalue = np.mean(evalue_scores) if evalue_scores else 1.0
            
            # Prepare top organisms
            top_organisms = []
            for organism, count in organism_counts.most_common(20):
                top_organisms.append({
                    'organism': organism,
                    'total_hits': count,
                    'unique_sequences': len(organism_sequences[organism]),
                    'hit_percentage': (count / total_hits * 100) if total_hits > 0 else 0
                })
            
            parsed_data = {
                'total_sequences': total_sequences,
                'sequences_classified': sequences_with_hits,
                'classification_rate': classification_rate,
                'total_hits': total_hits,
                'unique_organisms': len(organism_counts),
                'avg_identity': avg_identity,
                'avg_evalue': avg_evalue,
                'quality_distribution': quality_dist,
                'top_organisms': top_organisms,
                'identity_scores': identity_scores
            }
            
            self._quality_data = parsed_data
            
            print(f"  ✓ Parsed BLAST data: {sequences_with_hits}/{total_sequences} sequences classified "
                  f"({classification_rate:.1f}%)")
            
            return parsed_data
            
        except Exception as e:
            print(f"  ✗ Error parsing BLAST data: {e}")
            return None
    
    def _calculate_diversity_metrics(self, counts: np.ndarray) -> Dict[str, float]:
        """
        Calculate alpha diversity metrics.
        
        Args:
            counts: Array of read counts per taxon
            
        Returns:
            Dictionary of diversity metrics
        """
        counts = counts.astype(int)
        if counts.sum() == 0 or len(counts) == 0:
            return {
                'shannon': 0.0,
                'simpson': 0.0,
                'chao1': 0.0,
                'observed_species': 0,
                'evenness': 0.0,
                'dominance': 0.0
            }
        
        metrics = {}
        
        # Calculate using scikit-bio
        for metric_name in ['shannon', 'simpson', 'chao1', 'sobs']:
            try:
                result = alpha_diversity(metric_name, [counts], validate=False)
                metrics[metric_name if metric_name != 'sobs' else 'observed_species'] = float(result.iloc[0])
            except Exception:
                metrics[metric_name if metric_name != 'sobs' else 'observed_species'] = 0.0
        
        # Calculate evenness (Pielou's J)
        shannon = metrics.get('shannon', 0)
        sobs = metrics.get('observed_species', 0)
        metrics['evenness'] = (shannon / np.log(sobs)) if sobs > 1 else 0.0
        
        # Calculate dominance
        simpson = metrics.get('simpson', 0)
        metrics['dominance'] = 1 - simpson if isinstance(simpson, (int, float)) else 0.0
        
        return metrics
    
    def _generate_text_report(self, bracken_data: Optional[Dict], 
                             blast_data: Optional[Dict],
                             pipeline_config: Dict) -> List[str]:
        """Generate comprehensive text report."""
        content = []
        
        # Header
        content.extend(self._create_header(
            f"TAXONOMIC CLASSIFICATION REPORT - {pipeline_config['mode']} Analysis",
            pipeline_config['tool_info'],
            pipeline_config['database_info']
        ))
        
        # Executive Summary
        content.extend(self._create_section_header("EXECUTIVE SUMMARY", level=1))
        content.extend(self._generate_executive_summary(bracken_data, blast_data))
        
        # FASTQ Analysis Section
        if bracken_data:
            content.extend(self._create_section_header("FASTQ TAXONOMIC PROFILING", level=1))
            content.extend(self._generate_bracken_section(bracken_data))
            content.extend(self._generate_diversity_section(bracken_data))
        
        # FASTA Analysis Section
        if blast_data:
            content.extend(self._create_section_header("FASTA SEQUENCE CLASSIFICATION", level=1))
            content.extend(self._generate_blast_section(blast_data))
        
        # Quality Assessment
        content.extend(self._create_section_header("QUALITY ASSESSMENT", level=1))
        content.extend(self._generate_quality_section(bracken_data, blast_data))
        
        # Recommendations
        content.extend(self._create_section_header("RECOMMENDATIONS", level=1))
        content.extend(self._generate_recommendations(bracken_data, blast_data))
        
        # Footer
        content.extend(self._create_footer())
        
        return content
    
    def _generate_executive_summary(self, bracken_data: Optional[Dict], 
                                   blast_data: Optional[Dict]) -> List[str]:
        """Generate executive summary section."""
        summary_items = []
        
        if bracken_data:
            species_count = bracken_data['species_count']
            total_reads = bracken_data['total_classified_reads']
            top_taxon = bracken_data['top_taxa'][0] if bracken_data['top_taxa'] else None
            
            if top_taxon:
                top_abundance = top_taxon['abundance'] * 100
                diversity_metrics = bracken_data['diversity_metrics']
                shannon = diversity_metrics.get('shannon', 0)
                
                # Diversity interpretation
                if shannon < 1.0:
                    diversity_desc = "extremely low diversity with monodominant community structure"
                elif shannon < 2.0:
                    diversity_desc = "low diversity indicating selective environmental conditions"
                elif shannon < 3.0:
                    diversity_desc = "moderate diversity with balanced community composition"
                else:
                    diversity_desc = "high diversity suggesting complex ecosystem dynamics"
                
                summary_items.append(
                    f"Taxonomic profiling of metagenomic reads identified {species_count:,} distinct taxa "
                    f"from {total_reads:,} classified sequences. The microbial community exhibits {diversity_desc}. "
                    f"{top_taxon['name']} dominates the community at {top_abundance:.2f}% relative abundance."
                )
        
        if blast_data:
            classification_rate = blast_data['classification_rate']
            unique_organisms = blast_data['unique_organisms']
            avg_identity = blast_data['avg_identity']
            
            quality_assessment = "high" if classification_rate > 80 else "moderate" if classification_rate > 60 else "limited"
            
            summary_items.append(
                f"BLAST-based taxonomic assignment achieved {classification_rate:.1f}% classification success, "
                f"identifying {unique_organisms:,} unique organisms with {avg_identity:.1f}% average sequence identity. "
                f"Classification quality is assessed as {quality_assessment}."
            )
        
        if not summary_items:
            summary_items.append("No taxonomic classification data available for analysis.")
        
        return self._create_summary_box("Key Findings", summary_items)
    
    def _generate_bracken_section(self, bracken_data: Dict) -> List[str]:
        """Generate Bracken/Kraken analysis section."""
        lines = []
        
        # Overview statistics
        lines.extend(self._create_section_header("Classification Overview", level=2))
        lines.extend([
            f"Classification Method: {bracken_data['type']}",
            f"Total Classified Reads: {bracken_data['total_classified_reads']:,}",
            f"Taxa Detected (>0.01% abundance): {bracken_data['species_count']:,}",
            ""
        ])
        
        # Top taxa table
        lines.extend(self._create_section_header("Top 15 Most Abundant Taxa", level=2))
        
        headers = ["Rank", "Taxon", "Abundance", "Reads"]
        rows = []
        
        for i, taxon in enumerate(bracken_data['top_taxa'][:15], 1):
            rows.append([
                str(i),
                taxon['name'][:60],
                f"{taxon['abundance']*100:.3f}%",
                f"{int(taxon['reads']):,}"
            ])
        
        lines.extend(self._format_table(headers, rows, alignments=['left', 'left', 'right', 'right']))
        lines.append("")
        
        return lines
    
    def _generate_diversity_section(self, bracken_data: Dict) -> List[str]:
        """Generate diversity analysis section."""
        lines = self._create_section_header("Alpha Diversity Analysis", level=2)
        
        metrics = bracken_data['diversity_metrics']
        
        # Format metrics
        stats = {
            "Shannon Diversity Index (H')": f"{metrics['shannon']:.4f}",
            "Simpson's Index (1-D)": f"{metrics['simpson']:.4f}",
            "Observed Species": f"{int(metrics['observed_species']):,}",
            "Chao1 Richness Estimate": f"{metrics['chao1']:.1f}",
            "Pielou's Evenness (J)": f"{metrics['evenness']:.4f}",
            "Community Dominance": f"{metrics['dominance']:.4f}"
        }
        
        lines.extend(self._format_statistics(stats, "Diversity Metrics"))
        lines.append("")
        
        # Ecological interpretation
        lines.extend(self._create_section_header("Ecological Interpretation", level=3))
        
        shannon = metrics['shannon']
        evenness = metrics['evenness']
        
        if shannon < 1.0:
            interpretation = ("The community exhibits extremely low diversity (H' < 1.0), characteristic of "
                            "disturbed environments, selective culture conditions, or monoclonal dominance. "
                            "This pattern may indicate environmental stress or specialized ecological niches.")
        elif shannon < 2.0:
            interpretation = ("Moderate-low diversity (1.0 ≤ H' < 2.0) suggests limited taxonomic richness, "
                            "potentially due to selective environmental pressures or specialized metabolic requirements.")
        elif shannon < 3.0:
            interpretation = ("Moderate diversity (2.0 ≤ H' < 3.0) indicates a balanced microbial community "
                            "with reasonable taxonomic representation across multiple functional groups.")
        else:
            interpretation = ("High diversity (H' ≥ 3.0) reflects a complex, mature ecosystem with extensive "
                            "taxonomic richness and niche specialization, typical of stable environmental conditions.")
        
        lines.extend([interpretation, ""])
        
        if evenness < 0.3:
            lines.extend([
                "Low evenness (J < 0.3) indicates strong dominance by one or few taxa, suggesting "
                "competitive exclusion or environmental selection favoring specific organisms.",
                ""
            ])
        elif evenness > 0.7:
            lines.extend([
                "High evenness (J > 0.7) indicates relatively equal distribution of taxa, characteristic "
                "of stable ecosystems with balanced resource competition.",
                ""
            ])
        
        return lines
    
    def _generate_blast_section(self, blast_data: Dict) -> List[str]:
        """Generate BLAST analysis section."""
        lines = []
        
        # Classification statistics
        lines.extend(self._create_section_header("Classification Statistics", level=2))
        lines.extend([
            f"Total Sequences Analyzed: {blast_data['total_sequences']:,}",
            f"Successfully Classified: {blast_data['sequences_classified']:,} ({blast_data['classification_rate']:.1f}%)",
            f"Total BLAST Alignments: {blast_data['total_hits']:,}",
            f"Unique Organisms Identified: {blast_data['unique_organisms']:,}",
            f"Average Hits per Sequence: {blast_data['total_hits']/max(blast_data['sequences_classified'],1):.1f}",
            ""
        ])
        
        # Quality metrics
        lines.extend(self._create_section_header("Alignment Quality Metrics", level=2))
        
        identity_scores = blast_data.get('identity_scores', [])
        if identity_scores:
            lines.extend([
                f"Average Sequence Identity: {blast_data['avg_identity']:.2f}% ± {np.std(identity_scores):.2f}%",
                f"Identity Range: {min(identity_scores):.1f}% - {max(identity_scores):.1f}%",
                f"Median Identity: {np.median(identity_scores):.1f}%",
                f"Average E-value: {blast_data['avg_evalue']:.2e}",
                ""
            ])
        
        # Quality distribution
        quality_dist = blast_data['quality_distribution']
        total_hits = blast_data['total_hits']
        
        lines.extend(self._create_section_header("Quality Distribution", level=3))
        lines.extend([
            f"High Confidence (≥95% identity): {quality_dist['high']:,} alignments ({quality_dist['high']/max(total_hits,1)*100:.1f}%)",
            f"Medium Confidence (80-94% identity): {quality_dist['medium']:,} alignments ({quality_dist['medium']/max(total_hits,1)*100:.1f}%)",
            f"Low Confidence (<80% identity): {quality_dist['low']:,} alignments ({quality_dist['low']/max(total_hits,1)*100:.1f}%)",
            ""
        ])
        
        # Top organisms table
        lines.extend(self._create_section_header("Top 15 Organisms by Hit Count", level=2))
        
        headers = ["Rank", "Organism", "Total Hits", "Sequences", "Hit %"]
        rows = []
        
        for i, org in enumerate(blast_data['top_organisms'][:15], 1):
            rows.append([
                str(i),
                org['organism'][:50],
                f"{org['total_hits']:,}",
                f"{org['unique_sequences']:,}",
                f"{org['hit_percentage']:.2f}%"
            ])
        
        lines.extend(self._format_table(headers, rows, alignments=['left', 'left', 'right', 'right', 'right']))
        lines.append("")
        
        return lines
    
    def _generate_quality_section(self, bracken_data: Optional[Dict], 
                                  blast_data: Optional[Dict]) -> List[str]:
        """Generate quality assessment section."""
        lines = []
        
        if bracken_data:
            lines.extend(self._create_section_header("FASTQ Classification Quality", level=2))
            lines.extend([
                "Database: Kraken2 Standard-8 (16GB reference database)",
                "Confidence Threshold: 0.0 (default, all assignments reported)",
                f"Taxa Detected: {bracken_data['species_count']:,}",
                f"Classification Depth: Species-level resolution",
                "",
                "Assessment: Kraken2/Bracken provides rapid, accurate taxonomic profiling with high sensitivity "
                "for known organisms. Results are suitable for community composition analysis and diversity studies.",
                ""
            ])
        
        if blast_data:
            lines.extend(self._create_section_header("FASTA Classification Quality", level=2))
            
            classification_rate = blast_data['classification_rate']
            avg_identity = blast_data['avg_identity']
            high_quality_pct = blast_data['quality_distribution']['high'] / max(blast_data['total_hits'], 1) * 100
            
            # Quality interpretation
            if classification_rate > 80 and avg_identity > 90:
                quality = "Excellent"
                interpretation = ("High-quality sequences with comprehensive database coverage. "
                                "Results are highly reliable for taxonomic assignments.")
            elif classification_rate > 60 and avg_identity > 85:
                quality = "Good"
                interpretation = ("Acceptable classification rate with good sequence identity. "
                                "Results are suitable for most downstream analyses.")
            elif classification_rate > 40:
                quality = "Moderate"
                interpretation = ("Moderate classification success. Consider data quality review "
                                "and potential database limitations.")
            else:
                quality = "Limited"
                interpretation = ("Low classification rate may indicate novel sequences, poor data quality, "
                                "or limited database representation.")
            
            lines.extend([
                f"Classification Success Rate: {classification_rate:.1f}%",
                f"Average Sequence Identity: {avg_identity:.1f}%",
                f"High-Quality Alignments (>95% ID): {high_quality_pct:.1f}%",
                f"Overall Quality Assessment: {quality}",
                "",
                interpretation,
                ""
            ])
        
        # Combined assessment for dual-method
        if bracken_data and blast_data:
            lines.extend(self._create_section_header("Integrated Quality Assessment", level=2))
            lines.extend([
                "Dual-method validation provides robust taxonomic classification through complementary approaches:",
                "",
                "• Kraken2/Bracken: Rapid k-mer based profiling for community composition",
                "• BLAST: Alignment-based validation for individual sequence classification",
                "",
                "Cross-method concordance enhances confidence in taxonomic assignments and enables "
                "detection of potential classification discrepancies.",
                ""
            ])
        
        return lines
    
    def _generate_recommendations(self, bracken_data: Optional[Dict], 
                                  blast_data: Optional[Dict]) -> List[str]:
        """Generate recommendations section."""
        lines = []
        recommendations = []
        
        # Data-driven recommendations
        if bracken_data:
            diversity_metrics = bracken_data['diversity_metrics']
            shannon = diversity_metrics['shannon']
            
            if shannon < 1.0:
                recommendations.append(
                    "Extremely low diversity detected. Verify sample collection and storage protocols. "
                    "Consider contamination screening and validate dominant taxa through alternative methods."
                )
            
            if bracken_data['species_count'] < 10:
                recommendations.append(
                    "Limited taxonomic richness observed. Evaluate sequencing depth sufficiency and "
                    "consider increasing read coverage for rare taxa detection."
                )
        
        if blast_data:
            classification_rate = blast_data['classification_rate']
            avg_identity = blast_data['avg_identity']
            
            if classification_rate < 60:
                recommendations.append(
                    "Low classification rate suggests potential novel sequences or database gaps. "
                    "Consider: (1) Alternative reference databases, (2) De novo assembly for novel content, "
                    "(3) Phylogenetic placement for unclassified sequences."
                )
            
            if avg_identity < 85:
                recommendations.append(
                    "Lower average sequence identity may indicate distant taxonomic relationships or "
                    "sequence quality issues. Review quality filtering parameters and consider "
                    "specialized databases for your sample type."
                )
        
        # General best practices
        recommendations.extend([
            "Cross-reference taxonomic assignments with expected sample characteristics and metadata.",
            "Validate critical findings through orthogonal methods (qPCR, culture, amplicon sequencing).",
            "Compare community profiles with appropriate controls and reference datasets.",
            "Consider functional profiling to complement taxonomic characterization.",
            "Document all bioinformatics parameters for reproducibility."
        ])
        
        # Format recommendations
        for i, rec in enumerate(recommendations, 1):
            lines.append(f"{i}. {rec}")
            lines.append("")
        
        return lines
    
    def _generate_json_report(self, bracken_data: Optional[Dict], 
                             blast_data: Optional[Dict],
                             pipeline_config: Dict) -> Dict:
        """Generate structured JSON report."""
        json_data = {
            'summary': {
                'pipeline_mode': pipeline_config['mode'],
                'analysis_timestamp': self.timestamp.isoformat(),
                'total_taxa_detected': 0,
                'classification_methods': []
            },
            'fastq_analysis': None,
            'fasta_analysis': None,
            'quality_metrics': {}
        }
        
        if bracken_data:
            json_data['summary']['classification_methods'].append('Kraken2/Bracken')
            json_data['summary']['total_taxa_detected'] += bracken_data['species_count']
            
            json_data['fastq_analysis'] = {
                'method': bracken_data['type'],
                'total_classified_reads': bracken_data['total_classified_reads'],
                'species_count': bracken_data['species_count'],
                'diversity_metrics': bracken_data['diversity_metrics'],
                'top_taxa': bracken_data['top_taxa'][:20]
            }
        
        if blast_data:
            json_data['summary']['classification_methods'].append('BLAST')
            json_data['summary']['total_taxa_detected'] += blast_data['unique_organisms']
            
            json_data['fasta_analysis'] = {
                'total_sequences': blast_data['total_sequences'],
                'sequences_classified': blast_data['sequences_classified'],
                'classification_rate': round(blast_data['classification_rate'], 2),
                'average_identity': round(blast_data['avg_identity'], 2),
                'unique_organisms': blast_data['unique_organisms'],
                'quality_distribution': blast_data['quality_distribution'],
                'top_organisms': blast_data['top_organisms'][:20]
            }
            
            json_data['quality_metrics']['blast'] = {
                'classification_success_rate': blast_data['classification_rate'],
                'average_sequence_identity': blast_data['avg_identity'],
                'high_quality_alignments_percent': (
                    blast_data['quality_distribution']['high'] / 
                    max(blast_data['total_hits'], 1) * 100
                )
            }
        
        return json_data