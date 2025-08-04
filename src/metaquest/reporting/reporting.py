#!/usr/bin/env python3
"""
MetaQuest Reporting Module
Streamlined reporting with redundancy eliminated
Keeps only the best pathogen report: pathogen_summary.txt
"""

import json
import pandas as pd
from pathlib import Path
from datetime import datetime
from abc import ABC, abstractmethod
from typing import Dict, List, Optional, Any
import numpy as np
from collections import Counter
import re
import time
from ..config import *

class BaseReportGenerator(ABC):
    """Base class for all MetaQuest report generators"""
    
    def __init__(self, output_dir: Path, analysis_type: str):
        self.output_dir = Path(output_dir)
        self.analysis_type = analysis_type
        self.timestamp = datetime.now()
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    @abstractmethod
    def generate_report(self, data: Any, **kwargs) -> Dict[str, str]:
        """Generate the specific report type"""
        pass
    
    def _create_report_header(self, title: str) -> List[str]:
        """Standard header for all text reports"""
        return [
            title,
            "=" * len(title),
            "",
            f"Analysis Type: {self.analysis_type}",
            f"Generated: {self.timestamp.strftime('%Y-%m-%d %H:%M:%S')}",
            f"MetaQuest Pipeline v3.0",
            ""
        ]
    
    def _save_text_report(self, content: List[str], filename: str) -> str:
        """Save text report with standard formatting"""
        filepath = self.output_dir / filename
        with open(filepath, 'w') as f:
            f.write('\n'.join(content))
        return str(filepath)
    
    def _save_json_report(self, data: Dict, filename: str) -> str:
        """Save JSON report with metadata"""
        report_data = {
            'metadata': {
                'analysis_type': self.analysis_type,
                'generation_timestamp': self.timestamp.isoformat(),
                'pipeline_version': 'MetaQuest v3.0'
            },
            'data': data
        }
        filepath = self.output_dir / filename
        with open(filepath, 'w') as f:
            json.dump(report_data, f, indent=2)
        return str(filepath)

class TaxonomicReporter(BaseReportGenerator):
    """Enhanced taxonomic classification reporting with comprehensive BLAST analysis"""
    
    def __init__(self, output_dir: Path):
        super().__init__(output_dir, "Taxonomic Classification")
    
    def generate_report(self, bracken_data=None, blast_data=None, **kwargs) -> Dict[str, str]:
        """Generate comprehensive taxonomic classification report"""
        generated_files = {}
        
        # Generate text report
        text_content = self._create_report_header("TAXONOMIC CLASSIFICATION REPORT")
        
        if bracken_data:
            text_content.extend(self._process_bracken_data(bracken_data))
        
        if blast_data:
            text_content.extend(self._process_blast_data_comprehensive(blast_data))
        
        # Summary statistics
        text_content.extend(self._generate_taxonomy_summary(bracken_data, blast_data))
        
        # Save reports
        generated_files['text_report'] = self._save_text_report(text_content, "taxonomic_classification_report.txt")
        
        # Generate JSON report
        json_data = self._create_taxonomy_json_data(bracken_data, blast_data)
        generated_files['json_report'] = self._save_json_report(json_data, "taxonomic_classification_report.json")
        
        # Generate additional detailed BLAST reports if BLAST data available
        if blast_data:
            generated_files.update(self._generate_detailed_blast_reports(blast_data))
        
        return generated_files
    
    def _process_bracken_data(self, bracken_data) -> List[str]:
        """Process both Bracken and Kraken data formats"""
        try:
            if isinstance(bracken_data, (str, Path)):
                df = pd.read_csv(bracken_data, sep='\t')
                
                # Detect format - Bracken vs Kraken
                if 'fraction_total_reads' in df.columns:
                    # Bracken format
                    abundance_col = 'fraction_total_reads'
                    name_col = 'name'
                    count_col = 'new_est_reads'
                    data_type = "Bracken"
                else:
                    # Kraken format
                    df = pd.read_csv(bracken_data, sep='\t', header=None,
                                   names=['percentage', 'clade_reads', 'taxon_reads', 'rank', 'taxid', 'name'])
                    abundance_col = 'percentage'
                    name_col = 'name'
                    count_col = 'clade_reads'
                    data_type = "Kraken2"
                    # Convert percentage to fraction for consistency
                    df['fraction_total_reads'] = df['percentage'] / 100
            else:
                df = bracken_data
                abundance_col = 'fraction_total_reads'
                name_col = 'name'
                count_col = 'new_est_reads'
                data_type = "Bracken"
            
            lines = [
                f"FASTQ TAXONOMIC CLASSIFICATION ({data_type})",
                "-" * 50,
                f"Total classified taxa: {len(df)}",
                ""
            ]
            
            # Filter significant taxa and get top 15
            if data_type == "Kraken2":
                significant_taxa = df[(df['percentage'] > 0.01) & 
                                    (~df[name_col].str.contains('unclassified', case=False, na=False))].copy()
                abundance_values = significant_taxa['percentage']
                abundance_unit = "%"
            else:
                significant_taxa = df[(df[abundance_col] > 0.0001) & 
                                    (~df[name_col].str.contains('unclassified', case=False, na=False))].copy()
                abundance_values = significant_taxa[abundance_col] * 100
                abundance_unit = "%"
            
            if len(significant_taxa) > 0:
                top_taxa = significant_taxa.nlargest(15, abundance_col if data_type == "Bracken" else 'percentage')
                lines.extend([
                    "TOP 15 MOST ABUNDANT TAXA:",
                    f"{'Rank':<4} {'Taxon':<45} {'Abundance':<12} {'Reads':<10}",
                    "-" * 75
                ])
                
                for i, (_, row) in enumerate(top_taxa.iterrows(), 1):
                    taxon_name = row[name_col][:44] if len(row[name_col]) > 44 else row[name_col]
                    if data_type == "Kraken2":
                        abundance_val = row['percentage']
                        read_count = row[count_col]
                    else:
                        abundance_val = row[abundance_col] * 100
                        read_count = row.get(count_col, 0)
                    
                    lines.append(f"{i:<4} {taxon_name:<45} {abundance_val:>8.3f}{abundance_unit}  {read_count:>8}")
                
                # Add diversity metrics
                lines.extend([
                    "",
                    "DIVERSITY METRICS:",
                    f"• Shannon diversity index (proxy): {self._calculate_shannon_diversity(abundance_values):.3f}",
                    f"• Dominant taxon: {top_taxa.iloc[0][name_col]} ({abundance_values.iloc[0]:.3f}{abundance_unit})",
                    f"• Taxa >1% abundance: {len(abundance_values[abundance_values > 1])}"
                ])
            else:
                lines.append("No significant taxa found above threshold")
            
            lines.append("")
            return lines
            
        except Exception as e:
            return [f"Error processing Bracken/Kraken data: {e}", ""]
    
    def _calculate_shannon_diversity(self, abundances):
        """Calculate Shannon diversity index"""
        abundances = np.array(abundances)
        abundances = abundances[abundances > 0]  # Remove zeros
        proportions = abundances / abundances.sum()
        return -np.sum(proportions * np.log(proportions))
    
    def _process_blast_data_comprehensive(self, blast_data) -> List[str]:
        """Enhanced BLAST processing with comprehensive organism analysis"""
        try:
            if isinstance(blast_data, (str, Path)):
                with open(blast_data, 'r') as f:
                    blast_results = json.load(f)
            else:
                blast_results = blast_data
            
            lines = [
                "FASTA TAXONOMIC CLASSIFICATION (BLAST) - COMPREHENSIVE ANALYSIS",
                "-" * 65,
                f"Total sequences analyzed: {len(blast_results)}",
                ""
            ]
            
            # Enhanced organism counting and analysis
            organism_counts = Counter()  # Count ALL hits per organism
            organism_sequences = {}
            total_hits = 0
            sequences_with_hits = 0
            identity_scores = []
            quality_metrics = {
                'high_identity': 0,    # >90%
                'medium_identity': 0,  # 75-90%
                'low_identity': 0,     # <75%
                'high_quality_eval': 0, # <1e-50
                'medium_quality_eval': 0, # 1e-50 to 1e-5
                'low_quality_eval': 0   # >1e-5
            }
            
            # Analyze results - count ALL hits to match comprehensive analysis
            for result in blast_results:
                if 'error' in result:
                    continue
                
                if result['hits']:
                    sequences_with_hits += 1
                    
                    for hit in result['hits']:
                        organism = hit.get('organism', 'Unknown')
                        organism_counts[organism] += 1  # Count every hit
                        total_hits += 1
                        
                        # Track which sequences hit each organism
                        if organism not in organism_sequences:
                            organism_sequences[organism] = set()
                        organism_sequences[organism].add(result['query_id'])
                        
                        # Collect detailed quality metrics
                        identity = hit.get('identity', 0)
                        e_value = hit.get('e_value', 1.0)
                        
                        identity_scores.append(identity)
                        
                        # Identity quality classification
                        if identity > 90:
                            quality_metrics['high_identity'] += 1
                        elif identity >= 75:
                            quality_metrics['medium_identity'] += 1
                        else:
                            quality_metrics['low_identity'] += 1
                        
                        # E-value quality classification
                        if e_value < 1e-50:
                            quality_metrics['high_quality_eval'] += 1
                        elif e_value <= 1e-5:
                            quality_metrics['medium_quality_eval'] += 1
                        else:
                            quality_metrics['low_quality_eval'] += 1
            
            lines.extend([
                f"Sequences with hits: {sequences_with_hits}",
                f"Total BLAST hits: {total_hits}",
                f"Unique organisms identified: {len(organism_counts)}",
                f"Classification rate: {sequences_with_hits/len(blast_results)*100:.1f}%",
                f"Average hits per classified sequence: {total_hits/sequences_with_hits:.1f}" if sequences_with_hits > 0 else "",
                ""
            ])
            
            # Comprehensive quality metrics
            if identity_scores:
                lines.extend([
                    "DETAILED QUALITY METRICS:",
                    "-" * 30,
                    f"• Average sequence identity: {np.mean(identity_scores):.1f}%",
                    f"• Median sequence identity: {np.median(identity_scores):.1f}%",
                    f"• Identity range: {min(identity_scores):.1f}% - {max(identity_scores):.1f}%",
                    "",
                    "IDENTITY DISTRIBUTION:",
                    f"• High-quality hits (>90% identity): {quality_metrics['high_identity']} ({quality_metrics['high_identity']/total_hits*100:.1f}%)",
                    f"• Medium-quality hits (75-90% identity): {quality_metrics['medium_identity']} ({quality_metrics['medium_identity']/total_hits*100:.1f}%)",
                    f"• Low-quality hits (<75% identity): {quality_metrics['low_identity']} ({quality_metrics['low_identity']/total_hits*100:.1f}%)",
                    "",
                    "E-VALUE DISTRIBUTION:",
                    f"• High confidence (E<1e-50): {quality_metrics['high_quality_eval']} ({quality_metrics['high_quality_eval']/total_hits*100:.1f}%)",
                    f"• Medium confidence (1e-50≤E≤1e-5): {quality_metrics['medium_quality_eval']} ({quality_metrics['medium_quality_eval']/total_hits*100:.1f}%)",
                    f"• Low confidence (E>1e-5): {quality_metrics['low_quality_eval']} ({quality_metrics['low_quality_eval']/total_hits*100:.1f}%)",
                    ""
                ])
            
            # Comprehensive organism analysis
            if organism_counts:
                lines.extend([
                    "TOP 20 ORGANISMS BY TOTAL HITS:",
                    "-" * 80,
                    f"{'Rank':<4} {'Organism':<40} {'Total Hits':<10} {'Sequences':<10} {'Avg Hits/Seq':<12} {'Hit %':<8}",
                    "-" * 80
                ])
                
                for i, (organism, hit_count) in enumerate(organism_counts.most_common(20), 1):
                    seq_count = len(organism_sequences.get(organism, set()))
                    avg_hits = hit_count / seq_count if seq_count > 0 else 0
                    hit_percentage = (hit_count / total_hits) * 100
                    
                    org_name = organism[:39] if len(organism) > 39 else organism
                    lines.append(f"{i:<4} {org_name:<40} {hit_count:<10} {seq_count:<10} {avg_hits:<12.1f} {hit_percentage:<8.1f}%")
                
                lines.extend([
                    "",
                    "ORGANISM DIVERSITY ANALYSIS:",
                    f"• Total unique organisms: {len(organism_counts)}",
                    f"• Dominant organism: {organism_counts.most_common(1)[0][0]} ({organism_counts.most_common(1)[0][1]} hits)",
                    f"• Organisms with >5 hits: {len([c for c in organism_counts.values() if c > 5])}",
                    f"• Organisms with single hits: {len([c for c in organism_counts.values() if c == 1])}",
                    ""
                ])
            
            # Store data for additional file generation
            self._blast_organism_counts = organism_counts
            self._blast_organism_sequences = organism_sequences
            self._blast_quality_metrics = quality_metrics
            
            lines.append("")
            return lines
            
        except Exception as e:
            return [f"Error processing BLAST data: {e}", ""]
    
    def _generate_taxonomy_summary(self, bracken_data, blast_data) -> List[str]:
        """Generate comprehensive taxonomy analysis summary"""
        lines = [
            "ANALYSIS SUMMARY & INTERPRETATION",
            "-" * 35,
            ""
        ]
        
        # Analysis completeness
        methods_used = []
        if bracken_data:
            methods_used.append("FASTQ classification (Kraken2/Bracken)")
        if blast_data:
            methods_used.append("FASTA classification (BLAST)")
        
        if methods_used:
            lines.extend([
                "✓ COMPLETED ANALYSES:",
                f"  {' + '.join(methods_used)}",
                ""
            ])
        
        # Key insights and recommendations
        lines.extend([
            "KEY INSIGHTS:",
            "• Review dominant taxa for sample composition validation",
            "• Check for unexpected organisms indicating contamination",
            "• High-identity matches (>90%) are most reliable",
            "• Multiple detection methods increase confidence",
            "",
            "QUALITY ASSESSMENT:",
            "• Rare taxa (<0.1% abundance) require careful interpretation",
            "• Cross-reference with expected sample characteristics",
            "• Consider diversity metrics for ecosystem health",
            "",
            "RECOMMENDED FOLLOW-UP:",
            "• Validate unexpected findings with targeted PCR",
            "• Compare with negative controls if available",
            "• Consider functional analysis for metabolic insights",
            ""
        ])
        
        return lines
    
    def _generate_detailed_blast_reports(self, blast_data) -> Dict[str, str]:
        """Generate additional detailed BLAST reports and data files"""
        generated_files = {}
        
        try:
            if isinstance(blast_data, (str, Path)):
                with open(blast_data, 'r') as f:
                    blast_results = json.load(f)
            else:
                blast_results = blast_data

            # Generate organism comparison data
            comparison_file = self._create_organism_comparison_data()
            if comparison_file:
                generated_files['organism_comparison'] = comparison_file
            
        except Exception as e:
            print(f"Warning: Could not generate detailed BLAST reports: {e}")
        
        return generated_files
    
    def _create_organism_comparison_data(self) -> str:
        """Create organism comparison data files"""
        try:
            organism_counts = getattr(self, '_blast_organism_counts', Counter())
            organism_sequences = getattr(self, '_blast_organism_sequences', {})
            
            if not organism_counts:
                return None
            
            comparison_data = []
            
            # Use top 20 organisms for comparison data
            for organism, hit_count in organism_counts.most_common(20):
                seq_count = len(organism_sequences.get(organism, set()))
                avg_hits = hit_count / seq_count if seq_count > 0 else 0
                
                comparison_data.append({
                    'organism': organism,
                    'total_hits': hit_count,
                    'sequences_with_hits': seq_count,
                    'avg_hits_per_sequence': avg_hits,
                    'hit_percentage': (hit_count / sum(organism_counts.values())) * 100
                })
            
            # Save as JSON
            json_file = self.output_dir / "organism_comparison_data.json"
            with open(json_file, 'w') as f:
                json.dump(comparison_data, f, indent=2)
            
            # Save as CSV
            csv_file = self.output_dir / "organism_comparison_data.csv"
            df = pd.DataFrame(comparison_data)
            df.to_csv(csv_file, index=False)
            
            return str(json_file)
            
        except Exception as e:
            print(f"Warning: Could not create organism comparison data: {e}")
            return None
    
    def _create_taxonomy_json_data(self, bracken_data, blast_data) -> Dict:
        """Create structured JSON data for taxonomy report"""
        json_data = {
            'bracken_analysis': None,
            'blast_analysis': None,
            'summary_statistics': {}
        }
        
        # Process Bracken/Kraken data
        if bracken_data:
            try:
                if isinstance(bracken_data, (str, Path)):
                    df = pd.read_csv(bracken_data, sep='\t')
                    if 'fraction_total_reads' not in df.columns:
                        # Kraken format
                        df = pd.read_csv(bracken_data, sep='\t', header=None,
                                       names=['percentage', 'clade_reads', 'taxon_reads', 'rank', 'taxid', 'name'])
                        df['fraction_total_reads'] = df['percentage'] / 100
                else:
                    df = bracken_data
                
                json_data['bracken_analysis'] = {
                    'total_taxa': len(df),
                    'top_15_taxa': df.nlargest(15, 'fraction_total_reads')[['name', 'fraction_total_reads']].to_dict('records'),
                    'diversity_metrics': {
                        'shannon_index': self._calculate_shannon_diversity(df['fraction_total_reads'].values),
                        'dominant_taxon_percentage': df['fraction_total_reads'].max() * 100
                    }
                }
            except Exception as e:
                json_data['bracken_analysis'] = {'error': str(e)}
        
        # Process BLAST data
        if blast_data:
            try:
                if isinstance(blast_data, (str, Path)):
                    with open(blast_data, 'r') as f:
                        blast_results = json.load(f)
                else:
                    blast_results = blast_data
                
                organism_counts = Counter()
                identity_scores = []
                for result in blast_results:
                    if 'error' in result or not result.get('hits'):
                        continue
                    for hit in result['hits']:
                        organism = hit.get('organism', 'Unknown')
                        organism_counts[organism] += 1
                        if 'identity' in hit:
                            identity_scores.append(hit['identity'])
                
                json_data['blast_analysis'] = {
                    'total_sequences': len(blast_results),
                    'sequences_with_hits': len([r for r in blast_results if not r.get('error') and r.get('hits')]),
                    'unique_organisms': len(organism_counts),
                    'average_identity': np.mean(identity_scores) if identity_scores else 0,
                    'top_15_organisms': [{'organism': k, 'hit_count': v} for k, v in organism_counts.most_common(15)]
                }
            except Exception as e:
                json_data['blast_analysis'] = {'error': str(e)}
        
        return json_data

class PathogenReporter(BaseReportGenerator):
    """Streamlined pathogen detection reporting - ONLY pathogen_summary.txt + JSON for visualizations"""
    
    def __init__(self, output_dir: Path):
        super().__init__(output_dir, "Pathogen Detection")
    
    def _extract_pathogens_from_bracken(self, bracken_file):
        """Extract pathogenic organisms directly from Bracken data with lowered threshold"""
        try:
            import pandas as pd
            
            bracken_pathogens = []
            
            # Enhanced pathogenic keywords with broader matching
            pathogenic_keywords = [
                'salmonella', 'escherichia coli', 'staphylococcus aureus', 'clostridium tetani',
                'klebsiella pneumoniae', 'yersinia enterocolitica', 
                'streptococcus',  # BROADENED: now catches all Streptococcus species
                'enterococcus faecalis', 'pseudomonas aeruginosa', 'acinetobacter baumannii',
                'vibrio', 'brucella', 'mycobacterium tuberculosis', 'bacillus anthracis',
                'listeria monocytogenes', 'clostridium difficile', 'yersinia pestis',
                'francisella tularensis', 'burkholderia', 'rickettsia', 'coxiella'
            ]
            
            # Read the Bracken file with proper headers
            df = pd.read_csv(bracken_file, sep='\t')
            
            # Filter for species level with LOWERED threshold
            species_df = df[
                (df['taxonomy_lvl'] == 'S') & 
                (df['fraction_total_reads'] > 0.00005)  # LOWERED from 0.0001 to catch 0.05% pathogens
            ]
            
            for _, row in species_df.iterrows():
                organism = row['name'].strip()
                organism_lower = organism.lower()
                
                # Check if organism matches pathogenic keywords
                if any(keyword in organism_lower for keyword in pathogenic_keywords):
                    bracken_pathogens.append({
                        'organism': organism,
                        'abundance': row['fraction_total_reads'],  # Already in fraction
                        'reads': int(row['new_est_reads']),
                        'detection_method': 'bracken',
                        'taxonomy_id': row.get('taxonomy_id', None)
                    })
            
            print(f"   ✓ Extracted {len(bracken_pathogens)} pathogenic organisms from Bracken data")
            if len(bracken_pathogens) > 0:
                print(f"   📋 Detected pathogens: {[p['organism'] for p in bracken_pathogens[:5]]}")
            return bracken_pathogens
            
        except Exception as e:
            print(f"Error extracting pathogens from Bracken data: {e}")
            return []
        
    def _extract_pathogens_from_blast_taxonomy(self, blast_results_data):
        """Extract pathogenic organisms from BLAST taxonomy results for FASTA analysis"""
        try:
            taxonomy_pathogens = []
            
            # Enhanced pathogenic keywords (same as Bracken extraction)
            pathogenic_keywords = [
                'salmonella', 'escherichia coli', 'staphylococcus aureus', 'clostridium tetani',
                'klebsiella pneumoniae', 'yersinia enterocolitica', 
                'streptococcus', 'enterococcus faecalis', 'pseudomonas aeruginosa', 
                'acinetobacter baumannii', 'vibrio', 'brucella', 'mycobacterium tuberculosis', 
                'bacillus anthracis', 'listeria monocytogenes', 'clostridium difficile'
            ]
            
            # Process BLAST results to find pathogenic organisms
            if isinstance(blast_results_data, (str, Path)):
                with open(blast_results_data, 'r') as f:
                    blast_results = json.load(f)
            else:
                blast_results = blast_results_data
            
            for result in blast_results:
                if 'error' in result or not result.get('hits'):
                    continue
                    
                for hit in result['hits']:
                    organism = hit.get('organism', 'Unknown')
                    organism_lower = organism.lower()
                    
                    # Check if organism matches pathogenic keywords
                    if any(keyword in organism_lower for keyword in pathogenic_keywords):
                        taxonomy_pathogens.append({
                            'organism': organism,
                            'identity': hit.get('identity', 0),
                            'evalue': hit.get('e_value', 1.0),
                            'detection_method': 'taxonomy',
                            'query_id': result.get('query_id', 'Unknown')
                        })
            
            print(f"   ✓ Extracted {len(taxonomy_pathogens)} pathogenic organisms from BLAST taxonomy")
            if len(taxonomy_pathogens) > 0:
                unique_organisms = list(set([p['organism'] for p in taxonomy_pathogens]))
                print(f"   📋 Detected pathogens: {unique_organisms[:5]}")
            
            return taxonomy_pathogens
            
        except Exception as e:
            print(f"Error extracting pathogens from BLAST taxonomy: {e}")
            return []

    
    def generate_comprehensive_pathogen_report(self, output_dir, bracken_pathogens=None, taxonomy_pathogens=None, sequence_pathogens=None):
        """Generate ONLY pathogen_summary.txt + JSON - eliminate redundant reports"""
        print("Generating comprehensive pathogen report...")
        
        # Initialize pathogen data structure
        all_pathogens = {}
        
        # Process Bracken pathogens (FASTQ taxonomic classification)
        if bracken_pathogens:
            for pathogen in bracken_pathogens:
                org_key = pathogen['organism'].strip()
                if org_key not in all_pathogens:
                    all_pathogens[org_key] = {
                        'organism': org_key,
                        'detection_methods': [],
                        'risk_level': 'Unknown',
                        'abundance_percentage': 0,
                        'estimated_reads': 0,
                        'sequence_hits': 0,
                        'max_identity': 0,
                        'min_evalue': 1.0,
                        'taxonomy_id': None,
                        'detection_sources': []
                    }
                
                all_pathogens[org_key]['detection_methods'].append('bracken')
                all_pathogens[org_key]['detection_sources'].append('FASTQ taxonomic classification')
                all_pathogens[org_key]['abundance_percentage'] = max(all_pathogens[org_key]['abundance_percentage'], pathogen['abundance'])
                all_pathogens[org_key]['estimated_reads'] = max(all_pathogens[org_key]['estimated_reads'], pathogen['reads'])
                if 'taxonomy_id' in pathogen:
                    all_pathogens[org_key]['taxonomy_id'] = pathogen['taxonomy_id']
        
        # Process taxonomy BLAST pathogens
        if taxonomy_pathogens:
            for pathogen in taxonomy_pathogens:
                org_key = pathogen['organism'].strip()
                if org_key not in all_pathogens:
                    all_pathogens[org_key] = {
                        'organism': org_key,
                        'detection_methods': [],
                        'risk_level': 'Unknown',
                        'abundance_percentage': 0,
                        'estimated_reads': 0,
                        'sequence_hits': 0,
                        'max_identity': 0,
                        'min_evalue': 1.0,
                        'taxonomy_id': None,
                        'detection_sources': []
                    }
                
                if 'taxonomy' not in all_pathogens[org_key]['detection_methods']:
                    all_pathogens[org_key]['detection_methods'].append('taxonomy')
                    all_pathogens[org_key]['detection_sources'].append('BLAST taxonomy classification')
                all_pathogens[org_key]['sequence_hits'] += 1
                all_pathogens[org_key]['max_identity'] = max(all_pathogens[org_key]['max_identity'], pathogen['identity'])
                all_pathogens[org_key]['min_evalue'] = min(all_pathogens[org_key]['min_evalue'], pathogen['evalue'])
        
        # Process sequence-based pathogens (custom database)
        if sequence_pathogens:
            for pathogen in sequence_pathogens:
                org_key = pathogen['organism'].strip()
                if org_key not in all_pathogens:
                    all_pathogens[org_key] = {
                        'organism': org_key,
                        'detection_methods': [],
                        'risk_level': 'Unknown',
                        'abundance_percentage': 0,
                        'estimated_reads': 0,
                        'sequence_hits': 0,
                        'max_identity': 0,
                        'min_evalue': 1.0,
                        'taxonomy_id': None,
                        'detection_sources': []
                    }
                
                if 'sequence' not in all_pathogens[org_key]['detection_methods']:
                    all_pathogens[org_key]['detection_methods'].append('sequence')
                    all_pathogens[org_key]['detection_sources'].append('Custom pathogen database')
                all_pathogens[org_key]['sequence_hits'] += 1
                all_pathogens[org_key]['max_identity'] = max(all_pathogens[org_key]['max_identity'], pathogen['identity'])
                all_pathogens[org_key]['min_evalue'] = min(all_pathogens[org_key]['min_evalue'], pathogen['evalue'])
        
        # Enhanced risk assessment with comprehensive keywords
        high_risk_keywords = [
            'brucella', 'salmonella', 'escherichia coli', 'staphylococcus aureus',
            'clostridium difficile', 'bacillus anthracis', 'yersinia pestis',
            'francisella tularensis', 'mycobacterium tuberculosis', 'pseudomonas aeruginosa',
            'listeria monocytogenes', 'vibrio cholerae', 'coxiella burnetii',
            'burkholderia mallei', 'burkholderia pseudomallei', 'rickettsia prowazekii',
            'variola virus', 'marburg virus', 'ebola virus', 'lassa virus'
        ]
        
        medium_risk_keywords = [
            'klebsiella', 'acinetobacter', 'enterococcus', 'streptococcus pyogenes',
            'campylobacter', 'shigella', 'bacillus cereus', 'clostridium tetani',
            'mycobacterium', 'corynebacterium', 'proteus', 'citrobacter',
            'enterobacter', 'serratia', 'morganella', 'providencia'
        ]
        
        # Assign risk levels with enhanced multi-factor criteria
        for pathogen_data in all_pathogens.values():
            org_lower = pathogen_data['organism'].lower()
            base_risk = 'LOW'
            
            # Primary risk assessment based on organism type
            if any(keyword in org_lower for keyword in high_risk_keywords):
                base_risk = 'HIGH'
            elif any(keyword in org_lower for keyword in medium_risk_keywords):
                base_risk = 'MEDIUM'
            
            # Risk enhancement based on detection strength
            method_count = len(pathogen_data['detection_methods'])
            abundance = pathogen_data['abundance_percentage']
            identity = pathogen_data['max_identity']
            hits = pathogen_data['sequence_hits']
            
            # Multi-method detection increases confidence
            if method_count >= 2:
                if base_risk == 'MEDIUM':
                    base_risk = 'HIGH'
                elif base_risk == 'LOW':
                    base_risk = 'MEDIUM'
            
            # High abundance or quality increases risk
            if abundance > 0.1 or (identity > 95 and hits > 5):
                if base_risk == 'LOW':
                    base_risk = 'MEDIUM'
                elif base_risk == 'MEDIUM':
                    base_risk = 'HIGH'
            
            pathogen_data['risk_level'] = base_risk
            pathogen_data['confidence_score'] = self._calculate_confidence_score(pathogen_data)
        
        # Create enhanced summary statistics
        total_pathogens = len(all_pathogens)
        high_risk_count = sum(1 for p in all_pathogens.values() if p['risk_level'] == 'HIGH')
        medium_risk_count = sum(1 for p in all_pathogens.values() if p['risk_level'] == 'MEDIUM')
        low_risk_count = sum(1 for p in all_pathogens.values() if p['risk_level'] == 'LOW')
        
        # Determine overall risk with enhanced criteria
        if high_risk_count >= 3:
            overall_risk = 'CRITICAL'
        elif high_risk_count >= 1:
            overall_risk = 'HIGH'
        elif medium_risk_count >= 3:
            overall_risk = 'MEDIUM'
        elif total_pathogens == 0:
            overall_risk = 'LOW'
        else:
            overall_risk = 'LOW'
        
        # Sort pathogens by risk and evidence strength
        risk_order = {'HIGH': 3, 'MEDIUM': 2, 'LOW': 1, 'Unknown': 0}
        sorted_pathogens = sorted(all_pathogens.values(),
                                key=lambda x: (risk_order[x['risk_level']], 
                                              len(x['detection_methods']),
                                              x['max_identity'],
                                              x['abundance_percentage']),
                                reverse=True)
        
        # Generate comprehensive pathogen report
        pathogen_report = {
            'analysis_summary': {
                'total_pathogens_detected': total_pathogens,
                'high_risk_pathogens': high_risk_count,
                'medium_risk_pathogens': medium_risk_count,
                'low_risk_pathogens': low_risk_count,
                'overall_risk_assessment': overall_risk,
                'detection_methods_used': ['bracken_taxonomic', 'taxonomy_blast', 'sequence_database_search'],
                'analysis_timestamp': self.timestamp.isoformat(),
                'pipeline_version': 'MetaQuest v3.0'
            },
            'pathogen_detections': {}
        }
        
        # Add detailed pathogen information
        for pathogen in sorted_pathogens:
            pathogen_report['pathogen_detections'][pathogen['organism']] = {
                'risk_level': pathogen['risk_level'],
                'detection_methods': pathogen['detection_methods'],
                'detection_sources': pathogen['detection_sources'],
                'abundance_percentage': round(pathogen['abundance_percentage'] * 100, 3) if pathogen['abundance_percentage'] > 0 else 0,
                'estimated_reads': int(pathogen['estimated_reads']),
                'sequence_identity': round(pathogen['max_identity'], 1) if pathogen['max_identity'] > 0 else 'N/A',
                'sequence_hits': pathogen['sequence_hits'],
                'min_evalue': pathogen['min_evalue'] if pathogen['min_evalue'] < 1.0 else 'N/A',
                'taxonomy_id': pathogen['taxonomy_id'],
                'confidence_score': pathogen['confidence_score']
            }
        
        # Save ONLY JSON report for visualizations
        report_file = output_dir / "pathogen_detection_report.json"
        with open(report_file, 'w') as f:
            json.dump(pathogen_report, f, indent=2)
        
        # Create ONLY the comprehensive pathogen summary (BEST REPORT)
        self._create_enhanced_pathogen_summary(output_dir, sorted_pathogens, overall_risk, 
                                             high_risk_count, medium_risk_count, total_pathogens)
        
        print(f"✓ Comprehensive pathogen report generated:")
        print(f"  • {total_pathogens} pathogens detected")
        print(f"  • {high_risk_count} high-risk pathogens")
        print(f"  • {medium_risk_count} medium-risk pathogens")
        print(f"  • Overall risk: {overall_risk}")
        
        return pathogen_report
    
    def _calculate_confidence_score(self, pathogen_data):
        """Calculate confidence score based on multiple detection criteria"""
        score = 0
        
        # Method diversity (30% weight)
        method_bonus = min(len(pathogen_data['detection_methods']) * 0.15, 0.3)
        score += method_bonus
        
        # Sequence quality (40% weight)
        if pathogen_data['max_identity'] > 0:
            identity_score = (pathogen_data['max_identity'] / 100) * 0.4
            score += identity_score
        
        # Hit abundance (20% weight)
        if pathogen_data['sequence_hits'] > 0:
            hit_score = min(pathogen_data['sequence_hits'] / 10, 0.2)
            score += hit_score
        
        # Read support (10% weight)
        if pathogen_data['estimated_reads'] > 0:
            read_score = min(np.log10(pathogen_data['estimated_reads']) / 50, 0.1)
            score += read_score
        
        return round(min(score, 1.0), 3)
    
    def _create_enhanced_pathogen_summary(self, output_dir, sorted_pathogens, overall_risk, 
                                        high_risk_count, medium_risk_count, total_pathogens):
        """Create THE BEST comprehensive pathogen summary - KEEP THIS ONE"""
        summary_file = output_dir / "pathogen_summary.txt"
        
        with open(summary_file, 'w') as f:
            f.write("COMPREHENSIVE PATHOGEN DETECTION SUMMARY\n")
            f.write("=" * 60 + "\n\n")
            f.write(f"Overall Risk Assessment: {overall_risk}\n")
            f.write(f"Total Pathogens Detected: {total_pathogens}\n")
            f.write(f"High Risk Pathogens: {high_risk_count}\n")
            f.write(f"Medium Risk Pathogens: {medium_risk_count}\n")
            f.write(f"Analysis Date: {self.timestamp.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Pipeline Version: MetaQuest v3.0\n\n")
            
            # Detection methodology
            f.write("DETECTION METHODOLOGY:\n")
            f.write("• Multi-source pathogen screening approach\n")
            f.write("• FASTQ taxonomic classification (Kraken2/Bracken)\n")
            f.write("• FASTA BLAST taxonomy classification\n")
            f.write("• Custom pathogen database screening\n")
            f.write("• Enhanced risk stratification algorithms\n\n")
            
            # High risk pathogens section
            if high_risk_count > 0:
                f.write("🚨 HIGH RISK PATHOGENS DETECTED:\n")
                f.write("-" * 50 + "\n")
                for pathogen in sorted_pathogens:
                    if pathogen['risk_level'] == 'HIGH':
                        f.write(f"• {pathogen['organism']}\n")
                        f.write(f"  Risk Level: HIGH\n")
                        f.write(f"  Detection Methods: {', '.join(pathogen['detection_methods'])}\n")
                        f.write(f"  Sources: {'; '.join(pathogen['detection_sources'])}\n")
                        if pathogen['abundance_percentage'] > 0:
                            f.write(f"  Abundance: {pathogen['abundance_percentage']*100:.3f}%\n")
                        if pathogen['max_identity'] > 0:
                            f.write(f"  Max Identity: {pathogen['max_identity']:.1f}%\n")
                        f.write(f"  Confidence Score: {pathogen['confidence_score']:.3f}\n")
                        f.write(f"  Clinical Priority: IMMEDIATE ATTENTION REQUIRED\n")
                        f.write("\n")
            
            # Medium risk pathogens section
            if medium_risk_count > 0:
                f.write("⚠️ MEDIUM RISK PATHOGENS:\n")
                f.write("-" * 35 + "\n")
                for pathogen in sorted_pathogens:
                    if pathogen['risk_level'] == 'MEDIUM':
                        f.write(f"• {pathogen['organism']}\n")
                        f.write(f"  Detection Methods: {', '.join(pathogen['detection_methods'])}\n")
                        if pathogen['abundance_percentage'] > 0:
                            f.write(f"  Abundance: {pathogen['abundance_percentage']*100:.3f}%\n")
                        if pathogen['max_identity'] > 0:
                            f.write(f"  Max Identity: {pathogen['max_identity']:.1f}%\n")
                        f.write(f"  Confidence: {pathogen['confidence_score']:.3f}\n")
                        f.write("\n")
            
            # Clinical recommendations
            f.write("CLINICAL RECOMMENDATIONS:\n")
            f.write("-" * 30 + "\n")
            if overall_risk == 'CRITICAL':
                f.write("🚨 CRITICAL RISK - EMERGENCY PROTOCOLS REQUIRED\n")
                f.write("• Multiple high-risk pathogens detected across sources\n")
                f.write("• Implement immediate containment measures\n")
                f.write("• Consider broad-spectrum antimicrobial therapy\n")
                f.write("• Emergency infectious disease consultation\n")
                f.write("• Enhanced patient monitoring protocols\n")
            elif overall_risk == 'HIGH':
                f.write("🔴 HIGH RISK - URGENT CLINICAL EVALUATION\n")
                f.write("• Significant pathogenic burden identified\n")
                f.write("• Targeted antimicrobial therapy recommended\n")
                f.write("• Close patient monitoring required\n")
                f.write("• Consider isolation precautions\n")
                f.write("• Follow-up diagnostic confirmation advised\n")
            elif overall_risk == 'MEDIUM':
                f.write("🟡 MEDIUM RISK - ENHANCED MONITORING\n")
                f.write("• Some pathogenic indicators detected\n")
                f.write("• Standard monitoring protocols with vigilance\n")
                f.write("• Consider targeted diagnostics if symptomatic\n")
                f.write("• Implement standard infection control measures\n")
                f.write("• Regular reassessment recommended\n")
            else:
                f.write("🟢 LOW RISK - ROUTINE MONITORING\n")
                f.write("• No significant pathogenic signatures detected\n")
                f.write("• Continue standard care protocols\n")
                f.write("• Routine infection control measures sufficient\n")
                f.write("• Periodic reassessment as clinically indicated\n")
            
            f.write("\nDISCLAIMER:\n")
            f.write("This analysis is for research and clinical decision support.\n")
            f.write("Results should be interpreted in clinical context.\n")
            f.write("Confirmatory testing may be required for clinical decisions.\n")
    
    def generate_report(self, traditional_results=None, ml_results=None, ml_summary=None, **kwargs) -> Dict[str, str]:
        """Generate ONLY JSON report for ML integration - NO redundant text reports"""
        generated_files = {}
        
        # KEEP ONLY the JSON report for visualizations and ML integration
        if ml_results and ml_summary:
            json_data = self._create_pathogen_json_data(traditional_results, ml_results, ml_summary)
            generated_files['json_report'] = self._save_json_report(json_data, "pathogen_detection_report.json")
        
        return generated_files
    
    def generate_blast_ml_integrated_report(self, output_dir, blast_taxonomy_pathogens=None, ml_results=None, ml_summary=None):
        """Generate properly separated BLAST + ML pathogen report"""
        print("Generating integrated BLAST+ML pathogen report...")
        
        # SECTION 1: Process BLAST taxonomy pathogens (separate)
        blast_pathogens = {}
        if blast_taxonomy_pathogens:
            for pathogen in blast_taxonomy_pathogens:
                org_key = pathogen['organism'].strip()
                if org_key not in blast_pathogens:
                    blast_pathogens[org_key] = {
                        'organism': org_key,
                        'blast_identity': 0,
                        'blast_evalue': 1.0,
                        'blast_hits': 0,
                        'risk_level': 'Unknown'
                    }
                
                blast_pathogens[org_key]['blast_identity'] = max(blast_pathogens[org_key]['blast_identity'], pathogen['identity'])
                blast_pathogens[org_key]['blast_evalue'] = min(blast_pathogens[org_key]['blast_evalue'], pathogen['evalue'])
                blast_pathogens[org_key]['blast_hits'] += 1

        # SECTION 2: Process ML predictions (separate)
        ml_data = {
            'total_proteins': 0,
            'pathogenic_predictions': 0,
            'high_confidence_predictions': 0,
            'average_confidence': 0,
            'individual_predictions': []
        }
        
        if ml_results and ml_summary:
            ml_data['total_proteins'] = ml_summary.get('total_sequences_analyzed', 0)
            ml_data['pathogenic_predictions'] = ml_summary.get('pathogenic_predictions', 0)
            ml_data['high_confidence_predictions'] = ml_summary.get('high_confidence_predictions', 0)
            ml_data['average_confidence'] = ml_summary.get('average_confidence', 0)
            
            # Add individual protein predictions
            if ml_results:
                for prediction in ml_results:
                    ml_data['individual_predictions'].append({
                        'protein_id': prediction.get('sequence_id', 'Unknown'),
                        'confidence': prediction.get('confidence', 0),
                        'is_high_confidence': prediction.get('high_confidence', False)
                    })

        # SECTION 3: Overall risk assessment (combined logic)
        blast_high_risk = len([p for p in blast_pathogens.values() if self._assess_blast_risk(p) == 'HIGH'])
        ml_pathogenic_ratio = ml_data['pathogenic_predictions'] / max(ml_data['total_proteins'], 1)
        
        if blast_high_risk >= 2 and ml_pathogenic_ratio > 0.8:
            overall_risk = 'CRITICAL'
        elif blast_high_risk >= 1 and ml_pathogenic_ratio > 0.6:
            overall_risk = 'HIGH'
        else:
            overall_risk = 'MEDIUM'
        
        # Generate corrected JSON report
        integrated_report = {
            'analysis_summary': {
                'analysis_type': 'BLAST_ML_Separated_Integration',
                'overall_risk_assessment': overall_risk,
                'analysis_timestamp': self.timestamp.isoformat()
            },
            'blast_taxonomy_section': {
                'pathogenic_organisms': len(blast_pathogens),
                'detections': blast_pathogens
            },
            'ml_prediction_section': ml_data,
            'integrated_assessment': {
                'risk_justification': f"BLAST: {blast_high_risk} high-risk species, ML: {ml_pathogenic_ratio:.1%} pathogenic proteins"
            }
        }
        
        # Save corrected integrated report
        report_file = output_dir / "blast_ml_integrated_pathogen_report.json"
        with open(report_file, 'w') as f:
            json.dump(integrated_report, f, indent=2)
        
        # Create corrected text summary
        self._create_separated_blast_ml_summary(output_dir, blast_pathogens, ml_data, overall_risk)
        
        print(f"✓ Integrated BLAST+ML pathogen report generated:")
        print(f"  • {len(blast_pathogens)} pathogenic organisms detected (BLAST taxonomy)")
        print(f"  • {ml_data['pathogenic_predictions']} pathogenic proteins predicted (ML)")
        print(f"  • {blast_high_risk} high-risk organisms (BLAST)")
        print(f"  • Overall risk: {overall_risk}")
        
        return integrated_report

    
    def _assess_blast_risk(self, pathogen_data):
        """Assess risk level for BLAST-detected pathogens"""
        
        # High-risk pathogenic keywords
        high_risk_keywords = [
            'brucella', 'salmonella', 'escherichia coli', 'staphylococcus aureus',
            'clostridium difficile', 'bacillus anthracis', 'yersinia pestis',
            'francisella tularensis', 'mycobacterium tuberculosis', 'pseudomonas aeruginosa',
            'listeria monocytogenes', 'vibrio cholerae', 'coxiella burnetii'
        ]
        
        # Medium-risk pathogenic keywords
        medium_risk_keywords = [
            'klebsiella', 'acinetobacter', 'enterococcus', 'streptococcus',
            'campylobacter', 'shigella', 'bacillus cereus', 'clostridium tetani',
            'mycobacterium', 'corynebacterium', 'proteus', 'citrobacter'
        ]
        
        organism_lower = pathogen_data['organism'].lower()
        base_risk = 'LOW'
        
        # Primary risk assessment based on organism type
        if any(keyword in organism_lower for keyword in high_risk_keywords):
            base_risk = 'HIGH'
        elif any(keyword in organism_lower for keyword in medium_risk_keywords):
            base_risk = 'MEDIUM'
        
        # Risk enhancement based on BLAST quality metrics
        blast_identity = pathogen_data.get('blast_identity', 0)
        blast_evalue = pathogen_data.get('blast_evalue', 1.0)
        blast_hits = pathogen_data.get('blast_hits', 0)
        
        # Enhance risk based on high-quality BLAST matches
        if blast_identity > 98 and blast_evalue < 1e-100 and blast_hits >= 3:
            if base_risk == 'MEDIUM':
                base_risk = 'HIGH'
            elif base_risk == 'LOW':
                base_risk = 'MEDIUM'
        elif blast_identity > 95 and blast_evalue < 1e-50:
            if base_risk == 'LOW':
                base_risk = 'MEDIUM'
        
        # Store the risk level in the pathogen data
        pathogen_data['risk_level'] = base_risk
        
        return base_risk

    def _create_separated_blast_ml_summary(self, output_dir, blast_pathogens, ml_data, overall_risk):
        """Create separated BLAST+ML pathogen summary report"""
        summary_file = output_dir / "blast_ml_pathogen_summary.txt"
        
        with open(summary_file, 'w') as f:
            f.write("INTEGRATED BLAST+ML PATHOGEN DETECTION SUMMARY\n")
            f.write("=" * 65 + "\n\n")
            f.write(f"Overall Risk Assessment: {overall_risk}\n")
            f.write(f"Analysis Date: {self.timestamp.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Pipeline Version: MetaQuest v3.0\n\n")
            
            # Detection methodology
            f.write("INTEGRATED DETECTION METHODOLOGY:\n")
            f.write("• BLAST taxonomic classification against NCBI database\n")
            f.write("• Machine learning pathogenicity prediction\n")
            f.write("• Dual-method validation for enhanced confidence\n")
            f.write("• Clinical risk stratification with separated analysis\n\n")
            
            # SECTION 1: BLAST TAXONOMIC CLASSIFICATION
            f.write("SECTION 1: BLAST TAXONOMIC CLASSIFICATION\n")
            f.write("=" * 50 + "\n")
            f.write(f"• Pathogenic organisms detected: {len(blast_pathogens)}\n")
            
            if blast_pathogens:
                # Calculate total BLAST hits
                total_blast_hits = sum(p.get('blast_hits', 0) for p in blast_pathogens.values())
                avg_identity = sum(p.get('blast_identity', 0) for p in blast_pathogens.values()) / len(blast_pathogens)
                
                f.write(f"• Total BLAST hits: {total_blast_hits}\n")
                f.write(f"• Average identity: {avg_identity:.1f}%\n\n")
                
                # Group by risk level for BLAST pathogens
                high_risk_blast = [p for p in blast_pathogens.values() if p.get('risk_level') == 'HIGH']
                medium_risk_blast = [p for p in blast_pathogens.values() if p.get('risk_level') == 'MEDIUM']
                
                if high_risk_blast:
                    f.write("🚨 HIGH RISK PATHOGENS (BLAST):\n")
                    f.write("-" * 40 + "\n")
                    for pathogen in sorted(high_risk_blast, key=lambda x: x.get('blast_hits', 0), reverse=True):
                        f.write(f"• {pathogen['organism']}\n")
                        f.write(f"  - BLAST Identity: {pathogen.get('blast_identity', 0):.1f}%\n")
                        f.write(f"  - BLAST Hits: {pathogen.get('blast_hits', 0)}\n")
                        f.write(f"  - E-value: {pathogen.get('blast_evalue', 1.0):.2e}\n")
                        f.write(f"  - Risk Level: HIGH\n")
                        f.write(f"  - Clinical Priority: IMMEDIATE\n\n")
                
                if medium_risk_blast:
                    f.write("⚠️ MEDIUM RISK PATHOGENS (BLAST):\n")
                    f.write("-" * 40 + "\n")
                    for pathogen in sorted(medium_risk_blast, key=lambda x: x.get('blast_hits', 0), reverse=True):
                        f.write(f"• {pathogen['organism']}\n")
                        f.write(f"  - BLAST Identity: {pathogen.get('blast_identity', 0):.1f}%\n")
                        f.write(f"  - BLAST Hits: {pathogen.get('blast_hits', 0)}\n")
                        f.write(f"  - Risk Level: MEDIUM\n\n")
            
            # SECTION 2: MACHINE LEARNING PREDICTIONS
            f.write("SECTION 2: MACHINE LEARNING PREDICTIONS\n")
            f.write("=" * 50 + "\n")
            f.write(f"• Total proteins analyzed: {ml_data.get('total_proteins', 0)}\n")
            f.write(f"• Pathogenic predictions: {ml_data.get('pathogenic_predictions', 0)} ({ml_data.get('pathogenic_predictions', 0)/max(ml_data.get('total_proteins', 1), 1)*100:.1f}%)\n")
            f.write(f"• High-confidence predictions: {ml_data.get('high_confidence_predictions', 0)} ({ml_data.get('high_confidence_predictions', 0)/max(ml_data.get('pathogenic_predictions', 1), 1)*100:.1f}%)\n")
            f.write(f"• Average ML confidence: {ml_data.get('average_confidence', 0):.1%}\n\n")
            
            if ml_data.get('individual_predictions'):
                f.write("🤖 ML PATHOGENIC PROTEIN PREDICTIONS:\n")
                f.write("-" * 40 + "\n")
                # Sort by confidence
                sorted_predictions = sorted(ml_data['individual_predictions'], 
                                        key=lambda x: x.get('confidence', 0), reverse=True)
                
                for pred in sorted_predictions:
                    confidence = pred.get('confidence', 0)
                    confidence_level = "High" if pred.get('is_high_confidence', False) else "Medium"
                    f.write(f"• {pred.get('protein_id', 'Unknown')}: {confidence:.1%} confidence ({confidence_level})\n")
            
            # SECTION 3: INTEGRATED RISK ASSESSMENT
            f.write(f"\nSECTION 3: INTEGRATED RISK ASSESSMENT\n")
            f.write("=" * 50 + "\n")
            f.write(f"Overall Risk: {overall_risk}\n\n")
            
            # Risk justification
            blast_high_count = len([p for p in blast_pathogens.values() if p.get('risk_level') == 'HIGH'])
            ml_pathogenic_ratio = ml_data.get('pathogenic_predictions', 0) / max(ml_data.get('total_proteins', 1), 1)
            
            f.write("Risk Justification:\n")
            f.write(f"• BLAST: {blast_high_count} high-risk pathogenic species identified\n")
            f.write(f"• ML: {ml_pathogenic_ratio:.1%} of proteins predicted as pathogenic\n")
            f.write(f"• Confidence: Both methods show high reliability\n")
            f.write(f"• Clinical Impact: Multiple dangerous pathogens present\n\n")
            
            # Clinical recommendations based on integrated findings
            f.write("CLINICAL RECOMMENDATIONS:\n")
            f.write("-" * 30 + "\n")
            if overall_risk == 'CRITICAL':
                f.write("🚨 CRITICAL RISK - DUAL METHOD VALIDATION\n")
                f.write("• Multiple high-risk pathogens confirmed by both BLAST and ML\n")
                f.write("• High-confidence pathogenic signatures detected\n")
                f.write("• Implement immediate containment protocols\n")
                f.write("• Emergency infectious disease consultation required\n")
            elif overall_risk == 'HIGH':
                f.write("🔴 HIGH RISK - VALIDATED PATHOGENIC PRESENCE\n")
                f.write("• Significant pathogenic organisms detected\n")
                f.write("• ML predictions support BLAST taxonomy findings\n")
                f.write("• Enhanced monitoring and targeted interventions recommended\n")
                f.write("• Consider confirmatory diagnostics\n")
            else:
                f.write("🟡 STANDARD RISK - CONTINUE MONITORING\n")
                f.write("• Some pathogenic indicators detected\n")
                f.write("• Integrated analysis provides moderate confidence\n")
                f.write("• Standard protocols with enhanced vigilance\n")
            
            f.write("\nINTEGRATED ANALYSIS ADVANTAGES:\n")
            f.write("• BLAST provides taxonomic identification with sequence similarity\n")
            f.write("• ML provides pathogenicity assessment with confidence scoring\n")
            f.write("• Separated analysis prevents cross-contamination of evidence\n")
            f.write("• Dual validation increases clinical confidence in findings\n")
            
            f.write("\nDISCLAIMER:\n")
            f.write("This integrated analysis combines BLAST taxonomy and ML predictions.\n")
            f.write("Results should be interpreted in clinical context.\n")
            f.write("Confirmatory testing may be required for clinical decisions.\n")
        
        print(f"✓ Separated BLAST+ML pathogen summary created: {summary_file}")


    def _calculate_blast_ml_confidence(self, pathogen_data, ml_confidence):
        """Calculate confidence score for BLAST+ML integrated detection"""
        score = 0
        
        # BLAST quality (50% weight)
        blast_score = (pathogen_data['blast_identity'] / 100) * 0.5
        score += blast_score
        
        # ML confidence (30% weight)
        if ml_confidence > 0:
            ml_score = ml_confidence * 0.3
            score += ml_score
        
        # Multi-method detection bonus (20% weight)
        method_bonus = len(pathogen_data['detection_methods']) * 0.1
        score += min(method_bonus, 0.2)
        
        return round(min(score, 1.0), 3)
    
    def _create_pathogen_json_data(self, traditional_results, ml_results, ml_summary) -> Dict:
        """Create structured JSON data for pathogen report"""
        return {
            'traditional_analysis': traditional_results if isinstance(traditional_results, dict) else None,
            'ml_analysis': {
                'summary': ml_summary,
                'predictions': ml_results
            } if ml_results and ml_summary else None,
            'risk_assessment': {
                'timestamp': datetime.now().isoformat(),
                'methods_used': ['traditional_detection', 'ml_prediction'] if ml_results else ['traditional_detection'],
                'confidence_metrics': {
                    'traditional_confidence': 'database_based',
                    'ml_confidence': ml_summary.get('average_confidence', 0) if ml_summary else None
                }
            }
        }

class FunctionalReporter(BaseReportGenerator):
    """Functional annotation reporting with quality assessment"""
    
    def __init__(self, output_dir: Path):
        super().__init__(output_dir, "Functional Annotation")
    
    def generate_report(self, prokka_results=None, swissprot_results=None, **kwargs) -> Dict[str, str]:
        """Generate comprehensive functional annotation report"""
        generated_files = {}
        
        # Generate text report
        text_content = self._create_report_header("FUNCTIONAL ANNOTATION REPORT")
        
        if prokka_results:
            text_content.extend(self._process_prokka_results(prokka_results))
        
        if swissprot_results:
            text_content.extend(self._process_swissprot_results(swissprot_results))
        
        text_content.extend(self._generate_functional_insights(prokka_results, swissprot_results))
        
        # Save reports
        generated_files['text_report'] = self._save_text_report(text_content, "functional_annotation_report.txt")
        
        # Generate JSON report
        json_data = self._create_functional_json_data(prokka_results, swissprot_results)
        generated_files['json_report'] = self._save_json_report(json_data, "functional_annotation_report.json")
        
        return generated_files
    
    def _process_prokka_results(self, prokka_dir) -> List[str]:
        """Process Prokka results with enhanced quality assessment"""
        lines = [
            "GENE PREDICTION AND ANNOTATION (Prokka)",
            "-" * 42,
            ""
        ]
        
        try:
            prokka_path = Path(prokka_dir)
            
            # Count file types
            faa_files = list(prokka_path.glob("*.faa"))
            ffn_files = list(prokka_path.glob("*.ffn"))
            gff_files = list(prokka_path.glob("*.gff"))
            
            lines.extend([
                f"PROKKA OUTPUT SUMMARY:",
                f"• Analysis directory: {prokka_path.name}",
                f"• Protein sequences (*.faa): {len(faa_files)} files",
                f"• Gene sequences (*.ffn): {len(ffn_files)} files",
                f"• Annotation files (*.gff): {len(gff_files)} files",
                ""
            ])
            
            # Analyze protein sequences in detail
            if faa_files:
                protein_file = faa_files[0]
                try:
                    from Bio import SeqIO
                    proteins = list(SeqIO.parse(protein_file, "fasta"))
                    
                    if proteins:
                        lengths = [len(seq.seq) for seq in proteins]
                        
                        lines.extend([
                            f"PROTEIN SEQUENCE ANALYSIS:",
                            f"• Total proteins predicted: {len(proteins)}",
                            f"• Source file: {protein_file.name}",
                            "",
                            f"LENGTH DISTRIBUTION ANALYSIS:",
                            f"• Average length: {np.mean(lengths):.1f} amino acids",
                            f"• Median length: {np.median(lengths):.1f} amino acids",
                            f"• Shortest protein: {min(lengths)} amino acids",
                            f"• Longest protein: {max(lengths)} amino acids",
                            f"• Standard deviation: {np.std(lengths):.1f} amino acids",
                            "",
                            f"QUALITY ASSESSMENT:",
                        ])
                        
                        # Quality metrics based on length distribution
                        very_short = len([l for l in lengths if l < 50])
                        short = len([l for l in lengths if 50 <= l < 100])
                        medium = len([l for l in lengths if 100 <= l < 300])
                        long_proteins = len([l for l in lengths if l >= 300])
                        
                        lines.extend([
                            f"• Very short (<50 aa): {very_short} ({very_short/len(proteins)*100:.1f}%)",
                            f"• Short (50-99 aa): {short} ({short/len(proteins)*100:.1f}%)",
                            f"• Medium (100-299 aa): {medium} ({medium/len(proteins)*100:.1f}%)",
                            f"• Long (≥300 aa): {long_proteins} ({long_proteins/len(proteins)*100:.1f}%)",
                            ""
                        ])
                        
                        # Quality interpretation
                        if very_short / len(proteins) > 0.4:
                            lines.append("⚠️ HIGH PROPORTION OF VERY SHORT PROTEINS - May indicate fragmented assembly")
                        elif np.mean(lengths) < 150:
                            lines.append("⚠️ SHORT AVERAGE PROTEIN LENGTH - Consider assembly quality")
                        else:
                            lines.append("✓ PROTEIN LENGTH DISTRIBUTION APPEARS NORMAL")
                        
                        lines.append("")
                    else:
                        lines.append("⚠️ No proteins found in FASTA file")
                        
                except Exception as e:
                    lines.append(f"⚠️ Could not analyze protein sequences: {e}")
            
            return lines
            
        except Exception as e:
            return [f"Error processing Prokka results: {e}", ""]
    
    def _process_swissprot_results(self, swissprot_file) -> List[str]:
        """Process SwissProt annotation with detailed functional analysis"""
        lines = [
            "PROTEIN FUNCTIONAL ANNOTATION (SwissProt)",
            "-" * 45,
            ""
        ]
        
        try:
            if not Path(swissprot_file).exists():
                return ["⚠️ SwissProt annotation file not found", ""]
            
            df = pd.read_csv(swissprot_file, sep='\t')
            
            if df.empty:
                return ["⚠️ No SwissProt annotations found", ""]
            
            lines.extend([
                f"ANNOTATION SUMMARY:",
                f"• Total protein matches: {len(df)}",
                f"• Unique proteins annotated: {df['qseqid'].nunique()}",
                f"• Average sequence identity: {df['pident'].mean():.1f}%",
                f"• Average alignment length: {df['length'].mean():.1f} amino acids",
                f"• Best E-value: {df['evalue'].min():.2e}",
                ""
            ])
            
            # Quality distribution
            high_quality = len(df[df['pident'] >= 90])
            medium_quality = len(df[(df['pident'] >= 75) & (df['pident'] < 90)])
            low_quality = len(df[df['pident'] < 75])
            
            lines.extend([
                f"ANNOTATION QUALITY DISTRIBUTION:",
                f"• High quality (≥90% identity): {high_quality} ({high_quality/len(df)*100:.1f}%)",
                f"• Medium quality (75-89% identity): {medium_quality} ({medium_quality/len(df)*100:.1f}%)",
                f"• Low quality (<75% identity): {low_quality} ({low_quality/len(df)*100:.1f}%)",
                ""
            ])
            
            # Enhanced functional category analysis
            if 'stitle' in df.columns:
                descriptions = df['stitle'].str.lower()
                
                # Enhanced functional categories
                function_categories = {
                    'Metabolism': ['kinase', 'synthase', 'dehydrogenase', 'reductase', 'oxidase', 'transferase'],
                    'Transport': ['transporter', 'permease', 'channel', 'pump', 'efflux'],
                    'Regulation': ['regulator', 'repressor', 'activator', 'sensor', 'response'],
                    'Cell Structure': ['structural', 'ribosomal', 'membrane', 'cell wall', 'cytoskeleton'],
                    'DNA/RNA Processing': ['polymerase', 'helicase', 'ligase', 'nuclease', 'transcription'],
                    'Protein Processing': ['protease', 'peptidase', 'chaperone', 'folding', 'modification'],
                    'Stress Response': ['stress', 'heat shock', 'cold shock', 'oxidative', 'survival'],
                    'Virulence/Pathogenesis': ['virulence', 'toxin', 'pathogen', 'adhesin', 'invasion']
                }
                
                category_counts = {}
                for category, keywords in function_categories.items():
                    count = sum(descriptions.str.contains(keyword, na=False).sum() for keyword in keywords)
                    if count > 0:
                        category_counts[category] = count
                
                if category_counts:
                    lines.extend([
                        f"FUNCTIONAL CATEGORY DISTRIBUTION:",
                        f"{'Category':<25} {'Proteins':<10} {'Percentage':<10}",
                        "-" * 47
                    ])
                    
                    for category, count in sorted(category_counts.items(), key=lambda x: x[1], reverse=True):
                        percentage = (count / len(df)) * 100
                        lines.append(f"{category:<25} {count:<10} {percentage:<10.1f}%")
                    
                    lines.extend([
                        "",
                        f"• Total categorized functions: {sum(category_counts.values())}",
                        f"• Proteins with multiple functions may be counted in multiple categories",
                        ""
                    ])
            
            return lines
            
        except Exception as e:
            return [f"Error processing SwissProt results: {e}", ""]
    
    def _generate_functional_insights(self, prokka_results, swissprot_results) -> List[str]:
        """Generate functional analysis insights and recommendations"""
        lines = [
            "FUNCTIONAL ANALYSIS INSIGHTS & RECOMMENDATIONS",
            "-" * 48,
            ""
        ]
        
        # Analysis completeness
        analysis_components = []
        if prokka_results:
            analysis_components.append("Gene prediction (Prokka)")
        if swissprot_results:
            analysis_components.append("Protein annotation (SwissProt)")
        
        if analysis_components:
            lines.extend([
                f"COMPLETED ANALYSES:",
                f"• {' + '.join(analysis_components)}",
                ""
            ])
        
        # Quality assessment
        try:
            if prokka_results and swissprot_results:
                prokka_path = Path(prokka_results)
                faa_files = list(prokka_path.glob("*.faa"))
                
                if faa_files:
                    from Bio import SeqIO
                    proteins = list(SeqIO.parse(faa_files[0], "fasta"))
                    total_proteins = len(proteins)
                    
                    # Read SwissProt results
                    df = pd.read_csv(swissprot_results, sep='\t')
                    annotated_proteins = df['qseqid'].nunique()
                    
                    annotation_rate = (annotated_proteins / total_proteins) * 100 if total_proteins > 0 else 0
                    
                    lines.extend([
                        f"ANNOTATION COVERAGE ANALYSIS:",
                        f"• Predicted proteins: {total_proteins}",
                        f"• Successfully annotated: {annotated_proteins}",
                        f"• Annotation rate: {annotation_rate:.1f}%",
                        ""
                    ])
                    
                    # Coverage interpretation
                    if annotation_rate > 80:
                        lines.append("✓ EXCELLENT ANNOTATION COVERAGE - High-quality functional analysis")
                    elif annotation_rate > 60:
                        lines.append("✓ GOOD ANNOTATION COVERAGE - Reliable functional insights available")
                    elif annotation_rate > 40:
                        lines.append("⚠️ MODERATE ANNOTATION COVERAGE - Some functional gaps present")
                    else:
                        lines.append("⚠️ LOW ANNOTATION COVERAGE - May indicate novel/uncharacterized proteins")
                    
        except Exception as e:
            lines.append(f"⚠️ Could not calculate annotation coverage: {e}")
        
        lines.extend([
            "",
            "RECOMMENDATIONS FOR FURTHER ANALYSIS:",
            "• Cross-reference functional categories with sample type expectations",
            "• Investigate unannotated proteins for novel functions",
            "• Consider specialized databases for domain-specific annotations",
            "• Analyze functional diversity for ecosystem characterization",
            "• Compare functional profiles with reference samples if available",
            "",
            "QUALITY CONSIDERATIONS:",
            "• High-identity matches (>90%) are most reliable for functional assignment",
            "• Short proteins (<100 aa) may have limited functional annotation",
            "• Multiple database searches can improve annotation coverage",
            "• Consider experimental validation for critical findings",
            ""
        ])
        
        return lines
    
    def _create_functional_json_data(self, prokka_results, swissprot_results) -> Dict:
        """Create structured JSON data for functional report"""
        json_data = {
            'prokka_analysis': None,
            'swissprot_analysis': None,
            'coverage_metrics': None
        }
        
        # Process Prokka data
        if prokka_results:
            try:
                prokka_path = Path(prokka_results)
                faa_files = list(prokka_path.glob("*.faa"))
                
                if faa_files:
                    from Bio import SeqIO
                    proteins = list(SeqIO.parse(faa_files[0], "fasta"))
                    lengths = [len(seq.seq) for seq in proteins]
                    
                    json_data['prokka_analysis'] = {
                        'total_proteins': len(proteins),
                        'length_statistics': {
                            'mean': np.mean(lengths),
                            'median': np.median(lengths),
                            'min': min(lengths) if lengths else 0,
                            'max': max(lengths) if lengths else 0,
                            'std': np.std(lengths)
                        },
                        'length_distribution': {
                            'very_short': len([l for l in lengths if l < 50]),
                            'short': len([l for l in lengths if 50 <= l < 100]),
                            'medium': len([l for l in lengths if 100 <= l < 300]),
                            'long': len([l for l in lengths if l >= 300])
                        }
                    }
            except Exception as e:
                json_data['prokka_analysis'] = {'error': str(e)}
        
        # Process SwissProt data
        if swissprot_results and Path(swissprot_results).exists():
            try:
                df = pd.read_csv(swissprot_results, sep='\t')
                
                json_data['swissprot_analysis'] = {
                    'total_matches': len(df),
                    'unique_proteins': df['qseqid'].nunique(),
                    'identity_statistics': {
                        'mean': df['pident'].mean(),
                        'median': df['pident'].median(),
                        'min': df['pident'].min(),
                        'max': df['pident'].max()
                    },
                    'quality_distribution': {
                        'high_quality': len(df[df['pident'] >= 90]),
                        'medium_quality': len(df[(df['pident'] >= 75) & (df['pident'] < 90)]),
                        'low_quality': len(df[df['pident'] < 75])
                    }
                }
                
                # Calculate annotation coverage if both datasets available
                if json_data['prokka_analysis'] and not isinstance(json_data['prokka_analysis'], dict) or 'error' not in json_data['prokka_analysis']:
                    total_proteins = json_data['prokka_analysis']['total_proteins']
                    annotated_proteins = df['qseqid'].nunique()
                    
                    json_data['coverage_metrics'] = {
                        'annotation_rate': (annotated_proteins / total_proteins) * 100 if total_proteins > 0 else 0,
                        'total_predicted': total_proteins,
                        'total_annotated': annotated_proteins,
                        'unannotated': total_proteins - annotated_proteins
                    }
                    
            except Exception as e:
                json_data['swissprot_analysis'] = {'error': str(e)}
        
        return json_data

# Main reporting functions for pipeline integration - STREAMLINED
def generate_taxonomic_report(output_dir, bracken_data=None, blast_data=None):
    """Generate comprehensive taxonomic classification report"""
    reporter = TaxonomicReporter(Path(output_dir))
    return reporter.generate_report(bracken_data=bracken_data, blast_data=blast_data)

def generate_pathogen_report(output_dir, traditional_results=None, ml_results=None, ml_summary=None):
    """Generate ONLY JSON pathogen report for ML integration - NO redundant text"""
    reporter = PathogenReporter(Path(output_dir))
    return reporter.generate_report(traditional_results=traditional_results, ml_results=ml_results, ml_summary=ml_summary)

def generate_functional_report(output_dir, prokka_results=None, swissprot_results=None):
    """Generate comprehensive functional annotation report"""
    reporter = FunctionalReporter(Path(output_dir))
    return reporter.generate_report(prokka_results=prokka_results, swissprot_results=swissprot_results)

# THE ONLY COMPREHENSIVE PATHOGEN FUNCTION - replaces all redundant versions
def generate_pathogen_summary_report(output_dir, bracken_pathogens=None, taxonomy_pathogens=None, sequence_pathogens=None):
    """
    THE DEFINITIVE pathogen report function - ONLY ONE NEEDED
    Generates: pathogen_summary.txt (BEST) + pathogen_detection_report.json (for visualizations)
    """
    reporter = PathogenReporter(Path(output_dir))
    return reporter.generate_comprehensive_pathogen_report(output_dir, bracken_pathogens, taxonomy_pathogens, sequence_pathogens)
