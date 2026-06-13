"""
MetaQuest Pathogen Risk Reporter - Professional Edition v5.0.0
================================================================
Publication-ready pathogen risk assessment with integrated evidence.

Author: MetaQuest Development Team
License: MIT
"""

import pandas as pd
import json
import re
from pathlib import Path
from typing import Dict, List, Optional
from collections import defaultdict

from .base_reporter import BaseReporter
from .validation_engine import ValidationEngine

from ..io.output_formatter import get_formatter
formatter = get_formatter()


class PathogenRiskReporter(BaseReporter):
    """
    Generate comprehensive pathogen risk reports integrating:
      1. Taxonomic risk (Bracken + pathogen organisms)
      2. Functional risk (annotations + pathogen DB hits)
      3. ML prediction risk (MetaQuest predictions)
    
    SAMPLE-AGNOSTIC: Works for clinical, environmental, food, or any sample type.
    """
    
    def __init__(self, output_dir: Path):
        super().__init__(output_dir)
        self.section_name = "Pathogen Risk Assessment"
        
        # Risk level indicators
        self.risk_emojis = {
            'High': '[HIGH]',
            'Moderate': '[MOD]',
            'Low': '[LOW]'
        }
    
    # ========================================================================
    # MAIN REPORT GENERATION
    # ========================================================================
    
    def generate_report(self,
                       integrated_risk: Dict,
                       taxonomic_risk: Dict,
                       functional_risk: Dict,
                       ml_risk: Dict,
                       pathogen_hits_file: Optional[Path] = None,
                       ml_predictions_file: Optional[Path] = None,
                       validator: Optional[ValidationEngine] = None) -> str:
        """Generate professional pathogen risk report."""
        
        report_parts = []
        
        # Header
        report_parts.append(self._generate_header())
        
        # 1. Executive Risk Summary
        report_parts.append(self._generate_risk_summary(integrated_risk))
        
        # 2. Section 1: Taxonomic Pathogenicity
        report_parts.append(self._generate_section1_taxonomic(taxonomic_risk))
        
        # 3. Section 2: Functional Pathogenicity Markers
        if validator:
            report_parts.append(self._generate_section2_functional(
                functional_risk, pathogen_hits_file, validator
            ))
        else:
            report_parts.append(self._generate_section2_functional(
                functional_risk, pathogen_hits_file
            ))
        
        # 4. Section 3: Machine Learning Predictions
        report_parts.append(self._generate_section3_ml(ml_risk, ml_predictions_file))
        
        # 5. Section 4: Integrated Risk Assessment
        report_parts.append(self._generate_section4_integrated(integrated_risk))
        
        # 6. Section 5: Clinical Interpretation
        report_parts.append(self._generate_section5_clinical(integrated_risk))
        
        return '\n\n'.join(report_parts)
    
    # ========================================================================
    # SECTION GENERATORS
    # ========================================================================
    
    def _generate_header(self) -> str:
        """Generate professional report header - ✓ USING BASE REPORTER."""
        
        # ✓ USING format_header()
        header = self.format_header("PATHOGEN RISK ASSESSMENT REPORT", width=80, style="double")
        subtitle = "MetaQuest Metagenomic Analysis Pipeline".center(80)
        separator = "═" * 80
        
        return f"{header}\n{subtitle}\n{separator}"
    
    def _generate_risk_summary(self, risk: Dict) -> str:
        """Generate executive risk summary - ✓ USING BASE REPORTER."""
        
        final_score = risk['final_score']
        risk_level = risk['risk_level']
        emoji = self.risk_emojis[risk_level]
        
        # Extract tier scores
        tier_scores = risk.get('tier_scores', {})
        taxonomic_score = tier_scores.get('taxonomic', 0)
        functional_score = tier_scores.get('functional', 0)
        ml_score = tier_scores.get('ml', 0)
        
        lines = []
        lines.append("")
        lines.append(self.format_section("EXECUTIVE SUMMARY", level=1))
        lines.append("")
        
        # ✓ USING format_key_value()
        summary_metrics = [
            ("Overall Risk Score:", f"{final_score:.1f}/100"),
            ("Risk Classification:", f"{emoji} {risk_level.upper()}"),
            ("", ""),
            ("Component Scores:", ""),
            ("  Taxonomic Risk:", f"{taxonomic_score:.1f}/100"),
            ("  Functional Risk:", f"{functional_score:.1f}/100"),
            ("  ML Prediction Risk:", f"{ml_score:.1f}/100"),
        ]
        lines.append(self.format_key_value(summary_metrics, key_width=30))
        lines.append("")
        
        # Risk modifiers
        multipliers = risk['multipliers_applied']
        if any(multipliers.values()):
            lines.append("Risk Modifiers Applied:")
            if multipliers.get('ml_pathogen_db'):
                lines.append("  • ML-Pathogen DB convergence detected (+20%)")
            if multipliers.get('virulence_transposase'):
                lines.append("  • Virulence genes near mobile elements (+15%)")
            lines.append("")
        
        return '\n'.join(lines)
    
    def _generate_section1_taxonomic(self, risk: Dict) -> str:
        """Section 1: Taxonomic pathogenicity assessment - ✓ USING BASE REPORTER."""
        
        lines = []
        lines.append("")
        lines.append(self.format_section("SECTION 1: TAXONOMIC PATHOGENICITY ASSESSMENT", level=1))
        lines.append("")
        
        # ✓ USING format_key_value()
        tax_metrics = [
            ("Risk Score:", f"{risk['score']:.1f}/100"),
            ("Pathogens Detected:", f"{len(risk['pathogens_detected'])} species"),
        ]
        
        if risk.get('total_pathogen_abundance', 0) > 0:
            tax_metrics.append((
                "Total Pathogen Abundance:", 
                f"{risk['total_pathogen_abundance']*100:.2f}%"
            ))
        
        lines.append(self.format_key_value(tax_metrics, key_width=30))
        lines.append("")
        
        if not risk['pathogens_detected']:
            lines.extend([
                "Finding: No known pathogenic species detected at significant levels.",
                "Sample composition indicates predominantly commensal organisms.",
            ])
            return '\n'.join(lines)
        
        # 1.1 Detected Pathogenic Species
        lines.append("1.1 Detected Pathogenic Species")
        lines.append("")
        
        # ✓ USING format_table()
        headers = ["Species", "Abundance", "Risk"]
        rows = []
        
        for pathogen in risk['pathogens_detected'][:15]:
            species = pathogen['species']
            if len(species) > 45:
                species = species[:42] + "..."
            abundance = f"{pathogen['abundance']*100:.2f}%"
            risk_level = pathogen['risk_level']
            emoji = self.risk_emojis[risk_level]
            rows.append([species, abundance, f"{emoji} {risk_level}"])
        
        lines.append(self.format_table(headers, rows, alignments=['left', 'right', 'left']))
        
        if len(risk['pathogens_detected']) > 15:
            remaining = len(risk['pathogens_detected']) - 15
            lines.append(f"\n[{remaining} additional species not shown]")
        
        # 1.2 Dominant Pathogen
        if risk.get('dominant_pathogen'):
            dom = risk['dominant_pathogen']
            lines.append("")
            lines.append("1.2 Dominant Pathogen")
            lines.append("")
            
            dom_metrics = [
                ("Species:", dom['species']),
                ("Abundance:", f"{dom['abundance']*100:.2f}%"),
                ("Risk Level:", f"{self.risk_emojis[dom['risk_level']]} {dom['risk_level']}"),
            ]
            lines.append(self.format_key_value(dom_metrics, key_width=20, indent=0))
        
        return '\n'.join(lines)
    
    def _generate_section2_functional(self,
                                      risk: Dict,
                                      pathogen_hits_file: Optional[Path],
                                      validator: Optional[ValidationEngine] = None) -> str:
        """Section 2: Functional pathogenicity markers - SAMPLE AGNOSTIC."""
        
        lines = []
        lines.append("")
        lines.append(self.format_section("SECTION 2: FUNCTIONAL PATHOGENICITY MARKERS", level=1))
        lines.append("")
        lines.append(f"Risk Score: {risk['score']:.1f}/100")
        lines.append("")
        
        details = risk['details']
        
        # 2.1 Virulence Factors
        virulence_count = len(details['virulence_factors'])
        lines.append("2.1 Virulence Factors")
        lines.append("")
        lines.append(f"Total detected: {virulence_count} genes")
        
        # DENSITY-BASED ASSESSMENT (NO Z-SCORES)
        stat_analysis = details.get('statistical_analysis', {})
        if stat_analysis and 'virulence_pct' in stat_analysis:
            virulence_pct = stat_analysis.get('virulence_pct', 0)
            virulence_interpretation = stat_analysis.get('virulence_interpretation', 'N/A')
            lines.append(f"Percentage of proteome: {virulence_pct:.2f}%")
            lines.append(f"Density assessment: {virulence_interpretation}")
            lines.append("")
        else:
            lines.append("")
        
        if virulence_count > 0:
            # Categorize virulence genes
            adhesion = [g for g in details['virulence_factors']
                       if any(kw in g.get('description', '').lower()
                            for kw in ['fimbri', 'pili', 'adhesin'])]
            toxins = [g for g in details['virulence_factors']
                     if any(kw in g.get('description', '').lower()
                          for kw in ['toxin', 'hemolysin'])]
            secretion = [g for g in details['virulence_factors']
                        if 'secretion' in g.get('description', '').lower()]
            
            if adhesion:
                lines.append(f"Adhesion & colonization: {len(adhesion)} genes")
            if toxins:
                lines.append(f"Toxin genes: {len(toxins)} genes")
            if secretion:
                lines.append(f"Secretion systems: {len(secretion)} genes")
            
            # Critical findings
            critical_toxins = [t for t in toxins
                             if any(ct in t.get('description', '').lower()
                                  for ct in ['shiga', 'cholera', 'anthrax', 'botulinum'])]
            
            if critical_toxins:
                lines.append("")
                lines.append("CRITICAL FINDING:")
                lines.append(f"[HIGH] {len(critical_toxins)} high-impact toxin gene(s) detected")
                for toxin in critical_toxins[:3]:
                    desc = toxin.get('description', '')[:60]
                    lines.append(f"  • {desc}")
        else:
            lines.append("Status: None detected")
        
        lines.append("")
        
        # 2.2 Antimicrobial Resistance Genes
        amr_count = len(details['amr_genes'])
        lines.append("2.2 Antimicrobial Resistance Genes")
        lines.append("")
        lines.append(f"Total detected: {amr_count} genes")
        
        # DENSITY-BASED ASSESSMENT
        if stat_analysis and 'amr_pct' in stat_analysis:
            amr_pct = stat_analysis.get('amr_pct', 0)
            amr_interpretation = stat_analysis.get('amr_interpretation', 'N/A')
            lines.append(f"Percentage of proteome: {amr_pct:.2f}%")
            lines.append(f"Density assessment: {amr_interpretation}")
            lines.append("")
        else:
            lines.append("")
        
        if amr_count > 0:
            unique_amr = self._get_unique_genes(details['amr_genes'], limit=5)
            lines.append("Representative examples:")
            for amr in unique_amr:
                gene_id = amr['gene_id'][:30]
                identity = amr.get('identity', 0)
                lines.append(f"  • {gene_id} ({identity:.1f}% identity)")
        else:
            lines.append("Status: None detected")
        
        lines.append("")
        
        # 2.3 Mobile Genetic Elements
        transposon_count = len(details['transposases'])
        lines.append("2.3 Mobile Genetic Elements")
        lines.append("")
        lines.append(f"Transposases detected: {transposon_count}")
        lines.append("  (Unique transposase genes only - see Functional Report for total hits)")
        
        # DENSITY-BASED ASSESSMENT
        if stat_analysis and 'transposase_pct' in stat_analysis:
            trans_pct = stat_analysis.get('transposase_pct', 0)
            trans_interpretation = stat_analysis.get('transposase_interpretation', 'N/A')
            lines.append(f"Percentage of proteome: {trans_pct:.2f}%")
            lines.append(f"Density assessment: {trans_interpretation}")
            lines.append("")
        else:
            lines.append("")
        
        if transposon_count > 5:
            lines.append("WARNING: High mobile element activity detected")
            lines.append("Risk: Facilitation of horizontal gene transfer")
            lines.append("Impact: Resistance/virulence genes may spread between bacteria")
            
            # Transposase families
            families = self._categorize_transposases(details['transposases'])
            if families:
                lines.append("")
                lines.append("Transposase families detected:")
                for family, count in sorted(families.items(), key=lambda x: x[1], reverse=True)[:5]:
                    lines.append(f"  • {family}: {count} gene(s)")
        elif transposon_count > 0:
            lines.append(f"Status: Low activity ({transposon_count} detected)")
        else:
            lines.append("Status: None detected")
        
        # 2.4 Species-Specific Gene Arsenal (if validator available)
        if validator and pathogen_hits_file:
            lines.append("")
            lines.append("2.4 Species-Specific Virulence Arsenal")
            lines.append("")
            lines.extend(self._generate_species_arsenal(risk, pathogen_hits_file, validator))
        
        # 2.5 Pathogen Database Hits
        if pathogen_hits_file and pathogen_hits_file.exists():
            lines.append("")
            lines.append("2.5 Pathogen Database Matches")
            lines.append("")
            lines.append(self._format_pathogen_hits_professional(pathogen_hits_file, details))
        
        return '\n'.join(lines)
    
    def _generate_species_arsenal(self, risk: Dict, pathogen_hits_file: Path,
                                  validator: ValidationEngine) -> List[str]:
        """Generate species-specific gene arsenal."""
        lines = []
        
        # Get gene-to-species mapping
        functional_file = self.output_dir / 'functional_annotations.tsv'
        prokka_gff = self.output_dir / 'prokka_annotation' / 'sample.gff'
        kraken_classified = self.output_dir / 'kraken_classified.txt'
        
        gene_species_map = validator.link_genes_to_species(
            functional_file, pathogen_hits_file, prokka_gff, kraken_classified
        )
        
        # Group genes by species
        species_genes = defaultdict(lambda: {
            'virulence': [], 'amr': [], 'transposases': []
        })
        
        details = risk['details']
        for gene_list, category in [
            (details['virulence_factors'], 'virulence'),
            (details['amr_genes'], 'amr'),
            (details['transposases'], 'transposases')
        ]:
            for gene_info in gene_list:
                gene_id = gene_info.get('gene_id', gene_info.get('query_id', ''))
                if gene_id in gene_species_map:
                    species = gene_species_map[gene_id]['species']
                    species_genes[species][category].append(gene_info)
        
        # Display top 5 species by gene count
        sorted_species = sorted(
            species_genes.items(),
            key=lambda x: sum(len(v) for v in x[1].values()),
            reverse=True
        )
        
        for species, genes in sorted_species[:5]:
            virulence_count = len(genes['virulence'])
            amr_count = len(genes['amr'])
            transposon_count = len(genes['transposases'])
            
            if virulence_count == 0 and amr_count == 0:
                continue
            
            lines.append(f"Species: {species}")
            lines.append("")
            
            if virulence_count > 0:
                lines.append(f"  Virulence factors: {virulence_count} genes")
                
                # Subcategorize
                adhesion = [g for g in genes['virulence']
                           if any(kw in g.get('description', '').lower()
                                for kw in ['fimbri', 'pili', 'adhesin'])]
                toxins = [g for g in genes['virulence']
                         if any(kw in g.get('description', '').lower()
                              for kw in ['toxin', 'hemolysin'])]
                
                if adhesion:
                    lines.append(f"    - Adhesion/colonization: {len(adhesion)}")
                if toxins:
                    lines.append(f"    - Toxins: {len(toxins)}")
                    for toxin in toxins[:2]:
                        desc = toxin.get('description', '')[:50]
                        lines.append(f"      {desc}")
            
            if amr_count > 0:
                lines.append(f"  AMR genes: {amr_count} genes")
            
            if transposon_count > 5:
                lines.append(f"  Mobile elements: {transposon_count} transposases")
                lines.append("  [!] WARNING: High mobility - genes may transfer")
            
            lines.append("")
        
        return lines
    
    def _generate_section3_ml(self, risk: Dict, ml_file: Optional[Path]) -> str:
        """Section 3: Machine learning predictions."""
        
        pathogenic_pct = (risk['pathogenic_count'] / risk['total_sequences'] * 100) \
            if risk['total_sequences'] > 0 else 0
        
        lines = []
        lines.append("")
        lines.append(self.format_section("SECTION 3: MACHINE LEARNING PATHOGEN PREDICTIONS", level=1))
        lines.append("")
        lines.append(f"Risk Score: {risk['score']:.1f}/100")
        lines.append("")
        lines.append("3.1 Prediction Summary")
        lines.append("")
        lines.append(f"Total sequences analyzed: {risk['total_sequences']}")
        lines.append(f"Pathogenic predictions: {risk['pathogenic_count']} ({pathogenic_pct:.1f}%)")
        # high_confidence_percentage is already a percentage value (e.g. 72.8)
        hc_pct = risk.get('high_confidence_percentage',
                          risk['high_confidence_pathogenic'] / risk['total_sequences'] * 100
                          if risk['total_sequences'] > 0 else 0)
        lines.append(f"High-confidence pathogenic: {risk['high_confidence_pathogenic']} ({hc_pct:.1f}%)")
        lines.append(f"Non-pathogenic predictions: {risk['total_sequences'] - risk['pathogenic_count']}")
        lines.append("")
        
        if risk['high_confidence_pathogenic'] == 0:
            lines.append("Finding: No high-confidence pathogenic sequences detected.")
            lines.append("ML model assessment indicates low pathogenic potential.")
            return '\n'.join(lines)
        
        # 3.2 Confidence Distribution
        lines.append("3.2 Confidence Distribution")
        lines.append("")
        lines.append(self._generate_confidence_distribution(risk))
        
        # 3.3 High-Risk Predictions
        lines.append("")
        lines.append("3.3 High-Risk Pathogenic Predictions (>85% confidence)")
        lines.append("")
        
        high_risk = [p for p in risk['details'] if p['confidence'] > 0.85]
        
        if high_risk:
            lines.append(f"{'Sequence ID':<25} {'Confidence':>12} {'Path.Prob':>12} {'Length':>10}")
            lines.append("-" * 65)
            
            for pred in high_risk[:10]:
                seq_id = pred['sequence_id'][:23]
                conf = f"{pred['confidence']*100:.1f}%"
                prob = f"{pred['pathogenic_probability']*100:.1f}%"
                length = str(pred['length'])
                lines.append(f"{seq_id:<25} {conf:>12} {prob:>12} {length:>10}")
            
            if len(high_risk) > 10:
                lines.append(f"\n[{len(high_risk) - 10} additional predictions not shown]")
        else:
            lines.append("No predictions exceed 85% confidence threshold")
        
        return '\n'.join(lines)
    
    def _generate_section4_integrated(self, risk: Dict) -> str:
        """Section 4: Integrated risk assessment - ✓ USING BASE REPORTER."""
        
        lines = []
        lines.append("")
        lines.append(self.format_section("SECTION 4: INTEGRATED RISK ASSESSMENT", level=1))
        lines.append("")
        
        # 4.1 Risk Score Synthesis
        lines.append("4.1 Risk Score Synthesis")
        lines.append("")
        lines.append(f"Base risk score: {risk['base_score']:.1f}/100")
        lines.append(f"Final risk score: {risk['final_score']:.1f}/100")
        lines.append(f"Risk classification: {self.risk_emojis[risk['risk_level']]} " +
                    f"{risk['risk_level'].upper()}")
        
        # Score calculation transparency
        lines.append("")
        lines.append("Score Calculation Method:")
        lines.append(f"  Weight strategy: {risk.get('weight_strategy', 'standard')}")
        lines.append("")
        
        tier_scores = risk.get('tier_scores', {})
        weights_used = risk.get('weights_used', {})
        
        if tier_scores and weights_used:
            lines.append("Component Contributions:")
            
            taxonomic_contrib = tier_scores.get('taxonomic', 0) * weights_used.get('taxonomic', 0)
            functional_contrib = tier_scores.get('functional', 0) * weights_used.get('functional', 0)
            ml_contrib = tier_scores.get('ml', 0) * weights_used.get('ml', 0)
            
            lines.append(f"  Taxonomic ({weights_used.get('taxonomic', 0)*100:.0f}%): " +
                        f"{tier_scores.get('taxonomic', 0):.1f} × {weights_used.get('taxonomic', 0):.2f} = {taxonomic_contrib:.1f} points")
            lines.append(f"  Functional ({weights_used.get('functional', 0)*100:.0f}%): " +
                        f"{tier_scores.get('functional', 0):.1f} × {weights_used.get('functional', 0):.2f} = {functional_contrib:.1f} points")
            lines.append(f"  ML Prediction ({weights_used.get('ml', 0)*100:.0f}%): " +
                        f"{tier_scores.get('ml', 0):.1f} × {weights_used.get('ml', 0):.2f} = {ml_contrib:.1f} points")
            lines.append("")
        
        multipliers = risk.get('multipliers_applied', {})
        if multipliers:
            lines.append("Risk Multipliers Applied:")
            if multipliers.get('ml_pathogen_db'):
                lines.append(f"  × {multipliers['ml_pathogen_db']:.2f} (ML-PathogenDB convergence)")
            if multipliers.get('virulence_transposase'):
                lines.append(f"  × {multipliers['virulence_transposase']:.2f} (Virulence near mobile elements)")
            lines.append("")
        
        # 4.2 Convergent Risk Factors
        lines.append("4.2 Convergent Risk Factors")
        lines.append("")
        
        convergent = risk.get('convergent_risks', {})
        ml_pathogen_matches = convergent.get('ml_pathogen_db_matches', [])
        
        if ml_pathogen_matches:
            lines.append("The following genes are flagged by multiple detection methods:")
            lines.append("")
            
            for match in ml_pathogen_matches[:5]:
                lines.append(f"Gene: {match['gene_id']}")
                lines.append(f"  ML confidence: {match['ml_confidence']*100:.1f}% pathogenic")
                lines.append(f"  Pathogen DB match: {match['pathogen_db_hit'][:60]}")
                lines.append("  Priority: High - validation recommended")
                lines.append("")
        else:
            lines.append("No convergent risks detected.")
            lines.append("Finding: ML predictions and pathogen database hits do not overlap.")
            lines.append("Interpretation: Independent evidence sources provide complementary information.")
            lines.append("")
        
        # 4.3 Computational Interpretation
        lines.append("4.3 Computational Interpretation")
        lines.append("")
        
        interpretation = risk.get('interpretation', 'No interpretation available')
        lines.append(self._wrap_text(interpretation, width=76, indent=0))
        
        return '\n'.join(lines)
    
    def _generate_section5_clinical(self, risk: Dict) -> str:
        """Section 5: Sample-agnostic interpretation and recommendations."""
        
        risk_level = risk['risk_level']
        
        guides = {
            'Low': {
                'meaning': 'Sample composition dominated by organisms with minimal pathogenic markers detected.',
                'action': 'Routine monitoring recommended. Findings suggest low pathogenic potential based on metagenomic signatures.',
                'context': 'Metagenomic profile indicates predominantly commensal or environmental organisms. Interpret findings within sample context and study objectives.'
            },
            'Moderate': {
                'meaning': 'Opportunistic or conditional pathogens present and/or moderate pathogenic marker content detected.',
                'action': 'Contextual evaluation recommended. Consider sample origin, environmental conditions, and study objectives.',
                'context': 'Results suggest presence of organisms with pathogenic potential. Relevance depends on sample type (clinical, environmental, food, etc.) and specific research questions. Targeted validation may be warranted.'
            },
            'High': {
                'meaning': 'Known pathogens at significant abundance and/or multiple pathogenicity markers detected.',
                'action': 'Detailed characterization recommended. Findings warrant investigation based on sample context and biosafety considerations.',
                'context': 'Results indicate presence of organisms with established virulence mechanisms. For clinical samples, correlation with symptoms is essential. For environmental/food samples, consider contamination sources and public health implications.'
            }
        }
        
        guide = guides[risk_level]
        
        lines = []
        lines.append("")
        lines.append(self.format_section("SECTION 5: MICROBIOLOGICAL INTERPRETATION", level=1))
        lines.append("")
        lines.append(f"Risk Classification: {self.risk_emojis[risk_level]} {risk_level.upper()} " +
                    f"({risk['final_score']:.0f}/100)")
        lines.append("")
        lines.append("5.1 Interpretation")
        lines.append("")
        
        lines.append(self._wrap_text(guide['meaning'], width=76))
        
        lines.append("")
        lines.append("5.2 Recommended Action")
        lines.append("")
        
        lines.append(self._wrap_text(guide['action'], width=76))
        
        lines.append("")
        lines.append("5.3 Contextual Guidance")
        lines.append("")
        
        lines.append(self._wrap_text(guide['context'], width=76))
        
        # 5.4 Additional Considerations
        lines.append("")
        lines.append("5.4 Additional Considerations")
        lines.append("")
        
        considerations = []
        
        if risk.get('tier_scores', {}).get('taxonomic', 0) > 50:
            considerations.append(
                "Known pathogenic species detected - evaluate within sample context and study objectives"
            )
        
        if risk.get('tier_scores', {}).get('functional', 0) > 50:
            considerations.append(
                "Significant mobile genetic elements present - potential for horizontal gene transfer"
            )
        
        if risk.get('tier_scores', {}).get('ml', 0) > 40:
            considerations.append(
                "Machine learning predictions indicate pathogenic signatures - validation with orthogonal methods recommended"
            )
        

        convergent = risk.get('convergent_risks', {})
        if convergent.get('ml_pathogen_db_matches'):
            considerations.append(
                "Multiple detection methods converge on specific genes - high-priority findings"
            )
        
        if not considerations:
            considerations.append(
                "Standard clinical correlation and follow-up protocols apply"
            )
        
        for consideration in considerations:
            lines.append(f"  • {consideration}")
        
        # 5.5 Important Note
        lines.append("")
        lines.append("5.5 Important Note")
        lines.append("")
        lines.append("This analysis provides computational predictions based on metagenomic sequencing.")
        lines.append("Results should be interpreted within the specific context of sample origin,")
        lines.append("study objectives, and research questions. Validation with complementary methods")
        lines.append("(culture, PCR, microscopy) is recommended for definitive characterization.")
        lines.append("")
        lines.append("=" * 80)
        
        return '\n'.join(lines)
    
    # ========================================================================
    # HELPER METHODS (ALL PRESERVED FROM ORIGINAL)
    # ========================================================================
    
    def _get_unique_genes(self, gene_list: List[Dict], limit: int = 5) -> List[Dict]:
        """Get unique genes (no duplicates) from list."""
        unique = []
        seen_ids = set()
        
        for gene in gene_list:
            gene_id = gene.get('gene_id', gene.get('query_id', ''))
            if gene_id not in seen_ids:
                unique.append(gene)
                seen_ids.add(gene_id)
                if len(unique) >= limit:
                    break
        
        return unique
    
    def _categorize_transposases(self, transposases: List[Dict]) -> Dict[str, int]:
        """Categorize transposases by family."""
        families = {}
        
        for t in transposases:
            desc = t.get('description', '')
            if 'IS66' in desc:
                families['IS66'] = families.get('IS66', 0) + 1
            elif 'IS30' in desc:
                families['IS30'] = families.get('IS30', 0) + 1
            elif 'IS982' in desc:
                families['IS982'] = families.get('IS982', 0) + 1
            elif 'IS256' in desc:
                families['IS256'] = families.get('IS256', 0) + 1
            elif 'IS3' in desc:
                families['IS3'] = families.get('IS3', 0) + 1
            else:
                families['Other'] = families.get('Other', 0) + 1
        
        return families
    
    def _format_pathogen_hits_professional(self, pathogen_file: Path, details: Dict) -> str:
        """Format pathogen database hits as professional table."""
        
        # First try to get hits from details
        all_hits = (
            details.get('virulence_factors', []) +
            details.get('amr_genes', []) +
            details.get('secretion_systems', []) +
            details.get('other_pathogen_markers', [])
        )
        
        if all_hits:
            unique_hits = self._get_unique_genes(all_hits, limit=15)
            lines = [
                f"{'Query Gene':<20} {'Subject ID':<20} {'Identity':>10} {'Organism':<30}",
                "-" * 85,
            ]
            
            for hit in unique_hits:
                query = (hit.get('gene_id') or '')[:18]
                subject = (hit.get('subject_id') or 'Unknown')[:18]
                identity = f"{hit.get('identity', 0):.1f}%"
                organism = (hit.get('organism') or 'Unknown')[:28]
                lines.append(f"{query:<20} {subject:<20} {identity:>10} {organism:<30}")
            
            if len(all_hits) > 15:
                lines.append(f"\n[{len(all_hits) - 15} additional matches - see pathogen_results.tsv]")
            
            return '\n'.join(lines)
        
        # Fallback: Try to load from file
        if not pathogen_file or not pathogen_file.exists():
            return (
                "No pathogen database hits file available.\n"
                "Note: Pathogen database search was not performed or file is missing."
            )
        
        try:
            df = pd.read_csv(pathogen_file, sep='\t', header=None)
            if df.empty:
                return (
                    "No significant pathogen database matches detected.\n"
                    "Note: Pathogen database search completed with no hits meeting\n"
                    "      significance criteria (E-value < 1e-5, identity > 80%)."
                )
            
            # Column names
            col_names = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch',
                        'gapopen', 'qstart', 'qend', 'sstart', 'send',
                        'evalue', 'bitscore', 'qlen', 'slen', 'stitle']
            
            if len(df.columns) >= len(col_names):
                df.columns = col_names
            elif len(df.columns) == 15:
                df.columns = col_names
            else:
                df.columns = col_names[:len(df.columns)]
            
            # Check if we have minimum required columns
            required_cols = ['qseqid', 'sseqid', 'pident']
            if not all(col in df.columns for col in required_cols):
                return (
                    "No significant pathogen database matches detected.\n"
                    "Note: Pathogen database search completed with no hits meeting\n"
                    "      significance criteria (E-value < 1e-5, identity > 80%)."
                )
            
            lines = [
                f"{'Query Gene':<20} {'Subject ID':<20} {'Identity':>10} {'Organism':<30}",
                "-" * 85,
            ]
            
            for _, row in df.head(15).iterrows():
                query = str(row['qseqid'])[:18]
                subject = str(row['sseqid'])[:18]
                identity = f"{row['pident']:.1f}%"
                
                # Extract organism from stitle
                title = str(row.get('stitle', ''))
                org_match = re.findall(r'\[([^\]]+)\]', title)
                organism = org_match[-1][:28] if org_match else 'Unknown'
                
                lines.append(
                    f"{query:<20} {subject:<20} {identity:>10} {organism:<30}"
                )
            
            if len(df) > 15:
                lines.append(f"\n[{len(df) - 15} additional matches - see pathogen_results.tsv]")
            
            return '\n'.join(lines)
            
        except pd.errors.EmptyDataError:
            return (
                "No significant pathogen database matches detected.\n"
                "Note: Pathogen database search was performed but yielded no hits above\n"
                "      significance thresholds (E-value < 1e-5, identity > 80%)."
            )
        except Exception as e:
            return (
                "No significant pathogen database matches detected.\n"
                "Note: Pathogen database search completed with no hits meeting\n"
                "      significance criteria (E-value < 1e-5, identity > 80%)."
            )
    
    def _generate_confidence_distribution(self, risk: Dict) -> str:
        """Generate confidence distribution visualization."""
        details = risk['details']
        if not details:
            return "No prediction data available"
        
        # Group by confidence ranges
        ranges = {
            '70-80%': 0,
            '80-90%': 0,
            '90-95%': 0,
            '95-100%': 0
        }
        
        for pred in details:
            conf = pred['confidence']
            if 0.7 <= conf < 0.8:
                ranges['70-80%'] += 1
            elif 0.8 <= conf < 0.9:
                ranges['80-90%'] += 1
            elif 0.9 <= conf < 0.95:
                ranges['90-95%'] += 1
            elif conf >= 0.95:
                ranges['95-100%'] += 1
        
        # Create bar chart
        max_count = max(ranges.values()) if ranges.values() else 1
        lines = []
        
        for range_label, count in ranges.items():
            bar_length = int((count / max_count) * 40) if max_count > 0 else 0
            bar = '█' * bar_length
            lines.append(f"  {range_label:>10} | {bar} {count}")
        
        return '\n'.join(lines)
    
    def _wrap_text(self, text: str, width: int = 76, indent: int = 0) -> str:
        """Wrap text to specified width with optional indent."""
        words = text.split()
        lines = []
        current_line = []
        current_length = 0
        
        for word in words:
            if current_length + len(word) + 1 <= width:
                current_line.append(word)
                current_length += len(word) + 1
            else:
                lines.append(' ' * indent + ' '.join(current_line))
                current_line = [word]
                current_length = len(word)
        
        if current_line:
            lines.append(' ' * indent + ' '.join(current_line))
        
        return '\n'.join(lines)
    
    # ========================================================================
    # VISUALIZATION & EXPORT METHODS (Original API preserved)
    # ========================================================================
    
    def generate_visualization_data(self,
                                    integrated_risk: Dict,
                                    taxonomic_risk: Dict,
                                    functional_risk: Dict,
                                    ml_risk: Dict) -> Dict:
        """Generate data structure for visualization module."""
        
        return {
            'risk_scores': {
                'final': integrated_risk['final_score'],
                'taxonomic': taxonomic_risk['score'],
                'functional': functional_risk['score'],
                'ml': ml_risk['score'],
                'risk_level': integrated_risk['risk_level']
            },
            'pathogen_counts': {
                'total': len(taxonomic_risk['pathogens_detected']),
                'high_risk': len([p for p in taxonomic_risk['pathogens_detected']
                                 if p['risk_level'] == 'High']),
                'moderate_risk': len([p for p in taxonomic_risk['pathogens_detected']
                                     if p['risk_level'] == 'Moderate']),
                'low_risk': len([p for p in taxonomic_risk['pathogens_detected']
                                if p['risk_level'] == 'Low'])
            },
            'functional_markers': {
                'virulence': len(functional_risk['details']['virulence_factors']),
                'amr': len(functional_risk['details']['amr_genes']),
                'transposases': len(functional_risk['details']['transposases']),
                'secretion': len(functional_risk['details']['secretion_systems'])
            },
            'ml_predictions': {
                'total': ml_risk['total_sequences'],
                'pathogenic': ml_risk['pathogenic_count'],
                'high_conf_pathogenic': ml_risk['high_confidence_pathogenic'],
                'confidence_distribution': self._get_confidence_distribution(ml_risk)
            },
            'convergent_risks': integrated_risk['convergent_risks'],
            'interpretation': integrated_risk['interpretation']
        }
    
    def _get_confidence_distribution(self, ml_risk: Dict) -> Dict:
        """Get confidence score distribution for visualization."""
        details = ml_risk['details']
        distribution = {
            '0.7-0.8': 0,
            '0.8-0.9': 0,
            '0.9-0.95': 0,
            '0.95-1.0': 0
        }
        
        for pred in details:
            conf = pred['confidence']
            if 0.7 <= conf < 0.8:
                distribution['0.7-0.8'] += 1
            elif 0.8 <= conf < 0.9:
                distribution['0.8-0.9'] += 1
            elif 0.9 <= conf < 0.95:
                distribution['0.9-0.95'] += 1
            elif conf >= 0.95:
                distribution['0.95-1.0'] += 1
        
        return distribution
    
    def export_risk_summary(self, integrated_risk: Dict, output_file: Path):
        """Export risk summary to JSON for dashboard."""
        summary = {
            'risk_assessment': {
                'final_score': integrated_risk['final_score'],
                'risk_level': integrated_risk['risk_level'],
                'base_score': integrated_risk['base_score'],
                'tier_scores': integrated_risk.get('tier_scores', {}),
                'multipliers_applied': integrated_risk['multipliers_applied'],
                'interpretation': integrated_risk['interpretation']
            },
            'timestamp': pd.Timestamp.now().isoformat()
        }
        
        with open(output_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        formatter.info(f"Risk summary exported to {output_file}")
