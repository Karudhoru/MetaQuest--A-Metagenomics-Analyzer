"""
MetaQuest Pathogen Risk Reporter
Three-tier integrated pathogen risk assessment with ML predictions
"""

import pandas as pd
import json
from pathlib import Path
from typing import Dict, List, Optional
from .base_reporter import BaseReporter


class PathogenRiskReporter(BaseReporter):
    """
    Generate comprehensive pathogen risk reports integrating:
    1. Taxonomic risk (Bracken + pathogen organisms)
    2. Functional risk (annotations + pathogen DB hits)
    3. ML prediction risk (MetaQuest predictions)
    """
    
    def __init__(self, output_dir: Path):
        super().__init__(output_dir)
        self.section_name = "Pathogen Risk Assessment"
        
        # Risk level styling
        self.risk_emojis = {
            'High': '🔴',
            'Moderate': '🟡',
            'Low': '🟢'
        }
        
        self.risk_gauges = {
            'High': '[========💉========]',
            'Moderate': '[=====⚠️=====]',
            'Low': '[====✓====]'
        }
    
    def generate_report(self,
                       integrated_risk: Dict,
                       taxonomic_risk: Dict,
                       functional_risk: Dict,
                       ml_risk: Dict,
                       pathogen_hits_file: Optional[Path] = None,
                       ml_predictions_file: Optional[Path] = None) -> str:
        """
        Generate complete three-tier pathogen risk report.
        
        Args:
            integrated_risk: Final integrated risk scores
            taxonomic_risk: Taxonomic risk assessment
            functional_risk: Functional risk assessment
            ml_risk: ML prediction risk assessment
            pathogen_hits_file: Path to pathogen DB hits
            ml_predictions_file: Path to ML predictions JSON
        
        Returns:
            Formatted report string
        """
        self.print_section_header("PATHOGEN RISK ASSESSMENT")
        
        report_parts = []
        
        # 1. Overall Risk Score Dashboard
        report_parts.append(self._generate_risk_dashboard(integrated_risk))
        
        # 2. Tier 1: Taxonomic Pathogenicity
        report_parts.append(self._generate_tier1_taxonomic(taxonomic_risk))
        
        # 3. Tier 2: Functional Pathogenicity Markers
        report_parts.append(self._generate_tier2_functional(functional_risk, pathogen_hits_file))
        
        # 4. Tier 3: ML Predictions
        report_parts.append(self._generate_tier3_ml(ml_risk, ml_predictions_file))
        
        # 5. Integrated Risk Assessment
        report_parts.append(self._generate_integrated_assessment(integrated_risk))
        
        # 6. Clinical Interpretation
        report_parts.append(self._generate_interpretation_guide(integrated_risk))
        
        return '\n\n'.join(report_parts)
    
    def _generate_risk_dashboard(self, risk: Dict) -> str:
        """Generate executive risk score dashboard."""
        final_score = risk['final_score']
        risk_level = risk['risk_level']
        emoji = self.risk_emojis[risk_level]
        gauge = self.risk_gauges[risk_level]
        
        lines = [
            "┌" + "─" * 68 + "┐",
            "│" + " " * 20 + "OVERALL PATHOGENICITY RISK" + " " * 22 + "│",
            "│" + " " * 68 + "│",
            f"│{' ' * 20}{gauge} {final_score:.0f}/100{' ' * (23 - len(str(int(final_score))))}│",
            f"│{' ' * 25}{emoji} {risk_level.upper()} RISK{' ' * (27 - len(risk_level))}│",
            "│" + " " * 68 + "│",
            "├" + "─" * 68 + "┤",
            f"│  Taxonomic Risk:      {risk['taxonomic_score']:.0f}/100  (Known pathogens present)" + " " * (16 - len(str(int(risk['taxonomic_score'])))) + "│",
            f"│  Functional Risk:     {risk['functional_score']:.0f}/100  (Mobile elements + markers)" + " " * (15 - len(str(int(risk['functional_score'])))) + "│",
            f"│  ML Prediction Risk:  {risk['ml_score']:.0f}/100  ({risk.get('ml_score', 0):.0f}% sequences flagged)" + " " * (10 - len(str(int(risk.get('ml_score', 0))))) + "│",
            "└" + "─" * 68 + "┘",
        ]
        
        return '\n'.join(lines)
    
    def _generate_tier1_taxonomic(self, risk: Dict) -> str:
        """Generate Tier 1: Taxonomic pathogenicity section."""
        lines = [
            "",
            "=" * 70,
            "🦠 TIER 1: TAXONOMIC PATHOGENICITY ASSESSMENT",
            "=" * 70,
            "",
            f"Risk Score: {risk['score']:.2f}/100",
            f"Pathogens Detected: {len(risk['pathogens_detected'])} species",
        ]
        
        if risk['total_pathogen_abundance'] > 0:
            lines.append(
                f"Total Pathogen Abundance: {risk['total_pathogen_abundance']*100:.2f}%"
            )
        
        if not risk['pathogens_detected']:
            lines.extend([
                "",
                "✅ No known pathogenic species detected at significant levels.",
                "   Sample dominated by commensal organisms.",
            ])
            return '\n'.join(lines)
        
        # Pathogen risk matrix table
        lines.extend([
            "",
            "PATHOGEN RISK MATRIX:",
            "",
            f"{'Species':<40} {'Abundance':>12} {'Risk':>10} {'Action':>10}",
            "-" * 75,
        ])
        
        for pathogen in risk['pathogens_detected'][:10]:
            species = pathogen['species'][:38] + '..' if len(pathogen['species']) > 40 else pathogen['species']
            abundance = f"{pathogen['abundance']*100:.2f}%"
            risk_level = pathogen['risk_level']
            emoji = self.risk_emojis[risk_level]
            
            # Determine action
            if risk_level == 'High':
                action = "Monitor"
            elif risk_level == 'Moderate':
                action = "Observe"
            else:
                action = "Routine"
            
            lines.append(
                f"{species:<40} {abundance:>12} {emoji} {risk_level:>6} {action:>10}"
            )
        
        if len(risk['pathogens_detected']) > 10:
            lines.append(f"\n... and {len(risk['pathogens_detected']) - 10} more pathogenic species")
        
        # Highlight dominant pathogen
        if risk['dominant_pathogen']:
            dom = risk['dominant_pathogen']
            lines.extend([
                "",
                f"⚠️  DOMINANT PATHOGEN:",
                f"   {dom['species']} - {dom['abundance']*100:.2f}% abundance",
                f"   Risk Level: {dom['risk_level']}",
            ])
        
        return '\n'.join(lines)
    
    def _generate_tier2_functional(self, risk: Dict, pathogen_hits_file: Optional[Path]) -> str:
        """Generate Tier 2: Functional pathogenicity markers section."""
        lines = [
            "",
            "=" * 70,
            "🧬 TIER 2: FUNCTIONAL PATHOGENICITY MARKERS",
            "=" * 70,
            "",
            f"Risk Score: {risk['score']:.2f}/100",
        ]
        
        details = risk['details']
        
        # Virulence factors
        virulence_count = len(details['virulence_factors'])

        lines.extend([
            "",
            "┌──────────────────────────────────────────────┐",
            "│ 🦠 VIRULENCE FACTORS                         │",
            "├──────────────────────────────────────────────┤",
        ])

        if virulence_count > 0:
            lines.append(f"│ ✅ Detected: {virulence_count} genes" + " " * (28 - len(str(virulence_count))) + "│")
            
            # *** FIX: Show diverse examples, not just first 3 ***
            unique_examples = []
            seen_organisms = set()
            
            for vf in details['virulence_factors']:
                organism = vf.get('organism', 'Unknown')
                if organism not in seen_organisms and len(unique_examples) < 3:
                    unique_examples.append(vf)
                    seen_organisms.add(organism)
            
            for vf in unique_examples:
                gene_id = vf['gene_id'][:25]
                lines.append(f"│   • {gene_id}" + " " * (43 - len(gene_id)) + "│")
        else:
            lines.append("│ ❌ Not detected" + " " * 30 + "│")
        
        lines.append("└──────────────────────────────────────────────┘")
        
        # AMR genes
        amr_count = len(details['amr_genes'])

        lines.extend([
            "",
            "┌──────────────────────────────────────────────┐",
            "│ 💊 ANTIMICROBIAL RESISTANCE GENES            │",
            "├──────────────────────────────────────────────┤",
        ])

        if amr_count > 0:
            lines.append(f"│ ✅ Detected: {amr_count} genes" + " " * (28 - len(str(amr_count))) + "│")
            
            # *** FIX: Show REAL examples, not duplicates ***
            unique_examples = []
            seen_genes = set()
            
            for amr in details['amr_genes']:
                gene_id = amr['gene_id']
                if gene_id not in seen_genes and len(unique_examples) < 3:
                    unique_examples.append(amr)
                    seen_genes.add(gene_id)
            
            for amr in unique_examples:
                gene_id = amr['gene_id'][:25]
                lines.append(f"│   • {gene_id}" + " " * (43 - len(gene_id)) + "│")
        else:
            lines.append("│ ❌ Not detected" + " " * 30 + "│")
        
        lines.append("└──────────────────────────────────────────────┘")
        
        # Mobile genetic elements (CRITICAL!)
        lines.extend([
            "",
            "┌──────────────────────────────────────────────┐",
            "│ 🔄 MOBILE GENETIC ELEMENTS                   │",
            "├──────────────────────────────────────────────┤",
        ])
        
        transposon_count = len(details['transposases'])
        
        if transposon_count > 5:
            lines.extend([
                f"│ ⚠️  HIGH ACTIVITY: {transposon_count} transposases" + " " * (20 - len(str(transposon_count))) + "│",
                "│                                              │",
                "│ Risk: These elements can mobilize            │",
                "│ resistance/virulence genes between bacteria  │",
            ])
        elif transposon_count > 0:
            lines.append(f"│ ✅ Detected: {transposon_count} transposases" + " " * (27 - len(str(transposon_count))) + "│")
        else:
            lines.append("│ ❌ Not detected" + " " * 30 + "│")
        
        lines.append("└──────────────────────────────────────────────┘")
        
        # Transposase families
        if transposon_count > 0:
            families = {}
            for t in details['transposases']:
                # Extract IS family from description
                desc = t['description']
                if 'IS66' in desc:
                    families['IS66'] = families.get('IS66', 0) + 1
                elif 'IS30' in desc:
                    families['IS30'] = families.get('IS30', 0) + 1
                elif 'IS982' in desc:
                    families['IS982'] = families.get('IS982', 0) + 1
                elif 'IS256' in desc:
                    families['IS256'] = families.get('IS256', 0) + 1
                else:
                    families['Other'] = families.get('Other', 0) + 1
            
            if families:
                lines.extend([
                    "",
                    "Transposase Families Detected:",
                ])
                for family, count in sorted(families.items(), key=lambda x: x[1], reverse=True):
                    lines.append(f"  ├─ {family}: {count} gene(s)")
        
        # Pathogenicity-associated genes table
        if pathogen_hits_file and pathogen_hits_file.exists():
            lines.extend([
                "",
                "PATHOGENICITY-ASSOCIATED GENES:",
                "",
                self._format_pathogen_hits_table(pathogen_hits_file, details),
            ])
        
        return '\n'.join(lines)
    
    def _generate_tier3_ml(self, risk: Dict, ml_file: Optional[Path]) -> str:
        """Generate Tier 3: ML predictions section."""
        # Calculate percentage safely
        pathogenic_pct = (risk['pathogenic_count'] / risk['total_sequences'] * 100) if risk['total_sequences'] > 0 else 0

        lines = [
            "",
            "=" * 70,
            "🤖 TIER 3: MACHINE LEARNING PATHOGEN PREDICTIONS",
            "=" * 70,
            "",
            f"Risk Score: {risk['score']:.2f}/100",
            "",
            "SUMMARY:",
            f"  Total Sequences Analyzed:     {risk['total_sequences']}",
            f"  Pathogenic Predictions:       {risk['pathogenic_count']} ({pathogenic_pct:.1f}%)",
            f"  High Confidence Pathogenic:   {risk['high_confidence_pathogenic']} ({risk['high_confidence_proportion']*100:.1f}%)",
            f"  Non-pathogenic:               {risk['total_sequences'] - risk['pathogenic_count']}"
        ]
        
        if risk['high_confidence_pathogenic'] == 0:
            lines.extend([
                "",
                "✅ No high-confidence pathogenic sequences detected.",
                "   ML model suggests low pathogenic potential.",
            ])
            return '\n'.join(lines)
        
        # Confidence distribution visualization
        lines.extend([
            "",
            "CONFIDENCE DISTRIBUTION:",
            "",
            self._generate_confidence_viz(risk),
        ])
        
        # High-risk predictions table
        lines.extend([
            "",
            "🚨 HIGH-RISK PATHOGENIC PREDICTIONS (Confidence >85%):",
            "",
        ])
        
        high_risk_preds = [p for p in risk['details'] if p['confidence'] > 0.85]
        
        if high_risk_preds:
            lines.append(f"{'Sequence ID':<20} {'Confidence':>12} {'Path. Prob':>12} {'Length':>10}")
            lines.append("-" * 60)
            
            for pred in high_risk_preds[:10]:
                seq_id = pred['sequence_id']
                conf = f"{pred['confidence']*100:.1f}%"
                prob = f"{pred['pathogenic_probability']*100:.1f}%"
                length = str(pred['length'])
                
                lines.append(f"{seq_id:<20} {conf:>12} {prob:>12} {length:>10}")
        else:
            lines.append("No predictions with confidence >85%")
        
        # Load full predictions for cross-reference
        if ml_file and ml_file.exists():
            lines.extend([
                "",
                "CROSS-REFERENCE WITH FUNCTIONAL ANNOTATIONS:",
                "",
                self._cross_reference_ml_functional(ml_file, risk),
            ])
        
        return '\n'.join(lines)
    
    def _generate_integrated_assessment(self, risk: Dict) -> str:
        """Generate integrated risk assessment with convergence analysis."""
        lines = [
            "",
            "=" * 70,
            "🔬 INTEGRATED RISK ASSESSMENT",
            "=" * 70,
            "",
            "RISK SCORE BREAKDOWN:",
            f"  Base Score:        {risk['base_score']:.2f}/100",
            f"  Final Score:       {risk['final_score']:.2f}/100",
            f"  Risk Level:        {self.risk_emojis[risk['risk_level']]} {risk['risk_level'].upper()}",
        ]
        
        # Show multipliers if applied
        multipliers = risk['multipliers_applied']
        if any(multipliers.values()):
            lines.extend([
                "",
                "MULTIPLIERS APPLIED:",
            ])
            if multipliers['ml_pathogen_db']:
                lines.append("  ✅ ML + Pathogen DB convergence (×1.2)")
            if multipliers['virulence_transposase']:
                lines.append("  ✅ Virulence near transposase (×1.15)")
        
        # Convergent risks
        convergent = risk['convergent_risks']
        if convergent['ml_pathogen_db_matches']:
            lines.extend([
                "",
                "⚠️  CONVERGENT RISKS (Multiple tiers flag same genes):",
                "",
            ])
            
            for match in convergent['ml_pathogen_db_matches'][:5]:
                lines.extend([
                    f"Gene: {match['gene_id']}",
                    f"  ML Confidence: {match['ml_confidence']*100:.1f}% pathogenic",
                    f"  Pathogen DB: {match['pathogen_db_hit'][:60]}",
                    f"  → HIGH PRIORITY for validation",
                    "",
                ])
        
        # AI interpretation
        lines.extend([
            "",
            "🤖 INTELLIGENT INTERPRETATION:",
            "",
            self._wrap_text(risk['interpretation'], width=68, indent=2),
        ])
        
        return '\n'.join(lines)
    
    def _generate_interpretation_guide(self, risk: Dict) -> str:
        """Generate clinician-friendly interpretation guide."""
        risk_level = risk['risk_level']
        
        guides = {
            'Low': {
                'meaning': 'Dominated by commensal organisms with no significant pathogenic markers',
                'action': 'Routine monitoring is sufficient',
                'icon': '🟢',
                'clinical': 'No specific interventions indicated based on metagenomic profile'
            },
            'Moderate': {
                'meaning': 'Opportunistic pathogens present or some pathogenic markers detected',
                'action': 'Clinical correlation and monitoring for symptoms recommended',
                'icon': '🟡',
                'clinical': 'Consider patient history and clinical presentation before intervention'
            },
            'High': {
                'meaning': 'Known pathogens abundant or multiple pathogenicity markers present',
                'action': 'Clinical correlation required - targeted investigation strongly recommended',
                'icon': '🔴',
                'clinical': 'Consider antimicrobial therapy if patient is symptomatic'
            }
        }
        
        guide = guides[risk_level]
        
        lines = [
            "",
            "=" * 70,
            "🩺 CLINICAL INTERPRETATION GUIDE",
            "=" * 70,
            "",
            f"{guide['icon']} {risk_level.upper()} RISK ({risk['final_score']:.0f}/100)",
            "",
            "WHAT THIS MEANS:",
            f"  {guide['meaning']}",
            "",
            "RECOMMENDED ACTION:",
            f"  {guide['action']}",
            "",
            "CLINICAL GUIDANCE:",
            f"  {guide['clinical']}",
            "",
            "ADDITIONAL CONSIDERATIONS:",
        ]
        
        # Add context-specific considerations
        if risk['taxonomic_score'] > 50:
            lines.append("  • Known pathogenic species detected - correlate with symptoms")
        
        if risk['functional_score'] > 50:
            lines.append("  • Significant mobile genetic elements - potential for gene transfer")
        
        if risk['ml_score'] > 40:
            lines.append("  • ML model flags concerning sequences - consider further validation")
        
        lines.extend([
            "",
            "NOTE: This is a computational prediction. Clinical judgment and",
            "      additional testing are essential for patient management.",
            "=" * 70,
        ])
        
        return '\n'.join(lines)
    
    def _format_pathogen_hits_table(self, pathogen_file: Path, details: Dict) -> str:
        """Format pathogen DB hits as table."""
        # Combine all pathogen categories
        all_hits = (
            details['virulence_factors'] + 
            details['amr_genes'] + 
            details['secretion_systems'] + 
            details['other_pathogen_markers']
        )
        
        if not all_hits:
            return "No pathogen database hits"
        
        lines = [
            f"{'Gene ID':<18} {'Category':>12} {'Identity':>10} {'Organism':<25}",
            "-" * 70,
        ]
        
        # Categorize each hit
        for hit in all_hits[:15]:
            gene_id = hit['gene_id'][:16]
            desc = hit['description'].lower()
            organism = hit['organism'][:23]
            identity = f"{hit['identity']:.1f}%"
            
            # Determine category
            if any(kw in desc for kw in ['toxin', 'virulence', 'hemolysin']):
                category = '🦠 Virulence'
            elif any(kw in desc for kw in ['resistance', 'beta-lactam']):
                category = '💊 AMR'
            elif 'secretion' in desc:
                category = '🔧 Secretion'
            else:
                category = '⚠️  Other'
            
            lines.append(
                f"{gene_id:<18} {category:>12} {identity:>10} {organism:<25}"
            )
        
        if len(all_hits) > 15:
            lines.append(f"\n... and {len(all_hits) - 15} more pathogen-associated genes")
        
        return '\n'.join(lines)
    
    def _cross_reference_ml_functional(self, ml_file: Path, risk: Dict) -> str:
        """Cross-reference ML predictions with functional annotations."""
        with open(ml_file) as f:
            ml_data = json.load(f)
            predictions = ml_data.get('predictions', [])
        
        # Focus on high-confidence pathogenic predictions
        high_conf = [p for p in predictions 
                    if p['prediction'] == 'Pathogenic' and p['high_confidence']]
        
        if not high_conf:
            return "No high-confidence pathogenic predictions to cross-reference"
        
        lines = [
            "Sequences with both ML pathogenic prediction AND pathogen DB hit:",
            ""
        ]
        
        # This would need actual cross-referencing with functional annotations
        # For now, show high-confidence predictions
        for pred in high_conf[:5]:
            lines.extend([
                f"  {pred['sequence_id']}:",
                f"    ML Confidence: {pred['confidence']*100:.1f}%",
                f"    Sequence Length: {pred['sequence_length']} aa",
                f"    → Recommend functional validation",
                "",
            ])
        
        return '\n'.join(lines)
    
    def _generate_confidence_viz(self, risk: Dict) -> str:
        """Generate ASCII visualization of confidence distribution."""
        # Simple histogram
        details = risk['details']
        
        if not details:
            return "No data to visualize"
        
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
        
        # Create bars
        max_count = max(ranges.values()) if ranges.values() else 1
        lines = []
        
        for range_label, count in ranges.items():
            bar_length = int((count / max_count) * 40) if max_count > 0 else 0
            bar = '█' * bar_length
            lines.append(f"  {range_label:>10} │{bar} {count}")
        
        return '\n'.join(lines)
    
    def _wrap_text(self, text: str, width: int = 70, indent: int = 0) -> str:
        """Wrap text to specified width with indent."""
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
    
    def generate_visualization_data(self,
                                   integrated_risk: Dict,
                                   taxonomic_risk: Dict,
                                   functional_risk: Dict,
                                   ml_risk: Dict) -> Dict:
        """
        Generate data structure for visualization module.
        
        Returns:
            Dict with data ready for pathogenic_visualizer.py
        """
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
                'tier_scores': {
                    'taxonomic': integrated_risk['taxonomic_score'],
                    'functional': integrated_risk['functional_score'],
                    'ml': integrated_risk['ml_score']
                },
                'multipliers_applied': integrated_risk['multipliers_applied'],
                'interpretation': integrated_risk['interpretation']
            },
            'timestamp': pd.Timestamp.now().isoformat()
        }
        
        with open(output_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        self.logger.info(f"Risk summary exported to {output_file}")