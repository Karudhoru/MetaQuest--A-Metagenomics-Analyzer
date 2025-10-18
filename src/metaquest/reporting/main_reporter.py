"""
MetaQuest Main Reporter - Enhanced Orchestrator
Coordinates all reporting modules with integrated risk scoring
"""

import json
from pathlib import Path
from typing import Dict, Optional
import logging

from .taxonomy_reporter import TaxonomyReporter
from .functional_reporter import FunctionalReporter
from .pathogen_risk_reporter import PathogenRiskReporter
from .risk_scoring import RiskScorer, calculate_all_risks
from .base_reporter import BaseReporter


class MainReporter(BaseReporter):
    """
    Main orchestrator for MetaQuest reporting system.
    Generates comprehensive reports with integrated risk assessment.
    """
    
    def __init__(self, output_dir: Path, view_mode: str = 'both'):
        """
        Initialize main reporter.
        
        Args:
            output_dir: Output directory for reports
            view_mode: 'clinician', 'researcher', or 'both'
        """
        super().__init__(output_dir)
        self.view_mode = view_mode
        
        # Initialize sub-reporters
        self.taxonomy_reporter = TaxonomyReporter(output_dir)
        self.functional_reporter = FunctionalReporter(output_dir)
        self.pathogen_reporter = PathogenRiskReporter(output_dir)
        
        # Risk scorer
        self.risk_scorer = None
    
    def generate_report(self,
                       bracken_file: Path,
                       sample_info_file: Path,
                       functional_annotations_file: Path,
                       pathogen_hits_file: Path,
                       ml_predictions_file: Path,
                       pathogen_organism_file: Optional[Path] = None) -> str:
        """
        Generate complete integrated report.
        
        Args:
            bracken_file: Bracken taxonomy report
            sample_info_file: Sample assembly statistics
            functional_annotations_file: Functional annotations (SwissProt)
            pathogen_hits_file: Pathogen database hits
            ml_predictions_file: ML predictions JSON
            pathogen_organism_file: Pathogen organism risk levels
        
        Returns:
            Complete formatted report
        """
        self.logger.info("Generating comprehensive MetaQuest report...")
        
        # Initialize risk scorer
        self.risk_scorer = RiskScorer(pathogen_organism_file)
        
        # Calculate all risk scores
        self.logger.info("Calculating risk scores...")
        risk_data = calculate_all_risks(
            bracken_file=bracken_file,
            functional_file=functional_annotations_file,
            pathogen_hits_file=pathogen_hits_file,
            ml_predictions_file=ml_predictions_file,
            pathogen_organism_file=pathogen_organism_file
        )
        
        # Build complete report
        report_sections = []
        
        # 1. Executive Summary
        report_sections.append(self._generate_executive_summary(
            bracken_file, sample_info_file, risk_data
        ))
        
        # 2. Taxonomy Report
        report_sections.append(self.generate_taxonomy_report(
            bracken_file, risk_data
        ))
        
        # 3. Functional Report
        report_sections.append(self.generate_functional_report(
            sample_info_file, functional_annotations_file, risk_data
        ))
        
        # 4. Pathogen Risk Report
        report_sections.append(self.generate_pathogen_report(
            risk_data, pathogen_hits_file, ml_predictions_file
        ))
        
        # 5. Summary and Recommendations
        report_sections.append(self._generate_final_summary(risk_data))
        
        # Combine all sections
        complete_report = '\n\n'.join(report_sections)
        
        # Save report
        report_file = self.output_dir / 'comprehensive_report.txt'
        with open(report_file, 'w') as f:
            f.write(complete_report)
        
        self.logger.info(f"Complete report saved to {report_file}")
        
        # Export data for visualizations
        self._export_visualization_data(risk_data, bracken_file, 
                                       sample_info_file, functional_annotations_file)
        
        return complete_report
    
    def generate_taxonomy_report(self,
                                bracken_file: Path,
                                risk_data: Dict) -> str:
        """
        Generate standalone taxonomy report with risk assessment.
        
        Args:
            bracken_file: Path to Bracken taxonomy TSV
            risk_data: Complete risk assessment data
        
        Returns:
            Formatted taxonomy report
        """
        self.logger.info("Generating taxonomy report...")
        
        return self.taxonomy_reporter.generate_report(
            bracken_file=bracken_file,
            taxonomic_risk=risk_data['taxonomic'],
            mode=self.view_mode
        )
    
    def generate_functional_report(self,
                                   sample_info_file: Path,
                                   functional_annotations_file: Path,
                                   risk_data: Dict) -> str:
        """
        Generate standalone functional annotation report.
        
        Args:
            sample_info_file: Path to sample.txt with assembly stats
            functional_annotations_file: Path to functional_annotations.tsv
            risk_data: Complete risk assessment data
        
        Returns:
            Formatted functional report
        """
        self.logger.info("Generating functional report...")
        
        return self.functional_reporter.generate_report(
            sample_info_file=sample_info_file,
            annotation_file=functional_annotations_file,
            functional_risk=risk_data['functional']
        )
    
    def generate_pathogen_report(self,
                                risk_data: Dict,
                                pathogen_hits_file: Path,
                                ml_predictions_file: Path) -> str:
        """
        Generate standalone pathogen risk assessment report.
        
        Args:
            risk_data: Complete risk assessment data
            pathogen_hits_file: Path to pathogen database hits
            ml_predictions_file: Path to ML predictions JSON
        
        Returns:
            Formatted pathogen risk report
        """
        self.logger.info("Generating pathogen risk report...")
        
        return self.pathogen_reporter.generate_report(
            integrated_risk=risk_data['integrated'],
            taxonomic_risk=risk_data['taxonomic'],
            functional_risk=risk_data['functional'],
            ml_risk=risk_data['ml'],
            pathogen_hits_file=pathogen_hits_file,
            ml_predictions_file=ml_predictions_file
        )
    
    def _generate_executive_summary(self,
                               bracken_file: Path,
                               sample_file: Path,
                               risk_data: Dict) -> str:
        """Generate executive summary card."""
        import pandas as pd
        
        # Load basic data
        bracken_df = pd.read_csv(bracken_file, sep='\t')
        
        # *** FIX: Properly handle Pandas Series check ***
        # Old buggy code: dominant = bracken_df.iloc[0] if not bracken_df.empty else None
        # Old buggy check: if dominant:
        
        # New fixed code with explicit None check
        if not bracken_df.empty:
            dominant = bracken_df.iloc[0]
            has_dominant = True
        else:
            dominant = None
            has_dominant = False
        
        # Parse sample info
        sample_stats = self.functional_reporter._parse_sample_info(sample_file)
        
        # Get risk info
        integrated = risk_data['integrated']
        taxonomic = risk_data['taxonomic']
        
        risk_emoji = {
            'High': '🔴',
            'Moderate': '🟡',
            'Low': '🟢'
        }[integrated['risk_level']]
        
        lines = [
            "╔" + "═" * 68 + "╗",
            "║" + " " * 20 + "METAQUEST ANALYSIS REPORT" + " " * 23 + "║",
            "╠" + "═" * 68 + "╣",
            f"║  Sample: {sample_stats['organism']:<55} ║",
            f"║  Status: {risk_emoji} {integrated['risk_level'].upper()} RISK    Quality: {'⭐' * 4:<20} ║",
            "╠" + "═" * 68 + "╣",
            f"║  🧬 Assembly: {sample_stats['gene']} genes | {sample_stats['bases']:,} bp | {sample_stats['contigs']} contigs" + " " * (15 - len(str(sample_stats['contigs']))) + "║",
        ]
        
        # *** FIX: Use explicit boolean check instead of Series truthiness ***
        if has_dominant and dominant is not None:
            dominant_name = str(dominant['name'])[:40] if len(str(dominant['name'])) > 40 else str(dominant['name'])
            dominant_fraction = float(dominant['fraction_total_reads']) * 100
            padding_length = max(0, 4 - len(f"{dominant_fraction:.1f}"))
            
            lines.append(
                f"║  🦠 Dominant: {dominant_name:<40} ({dominant_fraction:.1f}%)" + " " * padding_length + "║"
            )
        
        pathogen_count = len(taxonomic['pathogens_detected'])
        if pathogen_count > 0:
            padding_length = max(0, 36 - len(str(pathogen_count)))
            lines.append(
                f"║  ⚠️  Pathogens: {pathogen_count} species detected" + " " * padding_length + "║"
            )
        
        ml_flagged = risk_data['ml']['high_confidence_pathogenic']
        ml_percentage = (ml_flagged / 40 * 100) if 40 > 0 else 0
        padding_length = max(0, 20 - len(str(ml_flagged)))
        lines.append(
            f"║  🧪 ML Risk: {ml_flagged}/40 sequences flagged ({ml_percentage:.0f}%)" + " " * padding_length + "║"
        )
        
        lines.append("╚" + "═" * 68 + "╝")
        
        return '\n'.join(lines)
    
    def _generate_final_summary(self, risk_data: Dict) -> str:
        """Generate final summary with actionable recommendations."""
        integrated = risk_data['integrated']
        
        lines = [
            "",
            "=" * 70,
            "📋 FINAL SUMMARY AND RECOMMENDATIONS",
            "=" * 70,
            "",
            f"OVERALL PATHOGENICITY RISK: {integrated['risk_level'].upper()} ({integrated['final_score']:.0f}/100)",
            "",
            "KEY FINDINGS:",
        ]
        
        # Taxonomic findings
        taxonomic = risk_data['taxonomic']
        if taxonomic['pathogens_detected']:
            lines.append(
                f"  • {len(taxonomic['pathogens_detected'])} pathogenic species detected"
            )
            if taxonomic['dominant_pathogen']:
                dom = taxonomic['dominant_pathogen']
                lines.append(
                    f"    → Dominant: {dom['species']} ({dom['abundance']*100:.2f}%)"
                )
        else:
            lines.append("  • No significant pathogenic species detected")
        
        # Functional findings
        functional = risk_data['functional']
        if functional['transposase_count'] > 5:
            lines.append(
                f"  • HIGH mobile genetic element activity ({functional['transposase_count']} transposases)"
            )
        
        if functional['virulence_count'] > 0:
            lines.append(
                f"  • {functional['virulence_count']} virulence-associated genes identified"
            )
        
        # ML findings
        ml = risk_data['ml']
        if ml['high_confidence_pathogenic'] > 0:
            lines.append(
                f"  • ML analysis flagged {ml['high_confidence_pathogenic']} sequences as pathogenic"
            )
        
        # Recommendations
        lines.extend([
            "",
            "ACTIONABLE RECOMMENDATIONS:",
        ])
        
        if integrated['risk_level'] == 'High':
            lines.extend([
                "  1. ⚠️  URGENT: Clinical correlation required",
                "  2. Consider targeted antimicrobial therapy if symptomatic",
                "  3. Implement infection control measures if in clinical setting",
                "  4. Further characterization of flagged genes recommended",
            ])
        elif integrated['risk_level'] == 'Moderate':
            lines.extend([
                "  1. Monitor patient/sample for clinical symptoms",
                "  2. Consider additional diagnostic testing if indicated",
                "  3. Track for changes in pathogen abundance over time",
                "  4. Correlate with patient history and presentation",
            ])
        else:
            lines.extend([
                "  1. Routine monitoring is sufficient",
                "  2. No immediate interventions indicated",
                "  3. Document baseline for future comparisons",
            ])
        
        # Data quality note
        lines.extend([
            "",
            "DATA QUALITY NOTES:",
            "  • This report integrates taxonomic, functional, and ML predictions",
            "  • Risk scores are computational predictions - not diagnostic results",
            "  • Clinical judgment is essential for patient management decisions",
            "  • Consider validation of high-risk findings with orthogonal methods",
            "",
            "=" * 70,
            "",
            f"Report generated by MetaQuest v4.0.0",
            f"For questions or support: metaquest-support@example.org",
            "",
        ])
        
        return '\n'.join(lines)
    
    def _export_visualization_data(self,
                                  risk_data: Dict,
                                  bracken_file: Path,
                                  sample_file: Path,
                                  annotation_file: Path):
        """Export all data needed for visualizations."""
        viz_dir = self.output_dir / 'visualization_data'
        viz_dir.mkdir(exist_ok=True)
        
        # Taxonomy visualization data
        taxonomy_viz = self.taxonomy_reporter.generate_visualization_data(
            bracken_file, risk_data['taxonomic']
        )
        with open(viz_dir / 'taxonomy_viz_data.json', 'w') as f:
            json.dump(taxonomy_viz, f, indent=2)
        
        # Functional visualization data
        functional_viz = self.functional_reporter.generate_visualization_data(
            sample_file, annotation_file
        )
        with open(viz_dir / 'functional_viz_data.json', 'w') as f:
            json.dump(functional_viz, f, indent=2)
        
        # Pathogen visualization data
        pathogen_viz = self.pathogen_reporter.generate_visualization_data(
            risk_data['integrated'],
            risk_data['taxonomic'],
            risk_data['functional'],
            risk_data['ml']
        )
        with open(viz_dir / 'pathogen_viz_data.json', 'w') as f:
            json.dump(pathogen_viz, f, indent=2)
        
        # Integrated risk data
        with open(viz_dir / 'integrated_risk_data.json', 'w') as f:
            json.dump(risk_data['integrated'], f, indent=2, default=str)
        
        self.logger.info(f"Visualization data exported to {viz_dir}")
    
    def export_tables(self,
                     bracken_file: Path,
                     annotation_file: Path,
                     risk_data: Dict):
        """Export all data tables to CSV format."""
        tables_dir = self.output_dir / 'tables'
        tables_dir.mkdir(exist_ok=True)
        
        # Taxonomy table
        self.taxonomy_reporter.export_table(
            bracken_file,
            tables_dir / 'taxonomy_complete.csv',
            risk_data['taxonomic']
        )
        
        # Functional annotations table
        self.functional_reporter.export_annotations_table(
            annotation_file,
            tables_dir / 'functional_annotations.csv'
        )
        
        # Risk summary
        self.pathogen_reporter.export_risk_summary(
            risk_data['integrated'],
            tables_dir / 'risk_summary.json'
        )
        
        self.logger.info(f"Data tables exported to {tables_dir}")


def generate_report(output_dir: Path,
                   bracken_file: Path,
                   sample_info_file: Path,
                   functional_annotations_file: Path,
                   pathogen_hits_file: Path,
                   ml_predictions_file: Path,
                   pathogen_organism_file: Optional[Path] = None,
                   view_mode: str = 'both') -> str:
    """
    Convenience function to generate MetaQuest report.
    
    Args:
        output_dir: Output directory
        bracken_file: Bracken taxonomy TSV
        sample_info_file: Sample statistics TXT
        functional_annotations_file: Functional annotations TSV
        pathogen_hits_file: Pathogen database hits TSV
        ml_predictions_file: ML predictions JSON
        pathogen_organism_file: Pathogen organism risk levels (optional)
        view_mode: 'clinician', 'researcher', or 'both'
    
    Returns:
        Complete report string
    """
    reporter = MainReporter(output_dir, view_mode=view_mode)
    
    return reporter.generate_report(
        bracken_file=bracken_file,
        sample_info_file=sample_info_file,
        functional_annotations_file=functional_annotations_file,
        pathogen_hits_file=pathogen_hits_file,
        ml_predictions_file=ml_predictions_file,
        pathogen_organism_file=pathogen_organism_file
    )


def generate_taxonomy_report(output_dir: Path,
                             bracken_file: Path,
                             pathogen_organism_file: Optional[Path] = None,
                             view_mode: str = 'both') -> str:
    """
    Generate standalone taxonomy report.
    
    Args:
        output_dir: Output directory
        bracken_file: Bracken taxonomy TSV
        pathogen_organism_file: Pathogen organism risk levels (optional)
        view_mode: 'clinician', 'researcher', or 'both'
    
    Returns:
        Taxonomy report string
    """
    reporter = MainReporter(output_dir, view_mode=view_mode)
    risk_scorer = RiskScorer(pathogen_organism_file)
    
    import pandas as pd
    bracken_df = pd.read_csv(bracken_file, sep='\t')
    taxonomic_risk = risk_scorer.calculate_taxonomic_risk(bracken_df)
    
    risk_data = {
        'taxonomic': taxonomic_risk,
        'functional': {'score': 0, 'details': {}},
        'ml': {'score': 0, 'details': []},
        'integrated': {'final_score': taxonomic_risk['score'], 'risk_level': 'Low'}
    }
    
    return reporter.generate_taxonomy_report(bracken_file, risk_data)


def generate_functional_report(output_dir: Path,
                               sample_info_file: Path,
                               functional_annotations_file: Path,
                               pathogen_hits_file: Optional[Path] = None) -> str:
    """
    Generate standalone functional annotation report.
    
    Args:
        output_dir: Output directory
        sample_info_file: Sample statistics TXT
        functional_annotations_file: Functional annotations TSV
        pathogen_hits_file: Optional pathogen database hits TSV
    
    Returns:
        Functional report string
    """
    reporter = MainReporter(output_dir)
    
    # Calculate functional risk if pathogen hits available
    if pathogen_hits_file and pathogen_hits_file.exists():
        import pandas as pd
        functional_df = pd.read_csv(functional_annotations_file, sep='\t', header=None)
        pathogen_df = pd.read_csv(pathogen_hits_file, sep='\t', header=None)
        
        ann_cols = ['query_id', 'subject_id', 'identity', 'length', 'mismatches', 
                    'gaps', 'q_start', 'q_end', 's_start', 's_end', 'evalue', 
                    'bitscore', 'description']
        functional_df.columns = ann_cols[:len(functional_df.columns)]
        pathogen_df.columns = ann_cols[:len(pathogen_df.columns)]
        
        risk_scorer = RiskScorer()
        functional_risk = risk_scorer.calculate_functional_risk(functional_df, pathogen_df)
    else:
        functional_risk = {'score': 0, 'details': {}}
    
    risk_data = {
        'taxonomic': {'score': 0, 'pathogens_detected': []},
        'functional': functional_risk,
        'ml': {'score': 0, 'details': []},
        'integrated': {'final_score': 0, 'risk_level': 'Low'}
    }
    
    return reporter.generate_functional_report(
        sample_info_file, functional_annotations_file, risk_data
    )


def generate_pathogen_report(output_dir: Path,
                             bracken_file: Path,
                             functional_annotations_file: Path,
                             pathogen_hits_file: Path,
                             ml_predictions_file: Path,
                             pathogen_organism_file: Optional[Path] = None) -> str:
    """
    Generate standalone pathogen risk assessment report.
    
    Args:
        output_dir: Output directory
        bracken_file: Bracken taxonomy TSV
        functional_annotations_file: Functional annotations TSV
        pathogen_hits_file: Pathogen database hits TSV
        ml_predictions_file: ML predictions JSON
        pathogen_organism_file: Pathogen organism risk levels (optional)
    
    Returns:
        Pathogen risk report string
    """
    reporter = MainReporter(output_dir)
    
    # Calculate all risks
    risk_data = calculate_all_risks(
        bracken_file=bracken_file,
        functional_file=functional_annotations_file,
        pathogen_hits_file=pathogen_hits_file,
        ml_predictions_file=ml_predictions_file,
        pathogen_organism_file=pathogen_organism_file
    )
    
    return reporter.generate_pathogen_report(
        risk_data, pathogen_hits_file, ml_predictions_file
    )   