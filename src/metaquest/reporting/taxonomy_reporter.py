"""
MetaQuest Taxonomy Reporter
Dual-mode reporting for clinicians and researchers with pathogen alerts
"""
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from .base_reporter import BaseReporter


class TaxonomyReporter(BaseReporter):
    """
    Generate taxonomy reports with pathogen risk assessment.
    Supports both clinician (simplified) and researcher (detailed) views.
    """
    

    def __init__(self, output_dir: Path):
        super().__init__(output_dir)
        self.section_name = "Taxonomy Report"
        
        # Risk level configurations - CASE SENSITIVE!
        self.risk_colors = {
            'High': '🔴',      # Must match exactly
            'Moderate': '🟡',  # Must match exactly
            'Low': '🟢',       # Must match exactly
            'Unknown': '⚪'
        }

        self.commensal_species = {
        'Streptococcus thermophilus',
        'Streptococcus salivarius',
        'Lactobacillus',
        'Bifidobacterium',
        'Lactococcus lactis'
    }
        
        # Pathogen organisms (loaded from risk scorer)
        self.pathogen_organisms = {}

    def _is_true_pathogen(self, species_name: str, abundance: float) -> bool:
        """
        Validate if species is a true pathogen, not a commensal.
        
        Args:
            species_name: Species name
            abundance: Relative abundance
        
        Returns:
            True if pathogen, False if commensal
        """
        species_lower = species_name.lower()
        
        # Check commensal list
        for commensal in self.commensal_species:
            if commensal.lower() in species_lower:
                # Exception: High abundance of normally commensal species
                if abundance > 0.10:  # >10% abundance is concerning
                    self.logger.warning(
                        f"{species_name} is normally commensal but at {abundance*100:.1f}% abundance"
                    )
                    return True  # Treat as opportunistic pathogen
                return False
        
        return True  # Not in commensal list, treat as potential pathogen
    
    def generate_report(self,
                       bracken_file: Path,
                       taxonomic_risk: Dict,
                       mode: str = 'both') -> str:
        """
        Generate complete taxonomy report.
        
        Args:
            bracken_file: Path to Bracken report TSV
            taxonomic_risk: Risk assessment from RiskScorer
            mode: 'clinician', 'researcher', or 'both'
        
        Returns:
            Formatted report string
        """
        self.print_section_header("TAXONOMY REPORT")
        
        # Load data
        bracken_df = pd.read_csv(bracken_file, sep='\t')
        
        report_parts = []
        
        # 1. Quick Overview (always shown)
        report_parts.append(self._generate_quick_overview(bracken_df, taxonomic_risk))
        
        # 2. Pathogen Alert Box (if pathogens detected)
        if taxonomic_risk['pathogens_detected']:
            report_parts.append(self._generate_pathogen_alerts(taxonomic_risk))
        
        # 3. Mode-specific content
        if mode in ['clinician', 'both']:
            report_parts.append(self._generate_clinician_view(bracken_df, taxonomic_risk))
        
        if mode in ['researcher', 'both']:
            report_parts.append(self._generate_researcher_view(bracken_df, taxonomic_risk))
        
        # 4. Diversity metrics
        report_parts.append(self._generate_diversity_metrics(bracken_df))
        
        return '\n\n'.join(report_parts)
    
    def _generate_quick_overview(self, bracken_df: pd.DataFrame, risk: Dict) -> str:
        """Generate executive summary box."""
        total_species = len(bracken_df)
        dominant = bracken_df.iloc[0] if not bracken_df.empty else None
        
        # BUG FIX: Add debugging for pathogen detection
        pathogen_count = len(risk['pathogens_detected'])
        
        # DEBUG: Log what we received
        self.logger.debug(f"Risk dict keys: {risk.keys()}")
        self.logger.debug(f"Pathogens detected: {pathogen_count}")
        self.logger.debug(f"Pathogen list: {risk['pathogens_detected']}")

        for pathogen in risk['pathogens_detected'][:3]:  # Check first 3
            if 'risk_level' not in pathogen:
                self.logger.warning(f"⚠️ Pathogen missing risk_level: {pathogen.get('species', 'Unknown')}")
            else:
                self.logger.debug(f"✓ Pathogen {pathogen['species']} has risk_level: {pathogen['risk_level']}")
        
        # BUG FIX: Check if pathogens_detected is actually populated
        if pathogen_count == 0:
            self.logger.warning("⚠️ pathogens_detected is empty - check pathogen_analysis.py")
            # Try to manually detect pathogens from bracken data as fallback
            pathogen_species = [
                'Clostridioides difficile', 'Clostridium difficile',
                'Staphylococcus aureus', 'Escherichia coli',
                'Klebsiella pneumoniae', 'Pseudomonas aeruginosa',
                'Streptococcus thermophilus', 'Streptococcus salivarius',
                'Campylobacter gracilis', 'Salmonella enterica'
            ]
            
            # Manual pathogen detection
            detected_pathogens = bracken_df[
                bracken_df['name'].isin(pathogen_species)
            ]
            
            if not detected_pathogens.empty:
                self.logger.warning(f"⚠️ FOUND {len(detected_pathogens)} PATHOGENS IN BRACKEN BUT NOT IN RISK DICT!")
                self.logger.warning(f"Pathogens found: {detected_pathogens['name'].tolist()}")
                # Override the count for display
                pathogen_count = len(detected_pathogens)
        
        lines = [
            "=" * 70,
            "📊 MICROBIAL COMPOSITION - QUICK OVERVIEW",
            "=" * 70,
        ]
        
        if dominant is not None:
            lines.extend([
                f"Dominant Species:       {dominant['name']}",
                f"  └─ Abundance:         {dominant['fraction_total_reads']*100:.2f}%",
                f"  └─ Reads:             {int(dominant['new_est_reads']):,}",
            ])
        
        lines.extend([
            f"Total Species Detected: {total_species}",
            f"Pathogenic Species:     {pathogen_count} {self.risk_colors.get('High' if pathogen_count > 0 else 'Low', '')}",
        ])
        
        if risk.get('total_pathogen_abundance', 0) > 0:
            lines.append(
                f"Total Pathogen Load:    {risk['total_pathogen_abundance']*100:.2f}%"
            )
        
        lines.append("=" * 70)
        
        return '\n'.join(lines)
    
    def _generate_pathogen_alerts(self, risk: Dict) -> str:
        """Generate pathogen alert box with clinical significance."""
        lines = [
            "",
            "⚠️  " + "=" * 66,
            "⚠️  CLINICALLY SIGNIFICANT ORGANISMS DETECTED",
            "⚠️  " + "=" * 66,
        ]
        
        true_pathogens = [
        p for p in risk['pathogens_detected']
        if self._is_true_pathogen(p['species'], p['abundance'])
        ]
        
        if not true_pathogens:
            lines.extend([
                "",
                "✅ No pathogenic species at clinically significant levels",
                "   Sample dominated by commensal organisms",
            ])
            return '\n'.join(lines)

        # Show top 5 pathogens by abundance
        top_pathogens = risk['pathogens_detected'][:5]
        
        for pathogen in top_pathogens:
            # FIX: Extract risk_level from pathogen dict
            risk_level = pathogen.get('risk_level', 'Low')  # Default to 'Low' if missing
            risk_icon = self.risk_colors.get(risk_level, '⚪')
            
            lines.extend([
                f"",
                f"{risk_icon} {pathogen['species']}",
                f"   Abundance:      {pathogen['abundance']*100:.2f}%",
                f"   Reads:          {int(pathogen['reads']):,}",
                f"   Risk Level:     {risk_level.upper()}",
                f"   Clinical Notes: {self._get_clinical_notes(pathogen['species'])}",
            ])
        
        if len(risk['pathogens_detected']) > 5:
            lines.append(f"\n... and {len(risk['pathogens_detected']) - 5} more pathogenic species")
        
        lines.append("=" * 70)
        
        return '\n'.join(lines)
    
    def _generate_clinician_view(self, bracken_df: pd.DataFrame, risk: Dict) -> str:
        """Generate simplified view for clinicians."""
        lines = [
            "",
            "👨‍⚕️ CLINICIAN VIEW - SIMPLIFIED SUMMARY",
            "=" * 70,
        ]
        
        # Categorize organisms
        pathogens = risk['pathogens_detected']
        pathogen_names = {p['species'] for p in pathogens}
        commensals = bracken_df[~bracken_df['name'].isin(pathogen_names)]
        
        # Show breakdown
        lines.extend([
            f"",
            f"Known Pathogens:        {len(pathogens)} species",
            f"Commensal Organisms:    {len(commensals)} species",
            f"",
            "KEY FINDINGS:",
        ])
        
        # Key findings
        if pathogens:
            lines.append(f"  • {len(pathogens)} pathogenic species require attention")
            
            high_risk = [p for p in pathogens if p['risk_level'] == 'High']
            if high_risk:
                lines.append(f"  • {len(high_risk)} HIGH-RISK pathogens detected:")
                for p in high_risk[:3]:
                    lines.append(f"      - {p['species']} ({p['abundance']*100:.2f}%)")
        else:
            lines.append("  • No known pathogenic species at significant levels")
        
        # Dominant commensal
        if not commensals.empty:
            dominant_commensal = commensals.iloc[0]
            lines.append(
                f"  • Dominant commensal: {dominant_commensal['name']} "
                f"({dominant_commensal['fraction_total_reads']*100:.2f}%)"
            )
        
        lines.extend([
            "",
            "CLINICAL RECOMMENDATIONS:",
            self._generate_clinical_recommendations(risk),
            "=" * 70,
        ])
        
        return '\n'.join(lines)
    
    def _generate_researcher_view(self, bracken_df: pd.DataFrame, risk: Dict) -> str:
        """Generate detailed view for researchers."""
        lines = [
            "",
            "🔬 RESEARCHER VIEW - DETAILED STATISTICS",
            "=" * 70,
        ]
        
        # Top 15 species table
        lines.extend([
            "",
            "TOP 15 SPECIES BY ABUNDANCE:",
            "",
            self._format_species_table(bracken_df.head(15), risk),
        ])
        
        # Taxonomic breakdown by genus
        lines.extend([
            "",
            "GENUS-LEVEL COMPOSITION:",
            "",
            self._generate_genus_breakdown(bracken_df),
        ])
        
        # Rare species summary
        rare_threshold = 0.01  # 1%
        rare_species = bracken_df[bracken_df['fraction_total_reads'] < rare_threshold]
        if not rare_species.empty:
            lines.extend([
                "",
                f"RARE SPECIES (<{rare_threshold*100}% abundance): {len(rare_species)} detected",
            ])
        
        lines.append("=" * 70)
        
        return '\n'.join(lines)
    
    def _generate_diversity_metrics(self, bracken_df: pd.DataFrame) -> str:
        """Calculate and display diversity metrics."""
        lines = [
            "",
            "📈 DIVERSITY METRICS",
            "=" * 70,
        ]
        
        # Calculate Shannon diversity
        abundances = bracken_df['fraction_total_reads'].values
        shannon = -sum(abundances * np.log(abundances + 1e-10))
        
        # Calculate Simpson diversity
        simpson = 1 - sum(abundances ** 2)
        
        # Species richness
        richness = len(bracken_df)
        
        lines.extend([
            f"Shannon Index:    {shannon:.3f}  (Higher = more diverse)",
            f"Simpson Index:    {simpson:.3f}  (Higher = more diverse)",
            f"Species Richness: {richness}",
            "",
            "Interpretation:",
        ])
        
        # Interpret diversity
        if shannon > 3.0:
            lines.append("  • HIGH diversity - complex microbial community")
        elif shannon > 2.0:
            lines.append("  • MODERATE diversity - balanced community")
        else:
            lines.append("  • LOW diversity - dominated by few species")
        
        lines.append("=" * 70)
        
        return '\n'.join(lines)
    
    def _format_species_table(self, df: pd.DataFrame, risk: Dict) -> str:
        """Format species data as table."""
        pathogen_dict = {p['species']: p for p in risk['pathogens_detected']}
        
        lines = [
            f"{'Species':<40} {'Abundance':>10} {'Reads':>12} {'Risk':>8}"
        ]
        lines.append("-" * 75)
        
        for _, row in df.iterrows():
            species = row['name']
            abundance = f"{row['fraction_total_reads']*100:.2f}%"
            reads = f"{int(row['new_est_reads']):,}"
            
            # Get risk level
            if species in pathogen_dict:
                # Extract risk_level from the pathogen dict
                pathogen_data = pathogen_dict[species]
                risk_level = pathogen_data.get('risk_level', 'Low')  # Use .get() with default
                risk_icon = self.risk_colors.get(risk_level, '⚪')
                risk_text = f"{risk_icon} {risk_level}"
            else:
                # Not a pathogen - mark as Low risk
                risk_text = f"{self.risk_colors['Low']} Low"
            
            # Truncate species name if too long
            species_display = species[:38] + '..' if len(species) > 40 else species
            
            lines.append(
                f"{species_display:<40} {abundance:>10} {reads:>12} {risk_text:>8}"
            )
        
        return '\n'.join(lines)
    
    def _generate_genus_breakdown(self, bracken_df: pd.DataFrame) -> str:
        """Group species by genus and show breakdown."""
        # Extract genus from species name
        bracken_df['genus'] = bracken_df['name'].str.split().str[0]
        
        # Group by genus
        genus_summary = bracken_df.groupby('genus').agg({
            'fraction_total_reads': 'sum',
            'new_est_reads': 'sum',
            'name': 'count'
        }).sort_values('fraction_total_reads', ascending=False).head(10)
        
        lines = [
            f"{'Genus':<25} {'Abundance':>12} {'Species':>10} {'Reads':>12}"
        ]
        lines.append("-" * 65)
        
        for genus, row in genus_summary.iterrows():
            abundance = f"{row['fraction_total_reads']*100:.2f}%"
            species_count = int(row['name'])
            reads = f"{int(row['new_est_reads']):,}"
            
            lines.append(
                f"{genus:<25} {abundance:>12} {species_count:>10} {reads:>12}"
            )
        
        return '\n'.join(lines)
    
    def _get_clinical_notes(self, species: str) -> str:
        """Return clinical significance notes for known pathogens."""
        clinical_notes = {
            'Clostridioides difficile': 'Antibiotic-associated diarrhea, pseudomembranous colitis',
            'Clostridium difficile': 'Antibiotic-associated diarrhea, pseudomembranous colitis',
            'Staphylococcus aureus': 'Skin/soft tissue infections, bacteremia, endocarditis',
            'Escherichia coli': 'UTIs, gastroenteritis, sepsis (strain-dependent)',
            'Klebsiella pneumoniae': 'Pneumonia, UTIs, bloodstream infections',
            'Pseudomonas aeruginosa': 'Opportunistic infections, nosocomial pneumonia',
            'Acinetobacter baumannii': 'Multidrug-resistant nosocomial infections',
            'Enterococcus faecium': 'Bloodstream and urinary tract infections',
            'Salmonella enterica': 'Gastroenteritis, typhoid fever',
            'Bacteroides fragilis': 'Anaerobic infections, intra-abdominal abscesses',
            'Streptococcus pneumoniae': 'Pneumonia, meningitis, otitis media',
            'Mycobacterium tuberculosis': 'Tuberculosis - requires immediate attention',
        }
        
        return clinical_notes.get(species, 'Opportunistic pathogen - monitor clinically')
    
    def _generate_clinical_recommendations(self, risk: Dict) -> str:
        """Generate clinical action recommendations."""
        pathogens = risk['pathogens_detected']
        
        if not pathogens:
            return "  • No specific interventions indicated\n  • Continue routine monitoring"
        
        recommendations = []
        
        # High-risk pathogens
        high_risk = [p for p in pathogens if p['risk_level'] == 'High']
        if high_risk:
            recommendations.append(
                "  • HIGH RISK: Clinical correlation required"
            )
            recommendations.append(
                "  • Consider targeted antimicrobial therapy if symptomatic"
            )
            
            # Specific recommendations
            for p in high_risk[:3]:
                if 'difficile' in p['species'].lower():
                    recommendations.append(
                        "  • C. difficile: Consider CDI testing, discontinue broad-spectrum antibiotics"
                    )
                elif 'aureus' in p['species'].lower():
                    recommendations.append(
                        "  • S. aureus: Screen for MRSA, consider isolation precautions"
                    )
        
        # Moderate risk
        moderate_risk = [p for p in pathogens if p['risk_level'] == 'Moderate']
        if moderate_risk:
            recommendations.append(
                "  • MODERATE RISK: Monitor for clinical symptoms"
            )
        
        # General
        recommendations.append(
            "  • Correlate with patient clinical presentation and history"
        )
        
        return '\n'.join(recommendations)
    
    def generate_visualization_data(self, bracken_file: Path, risk: Dict) -> Dict:
        """
        Generate data structure for visualization module.
        
        Returns:
            Dict with data ready for taxonomic_visualizer.py
        """
        bracken_df = pd.read_csv(bracken_file, sep='\t')
        
        # Top species for bar chart
        top_species = bracken_df.head(15).to_dict('records')
        
        # Genus-level for pie chart
        bracken_df['genus'] = bracken_df['name'].str.split().str[0]
        genus_data = bracken_df.groupby('genus').agg({
            'fraction_total_reads': 'sum'
        }).sort_values('fraction_total_reads', ascending=False).head(10)
        
        # Pathogen vs commensal
        pathogen_names = {p['species'] for p in risk['pathogens_detected']}
        pathogen_abundance = bracken_df[
            bracken_df['name'].isin(pathogen_names)
        ]['fraction_total_reads'].sum()
        commensal_abundance = bracken_df[
            ~bracken_df['name'].isin(pathogen_names)
        ]['fraction_total_reads'].sum()
        
        return {
            'top_species': top_species,
            'genus_composition': genus_data.to_dict(),
            'pathogen_vs_commensal': {
                'pathogenic': pathogen_abundance,
                'commensal': commensal_abundance
            },
            'pathogens': risk['pathogens_detected'],
            'total_species': len(bracken_df),
            'diversity': {
                'shannon': self._calculate_shannon(bracken_df),
                'simpson': self._calculate_simpson(bracken_df),
                'richness': len(bracken_df)
            }
        }
    
    def _calculate_shannon(self, df: pd.DataFrame) -> float:
        """Calculate Shannon diversity index."""
        import numpy as np
        abundances = df['fraction_total_reads'].values
        return float(-sum(abundances * np.log(abundances + 1e-10)))
    
    def _calculate_simpson(self, df: pd.DataFrame) -> float:
        """Calculate Simpson diversity index."""
        import numpy as np
        abundances = df['fraction_total_reads'].values
        return float(1 - sum(abundances ** 2))
    
    def export_table(self, bracken_file: Path, output_file: Path, risk: Dict):
        """Export complete taxonomy table to CSV."""
        bracken_df = pd.read_csv(bracken_file, sep='\t')
        
        # Add risk level column
        pathogen_dict = {p['species']: p['risk_level'] for p in risk['pathogens_detected']}
        bracken_df['risk_level'] = bracken_df['name'].map(pathogen_dict).fillna('Low')
        
        # Add clinical notes
        bracken_df['clinical_notes'] = bracken_df['name'].apply(self._get_clinical_notes)
        
        # Export
        bracken_df.to_csv(output_file, index=False)
        self.logger.info(f"Taxonomy table exported to {output_file}")


