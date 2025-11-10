"""
MetaQuest Risk Scoring Engine
Real-time pathogenicity risk calculation with three-tier integration
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional
from pathlib import Path
import json


class RiskScorer:
    """
    Dynamic risk scoring engine for pathogenicity assessment.
    
    Integrates three tiers:
    1. Taxonomic risk (from Bracken + pathogen organism list)
    2. Functional risk (from annotations + pathogen DB hits)
    3. ML prediction risk (from MetaQuest ML predictions)
    """
    
    def __init__(self, pathogen_organism_file: Optional[Path] = None):
        """
        Initialize risk scorer.
        
        Args:
            pathogen_organism_file: Path to organism risk level mapping
        """
        self.pathogen_organisms = self._load_pathogen_organisms(pathogen_organism_file)
        
        # Risk weights for final score calculation
        self.weights = {
            'taxonomic': 0.4,
            'functional': 0.3,
            'ml': 0.3
        }
        
        # Risk multipliers for special conditions
        self.multipliers = {
            'ml_and_pathogen_db': 1.2,
            'virulence_near_transposase': 1.15,
            'high_abundance_pathogen': 1.3
        }
        
        # Gene category risk scores
        self.gene_risk_scores = {
            'toxin': 100,
            'amr': 100,
            'secretion_system': 80,
            'adhesin': 60,
            'virulence': 80,
            'capsule': 70,
            'pathogen_general': 40,
            'transposase': 5  # Low direct risk, but enables HGT
        }
    
    def _load_pathogen_organisms(self, filepath: Optional[Path]) -> Dict[str, int]:
        """
        Load pathogen organism risk levels.
        
        Returns:
            Dict mapping organism name to risk score (0-100)
        """
        risk_mapping = {'high': 100, 'medium': 50, 'low': 25}
        organisms = {}
        
        if filepath and filepath.exists():
            with open(filepath, 'r') as f:
                for line in f:
                    if line.strip() and not line.startswith('#'):
                        parts = line.strip().split('\t')
                        if len(parts) >= 2:
                            organism = parts[0].strip()
                            risk_level = parts[1].strip().lower()
                            organisms[organism] = risk_mapping.get(risk_level, 0)
        
        return organisms
    
    def calculate_taxonomic_risk(self, bracken_df: pd.DataFrame) -> Dict:
        """
        Calculate taxonomic risk from Bracken results.
        
        Args:
            bracken_df: Bracken report DataFrame with 'name' and 'fraction_total_reads'
        
        Returns:
            Dict with risk score, pathogen list, and detailed breakdown
        """
        if bracken_df.empty:
            return {
                'score': 0,
                'pathogens_detected': [],
                'dominant_pathogen': None,
                'total_pathogen_abundance': 0,
                'details': []
            }
        
        risk_score = 0
        pathogens = []
        details = []
        
        for _, row in bracken_df.iterrows():
            species = row['name']
            abundance = row['fraction_total_reads']
            
            # Check if species is in pathogen list
            org_risk = self._match_organism(species)
            
            if org_risk > 0:
                contribution = abundance * org_risk
                risk_score += contribution
                
                pathogens.append({
                    'species': species,
                    'abundance': abundance,
                    'reads': row.get('new_est_reads', 0),
                    'risk_level': self._get_risk_label(org_risk),
                    'risk_score': org_risk,
                    'contribution': contribution
                })
                
                details.append({
                    'species': species,
                    'abundance_pct': abundance * 100,
                    'base_risk': org_risk,
                    'weighted_risk': contribution
                })
        
        # Sort pathogens by abundance
        pathogens.sort(key=lambda x: x['abundance'], reverse=True)
        
        # Calculate total pathogen abundance
        total_pathogen_abundance = sum(p['abundance'] for p in pathogens)
        
        # Determine dominant pathogen
        dominant = pathogens[0] if pathogens else None
        
        # Apply high abundance multiplier if needed
        if dominant and dominant['abundance'] > 0.05:  # >5% abundance
            risk_score *= self.multipliers['high_abundance_pathogen']
        
        # Cap at 100
        risk_score = min(risk_score, 100)
        
        return {
            'score': round(risk_score, 2),
            'pathogens_detected': pathogens,
            'dominant_pathogen': dominant,
            'total_pathogen_abundance': total_pathogen_abundance,
            'details': details
        }
    
    def calculate_functional_risk(self, 
                              functional_df: pd.DataFrame,
                              pathogen_hits_df: pd.DataFrame,
                              total_cds: int = None) -> Dict:
        """
        Calculate functional risk from annotations and pathogen DB hits.
        
        Args:
            functional_df: Functional annotations (SwissProt hits)
            pathogen_hits_df: Pathogen database hits
            total_cds: Total CDS count for normalization
        
        Returns:
            Dict with risk score and detailed breakdown
        """
        risk_score = 0
        details = {
            'virulence_factors': [],
            'amr_genes': [],
            'transposases': [],
            'secretion_systems': [],
            'other_pathogen_markers': []
        }
        
        # Count transposases from functional annotations
        transposase_genes = self._identify_transposases(functional_df)
        details['transposases'] = transposase_genes
        
        # Analyze pathogen DB hits
        if not pathogen_hits_df.empty:
            pathogen_analysis = self._analyze_pathogen_db_hits(pathogen_hits_df)
            
            details['virulence_factors'] = pathogen_analysis['virulence']
            details['amr_genes'] = pathogen_analysis['amr']
            details['secretion_systems'] = pathogen_analysis['secretion']
            details['other_pathogen_markers'] = pathogen_analysis['other']
        
        # *** FIX: Calculate risk as PERCENTAGE of CDS, not absolute count ***
        if total_cds and total_cds > 0:
            virulence_pct = (len(details['virulence_factors']) / total_cds) * 100
            amr_pct = (len(details['amr_genes']) / total_cds) * 100
            transposase_pct = (len(details['transposases']) / total_cds) * 100
            
            # Score components with diminishing returns
            # Typical gut microbiome: ~3% virulence, ~2% AMR, ~1% transposases
            virulence_score = min((virulence_pct / 5.0) * 40, 40)  # Max 40 points at 5%
            amr_score = min((amr_pct / 3.0) * 30, 30)              # Max 30 points at 3%
            transposase_score = min((transposase_pct / 2.0) * 30, 30)  # Max 30 points at 2%
            
            risk_score = virulence_score + amr_score + transposase_score
        else:
            # Fallback if total_cds not provided
            risk_score = min(
                len(details['virulence_factors']) * 0.05 +
                len(details['amr_genes']) * 0.04 +
                len(details['transposases']) * 0.01,
                100
            )
        
        # Already capped in calculation above
        risk_score = min(risk_score, 100)
        
        return {
            'score': round(risk_score, 2),
            'transposase_count': len(transposase_genes),
            'virulence_count': len(details['virulence_factors']),
            'amr_count': len(details['amr_genes']),
            'secretion_count': len(details['secretion_systems']),
            'details': details
        }
        
    def calculate_ml_risk(self, ml_predictions: List[Dict]) -> Dict:
        """
        Calculate risk from ML pathogen predictions.
        
        Args:
            ml_predictions: List of ML prediction dicts from JSON
        
        Returns:
            Dict with risk score and detailed breakdown
        """
    def calculate_ml_risk(self, ml_predictions: List[Dict]) -> Dict:
        if not ml_predictions:
            return {
                'score': 0,
                'total_sequences': 0,
                'pathogenic_count': 0,
                'high_confidence_pathogenic': 0,
                'high_confidence_proportion': 0,  # ADD THIS
                'details': [],
                'skipped': True  # ADD THIS FLAG
            }
            
        total = len(ml_predictions)
        pathogenic = [p for p in ml_predictions if p['prediction'] == 'Pathogenic']
        high_conf_pathogenic = [p for p in pathogenic if p['high_confidence']]
        
        # Risk is based on proportion of high-confidence pathogenic predictions
        risk_score = (len(high_conf_pathogenic) / total) * 100
        
        # Create detailed breakdown
        details = []
        for pred in high_conf_pathogenic:
            details.append({
                'sequence_id': pred['sequence_id'],
                'confidence': pred['confidence'],
                'pathogenic_probability': pred['pathogenic_probability'],
                'length': pred['sequence_length']
            })
        
        # Sort by confidence
        details.sort(key=lambda x: x['confidence'], reverse=True)
        
        return {
            'score': round(risk_score, 2),
            'total_sequences': total,
            'pathogenic_count': len(pathogenic),
            'high_confidence_pathogenic': len(high_conf_pathogenic),
            'high_confidence_proportion': len(high_conf_pathogenic) / total,
            'details': details
        }
    
    def calculate_integrated_risk(self,
                                  taxonomic_risk: Dict,
                                  functional_risk: Dict,
                                  ml_risk: Dict,
                                  functional_df: pd.DataFrame,
                                  pathogen_hits_df: pd.DataFrame,
                                  ml_predictions: List[Dict]) -> Dict:
        """
        Calculate final integrated risk score with cross-referencing.
        
        Returns:
            Dict with final risk score, tier scores, and convergent risks
        """
        # Base score from three tiers
        base_score = (
            taxonomic_risk['score'] * self.weights['taxonomic'] +
            functional_risk['score'] * self.weights['functional'] +
            ml_risk['score'] * self.weights['ml']
        )
        
        # Find convergent risks (genes that appear in multiple tiers)
        convergent = self._find_convergent_risks(
            functional_df, pathogen_hits_df, ml_predictions
        )
        
        # Apply multipliers
        final_score = base_score
        
        if convergent['ml_pathogen_db_matches']:
            final_score *= self.multipliers['ml_and_pathogen_db']
        
        if convergent['virulence_near_transposase']:
            final_score *= self.multipliers['virulence_near_transposase']
        
        # Cap at 100
        final_score = min(final_score, 100)
        
        # Determine risk level
        risk_level = self._get_risk_level(final_score)
        
        return {
            'final_score': round(final_score, 2),
            'risk_level': risk_level,
            'base_score': round(base_score, 2),
            'taxonomic_score': taxonomic_risk['score'],
            'functional_score': functional_risk['score'],
            'ml_score': ml_risk['score'],
            'multipliers_applied': {
                'ml_pathogen_db': convergent['ml_pathogen_db_matches'],
                'virulence_transposase': convergent['virulence_near_transposase']
            },
            'convergent_risks': convergent,
            'interpretation': self._generate_interpretation(
                final_score, taxonomic_risk, functional_risk, ml_risk
            )
        }
    
    def _match_organism(self, species_name: str) -> int:
        """Match species to pathogen organism list with fuzzy matching."""
        species_lower = species_name.lower()
        
        # Direct match (case-insensitive)
        for org_name, risk in self.pathogen_organisms.items():
            if org_name.lower() == species_lower:
                return risk
        
        # Try genus-level match
        genus = species_name.split()[0] if ' ' in species_name else species_name
        genus_lower = genus.lower()
        
        for org_name, risk in self.pathogen_organisms.items():
            org_lower = org_name.lower()
            # Check if pathogen DB entry starts with the genus
            if org_lower.startswith(genus_lower):
                return risk
            # Also check if the genus is anywhere in the name (for complex names)
            if genus_lower in org_lower:
                return risk
        
        # Partial match for known pathogens not in DB
        # Catch common pathogenic genera even without DB entry
        pathogenic_genera = {
            'clostridium': 100,
            'clostridioides': 100,
            'streptococcus': 75,
            'staphylococcus': 75,
            'campylobacter': 75,
            'salmonella': 100,
            'escherichia': 75,
            'shigella': 100,
            'listeria': 100,
            'mycobacterium': 100,
            'pseudomonas': 50,
            'klebsiella': 75,
            'acinetobacter': 75,
            'enterococcus': 50,
            'erysipelatoclostridium': 50,
            'lachnoclostridium': 25,
            'butyrivibrio': 0  # Commensal
        }
        
        genus_lower = genus.lower()
        if genus_lower in pathogenic_genera:
            return pathogenic_genera[genus_lower]
        
        return 0
    
    def _identify_transposases(self, functional_df: pd.DataFrame) -> List[Dict]:
        """Identify transposase genes from functional annotations."""
        if functional_df.empty:
            return []
        
        transposases = []
        keywords = ['transposase', 'IS66', 'IS30', 'IS982', 'IS256', 'IS element']
        
        for _, row in functional_df.iterrows():
            description = str(row.get('description', ''))
            if any(kw.lower() in description.lower() for kw in keywords):
                transposases.append({
                    'gene_id': row.get('query_id', row.iloc[0]),
                    'description': description,
                    'identity': row.get('identity', 0),
                    'evalue': row.get('evalue', 1)
                })
        
        return transposases
    
    def _analyze_pathogen_db_hits(self, pathogen_df: pd.DataFrame) -> Dict:
        """
        Analyze pathogen database hits and categorize.
        Works with raw Diamond output - no metadata needed.
        """
        categories = {
            'virulence': [],
            'amr': [],
            'secretion': [],
            'other': []
        }
        
        risk_contribution = 0
        
        # Simple keyword-based categorization from description
        virulence_keywords = ['toxin', 'virulence', 'hemolysin', 'adhesin', 'invasin', 
                             'fimbri', 'pili', 'flagell']
        amr_keywords = ['resistance', 'beta-lactamase', 'efflux', 'aminoglycoside',
                       'mec', 'van', 'tet', 'chloramphenicol']
        secretion_keywords = ['secretion system', 'type III', 'type IV', 'type VI',
                             'T3SS', 'T4SS', 'T6SS']
        
        # Get column index or name for description (last column in Diamond output)
        desc_col = 'description' if 'description' in pathogen_df.columns else pathogen_df.columns[-1]
        query_col = 'query_id' if 'query_id' in pathogen_df.columns else pathogen_df.columns[0]
        identity_col = 'identity' if 'identity' in pathogen_df.columns else pathogen_df.columns[2]
        
        for _, row in pathogen_df.iterrows():
            desc = str(row[desc_col]).lower()
            
            gene_info = {
                'gene_id': row[query_col],
                'description': row[desc_col],
                'organism': self._extract_organism(row[desc_col]),
                'identity': float(row[identity_col]) if pd.notna(row[identity_col]) else 0
            }
            
            # Categorize based on keywords in description
            if any(kw in desc for kw in virulence_keywords):
                categories['virulence'].append(gene_info)
                risk_contribution += self.gene_risk_scores['virulence']
            elif any(kw in desc for kw in amr_keywords):
                categories['amr'].append(gene_info)
                risk_contribution += self.gene_risk_scores['amr']
            elif any(kw in desc for kw in secretion_keywords):
                categories['secretion'].append(gene_info)
                risk_contribution += self.gene_risk_scores['secretion_system']
            else:
                # Any hit to pathogen DB is concerning
                categories['other'].append(gene_info)
                risk_contribution += self.gene_risk_scores['pathogen_general']
        
        return {
            'virulence': categories['virulence'],
            'amr': categories['amr'],
            'secretion': categories['secretion'],
            'other': categories['other'],
            'risk_contribution': risk_contribution
        }
    
    def _extract_organism(self, description: str) -> str:
        """Extract organism name from description [Organism name] format."""
        import re
        match = re.search(r'\[(.*?)\]$', description)
        return match.group(1) if match else 'Unknown'
    
    def _find_convergent_risks(self,
                               functional_df: pd.DataFrame,
                               pathogen_df: pd.DataFrame,
                               ml_predictions: List[Dict]) -> Dict:
        """Find genes that appear in multiple risk tiers."""
        convergent = {
            'ml_pathogen_db_matches': [],
            'virulence_near_transposase': False,
            'high_risk_clusters': []
        }
        
        # Get ML pathogenic predictions
        ml_pathogenic_ids = {
            p['sequence_id'] for p in ml_predictions 
            if p['prediction'] == 'Pathogenic' and p['high_confidence']
        }
        
        # Get pathogen DB hit IDs
        pathogen_hit_ids = set(pathogen_df['query_id'].unique()) if not pathogen_df.empty else set()
        
        # Find overlap
        overlap = ml_pathogenic_ids & pathogen_hit_ids
        if overlap:
            for gene_id in overlap:
                convergent['ml_pathogen_db_matches'].append({
                    'gene_id': gene_id,
                    'ml_confidence': next(
                        p['confidence'] for p in ml_predictions 
                        if p['sequence_id'] == gene_id
                    ),
                    'pathogen_db_hit': pathogen_df[
                        pathogen_df['query_id'] == gene_id
                    ].iloc[0]['description'] if not pathogen_df.empty else ''
                })
        
        # Check for virulence genes near transposases (genomic context)
        # This is simplified - would need genomic coordinates for full analysis
        if not pathogen_df.empty and not functional_df.empty:
            convergent['virulence_near_transposase'] = True  # Placeholder
        
        return convergent
    
    def _get_risk_level(self, score: float) -> str:
        """Convert numeric score to risk level label."""
        if score >= 61:
            return 'High'
        elif score >= 31:
            return 'Moderate'
        else:
            return 'Low'
    
    def _get_risk_label(self, score: int) -> str:
        """Convert organism risk score to label."""
        if score >= 75:
            return 'High'
        elif score >= 40:
            return 'Moderate'
        else:
            return 'Low'
    
    def _generate_interpretation(self,
                                final_score: float,
                                taxonomic: Dict,
                                functional: Dict,
                                ml: Dict) -> str:
        """Generate AI-style narrative interpretation."""
        risk_level = self._get_risk_level(final_score)
        
        # Build interpretation
        parts = []
        
        # Taxonomic context
        if taxonomic['pathogens_detected']:
            dominant = taxonomic['dominant_pathogen']
            parts.append(
                f"Taxonomic analysis detected {len(taxonomic['pathogens_detected'])} "
                f"known pathogenic species, with {dominant['species']} being most abundant "
                f"({dominant['abundance']*100:.2f}%)."
            )
        else:
            parts.append("No known pathogenic species detected at significant levels.")
        
        # Functional context
        if functional['transposase_count'] > 5:
            parts.append(
                f"High transposase content ({functional['transposase_count']} genes) "
                f"indicates active mobile genetic elements, suggesting potential for "
                f"horizontal gene transfer."
            )
        
        if functional['virulence_count'] > 0:
            parts.append(
                f"{functional['virulence_count']} virulence-associated genes identified "
                f"from pathogen database."
            )
        
        # ML context
        if ml['high_confidence_pathogenic'] > 0:
            parts.append(
                f"Machine learning analysis flagged {ml['high_confidence_pathogenic']} "
                f"sequences ({ml['high_confidence_proportion']*100:.1f}%) as potentially "
                f"pathogenic with high confidence."
            )
        
        # Final recommendation
        if risk_level == 'High':
            recommendation = (
                "Overall risk is HIGH. Clinical correlation and targeted investigation "
                "are strongly recommended."
            )
        elif risk_level == 'Moderate':
            recommendation = (
                "Overall risk is MODERATE. Monitor for clinical symptoms and consider "
                "further characterization if indicated."
            )
        else:
            recommendation = (
                "Overall risk is LOW. Sample appears to be dominated by commensal organisms. "
                "Routine monitoring is sufficient."
            )
        
        parts.append(recommendation)
        
        return ' '.join(parts)


def calculate_all_risks(bracken_file: Path,
                       functional_file: Path,
                       pathogen_hits_file: Path,
                       ml_predictions_file: Path,
                       pathogen_organism_file: Optional[Path] = None) -> Dict:
    """Calculate all risk scores from files."""
    scorer = RiskScorer(pathogen_organism_file)
    
    # Load data
    bracken_df = pd.read_csv(bracken_file, sep='\t')
    
    # *** FIX: Read files first to check column count ***
    functional_df = pd.read_csv(functional_file, sep='\t', header=None)
    pathogen_df = pd.read_csv(pathogen_hits_file, sep='\t', header=None)
    
    # Add column names - handle variable column counts
    ann_cols = ['query_id', 'subject_id', 'identity', 'length', 'mismatches', 
                'gaps', 'q_start', 'q_end', 's_start', 's_end', 'evalue', 
                'bitscore', 'description']
    
    # *** CRITICAL FIX: Only assign columns that exist ***
    if len(functional_df.columns) <= len(ann_cols):
        functional_df.columns = ann_cols[:len(functional_df.columns)]
    else:
        # More columns than expected - use positional names
        functional_df.columns = ann_cols + [f'extra_{i}' for i in range(len(functional_df.columns) - len(ann_cols))]
    
    if len(pathogen_df.columns) <= len(ann_cols):
        pathogen_df.columns = ann_cols[:len(pathogen_df.columns)]
    else:
        # More columns than expected - use positional names
        pathogen_df.columns = ann_cols + [f'extra_{i}' for i in range(len(pathogen_df.columns) - len(ann_cols))]
    
    # Rest of function continues...
    ml_predictions = []
    if ml_predictions_file and ml_predictions_file.exists():
        try:
            with open(ml_predictions_file) as f:
                ml_data = json.load(f)
                ml_predictions = ml_data.get('predictions', [])
        except Exception as e:
            print(f"Warning: Could not load ML predictions: {e}")
    
    # Calculate risks
    taxonomic = scorer.calculate_taxonomic_risk(bracken_df)
    total_cds = len(functional_df['query_id'].unique()) if not functional_df.empty else 0
    functional = scorer.calculate_functional_risk(functional_df, pathogen_df, total_cds)
    ml = scorer.calculate_ml_risk(ml_predictions)
    integrated = scorer.calculate_integrated_risk(
        taxonomic, functional, ml, functional_df, pathogen_df, ml_predictions
    )
    
    return {
        'taxonomic': taxonomic,
        'functional': functional,
        'ml': ml,
        'integrated': integrated
    }

