#!/usr/bin/env python3

"""
MetaQuest Protein Feature Extractor - Production Version
Environment-agnostic pathogen detection without ESM2/torch dependencies
Compatible with trained voting ensemble model
Complete implementation with all helper functions
"""

import pandas as pd
import numpy as np
import pickle
import joblib
import os
from collections import Counter
import json
from pathlib import Path
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

# Bio-related imports
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import molecular_weight
import re

# ML imports
from sklearn.base import BaseEstimator, TransformerMixin

class MetaQuestProteinFeatureExtractor(BaseEstimator, TransformerMixin):
    """
    Production MetaQuest protein feature extraction
    Compatible with trained voting ensemble model
    No ESM2/torch dependencies
    Environment-agnostic pathogen detection
    """
    
    def __init__(self, k_values=[2,3]):
        self.is_fitted = False
        self.k_values = k_values
        # Store precomputed lookup tables and patterns
        self.pathogen_motifs = {}
        self.virulence_patterns = {}
        self.conservation_cache = {}
        self.structural_cache = {}
        
    def fit(self, X=None, y=None):
        """Fit the feature extractor (for sklearn compatibility)"""
        print("🔧 Fitting MetaQuest Feature Extractor (Production)...")
        # Initialize lookup tables and patterns
        self._build_pathogen_motif_library()
        self._build_virulence_pattern_library()
        self._build_conservation_library()
        self._build_structural_pattern_library()
        self.is_fitted = True
        print("✅ Feature extractor fitted successfully!")
        return self
    
    def transform(self, sequences):
        """Transform sequences to features (sklearn compatible)"""
        if not self.is_fitted:
            print("⚠️ Feature extractor not fitted. Fitting now...")
            self.fit()
            
        if isinstance(sequences, pd.Series):
            sequences = sequences.tolist()
        elif isinstance(sequences, str):
            sequences = [sequences]
            
        feature_list = []
        print(f"🧬 Extracting features for {len(sequences)} sequences...")
        
        for seq in tqdm(sequences, desc="Processing sequences"):
            features = self.extract_features(seq)
            feature_list.append(features)
            
        return pd.DataFrame(feature_list)
    
    def fit_transform(self, X, y=None):
        """Fit and transform in one step"""
        return self.fit(X, y).transform(X)
    
    def get_feature_names_out(self, input_features=None):
        """Get output feature names for compatibility"""
        if not self.is_fitted:
            self.fit()
        
        # Generate feature names based on extraction method
        feature_names = []
        
        # Amino acid composition
        for aa in 'ACDEFGHIKLMNPQRSTVWY':
            feature_names.append(f'aa_percent_{aa}')
        
        # Traditional features
        traditional_features = [
            'length_log', 'length_sqrt', 'length_normalized',
            'hydrophobic_fraction', 'charged_fraction', 'polar_fraction', 
            'aromatic_fraction', 'small_fraction', 'tiny_fraction',
            'molecular_weight_norm', 'aromaticity', 'gravy', 'isoelectric_point',
            'charge_density', 'charge_distribution', 'net_charge',
            'helix_fraction', 'turn_fraction', 'sheet_fraction', 'structure_ratio'
        ]
        feature_names.extend(traditional_features)
        
        # Pathogen motifs
        for motif_class in ['type3_secretion', 'toxin_binding', 'cell_adhesion', 
                           'membrane_disruption', 'immune_evasion', 'host_invasion']:
            feature_names.append(f'{motif_class}_motifs')
        
        # Pathogen features
        pathogen_features = [
            'virulence_score', 'toxin_like_regions', 'host_binding_potential',
            'immune_recognition_sites', 'amr_potential', 'biofilm_potential'
        ]
        feature_names.extend(pathogen_features)
        
        # Evolutionary features
        evolutionary_features = [
            'conservation_mean', 'conservation_std', 'conservation_max', 
            'conservation_entropy', 'selection_pressure', 'codon_adaptation_index',
            'gc_content', 'hgt_signal', 'mobile_element_association', 
            'taxonomic_breadth', 'ancient_origin'
        ]
        feature_names.extend(evolutionary_features)
        
        # Structural features
        structural_features = [
            'disorder_fraction', 'disorder_longest_region', 'disorder_clusters',
            'tm_helices_count', 'tm_fraction', 'tm_topology_score',
            'signal_peptide_score', 'signal_peptide_type', 'coiled_coil_score',
            'beta_strand_content', 'turn_propensity'
        ]
        feature_names.extend(structural_features)
        
        # Length-optimized features
        length_features = [
            'n_terminal_hydrophobic', 'n_terminal_basic', 'n_terminal_signal_score',
            'c_terminal_charged', 'c_terminal_polar', 'core_complexity',
            'core_conservation', 'optimal_length_score'
        ]
        feature_names.extend(length_features)
        
        # Motif densities
        for motif in ['RGD', 'KKK', 'DDD', 'LLL', 'FFF']:
            feature_names.append(f'{motif}_density')
        
        # Advanced pathogen signatures
        advanced_features = [
            'sec_pathway_score', 'tat_pathway_score', 'adhesin_score',
            'invasin_score', 'protease_score', 'lipase_score',
            'heat_shock_score', 'oxidative_stress_score', 'iron_acquisition_score',
            'nutrient_scavenging_score', 'beta_lactamase_score', 'efflux_pump_score'
        ]
        feature_names.extend(advanced_features)
        
        return feature_names
    
    def _build_pathogen_motif_library(self):
        """Build comprehensive pathogen motif library"""
        self.pathogen_motifs = {
            'type3_secretion': {
                'patterns': ['PXXP', 'SxxE', 'ExE', 'YxxY'],
                'weights': [1.0, 1.0, 0.8, 1.2],
                'description': 'Type III secretion system motifs'
            },
            'toxin_binding': {
                'patterns': ['RGD', 'LDV', 'NGR', 'KGD'],
                'weights': [2.0, 1.5, 1.8, 1.5],
                'description': 'Toxin binding motifs'
            },
            'cell_adhesion': {
                'patterns': ['KXXK', 'RXXR', 'RGDX'],
                'weights': [2.0, 2.0, 2.5],
                'description': 'Cell adhesion motifs'
            },
            'membrane_disruption': {
                'patterns': ['LXXL', 'AXXXA', 'WXXW'],
                'weights': [2.0, 1.5, 2.5],
                'description': 'Membrane disruption motifs'
            },
            'immune_evasion': {
                'patterns': ['CxxC', 'HxxH', 'CxC'],
                'weights': [1.5, 1.2, 1.0],
                'description': 'Immune evasion motifs'
            },
            'host_invasion': {
                'patterns': ['NPXY', 'YxxΦ', 'FXXF'],
                'weights': [2.5, 2.0, 1.8],
                'description': 'Host cell invasion motifs'
            }
        }
    
    def _kmer_features(self, sequence: str) -> dict:
        """
        Calculate k-mer (dipeptide, tripeptide, etc.) frequencies.
        """
        features = {}
        seq_len = len(sequence)
        
        for k in self.k_values:
            if seq_len < k: continue
            
            kmers = Counter()
            for i in range(seq_len - k + 1):
                kmer = sequence[i:i+k]
                kmers[kmer] += 1
            
            total_kmers = sum(kmers.values())
            for kmer, count in kmers.items():
                features[f'kmer_{kmer}'] = count / total_kmers
                
        return features
    
    def _build_virulence_pattern_library(self):
        """Build virulence factor pattern library"""
        self.virulence_patterns = {
            'secretion_signals': {
                'patterns': ['SS', 'SP', 'LP', 'MKKL'],
                'weights': [1.0, 1.0, 0.8, 1.2]
            },
            'adhesion_motifs': {
                'patterns': ['RGD', 'LDV', 'RGDX'],
                'weights': [2.0, 1.5, 1.8]
            },
            'toxin_motifs': {
                'patterns': ['CxxC', 'HxxH', 'CxC'],
                'weights': [1.5, 1.2, 1.0]
            },
            'invasion_signals': {
                'patterns': ['NPXY', 'YxxΦ'],
                'weights': [2.5, 2.0]
            }
        }
    
    def _build_conservation_library(self):
        """Build conservation scoring library"""
        self.conservation_cache = {
            'highly_conserved': 'CWFYH',  # Structurally important
            'moderately_conserved': 'ILMV',  # Hydrophobic core
            'variable': 'GPST',  # Flexible regions
            'charged_conserved': 'DEKR'  # Functionally important
        }
    
    def _build_structural_pattern_library(self):
        """Build structural pattern library"""
        self.structural_cache = {
            'disorder_prone': 'PGQSRKEND',
            'helix_prone': 'ALEK',
            'sheet_prone': 'VTIFWY',
            'turn_prone': 'PGNDS',
            'hydrophobic_core': 'AILMFPWV',
            'surface_exposed': 'DEKRHNQST',
            'coiled_coil': 'ALIEKR',
            'tm_prone': 'AILMFPWV'
        }
    
    def extract_features(self, sequence):
        """
        Extract comprehensive pathogen detection features
        Main feature extraction method for compatibility with trained model
        """
        try:
            # Validate and clean protein sequence
            cleaned_sequence = self._validate_and_clean_sequence(sequence)
            if not cleaned_sequence:
                return self._get_default_features()
            
            features = {}
            
            # 1. Enhanced Traditional Features
            features.update(self._enhanced_traditional_features(cleaned_sequence))
            
            # 2. MetaQuest-Specific Pathogen Features
            features.update(self._metaquest_pathogen_features(cleaned_sequence))
            
            # 3. Advanced Evolutionary Features
            features.update(self._advanced_evolutionary_features(cleaned_sequence))
            
            # 4. Structural Prediction Features
            features.update(self._enhanced_structural_features(cleaned_sequence))
            
            # 5. Length-Optimized Features
            features.update(self._length_optimized_features(cleaned_sequence))
            
            # 6. Advanced Pathogen-Specific Features
            features.update(self._advanced_pathogen_signatures(cleaned_sequence))

            # 7. K-mer Frequency Features
            features.update(self._kmer_features(cleaned_sequence))
            
            return features
            
        except Exception as e:
            print(f"Warning: Error in feature extraction: {e}")
            return self._get_default_features()
    
    def _validate_and_clean_sequence(self, sequence):
        """Validate and clean protein sequence"""
        if not sequence or len(sequence) == 0:
            return None
        
        cleaned_sequence = ''.join([aa for aa in str(sequence).upper() if aa in 'ACDEFGHIKLMNPQRSTVWY'])
        
        if len(cleaned_sequence) < 50 or len(cleaned_sequence) > 2000:
            return None
            
        return cleaned_sequence
    
    def _enhanced_traditional_features(self, sequence):
        """Enhanced traditional features with error handling"""
        features = {}
        
        try:
            analyzer = ProteinAnalysis(sequence)
            seq_len = len(sequence)
            
            # Length features
            features['length_log'] = np.log(seq_len + 1)
            features['length_sqrt'] = np.sqrt(seq_len)
            features['length_normalized'] = (seq_len - 750) / 250
            
            # Amino acid composition
            aa_percent = analyzer.get_amino_acids_percent()
            for aa in 'ACDEFGHIKLMNPQRSTVWY':
                features[f'aa_percent_{aa}'] = aa_percent.get(aa, 0)
            
            # Amino acid groups
            features['hydrophobic_fraction'] = sum(aa_percent.get(aa, 0) for aa in 'AILMFPWV')
            features['charged_fraction'] = sum(aa_percent.get(aa, 0) for aa in 'DEKRH')
            features['polar_fraction'] = sum(aa_percent.get(aa, 0) for aa in 'NQST')
            features['aromatic_fraction'] = sum(aa_percent.get(aa, 0) for aa in 'FWY')
            features['small_fraction'] = sum(aa_percent.get(aa, 0) for aa in 'AGST')
            features['tiny_fraction'] = sum(aa_percent.get(aa, 0) for aa in 'AGS')
            
            # Physicochemical properties
            features['molecular_weight_norm'] = analyzer.molecular_weight() / seq_len
            features['aromaticity'] = analyzer.aromaticity()
            features['gravy'] = analyzer.gravy()
            features['isoelectric_point'] = analyzer.isoelectric_point()
            
            # Charge properties
            features['charge_density'] = self._calculate_enhanced_charge_density(sequence)
            features['charge_distribution'] = self._calculate_charge_distribution(sequence)
            features['net_charge'] = self._calculate_net_charge(sequence)
            
            # Secondary structure
            try:
                sec_struct = analyzer.secondary_structure_fraction()
                features['helix_fraction'] = sec_struct[0]
                features['turn_fraction'] = sec_struct[1]
                features['sheet_fraction'] = sec_struct[2]
                features['structure_ratio'] = sec_struct[0] / (sec_struct[2] + 0.001)
            except:
                features.update({
                    'helix_fraction': 0.0, 'turn_fraction': 0.0,
                    'sheet_fraction': 0.0, 'structure_ratio': 1.0
                })
            
        except Exception as e:
            print(f"Warning: Traditional features error: {e}")
            features.update(self._get_default_traditional_features())
        
        return features
    
    def _metaquest_pathogen_features(self, sequence):
        """MetaQuest pathogen detection features"""
        features = {}
        
        # Enhanced pathogenicity motifs using stored library
        for motif_class, motif_data in self.pathogen_motifs.items():
            motifs = motif_data['patterns']
            weights = motif_data['weights']
            total_score = 0
            for motif, weight in zip(motifs, weights):
                count = self._count_motif_pattern(sequence, motif)
                total_score += count * weight
            features[f'{motif_class}_motifs'] = (total_score / len(sequence) * 1000)
        
        # Virulence factors
        features['virulence_score'] = self._calculate_enhanced_virulence_score(sequence)
        features['toxin_like_regions'] = float(self._detect_enhanced_toxin_regions(sequence))
        
        # Host interaction features
        features['host_binding_potential'] = self._predict_enhanced_host_binding(sequence)
        features['immune_recognition_sites'] = self._predict_enhanced_immune_sites(sequence)
        
        # Additional pathogen features
        features['amr_potential'] = self._detect_amr_signatures(sequence)
        features['biofilm_potential'] = self._detect_biofilm_signatures(sequence)
        
        return features
    
    def _advanced_evolutionary_features(self, sequence):
        """Advanced evolutionary features"""
        features = {}
        
        # Conservation analysis
        conservation_scores = self._calculate_enhanced_conservation_scores(sequence)
        features['conservation_mean'] = np.mean(conservation_scores)
        features['conservation_std'] = np.std(conservation_scores)
        features['conservation_max'] = np.max(conservation_scores)
        features['conservation_entropy'] = self._calculate_conservation_entropy(conservation_scores)
        
        # Evolutionary pressure
        features['selection_pressure'] = self._calculate_enhanced_selection_pressure(sequence)
        features['codon_adaptation_index'] = self._calculate_enhanced_codon_adaptation(sequence)
        features['gc_content'] = self._calculate_gc_content_proxy(sequence)
        
        # HGT signatures
        features['hgt_signal'] = self._detect_enhanced_hgt_signals(sequence)
        features['mobile_element_association'] = self._detect_mobile_elements(sequence)
        
        # Phylogenetic features
        features['taxonomic_breadth'] = self._estimate_taxonomic_breadth(sequence)
        features['ancient_origin'] = self._detect_ancient_protein_signatures(sequence)
        
        return features
    
    def _enhanced_structural_features(self, sequence):
        """Enhanced structural prediction features"""
        features = {}
        
        # Disorder prediction
        disorder_scores = self._predict_enhanced_disorder_regions(sequence)
        features['disorder_fraction'] = np.mean(disorder_scores > 0.5)
        features['disorder_longest_region'] = self._longest_disordered_region(disorder_scores)
        features['disorder_clusters'] = self._count_disorder_clusters(disorder_scores)
        
        # Transmembrane prediction
        tm_regions = self._predict_enhanced_transmembrane_regions(sequence)
        features['tm_helices_count'] = len(tm_regions)
        features['tm_fraction'] = sum(end - start for start, end in tm_regions) / len(sequence) if tm_regions else 0
        features['tm_topology_score'] = self._calculate_tm_topology_score(tm_regions, len(sequence))
        
        # Signal peptides
        features['signal_peptide_score'] = self._predict_enhanced_signal_peptide(sequence)
        features['signal_peptide_type'] = self._classify_signal_peptide_type(sequence)
        
        # Other structural features
        features['coiled_coil_score'] = self._predict_enhanced_coiled_coil(sequence)
        features['beta_strand_content'] = self._predict_beta_strand_content(sequence)
        features['turn_propensity'] = self._calculate_turn_propensity(sequence)
        
        return features
    
    def _length_optimized_features(self, sequence):
        """Features optimized for sequence length"""
        features = {}
        seq_len = len(sequence)
        
        # Positional features
        n_terminal_20 = sequence[:20] if len(sequence) >= 20 else sequence
        c_terminal_20 = sequence[-20:] if len(sequence) >= 20 else sequence
        middle_region = sequence[seq_len//3:2*seq_len//3] if seq_len >= 60 else sequence
        
        # N-terminal features
        features['n_terminal_hydrophobic'] = sum(1 for aa in n_terminal_20 if aa in 'AILMFPWV') / len(n_terminal_20)
        features['n_terminal_basic'] = sum(1 for aa in n_terminal_20 if aa in 'KR') / len(n_terminal_20)
        features['n_terminal_signal_score'] = self._calculate_signal_score(n_terminal_20)
        
        # C-terminal features
        features['c_terminal_charged'] = sum(1 for aa in c_terminal_20 if aa in 'DEKRH') / len(c_terminal_20)
        features['c_terminal_polar'] = sum(1 for aa in c_terminal_20 if aa in 'NQST') / len(c_terminal_20)
        
        # Core region features
        features['core_complexity'] = self._calculate_sequence_complexity(middle_region)
        features['core_conservation'] = self._estimate_core_conservation(middle_region)
        
        # Motif densities
        for motif in ['RGD', 'KKK', 'DDD', 'LLL', 'FFF']:
            features[f'{motif}_density'] = sequence.count(motif) / seq_len * 1000
        
        # Length optimization score
        features['optimal_length_score'] = self._calculate_optimal_length_score(seq_len)
        
        return features
    
    def _advanced_pathogen_signatures(self, sequence):
        """Advanced pathogen-specific signatures"""
        features = {}
        
        # Secretion system signatures
        features['sec_pathway_score'] = self._calculate_sec_pathway_score(sequence)
        features['tat_pathway_score'] = self._calculate_tat_pathway_score(sequence)
        
        # Virulence factor classes
        features['adhesin_score'] = self._calculate_adhesin_score(sequence)
        features['invasin_score'] = self._calculate_invasin_score(sequence)
        features['protease_score'] = self._calculate_protease_score(sequence)
        features['lipase_score'] = self._calculate_lipase_score(sequence)
        
        # Stress response signatures
        features['heat_shock_score'] = self._calculate_heat_shock_score(sequence)
        features['oxidative_stress_score'] = self._calculate_oxidative_stress_score(sequence)
        
        # Metabolic pathogen signatures
        features['iron_acquisition_score'] = self._calculate_iron_acquisition_score(sequence)
        features['nutrient_scavenging_score'] = self._calculate_nutrient_scavenging_score(sequence)
        
        # Antibiotic resistance signatures
        features['beta_lactamase_score'] = self._calculate_beta_lactamase_score(sequence)
        features['efflux_pump_score'] = self._calculate_efflux_pump_score(sequence)
        
        return features
    
    # ===================== ALL HELPER METHODS =====================
    
    def _calculate_enhanced_charge_density(self, sequence):
        """Enhanced charge density calculation"""
        positive = sequence.count('K') + sequence.count('R') + sequence.count('H')
        negative = sequence.count('D') + sequence.count('E')
        return abs(positive - negative) / len(sequence)
    
    def _calculate_charge_distribution(self, sequence):
        """Calculate charge distribution pattern"""
        charges = []
        for aa in sequence:
            if aa in 'KRH':
                charges.append(1)
            elif aa in 'DE':
                charges.append(-1)
            else:
                charges.append(0)
        charge_changes = sum(1 for i in range(len(charges)-1) if charges[i] != charges[i+1])
        return charge_changes / len(sequence)
    
    def _calculate_net_charge(self, sequence):
        """Calculate net charge at physiological pH"""
        positive = sequence.count('K') + sequence.count('R') + sequence.count('H')
        negative = sequence.count('D') + sequence.count('E')
        return (positive - negative) / len(sequence)
    
    def _count_motif_pattern(self, sequence, pattern):
        """Count motif occurrences with wildcard support"""
        regex_pattern = pattern.replace('x', '.').replace('X', '.').replace('Φ', '[AILMFPWV]')
        try:
            return len(re.findall(regex_pattern, sequence))
        except:
            return 0
    
    def _calculate_enhanced_virulence_score(self, sequence):
        """Enhanced virulence factor scoring using stored patterns"""
        total_score = 0
        for category, pattern_data in self.virulence_patterns.items():
            motifs = pattern_data['patterns']
            weights = pattern_data['weights']
            for motif, weight in zip(motifs, weights):
                count = self._count_motif_pattern(sequence, motif)
                total_score += count * weight
        return total_score / len(sequence) * 1000
    
    def _detect_enhanced_toxin_regions(self, sequence):
        """Enhanced toxin region detection"""
        cysteine_density = sequence.count('C') / len(sequence)
        disulfide_potential = self._calculate_disulfide_potential(sequence)
        small_molecule_binding = self._detect_small_molecule_binding_sites(sequence)
        toxin_score = (cysteine_density * 2 + disulfide_potential + small_molecule_binding) / 4
        return toxin_score > 0.15
    
    def _predict_enhanced_host_binding(self, sequence):
        """Enhanced host binding prediction"""
        binding_residues = self.structural_cache['surface_exposed']
        binding_motifs = ['RGD', 'LDV', 'NGR', 'KGD']
        
        # Basic binding residue content
        binding_content = sum(1 for aa in sequence if aa in binding_residues) / len(sequence)
        
        # Specific binding motifs
        motif_score = sum(self._count_motif_pattern(sequence, motif) for motif in binding_motifs)
        motif_score = motif_score / len(sequence) * 1000
        
        # Surface accessibility approximation
        surface_score = self._estimate_surface_accessibility(sequence)
        
        return (binding_content * 0.4 + motif_score * 0.4 + surface_score * 0.2)
    
    def _predict_enhanced_immune_sites(self, sequence):
        """Enhanced immune recognition site prediction"""
        epitope_patterns = ['PXXP', 'GxxG', 'LXXL', 'KXXK', 'RXXR']
        immune_score = sum(self._count_motif_pattern(sequence, pattern) for pattern in epitope_patterns)
        
        # Add hydrophobic patches (common in epitopes)
        hydrophobic_patches = self._detect_hydrophobic_patches(sequence)
        
        # Add charged clusters
        charged_clusters = self._detect_charged_clusters(sequence)
        
        return immune_score + hydrophobic_patches * 0.5 + charged_clusters * 0.3
    
    def _detect_amr_signatures(self, sequence):
        """Detect antibiotic resistance signatures"""
        amr_score = 0
        
        # Look for β-lactamase signatures
        if 'SxxK' in sequence or 'SDN' in sequence:
            amr_score += 1
        
        # Look for efflux pump signatures
        if sequence.count('L') / len(sequence) > 0.15:  # High leucine content
            amr_score += 0.5
        
        return amr_score
    
    def _detect_biofilm_signatures(self, sequence):
        """Detect biofilm formation signatures"""
        # Simplified biofilm detection based on adhesion properties
        adhesion_motifs = ['RGD', 'NPQG', 'LPxTG']
        biofilm_score = sum(self._count_motif_pattern(sequence, motif) for motif in adhesion_motifs)
        
        # Add surface protein indicators
        if sequence.count('P') / len(sequence) > 0.1:  # High proline content
            biofilm_score += 0.5
        
        return biofilm_score
    
    def _calculate_enhanced_conservation_scores(self, sequence):
        """Calculate conservation scores using cached patterns"""
        conservation = []
        for aa in sequence:
            if aa in self.conservation_cache['highly_conserved']:
                conservation.append(np.random.beta(5, 2))
            elif aa in self.conservation_cache['variable']:
                conservation.append(np.random.beta(2, 5))
            else:
                conservation.append(np.random.beta(3, 3))
        return np.array(conservation)
    
    def _calculate_conservation_entropy(self, conservation_scores):
        """Calculate conservation entropy"""
        # Binned entropy calculation
        bins = np.histogram(conservation_scores, bins=5)[0]
        bins = bins[bins > 0]  # Remove empty bins
        probs = bins / bins.sum()
        return -np.sum(probs * np.log2(probs + 1e-10))
    
    def _calculate_enhanced_selection_pressure(self, sequence):
        """Calculate enhanced selection pressure"""
        # Multiple indicators of selection pressure
        hydrophobic = sum(1 for aa in sequence if aa in 'AILMFPWV')
        charged = sum(1 for aa in sequence if aa in 'DEKRH')
        aromatic = sum(1 for aa in sequence if aa in 'FWY')
        return (hydrophobic + charged + aromatic * 2) / len(sequence)
    
    def _calculate_enhanced_codon_adaptation(self, sequence):
        """Enhanced codon adaptation index"""
        rare_codons = sum(1 for aa in sequence if aa in 'WMYC')
        common_codons = sum(1 for aa in sequence if aa in 'ALSGTRE')
        if rare_codons + common_codons == 0:
            return 0.5
        return common_codons / (rare_codons + common_codons)
    
    def _calculate_gc_content_proxy(self, sequence):
        """Calculate GC content proxy from amino acids"""
        # GC-rich codons tend to encode certain amino acids
        gc_rich_aa = 'GAPRTW'
        return sum(1 for aa in sequence if aa in gc_rich_aa) / len(sequence)
    
    def _detect_enhanced_hgt_signals(self, sequence):
        """Detect horizontal gene transfer signals"""
        # Unusual amino acid composition as HGT indicator
        rare_aa_fraction = sum(1 for aa in sequence if aa in 'WMYC') / len(sequence)
        return min(rare_aa_fraction * 3, 1.0)
    
    def _detect_mobile_elements(self, sequence):
        """Detect mobile genetic element associations"""
        # Simplified mobile element detection
        mobile_signatures = ['DDE', 'HTH']  # Common in transposases
        return sum(self._count_motif_pattern(sequence, sig) for sig in mobile_signatures)
    
    def _estimate_taxonomic_breadth(self, sequence):
        """Estimate taxonomic distribution breadth"""
        # Based on conservation and composition
        conservation = self._calculate_enhanced_conservation_scores(sequence)
        breadth_score = np.mean(conservation) * (1 - np.std(conservation))
        return breadth_score
    
    def _detect_ancient_protein_signatures(self, sequence):
        """Detect ancient protein signatures"""
        ancient_motifs = ['WxxW', 'FxxF', 'YxxY', 'HxxH']
        return sum(self._count_motif_pattern(sequence, motif) for motif in ancient_motifs)
    
    def _predict_enhanced_disorder_regions(self, sequence):
        """Enhanced disorder prediction using cached patterns"""
        disorder_prone = self.structural_cache['disorder_prone']
        scores = np.array([1 if aa in disorder_prone else 0 for aa in sequence])
        
        # Smooth with window
        window = min(7, len(sequence))
        if window > 1:
            smoothed = np.convolve(scores, np.ones(window)/window, mode='same')
        else:
            smoothed = scores
        return smoothed
    
    def _longest_disordered_region(self, disorder_scores):
        """Find longest disordered region"""
        current_length = 0
        max_length = 0
        for score in disorder_scores:
            if score > 0.5:
                current_length += 1
                max_length = max(max_length, current_length)
            else:
                current_length = 0
        return max_length / len(disorder_scores) if len(disorder_scores) > 0 else 0
    
    def _count_disorder_clusters(self, disorder_scores):
        """Count disorder clusters"""
        clusters = 0
        in_cluster = False
        for score in disorder_scores:
            if score > 0.5 and not in_cluster:
                clusters += 1
                in_cluster = True
            elif score <= 0.5:
                in_cluster = False
        return clusters
    
    def _predict_enhanced_transmembrane_regions(self, sequence):
        """Enhanced transmembrane prediction using cached patterns"""
        hydrophobic = self.structural_cache['tm_prone']
        regions = []
        i = 0
        while i < len(sequence) - 19:
            window = sequence[i:i+20]
            hydrophobic_count = sum(1 for aa in window if aa in hydrophobic)
            if hydrophobic_count > 15:
                regions.append((i, i+20))
                i += 20
            else:
                i += 1
        return regions
    
    def _calculate_tm_topology_score(self, tm_regions, seq_length):
        """Calculate transmembrane topology score"""
        if not tm_regions:
            return 0
        # Score based on number and distribution of TM regions
        num_regions = len(tm_regions)
        coverage = sum(end - start for start, end in tm_regions) / seq_length
        return (num_regions * 0.1 + coverage * 0.9)
    
    def _predict_enhanced_signal_peptide(self, sequence):
        """Enhanced signal peptide prediction"""
        if len(sequence) < 20:
            return 0
        
        n_term = sequence[:25]  # Extended N-terminal region
        
        # Hydrophobic core region (positions 10-20)
        if len(n_term) >= 20:
            h_region = n_term[9:19]
            hydrophobic = sum(1 for aa in h_region if aa in 'AILMFPWV')
            h_score = hydrophobic / len(h_region)
        else:
            h_score = 0
        
        # N-terminal positive region
        n_region = n_term[:8]
        positive = sum(1 for aa in n_region if aa in 'KR')
        n_score = positive / len(n_region)
        
        # C-terminal cleavage site
        if len(n_term) >= 20:
            c_region = n_term[15:20]
            small_aa = sum(1 for aa in c_region if aa in 'AGSTC')
            c_score = small_aa / len(c_region)
        else:
            c_score = 0
        
        return (n_score * 0.3 + h_score * 0.5 + c_score * 0.2)
    
    def _classify_signal_peptide_type(self, sequence):
        """Classify signal peptide type"""
        sp_score = self._predict_enhanced_signal_peptide(sequence)
        if sp_score > 0.7:
            return 2  # Strong signal peptide
        elif sp_score > 0.4:
            return 1  # Weak signal peptide
        else:
            return 0  # No signal peptide
    
    def _predict_enhanced_coiled_coil(self, sequence):
        """Enhanced coiled coil prediction using cached patterns"""
        coiled_coil_residues = self.structural_cache['coiled_coil']
        
        # Look for heptad repeats pattern
        heptad_score = 0
        for i in range(0, len(sequence) - 7, 7):
            heptad = sequence[i:i+7]
            if len(heptad) == 7:
                # Positions 1 and 4 (a and d) should be hydrophobic
                if heptad[0] in 'AILMFPWV' and heptad[3] in 'AILMFPWV':
                    heptad_score += 1
        
        basic_score = sum(1 for aa in sequence if aa in coiled_coil_residues) / len(sequence)
        return (heptad_score / (len(sequence) // 7 + 1) + basic_score) / 2
    
    def _predict_beta_strand_content(self, sequence):
        """Predict beta strand content using cached patterns"""
        beta_prone = self.structural_cache['sheet_prone']
        beta_content = sum(1 for aa in sequence if aa in beta_prone) / len(sequence)
        return beta_content
    
    def _calculate_turn_propensity(self, sequence):
        """Calculate turn propensity using cached patterns"""
        turn_prone = self.structural_cache['turn_prone']
        turn_content = sum(1 for aa in sequence if aa in turn_prone) / len(sequence)
        return turn_content
    
    def _calculate_signal_score(self, n_terminal):
        """Calculate signal sequence score"""
        if len(n_terminal) < 10:
            return 0
        hydrophobic_score = sum(1 for aa in n_terminal if aa in 'AILMFPWV') / len(n_terminal)
        basic_score = sum(1 for aa in n_terminal[:5] if aa in 'KR') / 5
        return (hydrophobic_score + basic_score) / 2
    
    def _calculate_sequence_complexity(self, sequence):
        """Calculate sequence complexity"""
        if len(sequence) == 0:
            return 0
        
        # Shannon entropy-based complexity
        aa_counts = {}
        for aa in sequence:
            aa_counts[aa] = aa_counts.get(aa, 0) + 1
        
        entropy = 0
        for count in aa_counts.values():
            p = count / len(sequence)
            entropy -= p * np.log2(p)
        
        # Normalize by maximum possible entropy (log2(20))
        return entropy / np.log2(20)
    
    def _estimate_core_conservation(self, sequence):
        """Estimate conservation of core region using cached patterns"""
        conserved_residues = self.conservation_cache['highly_conserved']
        conservation_score = sum(1 for aa in sequence if aa in conserved_residues) / len(sequence)
        return conservation_score
    
    def _calculate_optimal_length_score(self, seq_len):
        """Calculate how close sequence is to optimal length"""
        optimal_center = 750  # Center of 500-1000 range
        distance_from_optimal = abs(seq_len - optimal_center) / 250
        return max(0, 1 - distance_from_optimal)
    
    def _calculate_disulfide_potential(self, sequence):
        """Calculate disulfide bond potential"""
        cys_positions = [i for i, aa in enumerate(sequence) if aa == 'C']
        if len(cys_positions) < 2:
            return 0.0
        
        optimal_spacings = 0
        for i in range(len(cys_positions) - 1):
            spacing = cys_positions[i+1] - cys_positions[i]
            if 8 <= spacing <= 25:  # Optimal spacing for disulfide bonds
                optimal_spacings += 1
        
        return optimal_spacings / len(cys_positions) if cys_positions else 0.0
    
    def _detect_small_molecule_binding_sites(self, sequence):
        """Detect small molecule binding sites"""
        binding_motifs = ['GxGxxG', 'DxxG', 'NxxD', 'HxxH']
        binding_score = sum(self._count_motif_pattern(sequence, motif) for motif in binding_motifs)
        return binding_score / len(sequence) * 100
    
    def _estimate_surface_accessibility(self, sequence):
        """Estimate surface accessibility using cached patterns"""
        surface_residues = self.structural_cache['surface_exposed']
        return sum(1 for aa in sequence if aa in surface_residues) / len(sequence)
    
    def _detect_hydrophobic_patches(self, sequence):
        """Detect hydrophobic patches using cached patterns"""
        hydrophobic = self.structural_cache['hydrophobic_core']
        patch_count = 0
        
        # Look for consecutive hydrophobic residues
        consecutive = 0
        for aa in sequence:
            if aa in hydrophobic:
                consecutive += 1
            else:
                if consecutive >= 4:  # Patch of 4+ hydrophobic residues
                    patch_count += 1
                consecutive = 0
        
        # Check final patch
        if consecutive >= 4:
            patch_count += 1
        
        return patch_count
    
    def _detect_charged_clusters(self, sequence):
        """Detect charged residue clusters"""
        charged = 'DEKRH'
        cluster_count = 0
        
        # Look for charged clusters
        consecutive = 0
        for aa in sequence:
            if aa in charged:
                consecutive += 1
            else:
                if consecutive >= 3:  # Cluster of 3+ charged residues
                    cluster_count += 1
                consecutive = 0
        
        # Check final cluster
        if consecutive >= 3:
            cluster_count += 1
        
        return cluster_count
    
    # === ADVANCED PATHOGEN SIGNATURE METHODS ===
    
    def _calculate_sec_pathway_score(self, sequence):
        """Calculate Sec pathway signal score"""
        sec_motifs = ['LXXL', 'AXXXA', 'MXXM']
        return sum(self._count_motif_pattern(sequence, motif) for motif in sec_motifs)
    
    def _calculate_tat_pathway_score(self, sequence):
        """Calculate Tat pathway signal score"""
        tat_motifs = ['RRxFLK', 'RRXXX']
        return sum(self._count_motif_pattern(sequence, motif) for motif in tat_motifs)
    
    def _calculate_adhesin_score(self, sequence):
        """Calculate adhesin protein score"""
        adhesin_motifs = ['RGD', 'LDV', 'YIGSR', 'REDV']
        return sum(self._count_motif_pattern(sequence, motif) for motif in adhesin_motifs)
    
    def _calculate_invasin_score(self, sequence):
        """Calculate invasin protein score"""
        invasin_motifs = ['NPXY', 'YxxΦ', 'FXXF']
        return sum(self._count_motif_pattern(sequence, motif) for motif in invasin_motifs)
    
    def _calculate_protease_score(self, sequence):
        """Calculate protease activity score"""
        protease_motifs = ['HExxH', 'DxxD', 'SxxS']
        catalytic_triads = sequence.count('HDS') + sequence.count('SHD')
        return sum(self._count_motif_pattern(sequence, motif) for motif in protease_motifs) + catalytic_triads
    
    def _calculate_lipase_score(self, sequence):
        """Calculate lipase activity score"""
        lipase_motifs = ['GxSxG', 'GDSAG']
        return sum(self._count_motif_pattern(sequence, motif) for motif in lipase_motifs)
    
    def _calculate_heat_shock_score(self, sequence):
        """Calculate heat shock protein signatures"""
        # Heat shock proteins are rich in certain amino acids
        hsp_residues = 'GPKR'  # Common in heat shock proteins
        return sum(1 for aa in sequence if aa in hsp_residues) / len(sequence)
    
    def _calculate_oxidative_stress_score(self, sequence):
        """Calculate oxidative stress response score"""
        # Proteins involved in oxidative stress response
        antioxidant_motifs = ['CxxC', 'CXXC', 'HxxH']
        return sum(self._count_motif_pattern(sequence, motif) for motif in antioxidant_motifs)
    
    def _calculate_iron_acquisition_score(self, sequence):
        """Calculate iron acquisition system score"""
        iron_motifs = ['HxxH', 'HxxxH', 'ExxE']
        return sum(self._count_motif_pattern(sequence, motif) for motif in iron_motifs)
    
    def _calculate_nutrient_scavenging_score(self, sequence):
        """Calculate nutrient scavenging score"""
        # Transport and binding proteins for nutrients
        transport_motifs = ['DxxG', 'NxxD', 'KxxK']
        return sum(self._count_motif_pattern(sequence, motif) for motif in transport_motifs)
    
    def _calculate_beta_lactamase_score(self, sequence):
        """Calculate beta-lactamase signature score"""
        blac_motifs = ['SxxK', 'SDN', 'KTG']
        return sum(self._count_motif_pattern(sequence, motif) for motif in blac_motifs)
    
    def _calculate_efflux_pump_score(self, sequence):
        """Calculate efflux pump signature score"""
        # Efflux pumps are typically membrane proteins with high hydrophobic content
        hydrophobic_content = sum(1 for aa in sequence if aa in 'AILMFPWV') / len(sequence)
        return 1.0 if hydrophobic_content > 0.6 else hydrophobic_content
    
    def _get_default_features(self):
        """Return default feature values for failed extractions"""
        features = {}
        
        # Default amino acid composition
        for aa in 'ACDEFGHIKLMNPQRSTVWY':
            features[f'aa_percent_{aa}'] = 0.05  # Equal distribution
        
        # Default traditional features
        features.update({
            'length_log': 6.0, 'length_sqrt': 25.0, 'length_normalized': 0.0,
            'hydrophobic_fraction': 0.4, 'charged_fraction': 0.2,
            'polar_fraction': 0.2, 'aromatic_fraction': 0.1,
            'small_fraction': 0.3, 'tiny_fraction': 0.2,
            'molecular_weight_norm': 110.0, 'aromaticity': 0.1,
            'gravy': 0.0, 'isoelectric_point': 7.0,
            'charge_density': 0.0, 'charge_distribution': 0.1,
            'net_charge': 0.0,
            'helix_fraction': 0.3, 'turn_fraction': 0.3,
            'sheet_fraction': 0.4, 'structure_ratio': 1.0
        })
        
        # Default pathogen features
        for motif_class in ['type3_secretion', 'toxin_binding', 'cell_adhesion', 
                           'membrane_disruption', 'immune_evasion', 'host_invasion']:
            features[f'{motif_class}_motifs'] = 0.0
        
        features.update({
            'virulence_score': 0.0, 'toxin_like_regions': 0.0,
            'host_binding_potential': 0.0, 'immune_recognition_sites': 0.0,
            'amr_potential': 0.0, 'biofilm_potential': 0.0
        })
        
        # Default evolutionary features
        features.update({
            'conservation_mean': 0.5, 'conservation_std': 0.2,
            'conservation_max': 0.8, 'conservation_entropy': 2.0,
            'selection_pressure': 0.5, 'codon_adaptation_index': 0.5,
            'gc_content': 0.5, 'hgt_signal': 0.1,
            'mobile_element_association': 0.0, 'taxonomic_breadth': 0.5,
            'ancient_origin': 0.0
        })
        
        # Default structural features
        features.update({
            'disorder_fraction': 0.3, 'disorder_longest_region': 0.1,
            'disorder_clusters': 2.0, 'tm_helices_count': 0.0,
            'tm_fraction': 0.0, 'tm_topology_score': 0.0,
            'signal_peptide_score': 0.0, 'signal_peptide_type': 0.0,
            'coiled_coil_score': 0.0, 'beta_strand_content': 0.2,
            'turn_propensity': 0.3
        })
        
        # Default length-optimized features
        features.update({
            'n_terminal_hydrophobic': 0.3, 'n_terminal_basic': 0.1,
            'n_terminal_signal_score': 0.0, 'c_terminal_charged': 0.2,
            'c_terminal_polar': 0.3, 'core_complexity': 0.5,
            'core_conservation': 0.3, 'optimal_length_score': 0.5
        })
        
        # Motif densities
        for motif in ['RGD', 'KKK', 'DDD', 'LLL', 'FFF']:
            features[f'{motif}_density'] = 0.0
        
        # Default advanced pathogen signatures
        features.update({
            'sec_pathway_score': 0.0, 'tat_pathway_score': 0.0,
            'adhesin_score': 0.0, 'invasin_score': 0.0,
            'protease_score': 0.0, 'lipase_score': 0.0,
            'heat_shock_score': 0.3, 'oxidative_stress_score': 0.0,
            'iron_acquisition_score': 0.0, 'nutrient_scavenging_score': 0.0,
            'beta_lactamase_score': 0.0, 'efflux_pump_score': 0.4
        })
        
        return features
    
    def _get_default_traditional_features(self):
        """Default traditional features for error cases"""
        features = {}
        for aa in 'ACDEFGHIKLMNPQRSTVWY':
            features[f'aa_percent_{aa}'] = 0.05
        
        features.update({
            'length_log': 6.0, 'length_sqrt': 25.0, 'length_normalized': 0.0,
            'hydrophobic_fraction': 0.4, 'charged_fraction': 0.2,
            'polar_fraction': 0.2, 'aromatic_fraction': 0.1,
            'small_fraction': 0.3, 'tiny_fraction': 0.2,
            'molecular_weight_norm': 110.0, 'aromaticity': 0.1,
            'gravy': 0.0, 'isoelectric_point': 7.0,
            'charge_density': 0.0, 'charge_distribution': 0.1,
            'net_charge': 0.0,
            'helix_fraction': 0.3, 'turn_fraction': 0.3,
            'sheet_fraction': 0.4, 'structure_ratio': 1.0
        })
        
        return features


# Utility functions for saving and loading
def save_feature_extractor(feature_extractor, output_path):
    """Save the trained feature extractor with metadata"""
    try:
        output_path = Path(output_path)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Save the feature extractor
        extractor_file = output_path / "metaquest_feature_extractor.pkl"
        joblib.dump(feature_extractor, extractor_file)
        
        # Save metadata
        metadata = {
            'version': '2.1.0_production_complete',
            'creation_date': pd.Timestamp.now().isoformat(),
            'esm2_enabled': False,
            'feature_count': len(feature_extractor.get_feature_names_out()),
            'feature_categories': [
                'amino_acid_composition', 'physicochemical_properties', 'structural_features',
                'pathogen_motifs', 'virulence_factors', 'evolutionary_signatures', 
                'length_optimized_features', 'advanced_pathogen_signatures'
            ],
            'production_ready': True,
            'voting_ensemble_compatible': True,
            'environment_agnostic': True
        }
        
        metadata_file = output_path / "extractor_metadata.json"
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        print(f"✅ Feature extractor saved to: {extractor_file}")
        print(f"📋 Metadata saved to: {metadata_file}")
        
        return {
            'feature_extractor': str(extractor_file),
            'metadata': str(metadata_file)
        }
        
    except Exception as e:
        print(f"❌ Error saving feature extractor: {e}")
        raise


def load_feature_extractor(extractor_path):
    """Load a saved feature extractor"""
    try:
        feature_extractor = joblib.load(extractor_path)
        print(f"✅ Feature extractor loaded from: {extractor_path}")
        return feature_extractor
    except Exception as e:
        print(f"❌ Error loading feature extractor: {e}")
        raise
