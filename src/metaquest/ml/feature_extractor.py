#!/usr/bin/env python3

"""
MetaQuest Protein Feature Extractor - Production Version v2.1
Updated to match training pipeline with enhanced feature extraction
Environment-agnostic pathogen detection without ESM2/torch dependencies
Compatible with trained voting ensemble model
"""

import pandas as pd
import numpy as np
import joblib
import re
from collections import Counter
import json
from pathlib import Path
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

# Bio-related imports
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import molecular_weight

# ML imports
from sklearn.base import BaseEstimator, TransformerMixin

# Local imports
from ..io.output_formatter import get_formatter

formatter = get_formatter()


class MetaQuestProteinFeatureExtractor(BaseEstimator, TransformerMixin):
    """
    Production MetaQuest protein feature extraction v2.1
    UPDATED: Matches training pipeline feature extraction exactly
    Compatible with trained voting ensemble model
    No ESM2/torch dependencies
    """
    
    def __init__(self):
        self.is_fitted = False
        
        # Amino acid groups (matching training pipeline)
        self.aa_groups = {
            'hydrophobic': 'AILMFPWV',
            'charged': 'DEKRH',
            'polar': 'NQST',
            'aromatic': 'FWY',
            'aliphatic': 'ILV',
            'small': 'AGST',
            'positive': 'KRH',
            'negative': 'DE',
            'tiny': 'AGS',
            'proline': 'P'
        }
        
        # Hydrophobicity scale (Kyte-Doolittle)
        self.kd_scale = {
            'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
            'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
            'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
            'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
        }
        
        # Compile pathogen-specific motifs
        self._compile_pathogen_motifs()
        
    def _compile_pathogen_motifs(self):
        """Compile comprehensive pathogen-specific regex patterns (from training)."""
        
        # Type III/IV/VI Secretion System signals
        self.secretion_motifs = [
            re.compile(r'[LMFI]RL..[LMFI]'),  # T3SS canonical
            re.compile(r'Y.C'),  # T3SS binding
            re.compile(r'L.RC.'),  # T4SS-like
            re.compile(r'RR[ILMV]'),  # Tat signal
            re.compile(r'[KR]{2,}'),  # Basic motif clusters
        ]
        
        # Toxin signatures
        self.toxin_motifs = [
            re.compile(r'C.{2,4}C.{3}C.{5,6}C'),  # Cysteine-rich
            re.compile(r'G.G..G'),  # GxGxxG
            re.compile(r'W.W'),  # Aromatic pairs
            re.compile(r'HE[ILMV]{2}H'),  # Metalloprotease
            re.compile(r'[DE]{3,}'),  # Acidic stretches
        ]
        
        # AMR-specific patterns
        self.amr_motifs = [
            re.compile(r'S[ILMV]{2}K'),  # Beta-lactamase
            re.compile(r'[ST]..K'),  # Serine active site
            re.compile(r'YGN'),  # Penicillin-binding
            re.compile(r'[KR].G[KR]'),  # Aminoglycoside binding
            re.compile(r'Q.{2}E'),  # Efflux pump
        ]
        
        # Adhesin motifs
        self.adhesion_motifs = [
            re.compile(r'RGD'),  # Integrin-binding
            re.compile(r'[ILMV]{3,}'),  # Hydrophobic stretches
            re.compile(r'[NQST]{4,}'),  # Polar-rich
        ]
        
        # Signal peptides
        self.signal_peptide = re.compile(r'^[MKRFLI]{1,3}[AILMFVWP]{5,20}[AGSTN]')
        
        # Transmembrane helix
        self.tm_helix = re.compile(r'[AILMFVWP]{18,25}')
        
        # Catalytic motifs
        self.catalytic_motifs = [
            re.compile(r'[HN].[DEST]'),  # Catalytic dyad
            re.compile(r'G.{1,2}[ST].{2,4}G'),  # Nucleotide binding
        ]
        
    def fit(self, X=None, y=None):
        """Fit the feature extractor (sklearn compatibility)"""
        print("🔧 Fitting MetaQuest Feature Extractor v2.1...")
        self.is_fitted = True
        print("Feature extractor fitted successfully!")
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
    
    def extract_features(self, sequence: str) -> dict:
        """
        Main feature extraction method matching training pipeline exactly.
        UPDATED: Now matches preprocessing.py feature extraction
        """
        if not isinstance(sequence, str) or len(sequence) < 20:
            return self._get_default_features()
        
        # Clean sequence (matching training)
        cleaned_seq = sequence.replace('X', 'A').replace('B', 'N').replace('Z', 'Q')
        cleaned_seq = cleaned_seq.replace('J', 'L').replace('U', 'C').replace('O', 'K')
        cleaned_seq = sequence.strip('*').upper()
        cleaned_seq = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', 'A', cleaned_seq)
        if len(cleaned_seq) < 20:
            return self._get_default_features()

        # Check quality
        ambiguous_count = sum(1 for aa in sequence if aa in 'XBZJUO')
        if ambiguous_count > len(sequence) * 0.3:
            return self._get_default_features()
            
        features = {}
        
        try:
            analyzer = ProteinAnalysis(cleaned_seq)
            
            # 1. Physicochemical features (matching training)
            features.update(self._physicochemical_features(cleaned_seq, analyzer))
            
            # 2. Compositional features (matching training)
            features.update(self._compositional_features(cleaned_seq, analyzer))
            
            # 3. Pathogen motif features (matching training)
            features.update(self._pathogen_motif_features(cleaned_seq))
            
            # 4. Sliding window features (matching training)
            features.update(self._sliding_window_features(cleaned_seq))
            
            # 5. Structural features (matching training)
            features.update(self._structural_features(cleaned_seq, analyzer))
            
            # 6. Complexity features (matching training)
            features.update(self._complexity_features(cleaned_seq))
            
            # 7. Positional features (matching training)
            features.update(self._positional_features(cleaned_seq))
            
        except Exception as e:
            self.formatter.debug(f"Feature extraction failed for a sequence: {e}")
            return self._get_default_features()
            
        return features
    
    def _physicochemical_features(self, sequence, analyzer):
        """Enhanced physicochemical properties (from training pipeline)."""
        seq_len = len(sequence)
        charge_ph7 = sum(1 for aa in sequence if aa in 'KR') - sum(1 for aa in sequence if aa in 'DE')
        
        features = {
            'length': seq_len,
            'molecular_weight': analyzer.molecular_weight() / 1000,
            'isoelectric_point': analyzer.isoelectric_point(),
            'gravy': analyzer.gravy(),
            'aromaticity': analyzer.aromaticity(),
            'instability_index': analyzer.instability_index(),
            'helix_fraction': analyzer.secondary_structure_fraction()[0],
            'turn_fraction': analyzer.secondary_structure_fraction()[1],
            'sheet_fraction': analyzer.secondary_structure_fraction()[2],
            'charge_ph7': charge_ph7 / seq_len,
            'abs_charge_ph7': abs(charge_ph7) / seq_len,
            'flexibility_index': analyzer.instability_index() / 100,
        }
        
        return features
    
    def _compositional_features(self, sequence, analyzer):
        """Comprehensive compositional analysis (from training pipeline)."""
        features = {}
        seq_len = len(sequence)
        aa_percent = analyzer.get_amino_acids_percent()
        
        # Amino acid group fractions
        for aa_group, residues in self.aa_groups.items():
            features[f'{aa_group}_fraction'] = sum(
                aa_percent.get(residue, 0) for residue in residues
            )
        
        # Key individual amino acids
        for aa in ['C', 'P', 'G', 'W']:
            features[f'{aa}_count'] = aa_percent.get(aa, 0)
        
        # Extended dipeptides (30 most informative)
        dipeptides = [sequence[i:i+2] for i in range(seq_len - 1)]
        dp_counts = Counter(dipeptides)
        
        informative_dp = [
            'LA', 'AG', 'GV', 'VA', 'LL', 'AA', 'GL', 'LV', 'GG', 'LG',
            'AL', 'AV', 'GA', 'SG', 'SA', 'VL', 'IG', 'IV', 'GI', 'VG',
            'RG', 'GR', 'KG', 'GK', 'DE', 'ED', 'RR', 'KK', 'WG', 'YG'
        ]
        
        for dp in informative_dp:
            features[f'dp_{dp}'] = dp_counts.get(dp, 0) / max(seq_len - 1, 1)
        
        # Tripeptides (pathogen-enriched)
        tripeptides = [sequence[i:i+3] for i in range(seq_len - 2)]
        tp_counts = Counter(tripeptides)
        
        key_tripeptides = ['RGD', 'YGN', 'GGG', 'AAA', 'LLL', 'RRR', 'DED', 'GNG']
        for tp in key_tripeptides:
            features[f'tp_{tp}'] = tp_counts.get(tp, 0) / max(seq_len - 2, 1)
        
        return features
    
    def _pathogen_motif_features(self, sequence):
        """Enhanced pathogen-specific motif detection (from training pipeline)."""
        seq_len = len(sequence)
        
        features = {
            'secretion_motif_count': sum(
                len(motif.findall(sequence)) for motif in self.secretion_motifs
            ),
            'secretion_motif_density': sum(
                len(motif.findall(sequence)) for motif in self.secretion_motifs
            ) * 1000 / seq_len,
            
            'toxin_motif_count': sum(
                len(motif.findall(sequence)) for motif in self.toxin_motifs
            ),
            'toxin_motif_density': sum(
                len(motif.findall(sequence)) for motif in self.toxin_motifs
            ) * 1000 / seq_len,
            
            'amr_motif_count': sum(
                len(motif.findall(sequence)) for motif in self.amr_motifs
            ),
            'amr_motif_density': sum(
                len(motif.findall(sequence)) for motif in self.amr_motifs
            ) * 1000 / seq_len,
            
            'adhesion_motif_count': sum(
                len(motif.findall(sequence)) for motif in self.adhesion_motifs
            ),
            
            'has_signal_peptide': 1 if self.signal_peptide.search(sequence[:30]) else 0,
            'tm_helix_count': len(self.tm_helix.findall(sequence)),
            
            'catalytic_motif_count': sum(
                len(motif.findall(sequence)) for motif in self.catalytic_motifs
            ),
        }
        
        return features
    
    def _sliding_window_features(self, sequence, window_size=20):
        """Localized property analysis (from training pipeline)."""
        features = {}
        seq_len = len(sequence)
        
        if seq_len < window_size:
            return {
                'hydrophobicity_max_local': 0,
                'hydrophobicity_std_local': 0,
                'charge_max_local': 0,
                'charge_min_local': 0,
                'hydrophobic_cluster_score': 0,
                'charge_cluster_score': 0,
            }
        
        window_hydrophobicities = []
        window_charges = []
        
        for i in range(seq_len - window_size + 1):
            window = sequence[i:i+window_size]
            
            hydro = sum(self.kd_scale.get(aa, 0) for aa in window) / window_size
            window_hydrophobicities.append(hydro)
            
            charge = sum(1 for aa in window if aa in 'KRH') - sum(1 for aa in window if aa in 'DE')
            window_charges.append(charge)
        
        features.update({
            'hydrophobicity_max_local': np.max(window_hydrophobicities),
            'hydrophobicity_std_local': np.std(window_hydrophobicities),
            'charge_max_local': np.max(window_charges),
            'charge_min_local': np.min(window_charges),
            'hydrophobic_cluster_score': np.max(window_hydrophobicities) * np.std(window_hydrophobicities),
            'charge_cluster_score': (np.max(np.abs(window_charges)) * np.std(window_charges)) if len(window_charges) > 0 else 0,
        })
        
        return features
    
    def _structural_features(self, sequence, analyzer):
        """Structural propensity features (from training pipeline)."""
        seq_len = len(sequence)
        beta_turn_residues = sum(1 for aa in sequence if aa in 'NPGS')
        helix_breaker = sum(1 for aa in sequence if aa in 'PG')
        cys_count = sequence.count('C')
        
        features = {
            'beta_turn_propensity': beta_turn_residues / seq_len,
            'helix_breaker_fraction': helix_breaker / seq_len,
            'cysteine_count': cys_count,
            'disulfide_potential': cys_count * (cys_count - 1) / 2 if cys_count > 1 else 0,
            'proline_content': sequence.count('P') / seq_len,
        }
        
        return features
    
    def _complexity_features(self, sequence):
        """Sequence complexity and disorder metrics (from training pipeline)."""
        seq_len = len(sequence)
        aa_counts = Counter(sequence)
        
        # Shannon entropy
        entropy = -sum(
            (count / seq_len) * np.log2(count / seq_len)
            for count in aa_counts.values()
        )
        
        # Max repeat length
        max_repeat = max((len(list(g)) for _, g in __import__('itertools').groupby(sequence)), default=0)
        
        # Disorder propensity
        disorder_aa = 'PQSAEKRGTD'
        disorder_score = sum(1 for aa in sequence if aa in disorder_aa) / seq_len
        
        features = {
            'shannon_entropy': entropy,
            'max_repeat_length': max_repeat,
            'repeat_fraction': max_repeat / seq_len,
            'disorder_propensity': disorder_score,
            'unique_aa_fraction': len(aa_counts) / 20,
        }
        
        return features
    
    def _positional_features(self, sequence):
        """Position-specific features (from training pipeline)."""
        seq_len = len(sequence)
        n_term_len = min(30, seq_len // 3)
        c_term_len = min(30, seq_len // 3)
        
        n_term = sequence[:n_term_len]
        c_term = sequence[-c_term_len:]
        
        features = {
            'n_term_hydrophobic': sum(1 for aa in n_term if aa in 'AILMFVWP') / n_term_len if n_term_len > 0 else 0,
            'n_term_positive': sum(1 for aa in n_term if aa in 'KR') / n_term_len if n_term_len > 0 else 0,
            'n_term_small': sum(1 for aa in n_term if aa in 'AGST') / n_term_len if n_term_len > 0 else 0,
            'c_term_hydrophobic': sum(1 for aa in c_term if aa in 'AILMFVWP') / c_term_len if c_term_len > 0 else 0,
            'c_term_charged': sum(1 for aa in c_term if aa in 'DEKRH') / c_term_len if c_term_len > 0 else 0,
            'starts_with_M': 1 if sequence[0] == 'M' else 0,
            'ends_with_polar': 1 if sequence[-1] in 'NQST' else 0,
        }
        
        return features
    
    def _get_default_features(self):
        """Zero-filled features for invalid sequences."""
        dummy_seq = "A" * 100
        dummy_features = self.extract_features(dummy_seq)
        return {key: 0.0 for key in dummy_features.keys()}
    
    def get_feature_names_out(self, input_features=None):
        """Get output feature names (sklearn compatibility)."""
        if not self.is_fitted:
            self.fit()
        
        # Generate feature names matching training pipeline order
        dummy_seq = "A" * 100
        features = self.extract_features(dummy_seq)
        return list(features.keys())


# Utility functions
def save_feature_extractor(feature_extractor, output_path):
    """Save the trained feature extractor with metadata"""
    try:
        output_path = Path(output_path)
        output_path.mkdir(parents=True, exist_ok=True)
        
        extractor_file = output_path / "metaquest_feature_extractor.pkl"
        joblib.dump(feature_extractor, extractor_file)
        
        metadata = {
            'version': '2.1.0_production_updated',
            'creation_date': pd.Timestamp.now().isoformat(),
            'matches_training_pipeline': True,
            'feature_count': len(feature_extractor.get_feature_names_out()),
            'production_ready': True,
            'voting_ensemble_compatible': True,
        }
        
        metadata_file = output_path / "extractor_metadata.json"
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        print(f"Feature extractor saved to: {extractor_file}")
        return {'feature_extractor': str(extractor_file), 'metadata': str(metadata_file)}
        
    except Exception as e:
        print(f"❌ Error saving feature extractor: {e}")
        raise


def load_feature_extractor(extractor_path):
    """Load a saved feature extractor"""
    try:
        feature_extractor = joblib.load(extractor_path)
        print(f"Feature extractor loaded from: {extractor_path}")
        return feature_extractor
    except Exception as e:
        print(f"❌ Error loading feature extractor: {e}")
        raise