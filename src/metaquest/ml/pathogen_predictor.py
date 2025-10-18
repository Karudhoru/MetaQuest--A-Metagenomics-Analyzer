#!/usr/bin/env python3

"""
MetaQuest Pathogen Predictor - Production Version v2.1
FIXED: Proper OutputFormatter usage and progress bars
Enhanced feature alignment with training pipeline
Integrated with trained voting ensemble model (GPU-accelerated training)
Compatible with Prokka output processing
No ESM2/torch dependencies
"""

import json
import logging
import joblib
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from datetime import datetime
from collections import Counter
import warnings
from Bio import SeqIO

# Import our updated feature extractor
from .feature_extractor import MetaQuestProteinFeatureExtractor
from ..config import MODEL_ARTIFACTS_DIR, ML_CONFIG, PATHOGEN_COLORS
from ..io.output_formatter import OutputFormatter

# Disable warnings
import os
os.environ['SKLEARN_ENABLE_METADATA_ROUTING'] = 'False'
warnings.filterwarnings('ignore')
logger = logging.getLogger(__name__)


class PathogenPredictor:
    """
    Production pathogen predictor v2.1 with robust feature alignment.
    FIXED: Uses OutputFormatter correctly with proper progress bars
    """
    def __init__(self, verbosity: str = 'standard'):
        """
        Initialize the predictor
        
        Args:
            verbosity: One of 'quiet', 'minimal', 'standard', 'debug'
        """
        self.verbosity = verbosity
        self.formatter = OutputFormatter(verbosity=verbosity)
        self.model_dir = Path(MODEL_ARTIFACTS_DIR)
        self.scaler = None
        self.model = None
        self.feature_selector = None
        self.all_feature_names = None  # Complete feature list from training
        self.selected_feature_names = None  # Features after selection
        
        self.formatter.operation("Initializing MetaQuest Feature Extractor v2.1")
        self.feature_extractor = MetaQuestProteinFeatureExtractor()
        self.feature_extractor.fit()
        self.formatter.success("Feature extractor ready")
        
        self.ml_config = ML_CONFIG.copy()
        self._load_model_artifacts()
    
    def _load_model_artifacts(self):
        """Load trained model artifacts with enhanced error handling."""
        try:
            self.formatter.operation("Loading trained model artifacts")
            
            # Load scaler
            self.scaler = joblib.load(self.model_dir / "scaler.pkl")
            self.formatter.info("Loaded scaler", indent=2)
            
            # Load model
            self.model = joblib.load(self.model_dir / "best_model.pkl")
            self.formatter.info(f"Loaded model: {type(self.model).__name__}", indent=2)
            
            # Load feature selector
            self.feature_selector = joblib.load(self.model_dir / "feature_selector.pkl")
            self.formatter.info("Loaded feature selector", indent=2)
            
            # Load ALL feature names (CRITICAL for alignment)
            all_feature_path = self.model_dir / "all_feature_names.json"
            if all_feature_path.exists():
                with open(all_feature_path, 'r') as f:
                    self.all_feature_names = json.load(f)
                self.formatter.info(f"Loaded all feature names ({len(self.all_feature_names)} features)", indent=2)
            else:
                raise FileNotFoundError(f"Critical file missing: {all_feature_path}")
            
            # Load selected feature names (for reference)
            selected_feature_path = self.model_dir / "selected_feature_names.json"
            if selected_feature_path.exists():
                with open(selected_feature_path, 'r') as f:
                    self.selected_feature_names = json.load(f)
                self.formatter.info(f"Loaded selected feature names ({len(self.selected_feature_names)} features)", indent=2)
            
            # Load metadata
            metadata_path = self.model_dir / "model_metadata.json"
            if metadata_path.exists():
                with open(metadata_path, 'r') as f:
                    self.metadata = json.load(f)
                training_date = self.metadata.get('training_date', 'unknown')
                self.formatter.info(f"Loaded model metadata (trained: {training_date})", indent=2)
            
            self.formatter.success("Complete ML pipeline loaded successfully", style='bold')
            
        except Exception as e:
            logger.error(f"Error loading model artifacts from {self.model_dir}: {e}")
            self.formatter.error(
                f"Failed to load required ML artifacts: {e}",
                solutions=[
                    "Verify model files exist in model_artifacts directory",
                    "Re-run training pipeline if models are missing",
                    "Check file permissions"
                ]
            )
            raise RuntimeError(f"Failed to load required ML artifacts: {e}")
    
    def check_capabilities(self) -> Dict[str, bool]:
        """Check and return current system capabilities."""
        capabilities = {
            'feature_extractor_available': True,
            'ml_pipeline_complete': all([self.scaler, self.model, self.feature_selector]),
            'all_feature_names_available': self.all_feature_names is not None,
            'selected_feature_names_available': self.selected_feature_names is not None,
            'voting_ensemble_loaded': 'voting' in str(type(self.model)).lower(),
            'stacking_ensemble_loaded': 'stacking' in str(type(self.model)).lower(),
            'model_type': type(self.model).__name__ if self.model else 'None',
            'feature_count': len(self.all_feature_names) if self.all_feature_names else 0,
            'selected_feature_count': len(self.selected_feature_names) if self.selected_feature_names else 0,
        }
        return capabilities
    
    def predict_from_prokka(self, prokka_dir: Path) -> Tuple[List[Dict], Dict]:
        """
        Main prediction function. Processes a Prokka output directory and returns
        both the detailed predictions and a summary.
        
        Returns:
            tuple: (results_list, summary_dict)
        """
        protein_files = list(prokka_dir.glob("*.faa"))
        
        if not protein_files:
            self.formatter.warning(f"No protein FASTA files (*.faa) found in {prokka_dir}")
            logger.warning(f"No protein FASTA files (*.faa) found in {prokka_dir}. Skipping ML prediction.")
            return [], {}
        
        protein_file = protein_files[0]
        self.formatter.operation(f"Processing Prokka proteins from: {protein_file.name}")
        sequences = list(SeqIO.parse(protein_file, "fasta"))
        total_sequences = len(sequences)
        self.formatter.info(f"Found {total_sequences:,} protein sequences to analyze")        

        if not sequences:
            self.formatter.warning(f"Protein file {protein_file.name} is empty")
            logger.warning(f"Protein file {protein_file} is empty. Skipping ML prediction.")
            return [], {}
        
        # 1. Extract features for ALL sequences first (this is the main loop)
        self.formatter.operation("Extracting features from protein sequences")
        feature_list = []
        
        # Use progress bar for feature extraction
        with self.formatter.progress_bar(total=total_sequences, desc="      Extracting Features", unit="proteins") as pbar:
            for rec in sequences:
                # The feature extractor is now silent and just returns data
                features = self.feature_extractor.extract_features(str(rec.seq))
                feature_list.append(features)
                pbar.update(1)
        
        self.formatter.success(f"Extracted features from {len(feature_list):,} proteins")
        feature_df = pd.DataFrame(feature_list)
        
        if feature_df.empty:
            self.formatter.warning("Feature extraction yielded no data")
            logger.warning("Feature extraction yielded no data. Aborting prediction.")
            return [], {}

        # 2. Perform all ML steps on the ENTIRE DataFrame at once
        self.formatter.operation("Preparing features for ML prediction")
        
        # Remove metadata columns before alignment
        metadata_cols = ['sequence_id', 'sequence_length', 'extraction_method', 
                       'feature_count', 'timestamp', 'batch_number']
        feature_df_clean = feature_df.drop(columns=[col for col in metadata_cols if col in feature_df.columns], errors='ignore')
        
        # Align features to match the training set (CRITICAL FIX)
        self.formatter.info("Aligning features to training set", indent=2)
        aligned_df = feature_df_clean.reindex(columns=self.all_feature_names, fill_value=0.0)
        self.formatter.info(f"Aligned: {aligned_df.shape[1]} features matching training", indent=2)

        # Apply transformations (on the whole matrix)
        self.formatter.info("Selecting features", indent=2)
        features_selected = self.feature_selector.transform(aligned_df)
        self.formatter.info(f"Selected: {features_selected.shape[1]} features", indent=2)
        
        self.formatter.info("Scaling features", indent=2)
        features_scaled = self.scaler.transform(features_selected)
        self.formatter.success(f"Features prepared: {features_scaled.shape}")

        # Make predictions (on the whole matrix)
        self.formatter.operation("Running ensemble model predictions")
        predictions = self.model.predict(features_scaled)
        probabilities = self.model.predict_proba(features_scaled)
        self.formatter.success(f"Generated predictions for {len(predictions)} sequences")

        # 3. Assemble results (this is a fast final loop)
        self.formatter.operation("Assembling prediction results")
        results = []
        confidence_threshold = self.ml_config.get('confidence_threshold', 0.7)
        
        for i, (record, pred, probs) in enumerate(zip(sequences, predictions, probabilities)):
            confidence = float(max(probs))
            pathogenic_prob = float(probs[1]) if len(probs) > 1 else (confidence if pred == 1 else (1 - confidence))
            
            result = {
                'sequence_id': record.id,
                'sequence_type': 'protein',
                'prediction': 'Pathogenic' if pred == 1 else 'Non-pathogenic',
                'confidence': confidence,
                'pathogenic_probability': pathogenic_prob,
                'high_confidence': confidence >= confidence_threshold,
                'confidence_threshold': confidence_threshold,
                'sequence_length': len(record.seq),
                'method': 'ML_Ensemble_MetaQuest_v2.1',
                'extraction_method': 'MetaQuest_v2.1_Production',
                'model_type': type(self.model).__name__,
                'features_used': aligned_df.shape[1],
                'features_selected': features_selected.shape[1],
                'features_scaled': features_scaled.shape[1],
                'timestamp': datetime.now().isoformat()
            }
            results.append(result)

        self.formatter.success(f"Assembled {len(results)} prediction results")

        # 4. Generate summary
        self.formatter.operation("Generating prediction summary")
        
        if not results:
            self.formatter.warning("No results to summarize")
            logger.warning("No results to summarize")
            return [], {}
        
        summary = self.get_prediction_summary(results)

        logger.info(f"Completed processing: {len(results)} predictions generated.")
        self.formatter.success(f"ML prediction complete: {len(results)} sequences analyzed")
        return results, summary
    
    def save_results(self, results: List[Dict], output_file: Path, format_type: str = "csv"):
        """Save prediction results with comprehensive metadata."""
        if not results:
            self.formatter.warning("No results to save")
            logger.warning("No results to save")
            return
        
        try:
            if format_type.lower() == "json":
                output_data = {
                    'metadata': {
                        'total_predictions': len(results),
                        'metaquest_version': '2.1.0',
                        'feature_alignment_method': 'training_feature_names',
                        'feature_count': len(self.all_feature_names) if self.all_feature_names else 0,
                        'selected_features': len(self.selected_feature_names) if self.selected_feature_names else 0,
                        'model_type': type(self.model).__name__ if self.model else 'Unknown',
                        'generation_timestamp': datetime.now().isoformat(),
                    },
                    'predictions': results
                }
                with open(output_file, 'w') as f:
                    json.dump(output_data, f, indent=2)
            else:
                df = pd.DataFrame(results)
                df.to_csv(output_file, index=False)
            
            self.formatter.success(f"Results saved: {output_file.name}")
            logger.info(f"Results saved to: {output_file}")
            
        except Exception as e:
            logger.error(f"Error saving results: {e}")
            self.formatter.error(f"Error saving results: {e}")
            raise
    
    def get_prediction_summary(self, results: List[Dict]) -> Dict:
        """Generate comprehensive prediction summary with enhanced statistics."""
        if not results:
            return {}
        
        # Validate that results is a list of dicts
        if not isinstance(results, list):
            logger.error(f"Invalid results type: expected list, got {type(results)}")
            self.formatter.error(f"Invalid results type: expected list, got {type(results)}")
            return {}
        
        if not all(isinstance(r, dict) for r in results):
            logger.error(f"Invalid results structure: not all items are dicts")
            self.formatter.error(f"Invalid results structure: not all items are dicts")
            return {}
        
        pathogenic_sequences = [r for r in results if r.get('prediction') == 'Pathogenic']
        high_confidence_predictions = [r for r in results if r.get('high_confidence', False)]
        
        # Calculate statistics
        confidence_scores = [r['confidence'] for r in results if 'confidence' in r]
        pathogenic_scores = [r.get('pathogenic_probability', 0) for r in pathogenic_sequences]
        
        summary = {
            'total_sequences_analyzed': len(results),
            'pathogenic_predictions': len(pathogenic_sequences),
            'pathogenic_percentage': (len(pathogenic_sequences) / len(results)) * 100 if results else 0,
            'non_pathogenic_predictions': len(results) - len(pathogenic_sequences),
            'non_pathogenic_percentage': ((len(results) - len(pathogenic_sequences)) / len(results)) * 100 if results else 0,
            
            'high_confidence_predictions': len(high_confidence_predictions),
            'high_confidence_percentage': (len(high_confidence_predictions) / len(results)) * 100 if results else 0,
            'low_confidence_predictions': len(results) - len(high_confidence_predictions),
            
            'confidence_threshold': self.ml_config.get('confidence_threshold', 0.7),
            'average_confidence': float(np.mean(confidence_scores)) if confidence_scores else 0,
            'median_confidence': float(np.median(confidence_scores)) if confidence_scores else 0,
            'min_confidence': float(np.min(confidence_scores)) if confidence_scores else 0,
            'max_confidence': float(np.max(confidence_scores)) if confidence_scores else 0,
            'std_confidence': float(np.std(confidence_scores)) if confidence_scores else 0,
            
            'average_pathogenic_score': float(np.mean(pathogenic_scores)) if pathogenic_scores else 0,
            'median_pathogenic_score': float(np.median(pathogenic_scores)) if pathogenic_scores else 0,
            
            'methods_distribution': dict(Counter(r.get('method', 'Unknown') for r in results)),
            'model_type': type(self.model).__name__ if self.model else 'Unknown',
            'ml_config_used': self.ml_config,
            'model_artifacts_dir': str(self.model_dir),
            
            'voting_ensemble_used': 'voting' in str(type(self.model)).lower(),
            'stacking_ensemble_used': 'stacking' in str(type(self.model)).lower(),
            'feature_alignment_used': 'training_feature_names' if self.all_feature_names else 'extractor_defaults',
            'features_total': len(self.all_feature_names) if self.all_feature_names else 0,
            'features_selected': len(self.selected_feature_names) if self.selected_feature_names else 0,
            
            'system_capabilities': self.check_capabilities(),
            'pathogen_colors': PATHOGEN_COLORS,
            
            'top_pathogenic_sequences': sorted(
                pathogenic_sequences, 
                key=lambda x: x['confidence'], 
                reverse=True
            )[:10],
            
            'top_high_confidence_sequences': sorted(
                high_confidence_predictions, 
                key=lambda x: x['confidence'], 
                reverse=True
            )[:10],
            
            'confidence_distribution': {
                'very_high': len([r for r in results if r.get('confidence', 0) >= 0.9]),
                'high': len([r for r in results if 0.7 <= r.get('confidence', 0) < 0.9]),
                'medium': len([r for r in results if 0.5 <= r.get('confidence', 0) < 0.7]),
                'low': len([r for r in results if r.get('confidence', 0) < 0.5]),
            }
        }
        
        return summary


# Main integration function for analysis pipeline
def run_ml_pathogen_prediction(prokka_dir, output_dir, verbosity='standard'):
    """
    Main integration function for MetaQuest analysis pipeline v2.1
    FIXED: Proper verbosity handling and error recovery
    Uses trained voting/stacking ensemble model for pathogen prediction
    
    Args:
        prokka_dir: Path to Prokka output directory
        output_dir: Path to output directory
        verbosity: Verbosity level ('quiet', 'minimal', 'standard', 'debug')
    
    Returns:
        tuple: (results_list, summary_dict) or ([], {}) on error
    """
    try:
        formatter = OutputFormatter(verbosity=verbosity)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize predictor
        formatter.section_header("ML PATHOGEN PREDICTION v2.1")
        formatter.info("Enhanced feature extraction with pathogen motifs", indent=1)
        
        predictor = PathogenPredictor(verbosity=verbosity)
        
        # Check capabilities
        capabilities = predictor.check_capabilities()
        formatter.operation("Verifying ML pipeline capabilities")
        formatter.info(f"ML Pipeline: {'✓ Complete' if capabilities['ml_pipeline_complete'] else '✗ Incomplete'}", indent=2)
        formatter.info(f"Model: {capabilities['model_type']}", indent=2)
        formatter.info(f"Ensemble: {'Voting' if capabilities['voting_ensemble_loaded'] else 'Stacking' if capabilities['stacking_ensemble_loaded'] else 'Standard'}", indent=2)
        formatter.info(f"Features: {capabilities['selected_feature_count']}/{capabilities['feature_count']} selected", indent=2)
        
        if not capabilities['ml_pipeline_complete']:
            formatter.warning("ML prediction failed: Incomplete model artifacts")
            formatter.info("Continuing with traditional pathogen detection", indent=1)
            return [], {}
        
        if not capabilities['all_feature_names_available']:
            formatter.warning("ML prediction failed: Feature alignment data missing")
            formatter.info("Continuing with traditional pathogen detection", indent=1)
            return [], {}
        
        formatter.success("ML pipeline ready")
        
        # Run predictions
        results, summary = predictor.predict_from_prokka(Path(prokka_dir))
        
        if results and isinstance(results, list):
            # Save results
            formatter.operation("Saving ML prediction results")
            predictor.save_results(results, output_dir / "ml_pathogen_predictions.csv", "csv")
            predictor.save_results(results, output_dir / "ml_pathogen_predictions.json", "json")
            
            # Save summary
            if summary:
                with open(output_dir / "ml_pathogen_summary.json", 'w') as f:
                    json.dump(summary, f, indent=2)
                formatter.success("Saved prediction summary")
            
            # Print summary
            formatter.section_header("ML Prediction Results")
            formatter.result({
                'Total Sequences': f"{summary.get('total_sequences_analyzed', 0):,}",
                'Pathogenic': f"{summary.get('pathogenic_predictions', 0):,} ({summary.get('pathogenic_percentage', 0):.1f}%)",
                'Non-pathogenic': f"{summary.get('non_pathogenic_predictions', 0):,} ({summary.get('non_pathogenic_percentage', 0):.1f}%)",
                'High Confidence': f"{summary.get('high_confidence_predictions', 0):,} ({summary.get('high_confidence_percentage', 0):.1f}%)"
            })
            
            formatter.success("ML pathogen prediction complete", style='bold')
            
            return results, summary
        else:
            formatter.warning("No predictions generated")
            return [], {}
            
    except Exception as e:
        formatter = OutputFormatter(verbosity=verbosity)
        formatter.warning(f"ML prediction failed: {e}")
        formatter.info("Continuing with traditional pathogen detection", indent=1)
        logger.error(f"ML prediction failed: {e}")
        return [], {}