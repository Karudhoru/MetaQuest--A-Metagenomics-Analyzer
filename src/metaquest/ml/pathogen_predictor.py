#!/usr/bin/env python3

"""
MetaQuest Pathogen Predictor - Production Version
Integrated with trained voting ensemble model
Compatible with Prokka output processing
No ESM2/torch dependencies
"""

import json
import logging
import pickle
import joblib
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Union
from datetime import datetime
from collections import Counter
import warnings
from Bio import SeqIO

# Import our feature extractor
from .feature_extractor import MetaQuestProteinFeatureExtractor
from ..config import MODEL_ARTIFACTS_DIR, ML_CONFIG, PATHOGEN_COLORS
# 🔧 DISABLE sklearn feature name validation globally
import sklearn
import warnings
warnings.filterwarnings('ignore')

# Disable feature name validation
import os
os.environ['SKLEARN_ENABLE_METADATA_ROUTING'] = 'False'

warnings.filterwarnings('ignore')
logger = logging.getLogger(__name__)

class PathogenPredictor:
    """
    Production pathogen predictor with robust feature alignment.
    """
    def __init__(self):
        self.model_dir = Path(MODEL_ARTIFACTS_DIR)
        self.scaler = None
        self.model = None
        self.feature_selector = None
        self.training_columns = None # This is the correct attribute for storing feature names

        print("🔧 Initializing MetaQuest Feature Extractor...")
        self.feature_extractor = MetaQuestProteinFeatureExtractor()
        self.feature_extractor.fit()
        print("✅ Feature extractor ready!")
        
        self.ml_config = ML_CONFIG.copy()
        self._load_model_artifacts()
    
    def _load_model_artifacts(self):
        """Load trained model artifacts."""
        try:
            print("📊 Loading trained model artifacts...")
            self.scaler = joblib.load(self.model_dir / "scaler.pkl")
            model_data = joblib.load(self.model_dir / "best_model.pkl")
            self.model = model_data['model']
            self.feature_selector = joblib.load(self.model_dir / "feature_selector.pkl")
            self.training_columns = joblib.load(self.model_dir / "all_feature_names.pkl")
            
            print(f"✅ Loaded scaler, model, and selector.")
            print(f"✅ Loaded training feature order ({len(self.training_columns)} features).")
            print("🎯 Complete ML pipeline loaded successfully!")
            
        except Exception as e:
            logger.error(f"❌ Error loading model artifacts from {self.model_dir}: {e}")
            raise RuntimeError(f"Failed to load required ML artifacts: {e}")
    
    
    def check_capabilities(self) -> Dict[str, bool]:
        """Check and return current system capabilities."""
        # The check now correctly verifies the 'training_columns' attribute.
        capabilities = {
            'feature_extractor_available': True,
            'ml_pipeline_complete': all([self.scaler, self.model, self.feature_selector]),
            'feature_names_available': self.training_columns is not None,
            'voting_ensemble_loaded': 'voting' in str(type(self.model)).lower(),
            'model_type': type(self.model).__name__ if self.model else 'None'
        }
        return capabilities
    
    def extract_protein_features(self, protein_seq: str, seq_id: str) -> Dict:
        """Extract features using MetaQuest Feature Extractor"""
        try:
            logger.debug(f"🧬 Extracting features for {seq_id}")
            
            # Extract features using our trained extractor
            features = self.feature_extractor.extract_features(protein_seq)
            
            # Add metadata
            features.update({
                'sequence_id': seq_id,
                'sequence_length': len(protein_seq),
                'extraction_method': 'MetaQuest_Production',
                'feature_count': len(features),
                'timestamp': datetime.now().isoformat()
            })
            
            logger.debug(f"✅ Extracted {len(features)} features for {seq_id}")
            return features
            
        except Exception as e:
            logger.error(f"❌ Feature extraction failed for {seq_id}: {e}")
            raise RuntimeError(f"Feature extraction failed for {seq_id}: {e}")
    
    def predict_from_prokka(self, prokka_dir: Path, batch_size: int = None) -> List[Dict]:
        """
        Main prediction function for processing Prokka output
        Compatible with both FASTA and FASTQ derived annotations
        """
        results = []
        
        # Use config-defined batch size if not provided
        if batch_size is None:
            batch_size = self.ml_config.get('batch_size', 100)
        
        # Look for protein FASTA file from Prokka
        protein_files = list(prokka_dir.glob("*.faa"))
        if not protein_files:
            raise FileNotFoundError(f"❌ No protein FASTA files found in {prokka_dir}")
        
        protein_file = protein_files[0]
        logger.info(f"🔬 Processing Prokka proteins from: {protein_file}")
        
        try:
            sequences = list(SeqIO.parse(protein_file, "fasta"))
            total_sequences = len(sequences)
            logger.info(f"📊 Found {total_sequences} protein sequences")
            
            # Process sequences in batches for efficiency
            if total_sequences > batch_size:
                logger.info(f"🚀 Using batch processing (batch size: {batch_size})")
                results = self._batch_predict_sequences(sequences, batch_size)
            else:
                logger.info("🔄 Using sequential processing")
                results = self._sequential_predict_sequences(sequences)
                
        except Exception as e:
            logger.error(f"❌ Error processing Prokka output: {e}")
            raise
        
        logger.info(f"✅ Completed processing: {len(results)} predictions generated")
        return results
    
    def _sequential_predict_sequences(self, sequences: List) -> List[Dict]:
        """Sequential processing for protein sequences"""
        results = []
        total = len(sequences)
        
        for i, record in enumerate(sequences):
            if i % 50 == 0:  # More frequent updates
                logger.info(f"🔄 Processing protein {i+1}/{total}")
            
            try:
                # Extract protein features
                features = self.extract_protein_features(str(record.seq), record.id)
                
                # Make prediction using trained model
                prediction_result = self._make_ml_prediction(features, record.id, "protein")
                
                if prediction_result:
                    results.append(prediction_result)
                    
            except Exception as e:
                logger.error(f"❌ Error processing sequence {record.id}: {e}")
                continue
        
        return results
    
    def _batch_predict_sequences(self, sequences: List, batch_size: int = 100) -> List[Dict]:
        """Batch processing for multiple sequences"""
        results = []
        total_batches = (len(sequences) - 1) // batch_size + 1
        
        for i in range(0, len(sequences), batch_size):
            batch = sequences[i:i+batch_size]
            batch_num = i // batch_size + 1
            logger.info(f"🚀 Processing batch {batch_num}/{total_batches} ({len(batch)} sequences)")
            
            try:
                # Extract sequences and IDs
                batch_sequences = [str(record.seq) for record in batch]
                batch_ids = [record.id for record in batch]
                
                # Batch feature extraction
                features_batch = self.feature_extractor.transform(batch_sequences)
                
                # Process each result in the batch
                for j, (seq_id, original_seq) in enumerate(zip(batch_ids, batch_sequences)):
                    try:
                        # Get features for this sequence
                        features = features_batch.iloc[j].to_dict()
                        
                        # Add metadata
                        features.update({
                            'sequence_id': seq_id,
                            'sequence_length': len(original_seq),
                            'extraction_method': 'MetaQuest_Production_Batch',
                            'batch_number': batch_num
                        })
                        
                        # Make prediction
                        prediction_result = self._make_ml_prediction(features, seq_id, "protein")
                        if prediction_result:
                            results.append(prediction_result)
                            
                    except Exception as seq_e:
                        logger.error(f"❌ Error processing sequence {seq_id} in batch: {seq_e}")
                        continue
                        
            except Exception as e:
                logger.error(f"❌ Error in batch {batch_num}: {e}")
                continue
        
        return results
    
    def _make_ml_prediction(self, features: Dict, seq_id: str, seq_type: str) -> Optional[Dict]:
        try:
            feature_df = pd.DataFrame([features])
            
            # --- THIS BLOCK REPLACES THE OLD, COMPLEX "PATCH" ---
            if self.training_columns is None:
                raise RuntimeError("Training column order not loaded. Cannot proceed with prediction.")

            # Reindex the new data to match the training data's column order exactly.
            # - Missing columns will be added and filled with 0.
            # - Extra columns not in the training data will be dropped.
            aligned_features = feature_df.reindex(columns=self.training_columns, fill_value=0.0)
            # --- END OF FIX ---

            # Convert to numpy array for the model
            features_numpy = aligned_features.values.astype(np.float64)
            
            # Apply feature selection and scaling (as before)
            features_selected = self.feature_selector.transform(features_numpy)
            features_scaled = self.scaler.transform(features_selected)
            
            # Make prediction (as before)
            prediction = self.model.predict(features_scaled)[0]
            probabilities = self.model.predict_proba(features_scaled)[0]
            confidence = float(max(probabilities))
            
            # Apply confidence threshold
            confidence_threshold = self.ml_config.get('confidence_threshold', 0.7)
            is_high_confidence = confidence >= confidence_threshold
            
            result = {
                'sequence_id': seq_id,
                'sequence_type': seq_type,
                'prediction': 'Pathogenic' if prediction == 1 else 'Non-pathogenic',
                'confidence': confidence,
                'pathogenic_probability': confidence,
                'high_confidence': is_high_confidence,
                'confidence_threshold': confidence_threshold,
                'sequence_length': features.get('sequence_length', 0),
                'method': 'ML_VotingEnsemble_MetaQuest',
                'extraction_method': features.get('extraction_method', 'MetaQuest_Production'),
                'model_type': type(self.model).__name__,
                'features_used': features_numpy.shape[1],
                'features_selected': features_selected.shape[1],
                'timestamp': datetime.now().isoformat()
            }
            
            print(f"✅ ML prediction for {seq_id}: {result['prediction']} (confidence: {confidence:.3f})")
            return result
            
        except Exception as e:
            print(f"❌ Detailed error for {seq_id}: {e}")
            import traceback
            traceback.print_exc()
            raise RuntimeError(f"ML prediction failed for {seq_id}: {e}")


    def save_results(self, results: List[Dict], output_file: Path, format_type: str = "csv"):
        """Save prediction results with comprehensive metadata."""
        if not results:
            return
        
        try:
            if format_type.lower() == "json":
                # --- FIX #1: Use the correct attribute name here ---
                output_data = {
                    'metadata': {
                        'total_predictions': len(results),
                        'feature_alignment_method': 'saved_feature_names' if self.training_columns is not None else 'extractor_defaults',
                        'generation_timestamp': datetime.now().isoformat(),
                    },
                    'predictions': results
                }
                with open(output_file, 'w') as f:
                    json.dump(output_data, f, indent=2)
            else:
                df = pd.DataFrame(results)
                df.to_csv(output_file, index=False)
            
            logger.info(f"✅ Results saved to: {output_file}")
            
        except Exception as e:
            logger.error(f"❌ Error saving results: {e}")
            raise

    def get_prediction_summary(self, results: List[Dict]) -> Dict:
        """Generate comprehensive prediction summary."""
        if not results:
            return {}
        
        pathogenic_sequences = [r for r in results if r['prediction'] == 'Pathogenic']
        high_confidence_predictions = [r for r in results if r.get('high_confidence', False)]
        
        # Calculate statistics
        confidence_scores = [r['confidence'] for r in results if 'confidence' in r]
        pathogenic_scores = [r.get('pathogenic_probability', 0) for r in pathogenic_sequences]

        summary = {
            'total_sequences_analyzed': len(results),
            'pathogenic_predictions': len(pathogenic_sequences),
            'pathogenic_percentage': (len(pathogenic_sequences) / len(results)) * 100 if results else 0,
            'high_confidence_predictions': len(high_confidence_predictions),
            'high_confidence_percentage': (len(high_confidence_predictions) / len(results)) * 100 if results else 0,
            'confidence_threshold': self.ml_config.get('confidence_threshold', 0.7),
            'average_confidence': np.mean(confidence_scores) if confidence_scores else 0,
            'median_confidence': np.median(confidence_scores) if confidence_scores else 0,
            'min_confidence': np.min(confidence_scores) if confidence_scores else 0,
            'max_confidence': np.max(confidence_scores) if confidence_scores else 0,
            'average_pathogenic_score': np.mean(pathogenic_scores) if pathogenic_scores else 0,
            'methods_distribution': Counter(r.get('method', 'Unknown') for r in results),
            'model_type': type(self.model).__name__ if self.model else 'Unknown',
            'ml_config_used': self.ml_config,
            'model_artifacts_dir': str(self.model_dir),
            'voting_ensemble_used': 'voting' in str(type(self.model)).lower(),
            'feature_alignment_used': 'saved_feature_names' if self.training_columns is not None else 'extractor_defaults',  # FIX #2: Use correct attribute
            'system_capabilities': self.check_capabilities(),
            'pathogen_colors': PATHOGEN_COLORS,
            'top_pathogenic_sequences': sorted(pathogenic_sequences, key=lambda x: x['confidence'], reverse=True)[:10],
            'top_high_confidence_sequences': sorted(high_confidence_predictions, key=lambda x: x['confidence'], reverse=True)[:10]
        }
        
        return summary

# Main integration function for analysis pipeline
def run_ml_pathogen_prediction(prokka_dir, output_dir, batch_size=None):
    """
    Main integration function for MetaQuest analysis pipeline
    Uses trained voting ensemble model for pathogen prediction
    🔧 UPDATED: Enhanced error handling and feature alignment
    """
    try:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize predictor with trained model
        print("🚀 Initializing MetaQuest Pathogen Predictor with trained voting ensemble...")
        predictor = PathogenPredictor()
        
        # Check capabilities
        capabilities = predictor.check_capabilities()
        if not capabilities['ml_pipeline_complete']:
            raise RuntimeError("❌ ML pipeline incomplete - missing model artifacts")
        
        print(f"✅ Model loaded: {capabilities['model_type']}")
        print(f"✅ Voting ensemble: {'Yes' if capabilities['voting_ensemble_loaded'] else 'No'}")
        print(f"✅ Feature names available: {'Yes' if capabilities['feature_names_available'] else 'No'}")
        
        # Warn if feature names not available
        if not capabilities['feature_names_available']:
            print("⚠️ WARNING: Feature name files not found!")
            print("   This may cause feature mismatch errors.")
            print("   Please ensure feature_names.pkl and all_feature_names.pkl are in model_artifacts/")
        
        # Run predictions
        print(f"🔬 Processing Prokka directory: {prokka_dir}")
        results = predictor.predict_from_prokka(Path(prokka_dir), batch_size)
        
        if results:
            # Generate summary
            summary = predictor.get_prediction_summary(results)
            
            # Save results
            predictor.save_results(results, output_dir / "ml_pathogen_predictions.csv", "csv")
            predictor.save_results(results, output_dir / "ml_pathogen_predictions.json", "json")
            
            # Save summary
            with open(output_dir / "ml_pathogen_summary.json", 'w') as f:
                json.dump(summary, f, indent=2)
            
            # Log results
            print(f"🎯 MetaQuest ML Pathogen Prediction Complete:")
            print(f" 📊 {summary['pathogenic_predictions']}/{summary['total_sequences_analyzed']} sequences predicted as pathogenic ({summary['pathogenic_percentage']:.1f}%)")
            print(f" ⭐ {summary['high_confidence_predictions']} high-confidence predictions")
            print(f" 🤖 Model: {summary['model_type']}")
            
            return results, summary
        else:
            print("⚠️ No predictions generated")
            return [], {}
            
    except Exception as e:
        print(f"❌ MetaQuest pathogen prediction failed: {e}")
        raise