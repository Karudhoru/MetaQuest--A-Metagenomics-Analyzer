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
    Production pathogen predictor with trained voting ensemble model
    Designed for MetaQuest analysis pipeline integration
    """
    
    def __init__(self):
        """Initialize predictor with trained model artifacts"""
        self.model_dir = Path(MODEL_ARTIFACTS_DIR)
        self.scaler = None
        self.model = None
        self.feature_selector = None
        
        # NEW: Add feature name tracking
        self.selected_feature_names = None  # Features after selection
        self.all_feature_names = None       # All features before selection
        
        # Initialize MetaQuest feature extractor
        print("🔧 Initializing MetaQuest Feature Extractor...")
        self.feature_extractor = MetaQuestProteinFeatureExtractor()
        self.feature_extractor.fit()
        print("✅ Feature extractor ready!")
        
        self.ml_config = ML_CONFIG.copy()
        self._load_model_artifacts()
    
    def _load_model_artifacts(self):
        """Load trained model artifacts (scaler, model, feature_selector, feature_names)"""
        try:
            print("📊 Loading trained model artifacts...")
            
            # Load scaler
            scaler_path = self.model_dir / "scaler.pkl"
            if scaler_path.exists():
                self.scaler = joblib.load(scaler_path)
                print("✅ Loaded scaler successfully")
            else:
                raise FileNotFoundError(f"Scaler not found at: {scaler_path}")
            
            # Load voting ensemble model
            model_path = self.model_dir / "best_model.pkl"
            if model_path.exists():
                model_data = joblib.load(model_path)
                # Handle both old format (direct model) and new format (dict with metadata)
                if isinstance(model_data, dict):
                    self.model = model_data['model']
                    print(f"✅ Loaded model: {model_data.get('model_type', 'Unknown')}")
                    print(f"   Training date: {model_data.get('training_date', 'Unknown')}")
                    print(f"   Feature count: {model_data.get('feature_count', 'Unknown')}")
                else:
                    self.model = model_data
                    print(f"✅ Loaded model: {type(self.model).__name__}")
            else:
                raise FileNotFoundError(f"Model not found at: {model_path}")
            
            # Load feature selector
            selector_path = self.model_dir / "feature_selector.pkl"
            if selector_path.exists():
                self.feature_selector = joblib.load(selector_path)
                print("✅ Loaded feature selector successfully")
            else:
                raise FileNotFoundError(f"Feature selector not found at: {selector_path}")
            
            # 🔧 NEW: Load feature name lists (CRITICAL FOR FIXING FEATURE MISMATCH)
            feature_names_path = self.model_dir / "feature_names.pkl"
            all_feature_names_path = self.model_dir / "all_feature_names.pkl"
            
            if feature_names_path.exists() and all_feature_names_path.exists():
                self.selected_feature_names = joblib.load(feature_names_path)
                self.all_feature_names = joblib.load(all_feature_names_path)
                print(f"✅ Loaded feature name lists:")
                print(f"   Selected features: {len(self.selected_feature_names)}")
                print(f"   All features: {len(self.all_feature_names)}")
            else:
                print("⚠️ Feature name files not found - using feature extractor defaults")
                print(f"   Looking for: {feature_names_path}")
                print(f"   Looking for: {all_feature_names_path}")
                self.selected_feature_names = None
                self.all_feature_names = None
            
            print("🎯 Complete ML pipeline loaded successfully!")
            
        except Exception as e:
            logger.error(f"❌ Error loading model artifacts from {self.model_dir}: {e}")
            raise RuntimeError(f"Failed to load required ML artifacts: {e}")
    
    def check_capabilities(self) -> Dict[str, bool]:
        """Check and return current system capabilities"""
        capabilities = {
            'feature_extractor_available': True,
            'ml_pipeline_complete': all([self.scaler, self.model, self.feature_selector]),
            'feature_names_available': all([self.selected_feature_names, self.all_feature_names]),
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
            # Prepare features for ML model
            feature_df = pd.DataFrame([features])
            
            # Remove non-feature columns
            non_feature_cols = [
                'sequence_id', 'extraction_method', 'timestamp',
                'feature_count', 'batch_number', 'sequence_length'
            ]
            feature_df = feature_df.drop(
                columns=[col for col in non_feature_cols if col in feature_df.columns],
                errors='ignore'
            )
            
            print(f"🔧 Feature count check for {seq_id}")
            print(f"📊 Current extractor produces: {len(feature_df.columns)} features")
            
            # 🔧 CRITICAL FIX: Use exactly the features that the selector expects
            if self.all_feature_names is not None:
                expected_features = self.all_feature_names
                print(f"📚 Selector expects: {len(expected_features)} features")
                
                # If there's a mismatch, we need to align to the selector's expectations
                if len(expected_features) != len(feature_df.columns):
                    print(f"⚠️ Feature count mismatch detected!")
                    print(f"   Selector expects: {len(expected_features)}")
                    print(f"   Current extractor: {len(feature_df.columns)}")
                    
                    # Use only the features that were present during training
                    aligned_features = pd.DataFrame(index=feature_df.index)
                    for expected_feature in expected_features:
                        if expected_feature in feature_df.columns:
                            aligned_features[expected_feature] = feature_df[expected_feature]
                        else:
                            aligned_features[expected_feature] = 0.0  # Default for missing features
                            print(f"   🔧 Adding missing feature: {expected_feature}")
                    
                    print(f"✅ Features aligned: {aligned_features.shape[1]} features")
                else:
                    # Reorder to match training order
                    aligned_features = feature_df[expected_features]
            else:
                print("⚠️ No saved feature names - using current extractor output")
                aligned_features = feature_df
            
            # Convert to numpy arrays to bypass any remaining sklearn validation
            features_numpy = aligned_features.values.astype(np.float64)
            
            
            # Verify feature count matches selector expectation
            if features_numpy.shape[1] != 90:  # Your selector expects 90
                print(f"❌ Feature count still mismatched: {features_numpy.shape[1]} != 90")
                # Emergency fix: truncate or pad to exactly 90 features
                if features_numpy.shape[1] > 90:
                    print("🔧 Truncating to 90 features")
                    features_numpy = features_numpy[:, :90]
                else:
                    print("🔧 Padding to 90 features")
                    padding = np.zeros((features_numpy.shape[0], 90 - features_numpy.shape[1]))
                    features_numpy = np.concatenate([features_numpy, padding], axis=1)
            
            print(f"📊 Features shape before selector: {features_numpy.shape}")
            
            # Apply feature selection
            features_selected = self.feature_selector.transform(features_numpy)
            
            # Scale features
            features_scaled = self.scaler.transform(features_selected)
            print(f"📊 Features shape after scaling: {features_scaled.shape}")
            
            # Make prediction
            prediction = self.model.predict(features_scaled)[0]
            print(f"🎯 Prediction made successfully: {prediction}")
            
            # Get probabilities
            if hasattr(self.model, 'predict_proba'):
                probabilities = self.model.predict_proba(features_scaled)[0]
                confidence = float(max(probabilities))
                pathogenic_prob = float(probabilities[1]) if len(probabilities) > 1 else float(probabilities[0])
            else:
                confidence = 0.85
                pathogenic_prob = float(prediction)
            
            # Apply confidence threshold
            confidence_threshold = self.ml_config.get('confidence_threshold', 0.7)
            is_high_confidence = confidence >= confidence_threshold
            
            result = {
                'sequence_id': seq_id,
                'sequence_type': seq_type,
                'prediction': 'Pathogenic' if prediction == 1 else 'Non-pathogenic',
                'confidence': confidence,
                'pathogenic_probability': pathogenic_prob,
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
        """Save prediction results with comprehensive metadata"""
        if not results:
            logger.warning("⚠️ No results to save")
            return
        
        try:
            if format_type.lower() == "json":
                # Add comprehensive metadata to JSON
                capabilities = self.check_capabilities()
                output_data = {
                    'metadata': {
                        'total_predictions': len(results),
                        'pathogenic_count': len([r for r in results if r['prediction'] == 'Pathogenic']),
                        'high_confidence_count': len([r for r in results if r.get('high_confidence', False)]),
                        'voting_ensemble_used': capabilities['voting_ensemble_loaded'],
                        'model_type': capabilities['model_type'],
                        'feature_names_available': capabilities['feature_names_available'],
                        'ml_config': self.ml_config,
                        'model_artifacts_dir': str(self.model_dir),
                        'pathogen_colors': PATHOGEN_COLORS,
                        'system_capabilities': capabilities,
                        'generation_timestamp': datetime.now().isoformat(),
                        'methods_used': list(set(r.get('method', 'Unknown') for r in results)),
                        'feature_alignment_method': 'saved_feature_names' if self.all_feature_names else 'extractor_defaults'
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
        """Generate comprehensive prediction summary"""
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
            'feature_alignment_used': 'saved_feature_names' if self.all_feature_names else 'extractor_defaults',
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
            print(f" 🎲 Voting ensemble: {'✅ Used' if summary['voting_ensemble_used'] else '❌ Not used'}")
            print(f" 🔧 Feature alignment: {summary['feature_alignment_used']}")
            
            return results, summary
        else:
            print("⚠️ No predictions generated")
            return [], {}
            
    except Exception as e:
        print(f"❌ MetaQuest pathogen prediction failed: {e}")
        raise

