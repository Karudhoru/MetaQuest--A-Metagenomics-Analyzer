#!/usr/bin/env python3
"""
MetaQuest Pathogen Predictor v3.0

Production predictor with:
- XGBoost + LightGBM ensemble (simple average)
- Isotonic probability calibration
- Uncertainty estimation (model disagreement)
- Per-prediction SHAP explanations

No GPU required. No external tool dependencies at inference.
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import joblib
import numpy as np
from Bio import SeqIO

from .feature_extractor import ProteinFeatureExtractor
from .explainability import explain_prediction, render_shap_text

logger = logging.getLogger(__name__)

# Default model artifacts location
_DEFAULT_MODEL_DIR = Path(__file__).parent / "model_artifacts"


class PathogenPredictor:
    """
    Production pathogen predictor v3.0.

    Predicts pathogenicity of protein sequences with calibrated
    probabilities and SHAP-based explanations.
    """

    def __init__(self, model_dir: Path = None):
        """
        Initialize predictor by loading model artifacts.

        Args:
            model_dir: Directory containing trained model files.
                       Defaults to src/metaquest/ml/model_artifacts/
        """
        self.model_dir = Path(model_dir) if model_dir else _DEFAULT_MODEL_DIR
        self.extractor = ProteinFeatureExtractor()

        self.xgb = None
        self.lgbm = None
        self.scaler = None
        self.calibrator = None
        self.shap_explainer = None
        self.feature_names = None

        self._load_artifacts()

    def _load_artifacts(self):
        """Load trained model artifacts."""
        try:
            # Required artifacts
            self.scaler = joblib.load(self.model_dir / "scaler.pkl")
            self.xgb = joblib.load(self.model_dir / "model_xgb.pkl")
            self.lgbm = joblib.load(self.model_dir / "model_lgbm.pkl")
            self.calibrator = joblib.load(self.model_dir / "calibrator.pkl")

            with open(self.model_dir / "feature_names.json") as f:
                self.feature_names = json.load(f)

            # Optional: SHAP explainer
            shap_path = self.model_dir / "shap_explainer.pkl"
            if shap_path.exists():
                self.shap_explainer = joblib.load(shap_path)

            logger.info(
                f"Loaded model: {len(self.feature_names)} features, "
                f"SHAP={'yes' if self.shap_explainer else 'no'}"
            )
        except FileNotFoundError as e:
            logger.warning(f"Model artifacts not found: {e}")
            logger.warning("ML predictions will be unavailable. Train a model first.")

    @property
    def is_ready(self) -> bool:
        """Check if model is loaded and ready for prediction."""
        return self.xgb is not None and self.lgbm is not None

    def predict_sequence(self, sequence: str) -> Dict:
        """
        Predict pathogenicity of a single protein sequence.

        Returns:
            Dict with: label, confidence, uncertainty, explanation
        """
        if not self.is_ready:
            return {
                "label": "unknown",
                "confidence": 0.0,
                "uncertainty": 1.0,
                "explanation": None,
            }

        # Extract features
        features = self.extractor.extract_single(sequence)

        # Align to training feature order
        feature_vector = np.array(
            [features.get(name, 0.0) for name in self.feature_names]
        ).reshape(1, -1)

        # Scale
        feature_scaled = self.scaler.transform(feature_vector)

        # Get probabilities from both models
        xgb_prob = self.xgb.predict_proba(feature_scaled)[0, 1]
        lgbm_prob = self.lgbm.predict_proba(feature_scaled)[0, 1]

        # Ensemble average
        raw_prob = (xgb_prob + lgbm_prob) / 2.0

        # Calibrate
        calibrated_prob = float(self.calibrator.predict(np.array([raw_prob]))[0])

        # Uncertainty (model disagreement)
        uncertainty = float(abs(xgb_prob - lgbm_prob))

        # Classification
        if uncertainty > 0.25:
            label = "inconclusive"
        elif calibrated_prob >= 0.5:
            label = "pathogenic"
        else:
            label = "non-pathogenic"

        # SHAP explanation
        explanation = None
        if self.shap_explainer is not None:
            try:
                shap_values = self.shap_explainer.shap_values(feature_scaled)[0]
                explanation = explain_prediction(
                    shap_values=shap_values,
                    feature_names=self.feature_names,
                    feature_values=feature_vector[0],
                    top_n=5,
                )
            except Exception as e:
                logger.debug(f"SHAP explanation failed: {e}")

        return {
            "label": label,
            "confidence": calibrated_prob,
            "raw_probability": float(raw_prob),
            "uncertainty": uncertainty,
            "explanation": explanation,
        }

    def predict_fasta(self, fasta_path: Path, min_length: int = 50) -> List[Dict]:
        """
        Predict pathogenicity for all proteins in a FASTA file.

        Args:
            fasta_path: Path to protein FASTA file.
            min_length: Minimum sequence length to predict.

        Returns:
            List of prediction dicts, one per protein.
        """
        if not self.is_ready:
            logger.warning("Model not loaded — returning empty predictions")
            return []

        predictions = []
        for rec in SeqIO.parse(fasta_path, "fasta"):
            seq = str(rec.seq).upper().replace("*", "")
            if len(seq) < min_length:
                continue

            result = self.predict_sequence(seq)
            result["protein_id"] = rec.id
            result["protein_length"] = len(seq)
            result["description"] = rec.description[:100]
            predictions.append(result)

        # Sort by confidence (highest risk first)
        predictions.sort(key=lambda x: x["confidence"], reverse=True)

        logger.info(
            f"Predicted {len(predictions)} proteins: "
            f"{sum(1 for p in predictions if p['label'] == 'pathogenic')} pathogenic, "
            f"{sum(1 for p in predictions if p['label'] == 'inconclusive')} inconclusive"
        )

        return predictions

    def model_info(self) -> Dict:
        """Return model metadata."""
        metrics_path = self.model_dir / "training_metrics.json"
        if metrics_path.exists():
            with open(metrics_path) as f:
                return json.load(f)
        return {"status": "no metadata available"}


def run_ml_pathogen_prediction(prokka_dir: Path, output_dir: Path):
    """
    Pipeline-compatible entry point for ML pathogen prediction.

    Scans Prokka output for protein FASTA files, runs predictions,
    writes results to JSON.

    Args:
        prokka_dir: Directory containing Prokka annotation output (.faa files)
        output_dir: Directory to write ml_pathogen_predictions.json

    Returns:
        (results_list, summary_dict)
    """
    prokka_dir = Path(prokka_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find protein FASTA files from Prokka
    faa_files = list(prokka_dir.glob("**/*.faa"))
    if not faa_files:
        logger.info("No protein FASTA (.faa) files found in %s", prokka_dir)
        return [], {"status": "no_proteins", "total": 0}

    predictor = PathogenPredictor()
    if not predictor.is_ready:
        logger.warning("ML model not loaded — skipping prediction")
        return [], {"status": "model_not_ready", "total": 0}

    all_predictions = []
    for faa_path in faa_files:
        predictions = predictor.predict_fasta(faa_path)
        all_predictions.extend(predictions)

    # Write JSON output
    output_file = output_dir / "ml_pathogen_predictions.json"
    with open(output_file, "w") as f:
        json.dump(all_predictions, f, indent=2, default=str)

    # Generate summary
    n_pathogenic = sum(1 for p in all_predictions if p["label"] == "pathogenic")
    n_inconclusive = sum(1 for p in all_predictions if p["label"] == "inconclusive")
    n_safe = sum(1 for p in all_predictions if p["label"] == "non-pathogenic")

    summary = {
        "status": "completed",
        "total": len(all_predictions),
        "pathogenic": n_pathogenic,
        "non_pathogenic": n_safe,
        "inconclusive": n_inconclusive,
        "high_confidence_pathogenic": sum(
            1 for p in all_predictions
            if p["label"] == "pathogenic" and p["confidence"] > 0.8
        ),
        "model_info": predictor.model_info(),
    }

    logger.info(
        "ML prediction complete: %d proteins — %d pathogenic, %d safe, %d inconclusive",
        len(all_predictions), n_pathogenic, n_safe, n_inconclusive,
    )

    return all_predictions, summary
