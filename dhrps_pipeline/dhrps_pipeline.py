#!/usr/bin/env python3
"""
DHRPs Pipeline: Multi-modal Feature Integration → Model Training → Ensemble Fusion → Evaluation
"""

import os
import sys
import logging
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import joblib

from sklearn.feature_selection import VarianceThreshold

# Import custom components
import onekey_algo.custom.components as okcomp
from onekey_algo.custom.components.comp1 import (
    normalize_df,
    create_clf_model,
    plot_feature_importance,
)
from onekey_algo.custom.components.metrics import analysis_pred_binary


class DHRPsPipeline:
    """Pipeline for multi-modal feature integration and ensemble modeling"""
    
    def __init__(self, config: Dict):
        """
        Initialize the DHRPs pipeline
        
        Args:
            config: Configuration dictionary with pipeline parameters
        """
        self.config = config
        self.setup_logging()
        self.setup_directories()
        self.set_seed()
        
    def setup_logging(self):
        """Configure logging for the pipeline"""
        logging.basicConfig(
            level=logging.INFO,
            format="[%(asctime)s] %(levelname)s - %(message)s",
            datefmt="%H:%M:%S",
        )
        self.logger = logging.getLogger("DHRPs")
        
    def setup_directories(self):
        """Create necessary directories for output"""
        base_dir = Path(self.config.get("base_dir", "DHRPs"))
        self.dirs = {
            "base": base_dir,
            "results": base_dir / "results",
            "img": base_dir / "img",
            "features": base_dir / "features",
            "models": base_dir / "models"
        }
        
        for d in self.dirs.values():
            d.mkdir(parents=True, exist_ok=True)
            
    def set_seed(self):
        """Set random seed for reproducibility"""
        seed = self.config.get("seed", 0)
        np.random.seed(seed)
        os.environ["PYTHONHASHSEED"] = str(seed)
        os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
        
        # Set plot parameters
        plt.rcParams["pdf.fonttype"] = 42
        plt.rcParams["ps.fonttype"] = 42
        plt.rcParams["figure.dpi"] = 300
        
    def load_data(self) -> pd.DataFrame:
        """Load and return the dataset"""
        data_path = Path(self.config["data_path"])
        if not data_path.exists():
            raise FileNotFoundError(f"Data file not found: {data_path}")
            
        data = pd.read_csv(data_path)
        return data
    
    def preprocess_data(self, data: pd.DataFrame) -> Tuple[pd.DataFrame, pd.Series]:
        """Preprocess and standardize the data"""
        if "ID" not in data.columns:
            raise KeyError("Input data missing 'ID' column")
            
        ids = data["ID"].copy()
        label_col = self.config["label_column"]
        split_col = self.config["split_column"]
        
        self.logger.info("Label distribution:\n%s", data[label_col].value_counts())
        
        # Separate external data
        external_data = data[data[split_col] == "EXTERNAL"].copy()
        main_data = data[data[split_col].isin(["train", "test", "validation"])].copy()
        
        # Normalize data
        norm_data = normalize_df(main_data, not_norm=[label_col, split_col]).dropna(axis=1)
        norm_external = normalize_df(external_data, not_norm=[label_col, split_col]).copy()
        
        # Combine datasets
        cols = norm_data.columns.tolist()
        norm_external = norm_external.reindex(columns=cols)
        combined_data = pd.concat([norm_data, norm_external], axis=0).fillna(0)
        
        return combined_data, ids
        
    def variance_filter(self, data: pd.DataFrame) -> pd.DataFrame:
        """Apply variance threshold filtering"""
        features = data.columns[2:]
        train_data = data[data[self.config["split_column"]] == "train"].copy()
        
        selector = VarianceThreshold(threshold=self.config["var_threshold"])
        selector.fit_transform(train_data[features])
        
        kept_features = features[selector.get_support()].tolist()
        self.logger.info("Variance filtering kept %d features", len(kept_features))
        
        return data[[self.config["label_column"], self.config["split_column"]] + kept_features]
        
    def split_datasets(self, data: pd.DataFrame, ids: pd.Series) -> Dict[str, Tuple]:
        """Split data into training and validation sets"""
        splits = {}
        
        for split_name in self.config["split_names"]:
            split_data = data[data[self.config["split_column"]] == split_name].reset_index(drop=True)
            split_ids = ids[split_data.index]
            
            X = split_data.drop(columns=[self.config["label_column"], self.config["split_column"]])
            y = split_data[self.config["label_column"]]
            
            splits[split_name] = (X, y, split_ids)
            
        return splits
        
    def train_models(self, X_train: pd.DataFrame, y_train: pd.DataFrame) -> Dict[str, object]:
        """Train multiple classification models"""
        models = create_clf_model(self.config["model_names"])
        models = {k: v for k, v in models.items() if k in self.config["model_names"]}
        
        for name, model in models.items():
            model.fit(X_train, y_train)
            
            # Save model
            model_path = self.dirs["models"] / f"{name}_{self.config['label_column']}.pkl"
            joblib.dump(model, model_path)
            
            # Try to plot feature importance
            try:
                plot_feature_importance(model, list(X_train.columns), 
                                      save_dir=str(self.dirs["img"]))
            except Exception as e:
                self.logger.warning("Could not plot feature importance for %s: %s", name, e)
                
        return models
        
    def predict(self, model: object, X: pd.DataFrame) -> np.ndarray:
        """Generate predictions from a model"""
        if X is None or len(X) == 0:
            return np.zeros((0, 2))
            
        if hasattr(model, "predict_proba"):
            proba = model.predict_proba(X)
            if proba.ndim == 1:
                proba = np.vstack([1 - proba, proba]).T
            elif proba.shape[1] == 1:
                proba = np.hstack([1 - proba, proba])
            return proba
            
        if hasattr(model, "decision_function"):
            scores = model.decision_function(X)
            scores = (scores - scores.min()) / (scores.max() - scores.min() + 1e-12)
            return np.vstack([1 - scores, scores]).T
            
        preds = model.predict(X)
        scores = preds.astype(float)
        return np.vstack([1 - scores, scores]).T
        
    def evaluate_models(self, models: Dict[str, object], splits: Dict[str, Tuple]) -> pd.DataFrame:
        """Evaluate model performance across all splits"""
        results = []
        
        for model_name, model in models.items():
            for split_name, (X, y, _) in splits.items():
                if X is None or X.empty:
                    continue
                    
                scores = self.predict(model, X)
                if scores.shape[0] == 0:
                    continue
                    
                # Calculate metrics
                metrics = analysis_pred_binary(y, scores[:, 1])
                acc, auc, ci, tpr, tnr, ppv, npv, precision, recall, f1, thres = metrics
                
                ci_str = f"{auc:.3f}({ci[0]:.3f} - {ci[1]:.3f})"
                
                results.append((
                    model_name, acc, ci_str, tpr, tnr, ppv, npv, 
                    precision, recall, f1, thres, f"{self.config['label_column']}-{split_name}"
                ))
                
        # Create results DataFrame
        columns = [
            "model", "accuracy", "auc_ci", "sensitivity", "specificity", 
            "ppv", "npv", "precision", "recall", "f1", "threshold", "dataset"
        ]
        
        return pd.DataFrame(results, columns=columns)
        
    def softmax(self, z: np.ndarray) -> np.ndarray:
        """Apply softmax function to input array"""
        t = np.exp(z - z.max(axis=1, keepdims=True))
        return t / (t.sum(axis=1, keepdims=True) + 1e-12)
        
    def weight_fusion(self, weights: np.ndarray, pred_scores: List) -> Tuple[List, List]:
        """Perform weighted fusion of model predictions"""
        fusion_scores_all, fusion_preds_all = [], []
        
        for weight, scores in zip(weights, pred_scores):
            fused_scores = tuple(
                np.mean([s * w for s, w in zip(ds_scores, weight)], axis=0)
                for ds_scores in zip(*scores)
            )
            
            # Normalize fused scores
            norm_fused_scores = []
            for arr in fused_scores:
                arr = arr.copy()
                arr /= (arr.sum(axis=1, keepdims=True) + 1e-12)
                norm_fused_scores.append(arr)
                
            fused_preds = tuple(np.argmax(s, axis=1) for s in norm_fused_scores)
            fusion_scores_all.append([tuple(norm_fused_scores)])
            fusion_preds_all.append([tuple(fused_preds)])
            
        return fusion_preds_all, fusion_scores_all
        
    def visualize_results(self, results: pd.DataFrame):
        """Create visualizations of model performance"""
        # Accuracy bar plot
        plt.figure(figsize=(12, 8))
        sns.barplot(x="model", y="accuracy", data=results, hue="dataset")
        plt.title("Model Accuracy by Dataset")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        plt.savefig(self.dirs["img"] / "model_accuracy_bar.pdf", 
                   bbox_inches="tight", format="pdf")
        plt.close()
        
    def run(self):
        """Execute the full DHRPs pipeline"""
        self.logger.info("Starting DHRPs pipeline")
        
        try:
            # Load data
            data = self.load_data()
            
            # Preprocess data
            processed_data, ids = self.preprocess_data(data)
            
            # Apply variance filtering
            filtered_data = self.variance_filter(processed_data)
            
            # Split data
            splits = self.split_datasets(filtered_data, ids)
            X_train, y_train, _ = splits["train"]
            
            # Train models
            models = self.train_models(X_train, y_train)
            
            # Collect predictions
            pred_scores = []
            for model_name in self.config["model_names"]:
                model = models[model_name]
                scores = []
                for split_name in self.config["split_names"]:
                    X, y, _ = splits[split_name]
                    if X is not None and not X.empty:
                        scores.append(self.predict(model, X))
                    else:
                        scores.append(np.zeros((0, 2)))
                pred_scores.append(tuple(scores))
            
            # Calculate weights based on test performance
            test_metrics = []
            for model_name, scores in zip(self.config["model_names"], pred_scores):
                test_scores = scores[1]  # Index 1 corresponds to test set
                if test_scores.shape[0] == 0:
                    acc = 0.0
                else:
                    test_y = splits["test"][1]
                    acc, _, _, _, _, _, _, _, _, _, _ = analysis_pred_binary(test_y, test_scores[:, 1])
                test_metrics.append(acc)
            
            # Apply softmax to get weights
            weights = self.softmax(np.array(test_metrics).reshape(1, -1))
            
            # Perform weighted fusion
            fusion_preds, fusion_scores = self.weight_fusion(weights, [pred_scores])
            
            # Evaluate fusion model
            fusion_results = []
            for split_name in self.config["split_names"]:
                split_idx = self.config["split_names"].index(split_name)
                X, y, _ = splits[split_name]
                
                if X is None or X.empty:
                    continue
                    
                fusion_score = fusion_scores[0][0][split_idx]
                if fusion_score.shape[0] == 0:
                    continue
                    
                metrics = analysis_pred_binary(y, fusion_score[:, 1])
                acc, auc, ci, tpr, tnr, ppv, npv, precision, recall, f1, thres = metrics
                
                ci_str = f"{auc:.3f}({ci[0]:.3f} - {ci[1]:.3f})"
                
                fusion_results.append((
                    "WeightFusion", acc, ci_str, tpr, tnr, ppv, npv, 
                    precision, recall, f1, thres, f"{self.config['label_column']}-{split_name}"
                ))
            
            # Create results DataFrame
            columns = [
                "model", "accuracy", "auc_ci", "sensitivity", "specificity", 
                "ppv", "npv", "precision", "recall", "f1", "threshold", "dataset"
            ]
            
            results_df = pd.DataFrame(fusion_results, columns=columns)
            
            # Save results
            results_path = self.dirs["results"] / "model_evaluation_results.csv"
            results_df.to_csv(results_path, index=False)
            self.logger.info("Results saved to %s", results_path)
            
            # Create visualizations
            self.visualize_results(results_df)
            
            self.logger.info("DHRPs pipeline completed successfully")
            return results_df
            
        except Exception as e:
            self.logger.error("Pipeline failed: %s", e)
            raise


def main():
    """Main function to run the DHRPs pipeline"""
    # Configuration
    config = {
        "seed": 0,
        "base_dir": "dhrps_pipeline",
        "data_path": "path/to/your/combined_data.csv",
        "label_column": "label",
        "split_column": "dataset_split",
        "split_names": ["train", "test", "validation", "external"],
        "var_threshold": 1.0,
        "model_names": [
            "SVM", "KNN", "RandomForest", "ExtraTrees", "XGBoost",
            "LightGBM", "NaiveBayes", "AdaBoost", "GradientBoosting",
            "LR", "MLP",
        ]
    }
    
    # Run pipeline
    pipeline = DHRPsPipeline(config)
    results = pipeline.run()
    
    return results


if __name__ == "__main__":
    try:
        results = main()
    except Exception as e:
        logging.error("Pipeline execution failed: %s", e)
        sys.exit(1)