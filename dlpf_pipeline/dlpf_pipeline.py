#!/usr/bin/env python3
"""
DLPF Pipeline: Data Preparation → Feature Engineering → Model Training → Evaluation → Visualization
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
from sklearn import linear_model
from sklearn.ensemble import RandomForestRegressor
from boruta import BorutaPy

# Import custom components
import onekey_algo.custom.components as okcomp
from onekey_algo.custom.components.comp1 import (
    normalize_df,
    draw_matrix,
    select_feature,
    create_clf_model,
    plot_feature_importance,
)
from onekey_algo.custom.components.metrics import analysis_pred_binary
from onekey_algo.custom.components import stats


class DLPFPipeline:
    """End-to-end pipeline for data processing, feature selection, and model training"""
    
    def __init__(self, config: Dict):
        """
        Initialize the DLPF pipeline
        
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
        self.logger = logging.getLogger("DLPF")
        
    def setup_directories(self):
        """Create necessary directories for output"""
        base_dir = Path(self.config.get("base_dir", "DLPF"))
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
        
    def load_data(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
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
        
        self.logger.info("Label distribution:\n%s", data[label_col].value_counts())
        
        # Separate CPTAC data
        cptac_data = data[data["split_column"] == "CPTAC"].copy()
        main_data = data[data["split_column"].isin(["train", "test", "validation"])].copy()
        
        # Normalize data
        norm_data = normalize_df(main_data, not_norm=[label_col, "split_column"]).dropna(axis=1)
        norm_cptac = normalize_df(cptac_data, not_norm=[label_col, "split_column"]).copy()
        
        # Combine datasets
        cols = norm_data.columns.tolist()
        norm_cptac = norm_cptac.reindex(columns=cols)
        combined_data = pd.concat([norm_data, norm_cptac], axis=0).fillna(0)
        
        return combined_data, ids
        
    def variance_filter(self, data: pd.DataFrame) -> pd.DataFrame:
        """Apply variance threshold filtering"""
        features = data.columns[2:]
        train_data = data[data["split_column"] == "train"].copy()
        
        selector = VarianceThreshold(threshold=self.config["var_threshold"])
        selector.fit_transform(train_data[features])
        
        kept_features = features[selector.get_support()].tolist()
        self.logger.info("Variance filtering kept %d features", len(kept_features))
        
        return data[["label", "split_column"] + kept_features]
        
    def correlation_filter(self, data: pd.DataFrame) -> List[str]:
        """Select features based on correlation"""
        train_data = data[data["split_column"] == "train"].copy()
        features = [c for c in train_data.columns if c not in [self.config["label_column"]]]
        
        corr_matrix = train_data[features].corr("spearman")
        
        # Visualize correlation matrix if not too large
        if len(features) < 100:
            plt.figure(figsize=(12, 10))
            draw_matrix(corr_matrix, annot=False, cmap="YlOrRd", cbar=False)
            plt.savefig(self.dirs["img"] / "feature_correlation.pdf", 
                       bbox_inches="tight", format="pdf")
            plt.close()
            
        selected = select_feature(
            corr_matrix, 
            threshold=self.config["corr_threshold"], 
            topn=self.config["corr_topn"],
            verbose=False
        )
        
        self.logger.info("Correlation filtering kept %d features", len(selected))
        return selected
        
    def boruta_selection(self, data: pd.DataFrame) -> List[str]:
        """Apply Boruta feature selection"""
        train_data = data[data["split_column"] == "train"].copy()
        features = train_data.columns[2:]
        
        if self.config["run_boruta"]:
            rf = RandomForestRegressor(n_jobs=-1, max_depth=5, 
                                     random_state=self.config["seed"])
            boruta = BorutaPy(
                rf,
                n_estimators=500,
                verbose=2,
                random_state=self.config["seed"],
                max_iter=100
            )
            
            boruta.fit(
                np.array(train_data[features]),
                np.array(train_data[self.config["label_column"]]).ravel()
            )
            
            # Get selected features
            selected = features[boruta.support_].tolist()
            
            # Optionally include tentative features
            if self.config.get("keep_tentative", False):
                tentative = features[boruta.support_weak_].tolist()
                selected = list(sorted(set(selected + tentative)))
                self.logger.info("Included %d tentative features", len(tentative))
                
            # Save selected features
            feature_file = self.dirs["features"] / "boruta_selected_features.csv"
            pd.Series(selected, name="feature").to_csv(feature_file, index=False)
        else:
            # Load pre-selected features
            feature_file = self.dirs["features"] / "boruta_selected_features.csv"
            if not feature_file.exists():
                raise FileNotFoundError(f"Boruta feature file not found: {feature_file}")
            selected = pd.read_csv(feature_file)["feature"].tolist()
            
        self.logger.info("Boruta selected %d features", len(selected))
        return selected
        
    def lasso_selection(self, X: pd.DataFrame, y: pd.DataFrame) -> Tuple[List[str], List[Tuple[str, float]]]:
        """Apply Lasso feature selection"""
        # Find optimal alpha
        alpha = okcomp.comp1.lasso_cv_coefs(
            X, y, 
            column_names=None, 
            alpha_logmin=self.config.get("lasso_alpha_min", -5)
        )
        plt.savefig(self.dirs["img"] / "lasso_feature_selection.pdf", 
                   bbox_inches="tight", format="pdf")
        plt.close()
        
        # Fit Lasso model
        selected = []
        coef_pairs = []
        feature_names = list(X.columns)
        
        clf = linear_model.Lasso(alpha=alpha, random_state=self.config["seed"])
        clf.fit(X, y[self.config["label_column"]])
        
        # Get features with significant coefficients
        for feature, coef in zip(feature_names, clf.coef_):
            if abs(coef) > self.config["lasso_coef_threshold"]:
                selected.append(feature)
                coef_pairs.append((feature, float(coef)))
                
        # Visualize feature coefficients
        coef_df = pd.DataFrame(coef_pairs, columns=["feature", "coefficient"])
        coef_df.plot(x="feature", y="coefficient", kind="barh", 
                    figsize=(8, max(4, len(coef_df) * 0.25)))
        plt.savefig(self.dirs["img"] / "lasso_feature_weights.pdf", 
                   bbox_inches="tight", format="pdf")
        plt.close()
        
        # Save coefficients
        coef_df.assign(abs_value=lambda d: d["coefficient"].abs()).to_csv(
            self.dirs["features"] / "lasso_feature_weights.csv", index=False)
            
        self.logger.info("Lasso selected %d features", len(selected))
        return selected, coef_pairs
        
    def split_datasets(self, data: pd.DataFrame, ids: pd.Series) -> Dict[str, Tuple]:
        """Split data into training and validation sets"""
        splits = {}
        
        for split_name in self.config["split_names"]:
            split_data = data[data["split_column"] == split_name].reset_index(drop=True)
            split_ids = ids[split_data.index]
            
            X = split_data.drop(columns=[self.config["label_column"], "split_column"])
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
        
    def predict(self, models: Dict[str, object], X: pd.DataFrame) -> np.ndarray:
        """Generate predictions from models"""
        if X is None or len(X) == 0:
            return np.zeros((0, 2))
            
        if hasattr(models, "predict_proba"):
            proba = models.predict_proba(X)
            if proba.ndim == 1:
                proba = np.vstack([1 - proba, proba]).T
            elif proba.shape[1] == 1:
                proba = np.hstack([1 - proba, proba])
            return proba
            
        if hasattr(models, "decision_function"):
            scores = models.decision_function(X)
            scores = (scores - scores.min()) / (scores.max() - scores.min() + 1e-12)
            return np.vstack([1 - scores, scores]).T
            
        preds = models.predict(X)
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
        
        # Accuracy line plot
        plt.figure(figsize=(12, 8))
        sns.lineplot(x="model", y="accuracy", data=results, hue="dataset", marker="o")
        plt.title("Model Accuracy by Dataset")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        plt.savefig(self.dirs["img"] / "model_accuracy_line.pdf", 
                   bbox_inches="tight", format="pdf")
        plt.close()
        
    def run(self):
        """Execute the full DLPF pipeline"""
        self.logger.info("Starting DLPF pipeline")
        
        try:
            # Load data
            data = self.load_data()
            
            # Preprocess data
            processed_data, ids = self.preprocess_data(data)
            
            # Apply variance filtering
            var_filtered = self.variance_filter(processed_data)
            
            # Apply correlation filtering
            corr_features = self.correlation_filter(var_filtered)
            corr_filtered = processed_data[["label", "split_column"] + corr_features]
            
            # Apply Boruta feature selection
            boruta_features = self.boruta_selection(corr_filtered)
            boruta_filtered = corr_filtered[["label", "split_column"] + boruta_features]
            
            # Split data
            splits = self.split_datasets(boruta_filtered, ids)
            X_train, y_train, _ = splits["train"]
            
            # Apply Lasso feature selection
            lasso_features, _ = self.lasso_selection(X_train, y_train)
            
            # Filter all splits with selected features
            for split_name in splits:
                X, y, ids = splits[split_name]
                if X is not None and not X.empty:
                    splits[split_name] = (X[lasso_features], y, ids)
            
            # Train models
            X_train_final, y_train_final, _ = splits["train"]
            models = self.train_models(X_train_final, y_train_final)
            
            # Evaluate models
            results = self.evaluate_models(models, splits)
            
            # Save results
            results_path = self.dirs["results"] / "model_evaluation_results.csv"
            results.to_csv(results_path, index=False)
            self.logger.info("Results saved to %s", results_path)
            
            # Create visualizations
            self.visualize_results(results)
            
            self.logger.info("DLPF pipeline completed successfully")
            return results, X_train_final
            
        except Exception as e:
            self.logger.error("Pipeline failed: %s", e)
            raise


def main():
    """Main function to run the DLPF pipeline"""
    # Configuration
    config = {
        "seed": 0,
        "base_dir": "dlpf_pipeline",
        "data_path": "path/to/your/data.csv",
        "label_column": "label",
        "split_column": "dataset_split",
        "split_names": ["train", "test", "validation", "external"],
        "var_threshold": 1.0,
        "corr_threshold": 0.8,
        "corr_topn": 32,
        "run_boruta": True,
        "keep_tentative": False,
        "lasso_coef_threshold": 1e-6,
        "lasso_alpha_min": -5,
        "model_names": [
            "SVM", "KNN", "RandomForest", "ExtraTrees", "XGBoost",
            "LightGBM", "NaiveBayes", "AdaBoost", "GradientBoosting",
            "LR", "MLP",
        ]
    }
    
    # Run pipeline
    pipeline = DLPFPipeline(config)
    results, features = pipeline.run()
    
    return results, features


if __name__ == "__main__":
    try:
        results, features = main()
    except Exception as e:
        logging.error("Pipeline execution failed: %s", e)
        sys.exit(1)