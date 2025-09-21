import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import onekey_algo.custom.components as okcomp
from onekey_algo.custom.components.metrics import analysis_pred_binary
from onekey_algo.custom.components import stats


def evaluate_models(data_path, model_names, task_column, subset_column, output_dir="results"):
    """
    Evaluate multiple ML models across different data subsets
    
    Args:
        data_path: Path to CSV file with prediction scores
        model_names: List of model names to evaluate
        task_column: Name of the target variable column
        subset_column: Name of the column indicating data subsets
        output_dir: Directory to save output files
    """
    # Setup
    os.makedirs(output_dir, exist_ok=True)
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['pdf.fonttype'] = 42
    
    # Load data
    data = pd.read_csv(data_path)
    data[model_names] = 1 - data[model_names]  # Invert scores
    
    # Initialize results storage
    all_metrics = {}
    all_hosmer = {}
    
    # Get unique subsets
    subsets = data[subset_column].unique()
    
    # Evaluate each subset
    for subset in subsets:
        subset_data = data[data[subset_column] == subset]
        
        # Prepare data
        gt = [np.array(subset_data[task_column]) for _ in model_names]
        pred_scores = [np.array(subset_data[mn]) for mn in model_names]
        
        # Generate ROC curves
        plt.figure()
        okcomp.comp1.draw_roc(gt, pred_scores, labels=model_names, title=f"{subset} Model AUC")
        plt.savefig(f'{output_dir}/roc_{subset}.pdf', bbox_inches='tight', format='pdf')
        plt.close()
        
        # Calculate metrics
        metrics = []
        for mname, y, score in zip(model_names, gt, pred_scores):
            acc, auc, ci, tpr, tnr, ppv, npv, precision, recall, f1, thres = analysis_pred_binary(y, score)
            ci_fmt = f"{ci[0]:.4f} - {ci[1]:.4f}"
            metrics.append((mname, acc, auc, ci_fmt, tpr, tnr, ppv, npv, 
                           precision, recall, f1, thres, subset))
            
        metrics_df = pd.DataFrame(metrics, columns=[
            'Model', 'Accuracy', 'AUC', '95% CI', 'Sensitivity', 'Specificity',
            'PPV', 'NPV', 'Precision', 'Recall', 'F1', 'Threshold', 'Cohort'
        ])
        
        all_metrics[subset] = metrics_df
        
        # Generate DCA plot
        plt.figure()
        okcomp.comp1.plot_DCA([subset_data[mn] for mn in model_names], subset_data[task_column],
                             title=f'{subset} Model DCA', labels=model_names, y_min=-0.15)
        plt.savefig(f'{output_dir}/dca_{subset}.pdf', bbox_inches='tight', format='pdf')
        plt.close()
        
        # Hosmer-Lemeshow test
        hosmer = [stats.hosmer_lemeshow_test(y_true, y_pred, bins=15)
                 for _, y_true, y_pred in zip(model_names, gt, pred_scores)]
        all_hosmer[subset] = pd.DataFrame([hosmer], columns=model_names)
    
    # Export results
    with pd.ExcelWriter(f"{output_dir}/evaluation_results.xlsx") as writer:
        for subset in all_metrics:
            all_metrics[subset].to_excel(writer, sheet_name=f"{subset}_metrics", index=False)
            all_hosmer[subset].to_excel(writer, sheet_name=f"{subset}_hosmer", index=False)
    
    print(f"Evaluation complete. Results saved to {output_dir}")


if __name__ == "__main__":
    # Configuration - replace with your actual values
    data_path = "path/to/your/data.csv"
    model_names = ["model_1", "model_2", "model_3", "model_4", "model_5"]
    task_column = "target"
    subset_column = "dataset_split"
    
    evaluate_models(data_path, model_names, task_column, subset_column)