# Model Evaluation Pipeline

## Overview
A streamlined pipeline for evaluating multiple machine learning models across different datasets.

## Features
- ROC curve analysis
- Performance metrics calculation
- Decision Curve Analysis (DCA)
- Hosmer-Lemeshow goodness-of-fit test
- Excel export of results

## Requirements
- pandas
- numpy
- matplotlib
- scikit-learn
- onekey_algo (custom package)

## Usage
1. Update the configuration in the main section:
   - data_path: Path to your CSV file
   - model_names: List of your model names
   - task_column: Name of your target variable column
   - subset_column: Name of your dataset split column

2. Run the script:
```bash
python model_evaluation.py


Input Format
CSV file should contain:

A column indicating data subsets

A column with ground truth labels

Columns with prediction scores for each model

Output
ROC curves (PDF)

DCA plots (PDF)

Excel file with metrics and statistical test results