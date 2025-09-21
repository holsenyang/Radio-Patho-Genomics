# DLPF Pipeline

## Overview
A comprehensive machine learning pipeline for Data preparation, Feature engineering, Model training, and Evaluation. The pipeline includes multiple feature selection methods and supports various classification algorithms.

## Features
- Data loading and preprocessing
- Multiple feature selection methods:
  - Variance threshold filtering
  - Correlation-based selection
  - Boruta feature selection
  - Lasso regularization
- Multiple model training and evaluation
- Comprehensive visualization of results
- Reproducible results through random seed control

## Requirements
### Python Packages
- pandas>=1.3.0
- numpy>=1.20.0
- matplotlib>=3.3.0
- seaborn>=0.11.0
- scikit-learn>=1.0.0
- boruta>=0.3.0
- joblib>=1.0.0
- onekey_algo (custom package)

### Installation
```bash
pip install pandas numpy matplotlib seaborn scikit-learn boruta-py joblib