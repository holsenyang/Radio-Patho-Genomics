# Multi-omics and Machine Learning Pipelines Collection

A comprehensive collection of bioinformatics and machine learning pipelines for multi-omics data analysis, model evaluation, and performance comparison.

## 📁 Project Structure

### 1. Multi-omics Data Analysis Pipeline

A comprehensive bioinformatics pipeline for multi-omics data analysis, including co-expression network construction, functional enrichment, and pathway activity analysis.

**Key Features:**
- Data preprocessing and quality control
- Co-expression network construction (WGCNA)
- Hub gene identification
- Module preservation analysis
- Functional enrichment analysis (GO, KEGG, Reactome, etc.)
- Pathway activity analysis (GSVA)

**Integrated Datasets:**
- Batch Data (internal cohort)
- TCGA Data 
- CPTAC Data

### 2. AI vs Doctors Performance Comparison

Python framework for comparing AI models against human doctors (pathologists) in binary classification tasks for medical diagnosis scenarios.

**Key Features:**
- Performance metrics calculation (AUC, accuracy, sensitivity, specificity, etc.)
- ROC curve visualization
- Statistical comparisons (DeLong test, McNemar test)
- Bootstrap confidence intervals
- Support for both probabilistic outputs and binary decisions

### 3. DHRPs Pipeline (Weighted Voting Ensemble)

Comprehensive machine learning pipeline for multi-modal feature integration, model training, ensemble fusion, and evaluation.

**Core Features:**
- Multi-modal feature integration
- Variance threshold filtering
- Multiple model training
- **Weighted voting ensemble strategy**
- Comprehensive evaluation and visualization

**Weighted Voting Ensemble Strategy:**
- Softmax weight calculation based on model accuracy
- Weighted average prediction probability
- Performance optimization and robustness improvement

### 4. DLPF Pipeline

Comprehensive machine learning pipeline for data preparation, feature engineering, model training, and evaluation.

**Feature Selection Methods:**
- Variance threshold filtering
- Correlation-based selection
- Boruta feature selection
- Lasso regularization

### 5. DLRF Pipeline

Another comprehensive machine learning pipeline focused on data preparation, feature engineering, and model training.

**Features:**
- Data loading and preprocessing
- Multiple feature selection methods
- Multi-model training and evaluation
- Results visualization

### 6. PATH Pipeline (HPF Pipeline)

Comprehensive machine learning pipeline supporting multiple feature selection methods and classification algorithms.

### 7. RAD Pipeline (HRF Pipeline)

Machine learning pipeline for data preparation, feature engineering, model training, and evaluation.

### 8. Model Evaluation Pipeline

Streamlined pipeline for evaluating multiple machine learning models across different datasets.

**Evaluation Features:**
- ROC curve analysis
- Performance metrics calculation
- Decision Curve Analysis (DCA)
- Hosmer-Lemeshow goodness-of-fit test
- Excel export of results

## 🛠️ Technical Stack

### Python Package Requirements
```bash
pip install pandas numpy matplotlib seaborn scikit-learn boruta-py joblib
