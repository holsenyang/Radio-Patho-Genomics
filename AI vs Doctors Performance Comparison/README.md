# AI vs Doctors Performance Comparison

This repository contains Python code for comparing the performance of AI models against human doctors (pathologists) in binary classification tasks, particularly in medical diagnosis scenarios.

## Overview

The code provides a comprehensive framework for:
- Calculating performance metrics (AUC, accuracy, sensitivity, specificity, etc.)
- Visualizing ROC curves and performance points
- Comparing AI models with human experts using statistical tests
- Generating confidence intervals using bootstrap methods
- Performing McNemar tests for sensitivity and specificity comparisons

## Features

1. **Performance Metrics Calculation**:
   - AUC with confidence intervals (bootstrap method)
   - Accuracy, sensitivity, specificity, PPV, NPV, F1-score
   - Optimal threshold selection using Youden's index

2. **Statistical Comparisons**:
   - DeLong test for AUC comparisons
   - McNemar test for sensitivity and specificity comparisons
   - Bootstrap confidence intervals for all metrics

3. **Visualization**:
   - ROC curves for probabilistic models
   - Scatter points for human evaluators
   - Error bars for average human performance

4. **Flexible Input**:
   - Supports both probabilistic outputs (AI models) and binary decisions (doctors)
   - Handles multiple human evaluators

##Requirements
Python 3.6+

scikit-learn

matplotlib

numpy

pandas

scipy