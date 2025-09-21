# DHRPs Pipeline

## Overview
A comprehensive machine learning pipeline for multi-modal feature integration, model training, ensemble fusion, and evaluation. This pipeline combines features from multiple modalities and uses weighted fusion to create an ensemble model.

## Features
- Multi-modal feature integration
- Variance threshold filtering
- Multiple model training
- Weighted ensemble fusion
- Comprehensive evaluation
- Visualization of results

## Requirements
### Python Packages
- pandas>=1.3.0
- numpy>=1.20.0
- matplotlib>=3.3.0
- seaborn>=0.11.0
- scikit-learn>=1.0.0
- joblib>=1.0.0
- onekey_algo (custom package)

### Installation
```bash
pip install pandas numpy matplotlib seaborn scikit-learn joblib


## Weighted Voting Ensemble Strategy

The DHRPs pipeline employs a weighted voting ensemble strategy that combines predictions from multiple base models. Unlike simple majority voting, this approach assigns different weights to each model based on their individual performance, giving more influence to better-performing models.

### Weight Calculation Process

#### Specific Weight Calculation (Formula):

To convert model accuracy into usable weights and perform normalization (so that the sum of all weights equals 1), we adopt the Softmax function. This function amplifies the differences between high-performance and low-performance models, ensuring that models with superior performance play a dominant role in the final decision.

- The specific calculation consists of two steps:

a) Obtain the accuracy of each model on the test set: Accuracy_i, where i denotes the i-th model.

b) Use the Softmax function to calculate the normalized weights.

                                  Weight_i=exp⁡(Accuracy_i )/(∑_(j=1)^11▒exp⁡(Accuracy_j )   )

#### Ensemble Prediction Process:

For a new sample, the final ensemble prediction probability is not simply a "one vote per model" approach. Instead, it is the weighted average of the prediction probabilities from all models. The "voting power" of each model in the final result is its Weight_i, which is calculated using the formula described above.

- The probability P_final of predicting the sample as high-grade IPA is calculated as:

                                     P_final=∑_(i=1)^11 (Weight_i×Prob_i )

Here, Prob_i is the probability predicted by the i-th model that the sample belongs to the high-grade IPA category.

### Advantages of Weighted Voting

- Performance Optimization: Better-performing models have greater influence on the final decision

- Robustness: Reduces the impact of poorly performing models

- Flexibility: Can incorporate diverse model types with different strengths

- Interpretability: Weight values provide insight into model contributions



