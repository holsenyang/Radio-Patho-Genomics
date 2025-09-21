# Classification Pipeline with TF-IDF Feature Extraction

## Overview
This code provides a complete pipeline for training a classification model and extracting TF-IDF features from the results. The pipeline includes model training with configurable parameters and post-processing of results to generate TF-IDF based features.

## Requirements
### Python Packages
- onekey_algo
- pandas
- Other standard Python libraries

### Data
The pipeline expects:
1. Training and validation data in CSV format
2. Labels file in text format
3. Result file from model training in tab-separated format

## Usage
1. Set the appropriate file paths in the main function
2. Configure model parameters as needed
3. Run the script to train the model and generate TF-IDF features

## Output
The pipeline generates:
- Trained model files
- TF-IDF features in CSV format (tfidf_features.csv)

## Configuration
Modify the model_params dictionary to adjust training parameters such as:
- Model architecture (model_name)
- Batch size and epochs
- Learning rate and optimizer
- GPU configuration