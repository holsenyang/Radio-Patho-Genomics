# Feature Extraction Pipeline

## Overview
This code provides a pipeline for extracting deep learning features from images using both pretrained models and custom trained models. The pipeline supports feature extraction from specific layers and processing of the resulting feature vectors.

## Requirements
### Python Packages
- onekey_algo
- pandas
- Other standard Python libraries

### Data
The pipeline expects image data in PNG or JPG format organized in a directory structure.

## Usage
1. Set the appropriate directory and file paths in the main function
2. Configure model parameters (pretrained model name or custom model path)
3. Specify the feature layer for extraction
4. Run the script to extract and process features

## Output
The pipeline generates:
- Feature vectors in CSV format for both pretrained and custom models
- Processed feature files with appropriate column names

## Configuration
Modify the following parameters in the main function:
- data_directory: Path to directory containing images
- pretrained_model_name: Name of pretrained model to use
- custom_model_path: Path to custom trained model
- feature_layer: Name of the layer from which to extract features