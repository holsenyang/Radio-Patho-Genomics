# Medical Imaging Analysis Pipeline

## Overview
This code provides a pipeline for processing 3D medical images and extracting radiomic features using conventional radiomics approaches. The pipeline includes directory setup, data loading, image-mask pairing, and feature extraction functionalities.

## Requirements
### Python Packages
- onekey_algo
- pandas
- pathlib
- Other standard Python libraries

### Data
The pipeline expects medical imaging data in NIfTI format with corresponding Region of Interest (ROI) masks. The data should be organized in a directory structure where images and masks can be paired based on naming conventions.

## Usage
1. Set the data directory path in the main function call
2. Optionally specify a parameter file for custom feature extraction settings
3. Run the script to process images and extract features

## Output
The pipeline generates:
- Processed images in the 'img' directory
- Analysis results in the 'results' directory
- Extracted radiomic features in CSV format in the 'features' directory

## Customization
To customize feature extraction parameters, create a YAML configuration file following the example structure and pass its path to the main function.