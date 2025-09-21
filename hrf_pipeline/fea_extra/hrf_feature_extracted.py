import os
import warnings
import pandas as pd
from pathlib import Path
from onekey_algo import OnekeyDS as okds
from onekey_algo.custom.components.Radiology import (
    diagnose_3d_image_mask_settings, 
    get_image_mask_from_dir,
    ConventionalRadiomics
)

# Suppress warnings
warnings.filterwarnings("ignore")

# Set environment variable
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

def setup_directories():
    """Create necessary directories"""
    for dir_name in ['img', 'results', 'features']:
        os.makedirs(dir_name, exist_ok=True)

def load_data(data_dir):
    """Load image and mask data from directory"""
    images, masks = get_image_mask_from_dir(data_dir, images='NII', masks='ROI')
    diagnose_3d_image_mask_settings(images, masks)
    print(f'Found {len(images)} samples.')
    return images, masks

def extract_features(images, masks, param_file=None, force_extract=False):
    """Extract radiomics features"""
    features_path = 'features/rad_features.csv'
    
    if not force_extract and os.path.exists(features_path):
        return pd.read_csv(features_path, header=0)
    
    radiomics = ConventionalRadiomics(param_file, correctMask=True)
    radiomics.extract(images, masks)
    feature_data = radiomics.get_label_data_frame(label=1)
    feature_data.columns = [c.replace('-', '_') for c in feature_data.columns]
    feature_data.to_csv(features_path, header=True, index=False)
    return feature_data

def main(data_directory, parameter_file=None):
    """Main pipeline function"""
    setup_directories()
    images, masks = load_data(data_directory)
    feature_data = extract_features(images, masks, parameter_file)
    return feature_data

if __name__ == "__main__":
    # Example usage
    data_dir = "path/to/data/directory"  # Replace with actual data path
    param_file = "path/to/parameter/file.yaml"  # Optional parameter file
    results = main(data_dir, param_file)
    print(results.head())