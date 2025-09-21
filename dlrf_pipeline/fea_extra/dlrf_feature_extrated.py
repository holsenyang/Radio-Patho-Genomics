import os
import pandas as pd
from functools import partial
from onekey_algo.custom.components.comp2 import (
    extract, 
    print_feature_hook, 
    reg_hook_on_module, 
    init_from_model, 
    init_from_onekey
)

def get_image_samples(directory_path, extensions=('.png', '.jpg')):
    """Get image samples from directory"""
    directory = os.path.expanduser(directory_path)
    samples = [os.path.join(directory, p) for p in os.listdir(directory) 
              if p.endswith(extensions)]
    return samples

def extract_features(samples, model_init_func, model_init_args, 
                    feature_name, output_file):
    """Extract features using specified model and hook"""
    # Initialize model
    model, transformer, device = model_init_func(**model_init_args)
    
    # Register hook and extract features
    with open(output_file, 'w') as outfile:
        hook = partial(print_feature_hook, fp=outfile)
        reg_hook_on_module(feature_name, model, hook)
        extract(samples, model, transformer, device, fp=outfile)
    
    return output_file

def process_feature_file(input_file, output_file, prefix='RES_'):
    """Process feature file by adding column names"""
    features = pd.read_csv(input_file, header=None)
    features.columns = ['ID'] + [prefix + str(col) for col in features.columns[1:]]
    features.to_csv(output_file, header=True, index=False)
    return features

def main():
    """Main function to run feature extraction pipeline"""
    # Configuration parameters
    data_directory = "path/to/data/directory"
    pretrained_model_name = "resnet18"
    custom_model_path = "path/to/custom/model"
    feature_layer = "avgpool"
    
    # Get image samples
    samples = get_image_samples(data_directory)
    
    # Option 1: Extract features using pretrained model
    pretrained_output = "pretrained_features.csv"
    extract_features(
        samples, 
        init_from_model, 
        {"model_name": pretrained_model_name},
        feature_layer,
        pretrained_output
    )
    pretrained_features = process_feature_file(pretrained_output, pretrained_output)
    
    # Option 2: Extract features using custom trained model
    custom_output = "custom_features.csv"
    extract_features(
        samples, 
        init_from_onekey, 
        {"model_path": custom_model_path},
        feature_layer,
        custom_output
    )
    custom_features = process_feature_file(custom_output, custom_output)
    
    return pretrained_features, custom_features

if __name__ == "__main__":
    pretrained_feats, custom_feats = main()
    print("Pretrained features shape:", pretrained_feats.shape)
    print("Custom features shape:", custom_feats.shape)