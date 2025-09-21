import os
import pandas as pd
from collections import namedtuple
from typing import Union, List
from onekey_algo.classification.run_classification import main as clf_main
from onekey_algo.custom.utils import key2

def train_model(train_file, valid_file, labels_file, data_pattern, model_params):
    """Train a classification model with given parameters"""
    # Set up parameters
    params = dict(
        train=train_file,
        valid=valid_file,
        labels_file=labels_file,
        data_pattern=data_pattern,
        **model_params
    )
    
    # Train model
    Args = namedtuple("Args", params)
    clf_main(Args(**params))

def process_results(result_file, group_column='group', corpus_columns='prob'):
    """Process model results and generate TF-IDF features"""
    # Read results
    log = pd.read_csv(result_file, sep='\t', names=['fname', 'prob', 'pred', 'gt'])
    log[[group_column]] = log[['fname']].applymap(lambda x: x.split('_')[0])
    
    # Generate TF-IDF features
    results = key2.key2tfidf(log, group_column=group_column, corpus_columns=corpus_columns)
    results.to_csv('tfidf_features.csv', header=True, index=True)
    
    return results

def main():
    """Main function to run the complete pipeline"""
    # Configuration parameters
    save_dir = "path/to/save/directory"
    train_file = "path/to/train/data.csv"
    valid_file = "path/to/validation/data.csv"
    labels_file = "path/to/labels/file.txt"
    data_pattern = os.path.join(save_dir, "normalized")
    result_file = "path/to/result/file.txt"
    
    # Model parameters
    model_params = {
        'j': 4,
        'max2use': None,
        'val_max2use': None,
        'batch_balance': False,
        'normalize_method': 'imagenet',
        'model_name': 'resnet18',
        'gpus': [0],
        'batch_size': 32,
        'epochs': 30,
        'init_lr': 0.01,
        'optimizer': 'sgd',
        'retrain': None,
        'model_root': '.',
        'iters_start': 0,
        'iters_verbose': 1,
        'save_per_epoch': False,
        'pretrained': True
    }
    
    # Run training
    train_model(train_file, valid_file, labels_file, data_pattern, model_params)
    
    # Process results
    tfidf_results = process_results(result_file)
    
    return tfidf_results

if __name__ == "__main__":
    results = main()
    print(results.head())