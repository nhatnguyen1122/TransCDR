#!/usr/bin/env python3
"""
Training script for TransCDR using adapted splitter (Clean Version).

Optimized for:
- Regression only
- Pre-trained drug models (no ESPF needed)
- Your GDSC1 IC50 dataset

Usage:
    # Single fold
    python train_with_splitter_clean.py --scenario cold_cell --fold 1

    # 10-fold CV
    python train_with_splitter_clean.py --scenario cold_cell --cv
"""

import argparse
import os
import sys

# Import adapted modules (clean versions)
from splitter import DataSplitter
from DataEncoding_adapted_clean import DataEncoding
from model import TransCDR


def get_omics_dimensions():
    """
    Get omics data dimensions from your files.

    Returns:
        dict: input_dim_rna, input_dim_genetic, input_dim_mrna
    """
    import pandas as pd
    import json

    # RNA expression
    rna_data = pd.read_csv('data/omics_data/gene_expression/Cell_line_RMA_clean.txt', index_col=0)
    input_dim_rna = rna_data.shape[0]  # number of genes

    # Gene mutation
    with open('data/omics_data/gene_mutation/mutation/list_of_gene_mut.txt', 'r') as f:
        gene_names = [line.strip() for line in f]
    input_dim_genetic = len(gene_names)

    # DNA methylation
    mrna = pd.read_csv('data/omics_data/dna_methylation/gene_cell_matrix_promoter_filter_na.csv', index_col=0)
    input_dim_mrna = mrna.shape[0]  # number of CpG sites

    return {
        'input_dim_rna': input_dim_rna,
        'input_dim_genetic': input_dim_genetic,
        'input_dim_mrna': input_dim_mrna
    }


def get_config(args, omics_dims):
    """Create configuration dictionary for the model."""
    config = {
        'model_type': 'regression',
        'omics': args.omics,
        'input_dim_drug': args.input_dim_drug,
        'input_dim_rna': omics_dims['input_dim_rna'],
        'input_dim_genetic': omics_dims['input_dim_genetic'],
        'input_dim_mrna': omics_dims['input_dim_mrna'],
        'KG': '',
        'lr': args.lr,
        'decay': 0,
        'BATCH_SIZE': args.batch_size,
        'train_epoch': args.train_epoch,
        'pre_train': args.pre_train,
        'screening': 'None',
        'fusion_type': args.fusion_type,
        'drug_encoder': args.drug_encoder,
        'drug_model': args.drug_model,
        'modeldir': args.modeldir,
        'seq_model': args.seq_model,
        'graph_model': args.graph_model,
        'external_dataset': 'None'
    }
    return config


def train_single_fold(data_path, scenario, fold_idx, args):
    """Train and evaluate a single fold."""
    print(f"\n{'='*70}")
    print(f"Training Fold {fold_idx}")
    print(f"{'='*70}")

    # Step 1: Create data splits
    print("\n[1/5] Creating data splits...")
    splitter = DataSplitter(
        scenarios=scenario,
        data_path=data_path,
        n_folds=10,
        random_state=2022
    )

    train_idx, val_idx, test_idx = splitter.split_dataset(
        fold_idx=fold_idx,
        num_folds=10
    )

    print(f"  Train: {len(train_idx)} samples")
    print(f"  Val:   {len(val_idx)} samples")
    print(f"  Test:  {len(test_idx)} samples")

    # Step 2: Load data
    print("\n[2/5] Loading dataset...")
    df = splitter.CDR

    train_df = df.iloc[train_idx].copy()
    val_df = df.iloc[val_idx].copy()
    test_df = df.iloc[test_idx].copy()

    # Step 3: Prepare data (create Label column)
    print("\n[3/5] Preparing data...")
    data_encoding = DataEncoding()
    train_set, test_set, val_set = data_encoding.encode(train_df, test_df, val_df)

    # Step 4: Get omics dimensions
    print("\n[4/5] Loading omics dimensions...")
    try:
        omics_dims = get_omics_dimensions()
        print(f"  RNA: {omics_dims['input_dim_rna']} genes")
        print(f"  Mutation: {omics_dims['input_dim_genetic']} genes")
        print(f"  Methylation: {omics_dims['input_dim_mrna']} sites")
    except FileNotFoundError as e:
        print(f"Error loading omics data: {e}")
        print("\nRequired omics files:")
        print("  - data/omics_data/gene_expression/Cell_line_RMA_clean.txt")
        print("  - data/omics_data/gene_mutation/mutation/list_of_gene_mut.txt")
        print("  - data/omics_data/dna_methylation/gene_cell_matrix_promoter_filter_na.csv")
        sys.exit(1)

    # Create config
    config = get_config(args, omics_dims)
    config['modeldir'] = os.path.join(args.modeldir, f'fold{fold_idx}')

    # Step 5: Train model
    print(f"\n[5/5] Training model...")
    print(f"  Model dir: {config['modeldir']}")
    print(f"  Drug model: {config['drug_model']} (pre-trained: {config['pre_train']})")
    print(f"  Omics: {config['omics']}")
    print(f"  Fusion: {config['fusion_type']}")
    print(f"  LR: {config['lr']}, Batch: {config['BATCH_SIZE']}, Epochs: {config['train_epoch']}")

    net = TransCDR(**config)
    net.train(train_drug=train_set, test_drug=test_set, val_drug=val_set)
    net.save_model()

    print(f"\n✓ Fold {fold_idx} training complete!")

    # Show test results
    metrics_file = os.path.join(config['modeldir'], 'test_markdowntable.txt')
    if os.path.exists(metrics_file):
        with open(metrics_file, 'r') as f:
            print(f"\nTest Results:")
            print(f.read())

    return config['modeldir']


def run_cv(data_path, scenario, args):
    """Run 10-fold cross-validation."""
    print("="*70)
    print("10-FOLD CROSS-VALIDATION")
    print("="*70)

    results = []

    for fold in range(1, 11):
        try:
            modeldir = train_single_fold(data_path, scenario, fold, args)
            results.append(modeldir)
        except Exception as e:
            print(f"\n✗ Error in fold {fold}: {e}")
            import traceback
            traceback.print_exc()
            continue

    print("\n" + "="*70)
    print("CROSS-VALIDATION COMPLETE")
    print("="*70)
    print(f"Successfully trained {len(results)}/10 folds")


def main():
    parser = argparse.ArgumentParser(description='Train TransCDR (Clean Version - No ESPF)')

    # Data
    parser.add_argument('--data_path', type=str,
                        default='./data/Drug Screening - IC50s/GDSC1_IC50_with_pubchem_and_smiles_filtered_by_mut_rna_mrna_Dec29.xlsx',
                        help='Path to IC50 dataset')
    parser.add_argument('--scenario', type=str, required=True,
                        choices=['warm_start', 'cold_drug', 'cold_cell', 'cold_scaffold', 'final'],
                        help='Splitting scenario')
    parser.add_argument('--fold', type=int, default=1,
                        help='Fold number (1-10)')
    parser.add_argument('--cv', action='store_true',
                        help='Run 10-fold cross-validation')

    # Model architecture
    parser.add_argument('--omics', type=str, default='expr + mutation + methylation',
                        help='Omics data to use')
    parser.add_argument('--drug_encoder', type=str, default='None',
                        help='Drug encoder (use None for pre-trained)')
    parser.add_argument('--drug_model', type=str, default='sequence + graph + FP',
                        help='Pre-trained drug model combination')
    parser.add_argument('--fusion_type', type=str, default='encoder',
                        choices=['encoder', 'decoder', 'concat'],
                        help='Fusion type')
    parser.add_argument('--input_dim_drug', type=int, default=2092,
                        help='Drug input dim (seq:768 + graph:300 + FP:1024 = 2092)')

    # Pre-trained models
    parser.add_argument('--pre_train', type=str, default='True',
                        help='Use pre-trained models (recommended)')
    parser.add_argument('--seq_model', type=str, default='seyonec/PubChem10M_SMILES_BPE_450k',
                        help='Pre-trained SMILES model')
    parser.add_argument('--graph_model', type=str, default='gin_supervised_masking',
                        help='Pre-trained graph model')

    # Training
    parser.add_argument('--lr', type=float, default=1e-5,
                        help='Learning rate')
    parser.add_argument('--batch_size', type=int, default=64,
                        help='Batch size')
    parser.add_argument('--train_epoch', type=int, default=100,
                        help='Training epochs')

    # Output
    parser.add_argument('--modeldir', type=str, default='./result/clean',
                        help='Model output directory')

    args = parser.parse_args()
    os.makedirs(args.modeldir, exist_ok=True)

    if args.cv:
        run_cv(args.data_path, args.scenario, args)
    else:
        train_single_fold(args.data_path, args.scenario, args.fold, args)


if __name__ == "__main__":
    main()
