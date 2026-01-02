# Splitter Usage Guide

## Overview

The `splitter.py` module provides data splitting functionality for drug response prediction tasks. It has been adapted for the GDSC1 IC50 dataset with the following specifications:

- **Regression-only** (classification code removed)
- **Column requirements**: DRUG_ID, COSMIC_ID, smiles, LN_IC50
- **Cell identifier handling**: Works with COSMIC_ID (verified 1-to-1 mapping with CELL_LINE_NAME and SANGER_MODEL_ID)

## Dataset Information

- **File**: `./data/Drug Screening - IC50s/GDSC1_IC50_with_pubchem_and_smiles_filtered_by_mut_rna_mrna_Dec29.xlsx`
- **Samples**: 168,449 drug-cell combinations
- **Drugs**: 238 unique compounds
- **Cell lines**: 835 unique cell lines
- **1-to-1 mapping verified** between COSMIC_ID ↔ CELL_LINE_NAME ↔ SANGER_MODEL_ID

## Splitting Scenarios

### 1. **warm_start**
Random K-fold split - both drugs and cells can appear in train and test sets.

```python
from splitter import DataSplitter

splitter = DataSplitter(
    scenarios="warm_start",
    data_path="./data/Drug Screening - IC50s/GDSC1_IC50_with_pubchem_and_smiles_filtered_by_mut_rna_mrna_Dec29.xlsx",
    n_folds=10,
    random_state=2022
)

train_idx, val_idx, test_idx = splitter.split_dataset(fold_idx=1, num_folds=10)
```

### 2. **cold_drug**
Split by drug - no drug appears in both train and test sets.

```python
splitter = DataSplitter(
    scenarios="cold_drug",
    data_path="./data/Drug Screening - IC50s/GDSC1_IC50_with_pubchem_and_smiles_filtered_by_mut_rna_mrna_Dec29.xlsx"
)

train_idx, val_idx, test_idx = splitter.split_dataset(fold_idx=1, num_folds=10)
```

### 3. **cold_cell**
Split by cell line - no cell appears in both train and test sets.

**Important**: Splits on `COSMIC_ID`. Since the 1-to-1 mapping is verified, this ensures no cell line appears in multiple splits under any identifier (COSMIC_ID, CELL_LINE_NAME, or SANGER_MODEL_ID).

```python
splitter = DataSplitter(
    scenarios="cold_cell",
    data_path="./data/Drug Screening - IC50s/GDSC1_IC50_with_pubchem_and_smiles_filtered_by_mut_rna_mrna_Dec29.xlsx"
)

train_idx, val_idx, test_idx = splitter.split_dataset(fold_idx=1, num_folds=10)
```

### 4. **cold_scaffold**
Split by molecular scaffold - no scaffold appears in both train and test sets.

```python
splitter = DataSplitter(
    scenarios="cold_scaffold",
    data_path="./data/Drug Screening - IC50s/GDSC1_IC50_with_pubchem_and_smiles_filtered_by_mut_rna_mrna_Dec29.xlsx"
)

train_idx, val_idx, test_idx = splitter.split_dataset(fold_idx=1, num_folds=10)
```

### 5. **final**
Stratified split by motif count (80/10/10 split). Uses `fold_idx=1` by default (ignored).

```python
splitter = DataSplitter(
    scenarios="final",
    data_path="./data/Drug Screening - IC50s/GDSC1_IC50_with_pubchem_and_smiles_filtered_by_mut_rna_mrna_Dec29.xlsx"
)

train_idx, val_idx, test_idx = splitter.split_dataset()
```

## Using Split Indices

The splitter returns positional indices that can be used with `pandas.DataFrame.iloc[]`:

```python
import pandas as pd

# Load your data
df = pd.read_excel("./data/Drug Screening - IC50s/GDSC1_IC50_with_pubchem_and_smiles_filtered_by_mut_rna_mrna_Dec29.xlsx")

# Get split indices
splitter = DataSplitter(scenarios="cold_drug", data_path="...")
train_idx, val_idx, test_idx = splitter.split_dataset(fold_idx=1)

# Create train/val/test DataFrames
train_df = df.iloc[train_idx]
val_df = df.iloc[val_idx]
test_df = df.iloc[test_idx]

# Access specific columns
train_smiles = train_df['smiles'].values
train_labels = train_df['LN_IC50'].values
```

## Integration with Model Training

Use the split indices with your `data_process_loader`:

```python
from model import data_process_loader, TransCDR
from torch.utils.data import DataLoader

# Get splits
splitter = DataSplitter(scenarios="cold_cell", data_path="...")
train_idx, val_idx, test_idx = splitter.split_dataset(fold_idx=1)

# Load full dataset
df = pd.read_excel("...")
df = shuffle(df, random_state=2022).reset_index(drop=True)

# Create data loaders
train_loader = DataLoader(
    data_process_loader(train_idx, df.loc[train_idx, 'LN_IC50'].values, df, **config),
    batch_size=32,
    shuffle=True
)

val_loader = DataLoader(
    data_process_loader(val_idx, df.loc[val_idx, 'LN_IC50'].values, df, **config),
    batch_size=32,
    shuffle=False
)

test_loader = DataLoader(
    data_process_loader(test_idx, df.loc[test_idx, 'LN_IC50'].values, df, **config),
    batch_size=32,
    shuffle=False
)
```

## 10-Fold Cross-Validation Example

```python
from splitter import DataSplitter
import pandas as pd

results = []

for fold in range(1, 11):
    print(f"\n{'='*50}")
    print(f"Fold {fold}/10")
    print(f"{'='*50}")

    splitter = DataSplitter(
        scenarios="cold_cell",
        data_path="./data/Drug Screening - IC50s/GDSC1_IC50_with_pubchem_and_smiles_filtered_by_mut_rna_mrna_Dec29.xlsx"
    )

    train_idx, val_idx, test_idx = splitter.split_dataset(fold_idx=fold, num_folds=10)

    # Train your model here
    # model = TransCDR(**config)
    # model.train(train_df, test_df, val_df)
    # metrics = model.test(...)

    # results.append(metrics)

# Aggregate results across folds
# mean_metrics = pd.DataFrame(results).mean()
```

## Key Changes from Original

1. **Removed classification support** - Only regression mode
2. **Updated data path** - Points to Dec29 version of the dataset
3. **Simplified initialization** - Removed `model_type`, `n_sampling`, `label_threshold` parameters
4. **Added verification** - Dataset structure is verified by `verify_dataset.py`
5. **Improved documentation** - Clear notes about cell identifier mapping

## Troubleshooting

### Error: "Missing required columns"

Ensure your Excel file has these columns:
- `DRUG_ID`
- `COSMIC_ID`
- `smiles`
- `LN_IC50`

### Cell Identifier Mismatch

If you modify the dataset and the 1-to-1 mapping breaks, run:

```bash
python verify_dataset.py
```

This will detect if COSMIC_ID, CELL_LINE_NAME, and SANGER_MODEL_ID no longer map 1-to-1, which could cause data leakage in cold_cell splits.

## Model.py Integration Notes

When using the split data with `model.py`, ensure your `data_process_loader` uses the correct column names:

- **RNA data**: Should index by `COSMIC_ID` (not `assay_name`)
- **Mutation data**: Should index by `SANGER_MODEL_ID` (not `COSMIC_ID`)
- **Methylation data**: Should index by `CELL_LINE_NAME` (not `cell_type`)

Based on your omics data loading code, you're already handling this correctly.
