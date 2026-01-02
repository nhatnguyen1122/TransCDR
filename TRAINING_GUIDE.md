# TransCDR Training Guide with Adapted Splitter

## Overview

This guide shows you how to train and evaluate the TransCDR model using your GDSC1 IC50 dataset with the adapted splitter.

## Prerequisites

### 1. Required Data Files

Your dataset file should be at:
```
./data/Drug Screening - IC50s/GDSC1_IC50_with_pubchem_and_smiles_filtered_by_mut_rna_mrna_Dec29.xlsx
```
✅ **Status**: Verified (168,449 samples, 238 drugs, 835 cell lines)

### 2. Omics Data Files

You need these omics data files:

```
data/omics_data/
├── gene_expression/
│   └── Cell_line_RMA_clean.txt          # RNA expression (genes x cells, with "DATA.{COSMIC_ID}" columns)
├── gene_mutation/
│   └── mutation/
│       ├── cell_gene_mapping.json        # {SANGER_MODEL_ID: [gene_indices]}
│       └── list_of_gene_mut.txt          # List of gene names (one per line)
└── dna_methylation/
    └── gene_cell_matrix_promoter_filter_na.csv  # Methylation (genes x cells, with "{CELL_LINE_NAME}_AVG.Beta" columns)
```

### 3. Additional Required Data

For drug encoding (Transformer encoder):
```
data/ESPF/
├── drug_codes_chembl_freq_1500.txt
└── subword_units_map_chembl_freq_1500.csv
```

## Setup Steps

### Step 1: Modify model.py

**CRITICAL**: You need to modify `model.py` to use your column names.

In the `data_process_loader` class (around line 238), replace the `__getitem__` method with:

```python
def __getitem__(self, index):
    'Generates one sample of data'
    index = self.list_IDs[index]
    y = self.labels[index]

    # MODIFIED: Use correct column names for your dataset
    # RNA expression - use COSMIC_ID (not assay_name)
    cosmic_id = str(self.drug_df.iloc[index]['COSMIC_ID'])
    v_rna = np.array(self.rna_data.loc[cosmic_id, :])

    # Gene mutation - use SANGER_MODEL_ID (not COSMIC_ID)
    sanger_id = self.drug_df.iloc[index]['SANGER_MODEL_ID']
    v_genetic = np.array(self.genetic.loc[sanger_id, :])

    # DNA methylation - use CELL_LINE_NAME (not cell_type)
    cell_name = self.drug_df.iloc[index]['CELL_LINE_NAME']
    v_mrna = np.array(self.mrna.loc[cell_name, :])

    # ... rest of the method stays the same
```

**Alternatively**, use the provided `model_adapted.py` which has this change already made.

### Step 2: Update omics data loading in model.py

In the `__init__` method of `data_process_loader` (around line 200-221), replace the omics loading section with:

```python
# RNA expression data - indexed by COSMIC_ID
self.rna_data = pd.read_csv('data/omics_data/gene_expression/Cell_line_RMA_clean.txt', index_col=0)
self.rna_data.columns = [col.replace("DATA.", "") for col in self.rna_data.columns]
self.rna_data = self.rna_data.T  # Transpose to cells x genes

# Gene mutation data - indexed by SANGER_MODEL_ID
with open('data/omics_data/gene_mutation/mutation/cell_gene_mapping.json', 'r') as f:
    cell_to_indices = json.load(f)
with open('data/omics_data/gene_mutation/mutation/list_of_gene_mut.txt', 'r') as f:
    gene_names = [line.strip() for line in f]

cells = list(cell_to_indices.keys())
binary_matrix = np.zeros((len(gene_names), len(cells)), dtype=int)
for col_idx, cell in enumerate(cells):
    indices = cell_to_indices[cell]
    binary_matrix[indices, col_idx] = 1
self.genetic = pd.DataFrame(binary_matrix, index=gene_names, columns=cells).T

# DNA methylation data - indexed by CELL_LINE_NAME
self.mrna = pd.read_csv('data/omics_data/dna_methylation/gene_cell_matrix_promoter_filter_na.csv', index_col=0)
self.mrna.columns = self.mrna.columns.str.replace("_AVG.Beta", "", regex=False)
self.mrna = self.mrna.T  # Transpose to cells x genes
```

### Step 3: Verify omics dimensions

Run this to check your omics data dimensions:

```bash
python model_adapted.py
```

This will print the dimensions you need to configure.

## Training

### Option 1: Train a Single Fold

Train fold 1 with cold_cell scenario:

```bash
python train_with_splitter.py \
    --scenario cold_cell \
    --fold 1 \
    --omics 'expr + mutation + methylation' \
    --drug_model 'sequence + graph + FP' \
    --fusion_type encoder \
    --input_dim_drug 2092 \
    --lr 1e-5 \
    --batch_size 64 \
    --train_epoch 100 \
    --pre_train True \
    --modeldir ./result/cold_cell
```

**Scenarios available:**
- `warm_start`: Random split
- `cold_drug`: New drugs in test set
- `cold_cell`: New cell lines in test set
- `cold_scaffold`: New scaffolds in test set
- `final`: Stratified 80/10/10 split

### Option 2: Run 10-Fold Cross-Validation

```bash
python train_with_splitter.py \
    --scenario cold_cell \
    --cv \
    --omics 'expr + mutation + methylation' \
    --drug_model 'sequence + graph + FP' \
    --fusion_type encoder \
    --input_dim_drug 2092 \
    --lr 1e-5 \
    --batch_size 64 \
    --train_epoch 100 \
    --pre_train True \
    --modeldir ./result/cold_cell_cv10
```

This will train all 10 folds sequentially.

### Quick Test (CPU, no pre-training)

For testing the pipeline without GPU/pre-trained models:

```bash
python train_with_splitter.py \
    --scenario warm_start \
    --fold 1 \
    --omics 'expr + mutation + methylation' \
    --drug_encoder Transformer \
    --drug_model None \
    --fusion_type encoder \
    --input_dim_drug 256 \
    --lr 1e-4 \
    --batch_size 32 \
    --train_epoch 5 \
    --pre_train False \
    --modeldir ./result/test
```

## Getting Results

### View Training Results

After training, check the model directory:

```bash
ls ./result/cold_cell/fold1/
```

You should see:
- `model.pt` - Trained model weights
- `test_markdowntable.txt` - Test set metrics
- `valid_markdowntable.txt` - Validation metrics per epoch
- `loss_curve.png` - Training loss curve
- `config.json` - Model configuration

### Extract Test RMSE

```bash
cat ./result/cold_cell/fold1/test_markdowntable.txt
```

Example output:
```
+------+------+---------------------+---------+----------+-----------+-------------------+
| MSE  | RMSE | Pearson Correlation | p-value | spearman | s_p-value | Concordance Index |
+------+------+---------------------+---------+----------+-----------+-------------------+
| 1.23 | 1.11 | 0.85                | 0.0     | 0.83     | 0.0       | 0.87              |
+------+------+---------------------+---------+----------+-----------+-------------------+
```

The **RMSE** is the second column.

### Aggregate Results from 10-Fold CV

After running 10-fold CV, extract results:

```python
import pandas as pd
import re

results = []
for fold in range(1, 11):
    with open(f'./result/cold_cell_cv10/fold{fold}/test_markdowntable.txt', 'r') as f:
        content = f.read()
        # Parse the table (you may need to adjust parsing logic)
        # Extract RMSE, Pearson, etc.
    results.append({'fold': fold, 'rmse': rmse, 'pearson': pearson, ...})

df = pd.DataFrame(results)
print("Mean metrics across 10 folds:")
print(df.mean())
print("\nStd metrics across 10 folds:")
print(df.std())
```

Or use the Step3_result.py script (may need adaptation):

```bash
python Step3_result.py --result_folder ./result/cold_cell_cv10
```

## Troubleshooting

### Error: Missing omics data files

```
FileNotFoundError: [Errno 2] No such file or directory: 'data/omics_data/...'
```

**Solution**: Ensure all omics data files are in place. Check with:
```bash
find data/omics_data -type f
```

### Error: KeyError when accessing cell identifiers

```
KeyError: 'COSMIC_ID' or 'SANGER_MODEL_ID' or 'CELL_LINE_NAME'
```

**Solution**:
1. Verify your IC50 dataset has these columns
2. Check that cell identifiers in IC50 data match those in omics data
3. Run verify_dataset.py to confirm data structure

### Error: Dimension mismatch

```
RuntimeError: size mismatch, m1: [32 x 18451], m2: [735 x 1024]
```

**Solution**: Update the config with correct omics dimensions from `model_adapted.py`.

### CUDA out of memory

**Solution**: Reduce batch size:
```bash
--batch_size 32  # or even 16
```

## Configuration Options

### Drug Encoders

- **None**: Use pre-trained sequence/graph/FP models (requires `--pre_train True`)
- **CNN**: 1D CNN on SMILES
- **RNN**: LSTM on SMILES
- **Transformer**: Transformer encoder on SMILES tokens
- **GCN**: Graph convolutional network
- **NeuralFP**: Neural fingerprint
- **AttentiveFP**: Attentive fingerprint

### Drug Models (for pre-training)

- **sequence**: BERT-based SMILES encoder (768-dim)
- **graph**: GIN pre-trained graph encoder (300-dim)
- **FP**: Morgan fingerprint (1024-dim)
- **sequence + graph**: Concatenated (1068-dim)
- **sequence + graph + FP**: All three (2092-dim)

### Fusion Types

- **concat**: Simple concatenation
- **encoder**: Transformer encoder fusion
- **decoder**: Cross-attention decoder fusion

## Example: Complete Training Pipeline

```bash
# 1. Verify data
python verify_dataset.py

# 2. Check omics dimensions
python model_adapted.py

# 3. Test on 1 fold
python train_with_splitter.py --scenario cold_cell --fold 1 --train_epoch 5

# 4. If test works, run full 10-fold CV
python train_with_splitter.py --scenario cold_cell --cv --train_epoch 100

# 5. Extract results
for fold in {1..10}; do
    echo "Fold $fold:"
    cat ./result/cold_cell_cv10/fold${fold}/test_markdowntable.txt
done
```

## Expected Runtime

On a modern GPU (e.g., RTX 3090):
- **Single fold**: ~30-60 minutes (depending on epochs)
- **10-fold CV**: ~5-10 hours

On CPU:
- **Single fold**: ~2-4 hours
- **10-fold CV**: ~20-40 hours

Early stopping (default: 5 epochs patience) will reduce training time.

## Summary

✅ **Splitter verified**: All scenarios working
✅ **Data loader adapted**: Uses your column names
✅ **Training script ready**: `train_with_splitter.py`
✅ **Documentation complete**: This guide

**Next steps:**
1. Set up your omics data files
2. Modify model.py or use model_adapted.py
3. Run a test training (1 fold, few epochs)
4. Run full 10-fold CV
5. Extract and report RMSE results
