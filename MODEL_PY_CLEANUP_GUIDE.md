# model.py Cleanup Guide

## Critical Changes Needed in model.py

Your `model.py` needs TWO types of changes:
1. **Adapt for your column names** (REQUIRED)
2. **Remove unused code** (OPTIONAL but recommended)

---

## REQUIRED: Adapt Column Names

### Change 1: Update `data_process_loader.__getitem__` (around line 238-250)

**FIND THIS:**
```python
def __getitem__(self, index):
    'Generates one sample of data'
    index = self.list_IDs[index]
    y = self.labels[index]

    if self.config['screening'] == 'TCGA':
        # TCGA code...
    else:
        if self.config['external_dataset'] == 'None':
            v_rna = np.array(self.rna_data.loc[self.drug_df.iloc[index]['assay_name'],:])
            v_genetic = np.array(self.genetic.loc[int(self.drug_df.iloc[index]['COSMIC_ID']),:])
            v_mrna = np.array(self.mrna.loc[self.drug_df.iloc[index]['cell_type'],:])
```

**REPLACE WITH:**
```python
def __getitem__(self, index):
    'Generates one sample of data'
    index = self.list_IDs[index]
    y = self.labels[index]

    # Use YOUR column names
    v_rna = np.array(self.rna_data.loc[str(self.drug_df.iloc[index]['COSMIC_ID']),:])
    v_genetic = np.array(self.genetic.loc[self.drug_df.iloc[index]['SANGER_MODEL_ID']],:])
    v_mrna = np.array(self.mrna.loc[self.drug_df.iloc[index]['CELL_LINE_NAME']],:])
```

**Note:** Remove all the TCGA/CCLE external dataset code - you only use GDSC.

### Change 2: Update omics data loading in `__init__` (around line 200-221)

**FIND THIS:**
```python
if self.config['external_dataset'] == 'None':
    self.rna_data = pd.read_csv('./data/GDSC/data_processed/RNA_n18451_1018_zscore.csv',index_col=0)
    self.rna_data.columns = [name.split('.')[0][1:] for name in self.rna_data.columns.values]
    self.rna_data = self.rna_data.T
    self.genetic = pd.read_csv('./data/GDSC/data_processed/Genetic_features_n969_735.txt',sep='\t',index_col=0)
    self.mrna = pd.read_csv('./data/GDSC/data_processed/mrna_n20617_1028_zscore.csv',index_col=0)
    self.mrna = self.mrna.T
```

**REPLACE WITH:**
```python
# RNA expression - indexed by COSMIC_ID
self.rna_data = pd.read_csv('data/omics_data/gene_expression/Cell_line_RMA_clean.txt', index_col=0)
self.rna_data.columns = [col.replace("DATA.", "") for col in self.rna_data.columns]
self.rna_data = self.rna_data.T  # Transpose to cells x genes

# Gene mutation - indexed by SANGER_MODEL_ID
import json
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

# DNA methylation - indexed by CELL_LINE_NAME
self.mrna = pd.read_csv('data/omics_data/dna_methylation/gene_cell_matrix_promoter_filter_na.csv', index_col=0)
self.mrna.columns = self.mrna.columns.str.replace("_AVG.Beta", "", regex=False)
self.mrna = self.mrna.T  # Transpose to cells x genes
```

**Note:** Remove all TCGA and CCLE dataset loading code.

---

## OPTIONAL: Remove Unused Code

These removals will simplify model.py:

### 1. Remove Classification Code

**In `TransCDR.test()` method (lines 482-520):**

Delete this entire section (you only need regression):
```python
if self.config['model_type'] == 'classification':
    m = torch.nn.Sigmoid()
    logits = torch.squeeze(m(score)).detach().cpu().numpy()
```

And this at the end:
```python
if self.config['model_type'] == 'classification':
    return roc_auc_score(y_label, y_pred), \
        average_precision_score(y_label, y_pred), \
        f1_score(y_label, outputs), \
        log_loss(y_label, outputs), \
        y_pred,y_label
```

**In `TransCDR.train()` method (lines 522-706):**

Delete classification training code:
```python
if self.config['model_type'] == 'classification':
    loss_fct = torch.nn.BCELoss()
    m = torch.nn.Sigmoid()
    n = torch.squeeze(m(score), 1)
    loss = loss_fct(n, label)
```

And classification validation:
```python
if self.config['model_type'] == 'classification':
    auc, auprc, f1, loss, logits,_ = self.test(validation_generator, self.model)
    ...
```

### 2. Remove ESPF/Transformer Encoding Code (if using pre-trained)

**Lines 130-162:** Remove if you don't use Transformer encoder
```python
# transformer
from subword_nmt.apply_bpe import BPE
import codecs
vocab_path = './data/ESPF/drug_codes_chembl_freq_1500.txt'
...
def get_transformer_feats(drug_df):
    ...
```

### 3. Remove External Dataset Code

**In `data_process_loader.__init__` (lines 191-221):**

Remove all this (you only use GDSC):
```python
if self.config['screening'] == 'TCGA':
    ...

if self.config['external_dataset'] == 'TCGA':
    ...

if self.config['external_dataset'] == 'CCLE':
    ...
```

**In `data_process_loader.__getitem__` (lines 238-272):**

Remove:
```python
if self.config['screening'] == 'TCGA':
    ...

if self.config['external_dataset'] == 'TCGA':
    ...

if self.config['external_dataset'] == 'CCLE':
    ...
```

### 4. Remove Unused Encoders (if using pre-trained)

If you're using `--pre_train True` (recommended), you don't need:
- CNN encoder (lines 103-117)
- RNN encoder
- Transformer encoder (lines 130-162)

But keep them if you want flexibility to try different encoders later.

---

## Summary of Changes

### REQUIRED (for your dataset to work):
1. ✅ Update column names in `__getitem__`: COSMIC_ID, SANGER_MODEL_ID, CELL_LINE_NAME
2. ✅ Update omics loading paths in `__init__`

### RECOMMENDED (cleanup):
3. Remove classification code
4. Remove ESPF/Transformer code (if using pre-trained)
5. Remove TCGA/CCLE external dataset code
6. Remove unused drug encoders

---

## Alternative: Use model_adapted.py as Reference

Instead of manually editing model.py, you can:
1. Open `model_adapted.py`
2. Copy the `data_process_loader` class
3. Replace the one in `model.py`

This gives you the adapted version with clean code.

---

## After Cleanup: What model.py Should Contain

**Keep:**
- `get_sequence_feats()`, `get_graph_feats()`, `get_FP_feats()` (for pre-trained)
- `data_process_loader` (adapted for your column names)
- `Classifier` class
- `TransCDR` class (regression mode only)

**Remove:**
- Classification logic
- ESPF/Transformer encoding (if using pre-trained)
- External dataset code (TCGA, CCLE)
- Unused drug encoders (CNN, RNN if using pre-trained)

---

## Testing After Cleanup

After modifying model.py, test with:
```bash
python train_with_splitter_clean.py --scenario warm_start --fold 1 --train_epoch 2
```

If it works for 2 epochs, you're good!
