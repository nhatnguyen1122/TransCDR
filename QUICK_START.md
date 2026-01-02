# Quick Start Guide: Training TransCDR with Your Dataset

## ✅ What's Ready

Your splitter is **fully working** and verified:
- ✅ All 5 scenarios tested (warm_start, cold_drug, cold_cell, cold_scaffold, final)
- ✅ 168,449 samples, 238 drugs, 835 cell lines
- ✅ No data leakage verified
- ✅ 1-to-1 cell identifier mapping confirmed

## 🚀 How to Train and Get RMSE

### Step 1: Set Up Omics Data

Place your omics data files at:

```
data/omics_data/
├── gene_expression/
│   └── Cell_line_RMA_clean.txt
├── gene_mutation/
│   └── mutation/
│       ├── cell_gene_mapping.json
│       └── list_of_gene_mut.txt
└── dna_methylation/
    └── gene_cell_matrix_promoter_filter_na.csv
```

**Note about ESPF:** You do **NOT** need ESPF files if you use pre-trained models (recommended).
ESPF is only required if using `--drug_encoder Transformer --pre_train False` (not recommended).

See `ESPF_EXPLAINED.md` for details.

### Step 2: Modify model.py

**Option A: Quick patch**

In `model.py`, find the `data_process_loader.__getitem__` method (around line 238) and replace these lines:

```python
# OLD:
v_rna = np.array(self.rna_data.loc[self.drug_df.iloc[index]['assay_name'],:])
v_genetic = np.array(self.genetic.loc[int(self.drug_df.iloc[index]['COSMIC_ID']),:])
v_mrna = np.array(self.mrna.loc[self.drug_df.iloc[index]['cell_type'],:])

# NEW:
v_rna = np.array(self.rna_data.loc[str(self.drug_df.iloc[index]['COSMIC_ID']),:])
v_genetic = np.array(self.genetic.loc[self.drug_df.iloc[index]['SANGER_MODEL_ID'],:])
v_mrna = np.array(self.mrna.loc[self.drug_df.iloc[index]['CELL_LINE_NAME'],:])
```

Also update the omics loading section in `__init__` (see TRAINING_GUIDE.md for details).

**Option B: Use adapted version**

Copy the data loader from `model_adapted.py` into `model.py`.

### Step 3: Train a Single Fold (Test)

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
    --train_epoch 5 \
    --pre_train True \
    --modeldir ./result/test
```

This will:
1. ✅ Split data using your adapted splitter
2. ✅ Load and encode drug features
3. ✅ Train the model for 5 epochs
4. ✅ Evaluate on test set
5. ✅ Save results to `./result/test/fold1/`

### Step 4: Get Test RMSE

```bash
cat ./result/test/fold1/test_markdowntable.txt
```

You'll see something like:
```
+------+------+---------------------+---------+----------+-----------+-------------------+
| MSE  | RMSE | Pearson Correlation | p-value | spearman | s_p-value | Concordance Index |
+------+------+---------------------+---------+----------+-----------+-------------------+
| 1.23 | 1.11 | 0.85                | 0.0     | 0.83     | 0.0       | 0.87              |
+------+------+---------------------+---------+----------+-----------+-------------------+
```

**RMSE is the second column** (1.11 in this example).

### Step 5: Run Full 10-Fold CV

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

This will train all 10 folds.

### Step 6: Collect Results from All Folds

```bash
# Simple: view all results
for fold in {1..10}; do
    echo "=== Fold $fold ==="
    cat ./result/cold_cell_cv10/fold${fold}/test_markdowntable.txt
done

# Or create a summary script
python -c "
import pandas as pd
import re

results = []
for fold in range(1, 11):
    try:
        with open(f'./result/cold_cell_cv10/fold{fold}/test_markdowntable.txt', 'r') as f:
            lines = f.readlines()
            # Parse the table (row with numbers is line 3)
            if len(lines) > 3:
                data_line = lines[3].strip()
                # Extract numbers from line
                numbers = [float(x.strip()) for x in data_line.split('|')[1:-1]]
                results.append({
                    'fold': fold,
                    'MSE': numbers[0],
                    'RMSE': numbers[1],
                    'Pearson': numbers[2],
                    'p-value': numbers[3],
                    'Spearman': numbers[4],
                    's_p-value': numbers[5],
                    'CI': numbers[6]
                })
    except Exception as e:
        print(f'Error reading fold {fold}: {e}')

df = pd.DataFrame(results)
print('\nMean Results Across 10 Folds:')
print(df.mean(numeric_only=True))
print('\nStd Results Across 10 Folds:')
print(df.std(numeric_only=True))
print(f'\nFinal Test RMSE: {df[\"RMSE\"].mean():.4f} ± {df[\"RMSE\"].std():.4f}')
"
```

## 📋 All Splitting Scenarios

Try different scenarios to evaluate model generalization:

### 1. Warm Start (Random Split)
```bash
python train_with_splitter.py --scenario warm_start --cv
```
Both drugs and cells can appear in train and test.

### 2. Cold Drug
```bash
python train_with_splitter.py --scenario cold_drug --cv
```
Test on **completely new drugs** not seen during training.

### 3. Cold Cell
```bash
python train_with_splitter.py --scenario cold_cell --cv
```
Test on **completely new cell lines** not seen during training.

### 4. Cold Scaffold
```bash
python train_with_splitter.py --scenario cold_scaffold --cv
```
Test on drugs with **completely new scaffolds**.

### 5. Final (80/10/10 Split)
```bash
python train_with_splitter.py --scenario final --fold 1
```
Stratified by motif count, single split (not CV).

## 🎯 Expected Metrics

Based on similar GDSC regression tasks:

| Scenario | Expected RMSE | Difficulty |
|----------|---------------|------------|
| Warm Start | 0.8 - 1.2 | Easy |
| Cold Drug | 1.2 - 1.8 | Hard |
| Cold Cell | 1.5 - 2.2 | Very Hard |
| Cold Scaffold | 1.3 - 2.0 | Hard |

Your actual results will depend on:
- Model configuration
- Hyperparameters
- Omics data quality
- Training epochs

## 🔧 Troubleshooting

### "FileNotFoundError: omics data"
→ Set up omics data files (see Step 1)

### "KeyError: COSMIC_ID/SANGER_MODEL_ID/CELL_LINE_NAME"
→ Modify model.py to use your column names (see Step 2)

### "CUDA out of memory"
→ Reduce batch size: `--batch_size 32` or `--batch_size 16`

### Training is slow
→ First try:
  - Smaller epochs for testing: `--train_epoch 5`
  - Without pre-training: `--pre_train False --drug_encoder Transformer`
  - Smaller batch: `--batch_size 32`

## 📊 What You Get

After training, each fold directory contains:

```
result/cold_cell_cv10/fold1/
├── model.pt                      # Trained model
├── test_markdowntable.txt        # ⭐ Test RMSE and metrics
├── valid_markdowntable.txt       # Validation metrics per epoch
├── loss_curve.png                # Training loss visualization
├── loss_curve_iter.pkl           # Loss history
├── logits.npy                    # Model predictions
└── config.json                   # Model configuration
```

## 🎓 Quick Example

```bash
# Complete example: Train fold 1 of cold_cell scenario
python train_with_splitter.py --scenario cold_cell --fold 1

# Get RMSE
cat ./result/splitter_cv/fold1/test_markdowntable.txt

# Expected output shows:
# MSE, RMSE, Pearson r, p-value, Spearman r, p-value, CI
```

## 📚 More Details

- Full training guide: `TRAINING_GUIDE.md`
- Splitter usage: `SPLITTER_USAGE.md`
- Adaptation summary: `ADAPTATION_SUMMARY.md`

---

**You're all set!** The splitter is working, the training pipeline is ready. Just set up your omics data files, modify model.py, and start training! 🚀
