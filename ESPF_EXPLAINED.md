# ESPF Requirements Guide

## What is ESPF?

**ESPF (Extended Sub-structure Fingerprint)** is a Byte-Pair Encoding (BPE) vocabulary for tokenizing SMILES strings into subword units. It consists of:
- `drug_codes_chembl_freq_1500.txt` - BPE merge operations
- `subword_units_map_chembl_freq_1500.csv` - Vocabulary mapping

Think of it like a dictionary that breaks SMILES strings into meaningful chemical fragments.

## Do You Need ESPF?

**Short Answer:** **NO** - if you use pre-trained models (recommended)

**Longer Answer:** It depends on your drug encoder configuration:

### ✅ You DON'T Need ESPF if you use:

1. **Pre-trained models** (`--pre_train True`, **RECOMMENDED**)
   ```bash
   --pre_train True \
   --drug_model 'sequence + graph + FP'
   ```
   This uses ChemBERTa, GIN, and Morgan fingerprints - **no ESPF needed**

2. **CNN encoder** (`--drug_encoder CNN`)
   - Uses simple character-level encoding
   - No external vocabulary needed

3. **RNN encoder** (`--drug_encoder RNN`)
   - Uses simple character-level encoding
   - No external vocabulary needed

4. **Graph encoders** (`--drug_encoder GCN/NeuralFP/AttentiveFP`)
   - Uses molecular graph structure
   - No tokenization needed

### ❌ You NEED ESPF if you use:

**Transformer encoder without pre-training:**
```bash
--pre_train False \
--drug_encoder Transformer
```

This is the ONLY case where ESPF is required.

## Recommendation for Your Setup

Since you mentioned:
> "we only use our dataset, no external one"

I recommend using **pre-trained models** (`--pre_train True`):

```bash
python train_with_splitter.py \
    --scenario cold_cell \
    --fold 1 \
    --pre_train True \
    --drug_model 'sequence + graph + FP' \
    --drug_encoder 'None' \
    --input_dim_drug 2092
```

**Benefits:**
- ✅ **No ESPF needed**
- ✅ **Better performance** (pre-trained on millions of molecules)
- ✅ **Faster training** (drug features pre-computed)
- ✅ Uses your IC50 dataset for drug response prediction only

**What gets pre-trained vs what's trained on your data:**
- **Pre-trained (external)**: Drug feature extractors
  - ChemBERTa: SMILES encoder (trained on PubChem)
  - GIN: Graph encoder (trained on ChEMBL)
  - Morgan FP: Classic fingerprint (no training)
- **Trained on your data**: Everything else
  - Omics encoders (RNA, mutation, methylation)
  - Fusion module
  - Final prediction layers

## If You Still Want ESPF

If you want to use Transformer encoder (not recommended), you need to either:

### Option 1: Get Original ESPF Files

Contact the original TransCDR authors or check the original repository for:
- `data/ESPF/drug_codes_chembl_freq_1500.txt`
- `data/ESPF/subword_units_map_chembl_freq_1500.csv`

### Option 2: Create Your Own BPE Vocabulary

```python
# Create BPE vocabulary from your SMILES
from subword_nmt.learn_bpe import learn_bpe
import pandas as pd

# Load your SMILES
df = pd.read_excel('your_ic50_data.xlsx')
smiles_list = df['smiles'].unique()

# Save to file
with open('smiles_corpus.txt', 'w') as f:
    for s in smiles_list:
        f.write(s + '\n')

# Learn BPE
# (Use subword-nmt library - but this requires some manual work)
```

But this is **not recommended** because:
- More complex setup
- Worse performance than pre-trained models
- Your vocabulary is smaller (238 drugs vs millions in ChEMBL/PubChem)

## Summary

### For Your Use Case:

**ANSWER: You do NOT need ESPF**

Just use pre-trained models:
```bash
--pre_train True \
--drug_model 'sequence + graph + FP'
```

This will:
- ✅ Work with your IC50 dataset only
- ✅ Not require ESPF
- ✅ Give better performance
- ✅ Train faster

The "external" part is just pre-trained drug encoders (like using pre-trained BERT for NLP). Your actual drug response prediction model is trained entirely on your dataset.

## Updated Training Command (No ESPF Needed)

```bash
python train_with_splitter.py \
    --scenario cold_cell \
    --fold 1 \
    --data_path './data/Drug Screening - IC50s/GDSC1_IC50_with_pubchem_and_smiles_filtered_by_mut_rna_mrna_Dec29.xlsx' \
    --omics 'expr + mutation + methylation' \
    --pre_train True \
    --drug_model 'sequence + graph + FP' \
    --drug_encoder 'None' \
    --fusion_type encoder \
    --input_dim_drug 2092 \
    --seq_model 'seyonec/PubChem10M_SMILES_BPE_450k' \
    --graph_model 'gin_supervised_masking' \
    --lr 1e-5 \
    --batch_size 64 \
    --train_epoch 100 \
    --modeldir ./result/cold_cell
```

**No ESPF files required!** ✅
