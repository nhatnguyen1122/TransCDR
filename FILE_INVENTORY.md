# File Inventory and Usage Analysis

## Files You NEED (Keep and Use)

### 1. Core Training Files

| File | Status | Purpose |
|------|--------|---------|
| `splitter.py` | ✅ Keep | Data splitting (all scenarios) |
| `DataEncoding_adapted.py` | ✅ Keep | Encode drugs for your dataset |
| `train_with_splitter.py` | ✅ Keep | Main training script |
| `model.py` | ✅ Keep (needs modification) | TransCDR model architecture |
| `model_helper.py` | ✅ Keep | Model components (attention, encoders, etc.) |
| `drug_bert_model.py` | ✅ Keep | Pre-trained drug model utilities |

### 2. Data Files

| File/Directory | Status | Purpose |
|----------------|--------|---------|
| `data/Drug Screening - IC50s/GDSC1_IC50_..._Dec29.xlsx` | ✅ Required | Your IC50 dataset |
| `data/omics_data/` | ✅ Required | RNA, mutation, methylation data |
| `data/ESPF/` | ❌ NOT needed | Only for Transformer encoder (we use pre-trained) |

### 3. Documentation (Optional but Recommended)

| File | Keep? | Purpose |
|------|-------|---------|
| `QUICK_START.md` | ✅ | Quick start guide |
| `TRAINING_GUIDE.md` | ✅ | Detailed training instructions |
| `ESPF_EXPLAINED.md` | ✅ | ESPF explanation |
| `SPLITTER_USAGE.md` | ✅ | Splitter documentation |
| `ADAPTATION_SUMMARY.md` | ✅ | Summary of changes |

### 4. Testing/Verification (Optional)

| File | Keep? | Purpose |
|------|-------|---------|
| `verify_dataset.py` | Optional | Dataset verification |
| `test_splitter.py` | Optional | Splitter testing |

## Files You DON'T NEED (Can Delete)

### Original Files (Replaced by Adapted Versions)

| File | Replace With | Reason |
|------|--------------|--------|
| `DataEncoding.py` | `DataEncoding_adapted.py` | Old version uses `lnIC50` not `LN_IC50` |
| `Step1_Data_split.py` | `splitter.py` | Old splitting method |
| `Step2_train_model.py` | `train_with_splitter.py` | Old training script |
| `Step3_result.py` | Manual result extraction | Not adapted for new structure |
| `Step4_Train_final_model.py` | `train_with_splitter.py --scenario final` | Replaced |
| `Step5.1_test_on_CCLE_data.py` | N/A | External dataset (you don't use) |
| `Step5.2_screening_drugs_for_TCGA_patients.py` | N/A | External dataset (you don't use) |
| `Step6_CDR_prediction.py` | Use trained model directly | Not needed initially |

### Scripts Directory (Old Shell Scripts)

| File | Keep? | Reason |
|------|-------|--------|
| `script/Step*.sh` | ❌ Delete | Replaced by `train_with_splitter.py` |

### Reference Files

| File | Status |
|------|--------|
| `model_adapted.py` | ❌ Delete (after using as reference) | Only needed to copy changes to model.py |

## Summary

### Keep (11 files):
1. `splitter.py`
2. `DataEncoding_adapted.py`
3. `train_with_splitter.py`
4. `model.py` (after modification)
5. `model_helper.py`
6. `drug_bert_model.py`
7. Documentation files (5 .md files)
8. Verification scripts (optional)

### Delete (13+ files):
1. `DataEncoding.py`
2. `Step1_Data_split.py`
3. `Step2_train_model.py`
4. `Step3_result.py`
5. `Step4_Train_final_model.py`
6. `Step5.1_test_on_CCLE_data.py`
7. `Step5.2_screening_drugs_for_TCGA_patients.py`
8. `Step6_CDR_prediction.py`
9. `model_adapted.py` (after copying to model.py)
10. All `script/*.sh` files

## Next Steps

1. Clean up code in files we're keeping
2. Remove unused imports
3. Remove classification code
4. Remove ESPF dependencies (we use pre-trained)
