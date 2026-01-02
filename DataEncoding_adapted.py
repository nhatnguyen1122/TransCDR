# -*- coding: utf-8 -*-
"""
Adapted DataEncoding for GDSC1 IC50 dataset with Dec29 column names.

Changes from original:
- Works with LN_IC50 column (not lnIC50)
- Works with your dataset structure
- Regression only
"""

import numpy as np
import pandas as pd
import codecs
from subword_nmt.apply_bpe import BPE


class DataEncoding:
    def __init__(self, **config):
        self.config = config

    def _drug2emb_encoder(self, smile):
        '''Get the token and mask of drugs for Transformer encoder'''

        vocab_path = './data/ESPF/drug_codes_chembl_freq_1500.txt'
        sub_csv = pd.read_csv('./data/ESPF/subword_units_map_chembl_freq_1500.csv')

        bpe_codes_drug = codecs.open(vocab_path)
        dbpe = BPE(bpe_codes_drug, merges=-1, separator='')

        idx2word_d = sub_csv['index'].values
        words2idx_d = dict(zip(idx2word_d, range(0, len(idx2word_d))))

        max_d = 50
        t1 = dbpe.process_line(smile).split()  # split
        try:
            i1 = np.asarray([words2idx_d[i] for i in t1])  # index
        except:
            i1 = np.array([0])

        l = len(i1)
        if l < max_d:
            i = np.pad(i1, (0, max_d - l), 'constant', constant_values=0)
            input_mask = ([1] * l) + ([0] * (max_d - l))  # mask: 1 for real tokens, 0 for padding
        else:
            i = i1[:max_d]
            input_mask = [1] * max_d

        return i, np.asarray(input_mask)

    def encode(self, traindata, testdata, valdata):
        """
        Encode drug SMILES and prepare Label column for regression.

        Args:
            traindata: Training DataFrame with 'smiles' and 'LN_IC50' columns
            testdata: Test DataFrame
            valdata: Validation DataFrame

        Returns:
            Encoded train, test, val DataFrames with 'Label' column
        """
        # Get unique SMILES from the data
        all_smiles = pd.concat([
            traindata['smiles'],
            testdata['smiles'],
            valdata['smiles']
        ]).unique()

        # Create encoding dictionary for all unique SMILES
        smile_encode = pd.Series(all_smiles).apply(self._drug2emb_encoder)
        uniq_smile_dict = dict(zip(all_smiles, smile_encode))

        # Encode each dataset
        traindata = traindata.copy()
        testdata = testdata.copy()
        valdata = valdata.copy()

        traindata['drug_encoding'] = [uniq_smile_dict[i] for i in traindata['smiles']]
        testdata['drug_encoding'] = [uniq_smile_dict[i] for i in testdata['smiles']]
        valdata['drug_encoding'] = [uniq_smile_dict[i] for i in valdata['smiles']]

        # Reset index to ensure sequential indexing for data loader
        traindata = traindata.reset_index(drop=True)
        testdata = testdata.reset_index(drop=True)
        valdata = valdata.reset_index(drop=True)

        # Create Label column from LN_IC50
        traindata['Label'] = traindata['LN_IC50']
        testdata['Label'] = testdata['LN_IC50']
        valdata['Label'] = valdata['LN_IC50']

        print(f"Encoded datasets:")
        print(f"  Train: {len(traindata)} samples")
        print(f"  Val:   {len(valdata)} samples")
        print(f"  Test:  {len(testdata)} samples")

        return traindata, testdata, valdata
