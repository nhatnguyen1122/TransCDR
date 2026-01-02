# -*- coding: utf-8 -*-
"""
Adapted DataEncoding for GDSC1 IC50 dataset (Regression Only, No ESPF).

This version is optimized for use with pre-trained drug models.
ESPF code removed since it's not needed with pre-trained models.

Changes from original:
- Works with LN_IC50 column (not lnIC50)
- Regression only
- No ESPF dependency (use pre-trained models instead)
"""

import pandas as pd


class DataEncoding:
    """
    Simple data encoder for regression tasks.

    When using pre-trained models (recommended), drug encoding happens
    in the data loader, not here. This class only prepares the Label column.
    """

    def __init__(self, **config):
        self.config = config

    def encode(self, traindata, testdata, valdata):
        """
        Prepare Label column for regression.

        Args:
            traindata: Training DataFrame with 'smiles' and 'LN_IC50' columns
            testdata: Test DataFrame
            valdata: Validation DataFrame

        Returns:
            Train, test, val DataFrames with 'Label' column
        """
        # Make copies to avoid modifying original data
        traindata = traindata.copy()
        testdata = testdata.copy()
        valdata = valdata.copy()

        # Reset index to ensure sequential indexing for data loader
        traindata = traindata.reset_index(drop=True)
        testdata = testdata.reset_index(drop=True)
        valdata = valdata.reset_index(drop=True)

        # Create Label column from LN_IC50 for regression
        traindata['Label'] = traindata['LN_IC50']
        testdata['Label'] = testdata['LN_IC50']
        valdata['Label'] = valdata['LN_IC50']

        print(f"Prepared datasets:")
        print(f"  Train: {len(traindata)} samples")
        print(f"  Val:   {len(valdata)} samples")
        print(f"  Test:  {len(testdata)} samples")

        return traindata, testdata, valdata
