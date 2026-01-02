"""
Adapted model.py data loader for GDSC1 IC50 dataset.

Key changes:
- Uses COSMIC_ID for RNA data (not assay_name)
- Uses SANGER_MODEL_ID for mutation data (not COSMIC_ID)
- Uses CELL_LINE_NAME for methylation data (not cell_type)

Usage:
  Replace the data_process_loader class in model.py with this version,
  OR import this as: from model_adapted import data_process_loader_adapted
"""

import os
import numpy as np
import pandas as pd
import json
import torch
from torch.utils import data
from tqdm import tqdm
from dgllife.utils import mol_to_bigraph, PretrainAtomFeaturizer, PretrainBondFeaturizer
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from dgllife.model import load_pretrained
from dgl.nn.pytorch.glob import AvgPooling
from transformers import AutoModel

# Import helper functions from original model.py
# (Assumes these functions exist in model.py)
try:
    from model import (
        get_sequence_feats, get_graph_feats, get_FP_feats,
        drug_2_embed, get_transformer_feats, get_GNN_feats
    )
except ImportError:
    print("Warning: Could not import helper functions from model.py")
    print("Make sure model.py is in the same directory")


device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


class data_process_loader_adapted(data.Dataset):
    """
    Adapted data loader for GDSC1 IC50 dataset.

    Column mapping:
    - RNA expression: COSMIC_ID (your data) instead of assay_name (original)
    - Gene mutation: SANGER_MODEL_ID (your data) instead of COSMIC_ID (original)
    - DNA methylation: CELL_LINE_NAME (your data) instead of cell_type (original)
    """

    def __init__(self, list_IDs, labels, drug_df, **config):
        'Initialization'
        self.config = config
        self.labels = labels
        self.list_IDs = list_IDs
        self.drug_df = drug_df

        # Load omics data based on your structure
        # RNA expression data - indexed by COSMIC_ID
        self.rna_data = pd.read_csv('data/omics_data/gene_expression/Cell_line_RMA_clean.txt', index_col=0)
        self.rna_data.columns = [col.replace("DATA.", "") for col in self.rna_data.columns]
        # Note: rna_data is genes x cells, need to transpose for cell x gene lookup
        self.rna_data = self.rna_data.T

        # Gene mutation data - indexed by SANGER_MODEL_ID
        with open('data/omics_data/gene_mutation/mutation/cell_gene_mapping.json', 'r') as f:
            cell_to_indices = json.load(f)
        with open('data/omics_data/gene_mutation/mutation/list_of_gene_mut.txt', 'r') as f:
            gene_names = [line.strip() for line in f]

        cells = list(cell_to_indices.keys())
        max_index = max(max(indices) if len(indices) > 0 else -1
                        for indices in cell_to_indices.values())
        assert max_index < len(gene_names), "Some gene indices exceed gene list length"

        # Create binary matrix: rows=genes, columns=cells
        binary_matrix = np.zeros((len(gene_names), len(cells)), dtype=int)
        for col_idx, cell in enumerate(cells):
            indices = cell_to_indices[cell]
            binary_matrix[indices, col_idx] = 1

        self.genetic = pd.DataFrame(binary_matrix, index=gene_names, columns=cells)
        # Transpose for cell x gene lookup
        self.genetic = self.genetic.T

        # DNA methylation data - indexed by CELL_LINE_NAME
        self.mrna = pd.read_csv('data/omics_data/dna_methylation/gene_cell_matrix_promoter_filter_na.csv', index_col=0)
        self.mrna.columns = self.mrna.columns.str.replace("_AVG.Beta", "", regex=False)
        # Transpose for cell x gene lookup
        self.mrna = self.mrna.T

        print(f"Loaded omics data:")
        print(f"  RNA: {self.rna_data.shape} (cells x genes)")
        print(f"  Mutation: {self.genetic.shape} (cells x genes)")
        print(f"  Methylation: {self.mrna.shape} (cells x genes)")

        # Prepare drug embeddings based on encoder type
        if self.config['pre_train']:
            self.embedded_drug1 = get_sequence_feats(drug_df, **self.config)
            self.embedded_drug2 = get_graph_feats(drug_df, **self.config)
            self.embedded_drug3 = get_FP_feats(drug_df, **self.config)
        else:
            if (self.config['drug_encoder'] == 'CNN') | (self.config['drug_encoder'] == 'RNN'):
                self.embedded_drug = drug_2_embed(drug_df)
            if self.config['drug_encoder'] == 'Transformer':
                self.embedded_drug, self.mask = get_transformer_feats(drug_df)
            if self.config['drug_encoder'] in ['GCN', 'NeuralFP', 'AttentiveFP']:
                self.embedded_drug = get_GNN_feats(drug_df, **self.config)

    def __len__(self):
        'Denotes the total number of samples'
        return len(self.list_IDs)

    def __getitem__(self, index):
        'Generates one sample of data'
        index = self.list_IDs[index]
        y = self.labels[index]

        # Get omics data using the correct column names
        # RNA expression - use COSMIC_ID
        cosmic_id = str(self.drug_df.iloc[index]['COSMIC_ID'])
        v_rna = np.array(self.rna_data.loc[cosmic_id, :])

        # Gene mutation - use SANGER_MODEL_ID
        sanger_id = self.drug_df.iloc[index]['SANGER_MODEL_ID']
        v_genetic = np.array(self.genetic.loc[sanger_id, :])

        # DNA methylation - use CELL_LINE_NAME
        cell_name = self.drug_df.iloc[index]['CELL_LINE_NAME']
        v_mrna = np.array(self.mrna.loc[cell_name, :])

        # Get drug features
        if self.config['pre_train']:
            v_d1 = self.embedded_drug1[self.drug_df.iloc[index]['smiles']]
            v_d2 = self.embedded_drug2[self.drug_df.iloc[index]['smiles']]
            v_d3 = self.embedded_drug3[self.drug_df.iloc[index]['smiles']]
            return v_d1, v_d2, v_d3, v_rna, v_genetic, v_mrna, y
        else:
            if self.config['drug_encoder'] in ['CNN', 'RNN', 'GCN', 'NeuralFP', 'AttentiveFP']:
                v_d = self.embedded_drug[self.drug_df.iloc[index]['smiles']]
                return v_d, v_rna, v_genetic, v_mrna, y
            if self.config['drug_encoder'] == 'Transformer':
                v_d = self.embedded_drug[self.drug_df.iloc[index]['smiles']]
                mask = self.mask[self.drug_df.iloc[index]['smiles']]
                return (v_d, mask), v_rna, v_genetic, v_mrna, y


def get_omics_dimensions():
    """
    Get the dimensions of omics data for your dataset.
    Use this to configure the model's input dimensions.

    Returns:
        dict: Dictionary with input dimensions for RNA, mutation, methylation
    """
    # RNA expression
    rna_data = pd.read_csv('data/omics_data/gene_expression/Cell_line_RMA_clean.txt', index_col=0)
    input_dim_rna = rna_data.shape[0]  # number of genes

    # Gene mutation
    with open('data/omics_data/gene_mutation/mutation/list_of_gene_mut.txt', 'r') as f:
        gene_names = [line.strip() for line in f]
    input_dim_genetic = len(gene_names)

    # DNA methylation
    mrna = pd.read_csv('data/omics_data/dna_methylation/gene_cell_matrix_promoter_filter_na.csv', index_col=0)
    input_dim_mrna = mrna.shape[0]  # number of methylation sites

    return {
        'input_dim_rna': input_dim_rna,
        'input_dim_genetic': input_dim_genetic,
        'input_dim_mrna': input_dim_mrna
    }


if __name__ == "__main__":
    # Test: check omics dimensions
    try:
        dims = get_omics_dimensions()
        print("Omics data dimensions:")
        for key, value in dims.items():
            print(f"  {key}: {value}")
    except FileNotFoundError as e:
        print(f"Error: {e}")
        print("\nPlease ensure your omics data files are at:")
        print("  - data/omics_data/gene_expression/Cell_line_RMA_clean.txt")
        print("  - data/omics_data/gene_mutation/mutation/cell_gene_mapping.json")
        print("  - data/omics_data/gene_mutation/mutation/list_of_gene_mut.txt")
        print("  - data/omics_data/dna_methylation/gene_cell_matrix_promoter_filter_na.csv")
