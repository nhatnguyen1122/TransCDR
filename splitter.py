import pandas as pd
from sklearn.utils import shuffle
from sklearn.model_selection import KFold, train_test_split
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

# -------------------------------------------------------------------------- #
#                           HELPER FUNCTIONS                                 #
# -------------------------------------------------------------------------- #

def get_mol(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        Chem.Kekulize(mol)
    except:
        pass
    return mol


def motif_decomp(mol):
    """Simplified motif decomposition - returns list of atom index lists."""
    from rdkit.Chem import BRICS

    if isinstance(mol, str):
        mol = Chem.MolFromSmiles(mol)
    if mol is None:
        return []

    n_atoms = mol.GetNumAtoms()
    if n_atoms == 1:
        return [[0]]

    cliques = []
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom().GetIdx()
        a2 = bond.GetEndAtom().GetIdx()
        cliques.append([a1, a2])

    res = list(BRICS.FindBRICSBonds(mol))
    if len(res) != 0:
        for bond in res:
            if [bond[0][0], bond[0][1]] in cliques:
                cliques.remove([bond[0][0], bond[0][1]])
            else:
                cliques.remove([bond[0][1], bond[0][0]])
            cliques.append([bond[0][0]])
            cliques.append([bond[0][1]])

    for c in range(len(cliques) - 1):
        if c >= len(cliques):
            break
        for k in range(c + 1, len(cliques)):
            if k >= len(cliques):
                break
            if len(set(cliques[c]) & set(cliques[k])) > 0:
                cliques[c] = list(set(cliques[c]) | set(cliques[k]))
                cliques[k] = []
        cliques = [c for c in cliques if len(c) > 0]
    cliques = [c for c in cliques if n_atoms > len(c) > 0]

    return cliques


def generate_scaffold(smiles, include_chirality=False):
    """
    Obtain scaffold from smiles
    """
    try:
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(
            smiles=smiles, includeChirality=include_chirality)
        return scaffold
    except:
        return smiles  # Return original smiles if scaffold generation fails


# -------------------------------------------------------------------------- #
#                               MAIN CLASS                                   #
# -------------------------------------------------------------------------- #

class DataSplitter:
    """
    Data splitter for drug response prediction (regression only).
    Supports: warm_start, final, cold_drug, cold_cell, cold_scaffold

    Returns positional indices compatible with MoleculeDataset.

    Dataset requirements:
    - Required columns: DRUG_ID, COSMIC_ID, smiles, LN_IC50
    - Optional columns: CELL_LINE_NAME, SANGER_MODEL_ID
    - Note: COSMIC_ID, CELL_LINE_NAME, and SANGER_MODEL_ID should have 1-to-1 mapping
          (verified by verify_dataset.py script)
    """

    def __init__(self,
                 scenarios: str,
                 data_path: str,
                 n_folds: int = 10,
                 random_state: int = 2022):
        """
        Args:
            scenarios: 'warm_start', 'final', 'cold_drug', 'cold_cell', 'cold_scaffold'
            data_path: Path to Excel file with drug response data
            n_folds: Number of folds for K-fold splitting
            random_state: Random seed for reproducibility
        """
        self.scenarios = scenarios
        self.n_folds = n_folds
        self.rs = random_state

        self._load_data(data_path)

    def _load_data(self, data_path):
        """Load data from Excel file."""
        # Load without index_col to ensure integer positional indices
        df = pd.read_excel(data_path)

        # Verify required columns exist
        required_cols = ['smiles', 'DRUG_ID', 'COSMIC_ID', 'LN_IC50']

        missing = [c for c in required_cols if c not in df.columns]
        if missing:
            raise ValueError(f"Missing required columns: {missing}")

        # Shuffle and reset index for proper positional indexing
        self.CDR = shuffle(df, random_state=self.rs).reset_index(drop=True)

        print(f"Loaded dataset: {len(self.CDR)} samples")
        print(f"  - Unique drugs: {self.CDR['DRUG_ID'].nunique()}")
        print(f"  - Unique cell lines (COSMIC_ID): {self.CDR['COSMIC_ID'].nunique()}")
        print(f"  - Unique SMILES: {self.CDR['smiles'].nunique()}")

    # ------------------------------------------------------------------ #
    #                         SPLIT FUNCTIONS                            #
    # ------------------------------------------------------------------ #

    def split_final(self, fold_idx=1, num_fold=10, return_train_y=False):
        """
        Stratified split by motif count (80/10/10).
        fold_idx is ignored for this method.
        """
        # Compute motif length for stratification
        unique_smiles = self.CDR['smiles'].unique()

        motif_map = {}
        for s in unique_smiles:
            mol = get_mol(s)
            if mol is not None:
                motif_map[s] = len(motif_decomp(mol))
            else:
                motif_map[s] = 0

        self.CDR['motif'] = self.CDR['smiles'].map(motif_map)
        labels = self.CDR['motif'].values

        # 90/10 split for train_val/test
        X_train_val, X_test = train_test_split(
            self.CDR, test_size=0.10, stratify=labels, random_state=42
        )

        # 80/10 split from train_val
        labels_train_val = X_train_val['motif'].values
        X_train, X_val = train_test_split(
            X_train_val, test_size=1/9, stratify=labels_train_val, random_state=42
        )

        train_idx = list(X_train.index)
        val_idx = list(X_val.index)
        test_idx = list(X_test.index)

        if return_train_y:
            return (X_train['LN_IC50'].values, X_test['smiles'].values), train_idx, val_idx, test_idx

        return train_idx, val_idx, test_idx

    def split_warm_start(self, fold_idx, num_fold=10, return_train_y=False):
        """Random K-fold split."""
        kf = KFold(n_splits=num_fold, shuffle=True, random_state=self.rs)

        for fold, (tr, te) in enumerate(kf.split(self.CDR), 1):
            if fold_idx == fold:
                tr_val_df = self.CDR.iloc[tr]
                te_df = self.CDR.iloc[te]

                val_ratio = 1 / (num_fold - 1)
                tr_df, val_df = train_test_split(
                    tr_val_df, test_size=val_ratio, random_state=self.rs
                )

                train_idx = list(tr_df.index)
                val_idx = list(val_df.index)
                test_idx = list(te_df.index)

                if return_train_y:
                    return (tr_df['LN_IC50'].values, te_df['smiles'].values), train_idx, val_idx, test_idx

                return train_idx, val_idx, test_idx

        raise ValueError(f"fold_idx {fold_idx} not found in {num_fold} folds")

    def split_cold_drug(self, fold_idx, num_fold=10, return_train_y=False):
        """Split by drug - no drug overlap between splits."""
        drugs = pd.DataFrame({"DRUG_ID": self.CDR.DRUG_ID.unique()})
        kf = KFold(n_splits=num_fold, shuffle=True, random_state=self.rs)

        for fold, (tr, te) in enumerate(kf.split(drugs), 1):
            if fold_idx == fold:
                tr_val_drugs = drugs.iloc[tr]
                te_drugs = drugs.iloc[te]

                tr_drugs, val_drugs = train_test_split(
                    tr_val_drugs, test_size=1/(num_fold-1), random_state=self.rs
                )

                tr_df = self.CDR[self.CDR.DRUG_ID.isin(tr_drugs.DRUG_ID)]
                val_df = self.CDR[self.CDR.DRUG_ID.isin(val_drugs.DRUG_ID)]
                te_df = self.CDR[self.CDR.DRUG_ID.isin(te_drugs.DRUG_ID)]

                train_idx = list(tr_df.index)
                val_idx = list(val_df.index)
                test_idx = list(te_df.index)

                if return_train_y:
                    return (tr_df['LN_IC50'].values, te_df['smiles'].values), train_idx, val_idx, test_idx

                return train_idx, val_idx, test_idx

        raise ValueError(f"fold_idx {fold_idx} not found in {num_fold} folds")

    def split_cold_cell(self, fold_idx, num_fold=10, return_train_y=False):
        """
        Split by cell line - no cell overlap between splits.

        Note: Splits on COSMIC_ID. Since COSMIC_ID has 1-to-1 mapping with
        CELL_LINE_NAME and SANGER_MODEL_ID (verified), this ensures no cell
        appears in multiple splits under any identifier.
        """
        cells = pd.DataFrame({"COSMIC_ID": self.CDR.COSMIC_ID.unique()})
        kf = KFold(n_splits=num_fold, shuffle=True, random_state=self.rs)

        for fold, (tr, te) in enumerate(kf.split(cells), 1):
            if fold_idx == fold:
                tr_val_cells = cells.iloc[tr]
                te_cells = cells.iloc[te]

                tr_cells, val_cells = train_test_split(
                    tr_val_cells, test_size=1/(num_fold-1), random_state=self.rs
                )

                tr_df = self.CDR[self.CDR.COSMIC_ID.isin(tr_cells.COSMIC_ID)]
                val_df = self.CDR[self.CDR.COSMIC_ID.isin(val_cells.COSMIC_ID)]
                te_df = self.CDR[self.CDR.COSMIC_ID.isin(te_cells.COSMIC_ID)]

                train_idx = list(tr_df.index)
                val_idx = list(val_df.index)
                test_idx = list(te_df.index)

                if return_train_y:
                    return (tr_df['LN_IC50'].values, te_df['smiles'].values), train_idx, val_idx, test_idx

                return train_idx, val_idx, test_idx

        raise ValueError(f"fold_idx {fold_idx} not found in {num_fold} folds")

    def split_cold_scaffold(self, fold_idx, num_fold=10, return_train_y=False):
        """Split by scaffold - no scaffold overlap between splits."""
        smiles_list = self.CDR['smiles'].unique()

        # Group SMILES by scaffold
        all_scaffolds = {}
        for smiles in smiles_list:
            scaffold = generate_scaffold(smiles, include_chirality=True)
            all_scaffolds.setdefault(scaffold, []).append(smiles)

        all_scaffolds = {k: sorted(v) for k, v in all_scaffolds.items()}
        scaffold_keys = list(all_scaffolds.keys())

        kf = KFold(n_splits=num_fold, random_state=self.rs, shuffle=True)

        for fold, (train_scaffold_idx, test_scaffold_idx) in enumerate(kf.split(scaffold_keys), 1):
            if fold == fold_idx:
                # Split train scaffolds into train and validation
                train_scaffold_idx, val_scaffold_idx = train_test_split(
                    train_scaffold_idx,
                    test_size=1/(num_fold-1),
                    random_state=self.rs
                )

                # Map scaffold indices to SMILES
                train_smiles = set()
                val_smiles = set()
                test_smiles = set()

                for idx in train_scaffold_idx:
                    train_smiles.update(all_scaffolds[scaffold_keys[idx]])
                for idx in val_scaffold_idx:
                    val_smiles.update(all_scaffolds[scaffold_keys[idx]])
                for idx in test_scaffold_idx:
                    test_smiles.update(all_scaffolds[scaffold_keys[idx]])

                tr_df = self.CDR[self.CDR['smiles'].isin(train_smiles)]
                val_df = self.CDR[self.CDR['smiles'].isin(val_smiles)]
                te_df = self.CDR[self.CDR['smiles'].isin(test_smiles)]

                train_idx = list(tr_df.index)
                val_idx = list(val_df.index)
                test_idx = list(te_df.index)

                if return_train_y:
                    return (tr_df['LN_IC50'].values, te_df['smiles'].values), train_idx, val_idx, test_idx

                return train_idx, val_idx, test_idx

        raise ValueError(f"fold_idx {fold_idx} not found in {num_fold} folds")

    def split_dataset(self, fold_idx=1, num_folds=10, return_train_y=False):
        """Main entry point - dispatches to appropriate split method."""
        if self.scenarios == 'warm_start':
            return self.split_warm_start(fold_idx, num_folds, return_train_y)
        elif self.scenarios == 'final':
            return self.split_final(fold_idx, num_folds, return_train_y)
        elif self.scenarios == 'cold_drug':
            return self.split_cold_drug(fold_idx, num_folds, return_train_y)
        elif self.scenarios == 'cold_cell':
            return self.split_cold_cell(fold_idx, num_folds, return_train_y)
        elif self.scenarios == 'cold_scaffold':
            return self.split_cold_scaffold(fold_idx, num_folds, return_train_y)
        else:
            raise ValueError(f"Invalid scenario: {self.scenarios}. "
                           f"Choose from: warm_start, final, cold_drug, cold_cell, cold_scaffold")


if __name__ == "__main__":
    # Example usage
    splitter = DataSplitter(
        scenarios="warm_start",
        data_path="./data/Drug Screening - IC50s/GDSC1_IC50_with_pubchem_and_smiles_filtered_by_mut_rna_mrna_Dec29.xlsx",
    )
    train_idx, val_idx, test_idx = splitter.split_dataset(fold_idx=1, num_folds=10)
    print(f"\nSplit results:")
    print(f"  Train: {len(train_idx)} samples")
    print(f"  Val:   {len(val_idx)} samples")
    print(f"  Test:  {len(test_idx)} samples")
