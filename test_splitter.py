#!/usr/bin/env python3
"""
Test script for splitter.py

This script tests all splitting scenarios to ensure they work correctly
with your dataset.
"""

import sys
from splitter import DataSplitter

# Data path
DATA_PATH = "./data/Drug Screening - IC50s/GDSC1_IC50_with_pubchem_and_smiles_filtered_by_mut_rna_mrna_Dec29.xlsx"

def test_scenario(scenario_name, fold_idx=1, num_folds=10):
    """Test a specific splitting scenario."""
    print(f"\n{'='*70}")
    print(f"Testing: {scenario_name}")
    print(f"{'='*70}")

    try:
        splitter = DataSplitter(
            scenarios=scenario_name,
            data_path=DATA_PATH,
            n_folds=num_folds,
            random_state=2022
        )

        train_idx, val_idx, test_idx = splitter.split_dataset(
            fold_idx=fold_idx,
            num_folds=num_folds
        )

        # Get the actual data to verify no overlap
        train_set = set(train_idx)
        val_set = set(val_idx)
        test_set = set(test_idx)

        # Verify no overlap
        assert len(train_set & val_set) == 0, "Train and validation sets overlap!"
        assert len(train_set & test_set) == 0, "Train and test sets overlap!"
        assert len(val_set & test_set) == 0, "Validation and test sets overlap!"

        total = len(train_idx) + len(val_idx) + len(test_idx)

        print(f"✓ Split successful!")
        print(f"  Train:      {len(train_idx):6d} samples ({len(train_idx)/total*100:5.2f}%)")
        print(f"  Validation: {len(val_idx):6d} samples ({len(val_idx)/total*100:5.2f}%)")
        print(f"  Test:       {len(test_idx):6d} samples ({len(test_idx)/total*100:5.2f}%)")
        print(f"  Total:      {total:6d} samples")
        print(f"✓ No overlap between splits")

        # Additional scenario-specific checks
        if scenario_name == "cold_drug":
            # Verify no drug overlap
            train_drugs = set(splitter.CDR.iloc[train_idx]['DRUG_ID'])
            test_drugs = set(splitter.CDR.iloc[test_idx]['DRUG_ID'])
            val_drugs = set(splitter.CDR.iloc[val_idx]['DRUG_ID'])

            assert len(train_drugs & test_drugs) == 0, "Drugs overlap between train and test!"
            assert len(train_drugs & val_drugs) == 0, "Drugs overlap between train and val!"
            assert len(val_drugs & test_drugs) == 0, "Drugs overlap between val and test!"

            print(f"✓ No drug overlap verified")
            print(f"  Train drugs: {len(train_drugs)}")
            print(f"  Val drugs:   {len(val_drugs)}")
            print(f"  Test drugs:  {len(test_drugs)}")

        elif scenario_name == "cold_cell":
            # Verify no cell overlap
            train_cells = set(splitter.CDR.iloc[train_idx]['COSMIC_ID'])
            test_cells = set(splitter.CDR.iloc[test_idx]['COSMIC_ID'])
            val_cells = set(splitter.CDR.iloc[val_idx]['COSMIC_ID'])

            assert len(train_cells & test_cells) == 0, "Cells overlap between train and test!"
            assert len(train_cells & val_cells) == 0, "Cells overlap between train and val!"
            assert len(val_cells & test_cells) == 0, "Cells overlap between val and test!"

            print(f"✓ No cell overlap verified")
            print(f"  Train cells: {len(train_cells)}")
            print(f"  Val cells:   {len(val_cells)}")
            print(f"  Test cells:  {len(test_cells)}")

        elif scenario_name == "cold_scaffold":
            # Verify no scaffold overlap
            from splitter import generate_scaffold

            def get_scaffolds(indices):
                smiles_list = splitter.CDR.iloc[indices]['smiles'].unique()
                return set(generate_scaffold(s, include_chirality=True) for s in smiles_list)

            train_scaffolds = get_scaffolds(train_idx)
            val_scaffolds = get_scaffolds(val_idx)
            test_scaffolds = get_scaffolds(test_idx)

            assert len(train_scaffolds & test_scaffolds) == 0, "Scaffolds overlap between train and test!"
            assert len(train_scaffolds & val_scaffolds) == 0, "Scaffolds overlap between train and val!"
            assert len(val_scaffolds & test_scaffolds) == 0, "Scaffolds overlap between val and test!"

            print(f"✓ No scaffold overlap verified")
            print(f"  Train scaffolds: {len(train_scaffolds)}")
            print(f"  Val scaffolds:   {len(val_scaffolds)}")
            print(f"  Test scaffolds:  {len(test_scaffolds)}")

        return True

    except FileNotFoundError:
        print(f"✗ Error: Data file not found at {DATA_PATH}")
        print(f"  Please ensure the file exists before running this test.")
        return False
    except Exception as e:
        print(f"✗ Error: {str(e)}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run all tests."""
    print("="*70)
    print("SPLITTER TEST SUITE")
    print("="*70)

    scenarios = [
        "warm_start",
        "cold_drug",
        "cold_cell",
        "cold_scaffold",
        "final"
    ]

    results = {}

    for scenario in scenarios:
        results[scenario] = test_scenario(scenario, fold_idx=1, num_folds=10)

    # Summary
    print(f"\n{'='*70}")
    print("TEST SUMMARY")
    print(f"{'='*70}")

    all_passed = all(results.values())

    for scenario, passed in results.items():
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"  {scenario:20s}: {status}")

    print(f"{'='*70}")

    if all_passed:
        print("✓ All tests passed!")
        return 0
    else:
        print("✗ Some tests failed. Please check the output above.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
