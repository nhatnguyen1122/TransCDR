#!/usr/bin/env python3
"""
Script to verify your dataset structure and identify required modifications
"""

import pandas as pd
import sys

# Update this path to your actual file
DATA_PATH = "./data/Drug Screening - IC50s/GDSC1_IC50_with_pubchem_and_smiles_filtered_by_mut_rna_mrna_Dec29.xlsx"

try:
    print("=" * 80)
    print("DATASET VERIFICATION FOR SPLITTER ADAPTATION")
    print("=" * 80)

    # Load the dataset
    df = pd.read_excel(DATA_PATH)

    print(f"\n1. Dataset Shape: {df.shape[0]} rows × {df.shape[1]} columns")

    print(f"\n2. Column Names:")
    for i, col in enumerate(df.columns, 1):
        print(f"   {i}. {col}")

    print(f"\n3. Required Columns Check:")
    required_cols = ['DRUG_ID', 'COSMIC_ID', 'smiles', 'LN_IC50',
                     'CELL_LINE_NAME', 'SANGER_MODEL_ID']
    for col in required_cols:
        status = "✓ FOUND" if col in df.columns else "✗ MISSING"
        print(f"   {col}: {status}")

    print(f"\n4. Sample Data (first 3 rows):")
    print(df.head(3).to_string())

    print(f"\n5. Missing Value Analysis:")
    missing_cols = ['COSMIC_ID', 'CELL_LINE_NAME', 'SANGER_MODEL_ID']
    for col in missing_cols:
        if col in df.columns:
            missing_count = df[col].isna().sum()
            missing_pct = (missing_count / len(df)) * 100
            print(f"   {col}: {missing_count} missing ({missing_pct:.2f}%)")

    print(f"\n6. Cell Identifier Uniqueness:")
    if 'COSMIC_ID' in df.columns:
        print(f"   Unique COSMIC_IDs: {df['COSMIC_ID'].nunique()}")
    if 'CELL_LINE_NAME' in df.columns:
        print(f"   Unique CELL_LINE_NAMEs: {df['CELL_LINE_NAME'].nunique()}")
    if 'SANGER_MODEL_ID' in df.columns:
        print(f"   Unique SANGER_MODEL_IDs: {df['SANGER_MODEL_ID'].nunique()}")

    print(f"\n7. Cell Identifier Relationship:")
    if all(col in df.columns for col in ['COSMIC_ID', 'CELL_LINE_NAME', 'SANGER_MODEL_ID']):
        # Check if these map 1-to-1
        mapping_df = df[['COSMIC_ID', 'CELL_LINE_NAME', 'SANGER_MODEL_ID']].drop_duplicates()
        print(f"   Unique combinations: {len(mapping_df)}")

        # Check if each COSMIC_ID has unique CELL_LINE_NAME and SANGER_MODEL_ID
        cosmic_to_cell = df.groupby('COSMIC_ID')['CELL_LINE_NAME'].nunique()
        cosmic_to_sanger = df.groupby('COSMIC_ID')['SANGER_MODEL_ID'].nunique()

        if (cosmic_to_cell == 1).all() and (cosmic_to_sanger == 1).all():
            print(f"   ✓ One-to-one mapping: Each COSMIC_ID maps to exactly one CELL_LINE_NAME and SANGER_MODEL_ID")
        else:
            print(f"   ✗ WARNING: COSMIC_ID does NOT map 1-to-1 to other identifiers!")
            print(f"     - COSMIC_IDs with multiple CELL_LINE_NAMEs: {(cosmic_to_cell > 1).sum()}")
            print(f"     - COSMIC_IDs with multiple SANGER_MODEL_IDs: {(cosmic_to_sanger > 1).sum()}")

    print(f"\n8. Drug Information:")
    if 'DRUG_ID' in df.columns:
        print(f"   Unique drugs: {df['DRUG_ID'].nunique()}")
        print(f"   Total drug-cell combinations: {len(df)}")

    print("\n" + "=" * 80)
    print("VERIFICATION COMPLETE")
    print("=" * 80)

except FileNotFoundError:
    print(f"ERROR: File not found at {DATA_PATH}")
    print("Please update the DATA_PATH variable in this script to point to your file.")
    sys.exit(1)
except Exception as e:
    print(f"ERROR: {str(e)}")
    sys.exit(1)
