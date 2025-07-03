#!/usr/bin/env python3
"""
CSV to SDF Converter for ChemBL Database
Converts a CSV file containing SMILES strings to SDF format.
"""

import pandas as pd
import sys
import os
from rdkit import Chem
from rdkit.Chem import AllChem
import argparse

def csv_to_sdf(csv_file, sdf_file, smiles_column="Smiles", id_column="Molecule ChEMBL ID",
               name_column="Molecule Name", mw_column="Molecular Weight", alogp_column="AlogP"):
    """
    Convert CSV file with SMILES to SDF format.

    Args:
        csv_file (str): Input CSV file path
        sdf_file (str): Output SDF file path
        smiles_column (str): Column name containing SMILES strings
        id_column (str): Column name for molecule ID
        name_column (str): Column name for molecule name
        mw_column (str): Column name for molecular weight
        alogp_column (str): Column name for AlogP
    """

    try:
        # Read CSV file
        print(f"Reading CSV file: {csv_file}")
        df = pd.read_csv(csv_file, sep=';', quotechar='"')

        print(f"Found {len(df)} rows in CSV file")
        print(f"Columns: {list(df.columns)}")

        # Check if required columns exist
        if smiles_column not in df.columns:
            print(f"Error: Column '{smiles_column}' not found in CSV file")
            print(f"Available columns: {list(df.columns)}")
            return False

        # Remove rows with empty SMILES
        df = df.dropna(subset=[smiles_column])
        df = df[df[smiles_column].str.strip() != '']

        print(f"After filtering empty SMILES: {len(df)} rows")

        # Create SDF writer
        writer = Chem.SDWriter(sdf_file)

        successful_conversions = 0
        failed_conversions = 0

        # Process each row
        for idx, row in df.iterrows():
            try:
                smiles = str(row[smiles_column]).strip()

                # Create molecule from SMILES
                mol = Chem.MolFromSmiles(smiles)

                if mol is None:
                    print(f"Warning: Could not parse SMILES at row {idx + 1}: {smiles}")
                    failed_conversions += 1
                    continue

                # Add 2D coordinates
                AllChem.Compute2DCoords(mol)

                # Set molecule properties
                if id_column in df.columns and pd.notna(row[id_column]):
                    mol.SetProp("ChEMBL_ID", str(row[id_column]))

                if name_column in df.columns and pd.notna(row[name_column]):
                    mol.SetProp("Name", str(row[name_column]))
                else:
                    # Use ChEMBL ID as name if name is not available
                    if id_column in df.columns and pd.notna(row[id_column]):
                        mol.SetProp("Name", str(row[id_column]))
                    else:
                        mol.SetProp("Name", f"Molecule_{idx + 1}")

                if mw_column in df.columns and pd.notna(row[mw_column]):
                    mol.SetProp("Molecular_Weight", str(row[mw_column]))

                if alogp_column in df.columns and pd.notna(row[alogp_column]):
                    mol.SetProp("AlogP", str(row[alogp_column]))

                # Add SMILES as property
                mol.SetProp("SMILES", smiles)

                # Add other available properties
                for col in df.columns:
                    if col not in [smiles_column, id_column, name_column, mw_column, alogp_column]:
                        if pd.notna(row[col]) and str(row[col]).strip() != '':
                            # Clean column name for SDF property
                            prop_name = col.replace(' ', '_').replace('#', 'Num_')
                            mol.SetProp(prop_name, str(row[col]))

                # Write molecule to SDF
                writer.write(mol)
                successful_conversions += 1

                if successful_conversions % 100 == 0:
                    print(f"Processed {successful_conversions} molecules...")

            except Exception as e:
                print(f"Error processing row {idx + 1}: {e}")
                failed_conversions += 1
                continue

        writer.close()

        print(f"\nConversion completed!")
        print(f"Successfully converted: {successful_conversions} molecules")
        print(f"Failed conversions: {failed_conversions}")
        print(f"Output SDF file: {sdf_file}")

        return True

    except Exception as e:
        print(f"Error: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description='Convert CSV file with SMILES to SDF format')
    parser.add_argument('input_csv', help='Input CSV file path')
    parser.add_argument('-o', '--output', help='Output SDF file path (default: input_name.sdf)')
    parser.add_argument('--smiles-col', default='Smiles', help='Column name containing SMILES (default: Smiles)')
    parser.add_argument('--id-col', default='Molecule ChEMBL ID', help='Column name for molecule ID (default: Molecule ChEMBL ID)')
    parser.add_argument('--name-col', default='Molecule Name', help='Column name for molecule name (default: Molecule Name)')
    parser.add_argument('--mw-col', default='Molecular Weight', help='Column name for molecular weight (default: Molecular Weight)')
    parser.add_argument('--alogp-col', default='AlogP', help='Column name for AlogP (default: AlogP)')

    args = parser.parse_args()

    # Determine output file name
    if args.output:
        sdf_file = args.output
    else:
        base_name = os.path.splitext(args.input_csv)[0]
        sdf_file = f"{base_name}.sdf"

    # Convert CSV to SDF
    success = csv_to_sdf(
        args.input_csv,
        sdf_file,
        smiles_column=args.smiles_col,
        id_column=args.id_col,
        name_column=args.name_col,
        mw_column=args.mw_col,
        alogp_column=args.alogp_col
    )

    if success:
        print(f"\nSuccessfully converted {args.input_csv} to {sdf_file}")
    else:
        print(f"\nFailed to convert {args.input_csv}")
        sys.exit(1)

if __name__ == "__main__":
    main()