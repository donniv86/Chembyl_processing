#!/usr/bin/env python3
"""
Simple script to convert jiang.csv to SDF format
"""

from csv_to_sdf_converter import csv_to_sdf

def main():
    input_file = "jiang.csv"
    output_file = "jiang.sdf"

    print("Converting jiang.csv to SDF format...")

    success = csv_to_sdf(
        csv_file=input_file,
        sdf_file=output_file,
        smiles_column="Smiles",
        id_column="Molecule ChEMBL ID",
        name_column="Molecule Name",
        mw_column="Molecular Weight",
        alogp_column="AlogP"
    )

    if success:
        print(f"\n✅ Successfully converted {input_file} to {output_file}")
        print(f"📁 Output file: {output_file}")
    else:
        print(f"\n❌ Failed to convert {input_file}")

if __name__ == "__main__":
    main()