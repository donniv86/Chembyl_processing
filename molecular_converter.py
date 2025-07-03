#!/usr/bin/env python3
"""
General Molecular Format Converter
Converts between various molecular file formats including CSV with SMILES, SDF, MOL, MOL2, etc.
"""

import pandas as pd
import sys
import os
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
from pathlib import Path

class MolecularConverter:
    """General molecular format converter supporting multiple input/output formats."""

    SUPPORTED_INPUT_FORMATS = ['csv', 'sdf', 'mol', 'mol2']
    SUPPORTED_OUTPUT_FORMATS = ['sdf', 'mol', 'mol2', 'pdb', 'xyz', 'csv']

    def __init__(self):
        self.molecules = []
        self.properties = []

    def read_csv(self, csv_file, smiles_column="SMILES", id_column=None,
                 name_column=None, separator=",", quotechar='"'):
        """Read molecules from CSV file containing SMILES."""
        try:
            print(f"Reading CSV file: {csv_file}")
            df = pd.read_csv(csv_file, sep=separator, quotechar=quotechar)

            print(f"Found {len(df)} rows in CSV file")
            print(f"Columns: {list(df.columns)}")

            if smiles_column not in df.columns:
                print(f"Error: Column '{smiles_column}' not found in CSV file")
                return False

            # Remove rows with empty SMILES
            df = df.dropna(subset=[smiles_column])
            df = df[df[smiles_column].str.strip() != '']

            print(f"After filtering empty SMILES: {len(df)} rows")

            successful = 0
            failed = 0

            for idx, row in df.iterrows():
                try:
                    smiles = str(row[smiles_column]).strip()
                    mol = Chem.MolFromSmiles(smiles)

                    if mol is None:
                        print(f"Warning: Could not parse SMILES at row {idx + 1}: {smiles}")
                        failed += 1
                        continue

                    # Add 2D coordinates
                    AllChem.Compute2DCoords(mol)

                    # Set molecule properties
                    props = {}
                    for col in df.columns:
                        if pd.notna(row[col]) and str(row[col]).strip() != '':
                            props[col] = str(row[col])

                    # Set name
                    if name_column and name_column in df.columns and pd.notna(row[name_column]):
                        mol.SetProp("Name", str(row[name_column]))
                    elif id_column and id_column in df.columns and pd.notna(row[id_column]):
                        mol.SetProp("Name", str(row[id_column]))
                    else:
                        mol.SetProp("Name", f"Molecule_{idx + 1}")

                    self.molecules.append(mol)
                    self.properties.append(props)
                    successful += 1

                except Exception as e:
                    print(f"Error processing row {idx + 1}: {e}")
                    failed += 1
                    continue

            print(f"Successfully loaded: {successful} molecules")
            print(f"Failed conversions: {failed}")
            return True

        except Exception as e:
            print(f"Error reading CSV file: {e}")
            return False

    def read_sdf(self, sdf_file):
        """Read molecules from SDF file."""
        try:
            print(f"Reading SDF file: {sdf_file}")
            suppl = Chem.SDMolSupplier(sdf_file)

            successful = 0
            failed = 0

            for mol in suppl:
                if mol is not None:
                    self.molecules.append(mol)
                    # Extract properties
                    props = {}
                    for prop_name in mol.GetPropNames():
                        props[prop_name] = mol.GetProp(prop_name)
                    self.properties.append(props)
                    successful += 1
                else:
                    failed += 1

            print(f"Successfully loaded: {successful} molecules")
            print(f"Failed conversions: {failed}")
            return True

        except Exception as e:
            print(f"Error reading SDF file: {e}")
            return False

    def read_mol(self, mol_file):
        """Read single molecule from MOL file."""
        try:
            print(f"Reading MOL file: {mol_file}")
            mol = Chem.MolFromMolFile(mol_file)

            if mol is not None:
                self.molecules.append(mol)
                props = {}
                for prop_name in mol.GetPropNames():
                    props[prop_name] = mol.GetProp(prop_name)
                self.properties.append(props)
                print("Successfully loaded 1 molecule")
                return True
            else:
                print("Error: Could not read MOL file")
                return False

        except Exception as e:
            print(f"Error reading MOL file: {e}")
            return False

    def read_mol2(self, mol2_file):
        """Read molecules from MOL2 file."""
        try:
            print(f"Reading MOL2 file: {mol2_file}")
            suppl = Chem.MolFromMol2File(mol2_file)

            if suppl is not None:
                if isinstance(suppl, list):
                    for mol in suppl:
                        if mol is not None:
                            self.molecules.append(mol)
                            self.properties.append({})
                else:
                    self.molecules.append(suppl)
                    self.properties.append({})

                print(f"Successfully loaded {len(self.molecules)} molecules")
                return True
            else:
                print("Error: Could not read MOL2 file")
                return False

        except Exception as e:
            print(f"Error reading MOL2 file: {e}")
            return False

    def write_sdf(self, output_file):
        """Write molecules to SDF file."""
        try:
            writer = Chem.SDWriter(output_file)
            for i, mol in enumerate(self.molecules):
                # Add properties back
                if i < len(self.properties):
                    for key, value in self.properties[i].items():
                        mol.SetProp(key, value)
                writer.write(mol)
            writer.close()
            print(f"Successfully wrote {len(self.molecules)} molecules to {output_file}")
            return True
        except Exception as e:
            print(f"Error writing SDF file: {e}")
            return False

    def write_mol(self, output_file):
        """Write molecules to individual MOL files."""
        try:
            base_name = Path(output_file).stem
            base_dir = Path(output_file).parent

            for i, mol in enumerate(self.molecules):
                mol_file = base_dir / f"{base_name}_{i+1}.mol"
                Chem.MolToMolFile(mol, str(mol_file))

            print(f"Successfully wrote {len(self.molecules)} molecules to individual MOL files")
            return True
        except Exception as e:
            print(f"Error writing MOL files: {e}")
            return False

    def write_mol2(self, output_file):
        """Write molecules to MOL2 file."""
        try:
            with open(output_file, 'w') as f:
                for i, mol in enumerate(self.molecules):
                    mol2_block = Chem.MolToMol2Block(mol)
                    f.write(mol2_block)
                    f.write('\n')

            print(f"Successfully wrote {len(self.molecules)} molecules to {output_file}")
            return True
        except Exception as e:
            print(f"Error writing MOL2 file: {e}")
            return False

    def write_pdb(self, output_file):
        """Write molecules to PDB file."""
        try:
            with open(output_file, 'w') as f:
                for i, mol in enumerate(self.molecules):
                    # Add 3D coordinates if not present
                    if mol.GetNumConformers() == 0:
                        AllChem.EmbedMolecule(mol, randomSeed=42)
                        AllChem.MMFFOptimizeMolecule(mol)

                    pdb_block = Chem.MolToPDBBlock(mol)
                    f.write(pdb_block)
                    f.write('\n')

            print(f"Successfully wrote {len(self.molecules)} molecules to {output_file}")
            return True
        except Exception as e:
            print(f"Error writing PDB file: {e}")
            return False

    def write_xyz(self, output_file):
        """Write molecules to XYZ file."""
        try:
            with open(output_file, 'w') as f:
                for i, mol in enumerate(self.molecules):
                    # Add 3D coordinates if not present
                    if mol.GetNumConformers() == 0:
                        AllChem.EmbedMolecule(mol, randomSeed=42)
                        AllChem.MMFFOptimizeMolecule(mol)

                    xyz_block = Chem.MolToXYZBlock(mol)
                    f.write(xyz_block)
                    f.write('\n')

            print(f"Successfully wrote {len(self.molecules)} molecules to {output_file}")
            return True
        except Exception as e:
            print(f"Error writing XYZ file: {e}")
            return False

    def write_csv(self, output_file, smiles_column="SMILES"):
        """Write molecules to CSV file with SMILES."""
        try:
            data = []
            for i, mol in enumerate(self.molecules):
                row = {smiles_column: Chem.MolToSmiles(mol)}

                # Add properties
                if i < len(self.properties):
                    row.update(self.properties[i])

                data.append(row)

            df = pd.DataFrame(data)
            df.to_csv(output_file, index=False)
            print(f"Successfully wrote {len(self.molecules)} molecules to {output_file}")
            return True
        except Exception as e:
            print(f"Error writing CSV file: {e}")
            return False

    def convert(self, input_file, output_file, input_format=None, output_format=None):
        """Main conversion method."""
        # Auto-detect formats if not specified
        if input_format is None:
            input_format = Path(input_file).suffix.lower().lstrip('.')

        if output_format is None:
            output_format = Path(output_file).suffix.lower().lstrip('.')

        # Validate formats
        if input_format not in self.SUPPORTED_INPUT_FORMATS:
            print(f"Error: Unsupported input format '{input_format}'")
            print(f"Supported formats: {self.SUPPORTED_INPUT_FORMATS}")
            return False

        if output_format not in self.SUPPORTED_OUTPUT_FORMATS:
            print(f"Error: Unsupported output format '{output_format}'")
            print(f"Supported formats: {self.SUPPORTED_OUTPUT_FORMATS}")
            return False

        # Read input file
        if input_format == 'csv':
            success = self.read_csv(input_file)
        elif input_format == 'sdf':
            success = self.read_sdf(input_file)
        elif input_format == 'mol':
            success = self.read_mol(input_file)
        elif input_format == 'mol2':
            success = self.read_mol2(input_file)

        if not success:
            return False

        # Write output file
        if output_format == 'sdf':
            return self.write_sdf(output_file)
        elif output_format == 'mol':
            return self.write_mol(output_file)
        elif output_format == 'mol2':
            return self.write_mol2(output_file)
        elif output_format == 'pdb':
            return self.write_pdb(output_file)
        elif output_format == 'xyz':
            return self.write_xyz(output_file)
        elif output_format == 'csv':
            return self.write_csv(output_file)

def main():
    parser = argparse.ArgumentParser(
        description='Convert between molecular file formats',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Convert CSV to SDF
  python molecular_converter.py input.csv output.sdf

  # Convert SDF to MOL2
  python molecular_converter.py input.sdf output.mol2

  # Convert CSV to individual MOL files
  python molecular_converter.py input.csv output.mol

  # Convert SDF to CSV with SMILES
  python molecular_converter.py input.sdf output.csv

  # Specify formats explicitly
  python molecular_converter.py input.csv output.mol2 --input-format csv --output-format mol2
        """
    )

    parser.add_argument('input_file', help='Input file path')
    parser.add_argument('output_file', help='Output file path')
    parser.add_argument('--input-format', choices=['csv', 'sdf', 'mol', 'mol2'],
                       help='Input file format (auto-detected from extension if not specified)')
    parser.add_argument('--output-format', choices=['sdf', 'mol', 'mol2', 'pdb', 'xyz', 'csv'],
                       help='Output file format (auto-detected from extension if not specified)')

    # CSV-specific options
    parser.add_argument('--smiles-col', default='SMILES', help='SMILES column name for CSV input')
    parser.add_argument('--id-col', help='ID column name for CSV input')
    parser.add_argument('--name-col', help='Name column name for CSV input')
    parser.add_argument('--separator', default=',', help='CSV separator (default: ,)')
    parser.add_argument('--quotechar', default='"', help='CSV quote character (default: ")')

    args = parser.parse_args()

    # Check if input file exists
    if not os.path.exists(args.input_file):
        print(f"Error: Input file '{args.input_file}' does not exist")
        sys.exit(1)

    # Create converter and perform conversion
    converter = MolecularConverter()

    # Set CSV options if needed
    if args.input_format == 'csv' or (args.input_format is None and args.input_file.lower().endswith('.csv')):
        converter.read_csv = lambda f: converter.read_csv(
            f, smiles_column=args.smiles_col, id_column=args.id_col,
            name_column=args.name_col, separator=args.separator, quotechar=args.quotechar
        )

    success = converter.convert(
        args.input_file,
        args.output_file,
        args.input_format,
        args.output_format
    )

    if success:
        print(f"\n✅ Successfully converted {args.input_file} to {args.output_file}")
    else:
        print(f"\n❌ Failed to convert {args.input_file}")
        sys.exit(1)

if __name__ == "__main__":
    main()