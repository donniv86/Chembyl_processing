#!/usr/bin/env python3
"""
General Molecular Format Converter
Converts between various molecular file formats including CSV with SMILES, SDF, MOL, MOL2, etc.
Enhanced with robust error handling and encoding detection.
"""

import pandas as pd
import sys
import os
import argparse
import chardet
from rdkit import Chem
from rdkit.Chem import AllChem
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

class MolecularConverter:
    """General molecular format converter supporting multiple input/output formats."""

    SUPPORTED_INPUT_FORMATS = ['csv', 'tsv', 'xls', 'xlsx', 'sdf', 'mol', 'mol2']
    SUPPORTED_OUTPUT_FORMATS = ['sdf', 'mol', 'mol2', 'pdb', 'xyz', 'csv']

    # Common encodings to try
    ENCODINGS_TO_TRY = ['utf-8', 'cp1252', 'latin-1', 'iso-8859-1', 'utf-8-sig']

    SMILES_CANDIDATES = [
        'smiles', 'canonical_smiles', 'isomeric_smiles', 'structure', 'structure_smiles', 'mol_smiles', 'smile'
    ]

    def __init__(self):
        self.molecules = []
        self.properties = []

    def detect_encoding(self, file_path):
        """Detect the encoding of a file."""
        try:
            with open(file_path, 'rb') as f:
                raw_data = f.read(10000)  # Read first 10KB for detection
                result = chardet.detect(raw_data)
                detected_encoding = result['encoding']
                confidence = result['confidence']

                print(f"Detected encoding: {detected_encoding} (confidence: {confidence:.2f})")
                return detected_encoding
        except Exception as e:
            print(f"Warning: Could not detect encoding: {e}")
            return 'utf-8'

    def try_read_csv(self, csv_file, **kwargs):
        """Try to read CSV file with different encodings and format detection."""
        encodings_to_try = [kwargs.get('encoding')] if kwargs.get('encoding') else self.ENCODINGS_TO_TRY

        # Add detected encoding to the list
        detected_encoding = self.detect_encoding(csv_file)
        if detected_encoding and detected_encoding not in encodings_to_try:
            encodings_to_try.insert(0, detected_encoding)

        for encoding in encodings_to_try:
            try:
                delim = kwargs.get('sep', ',')
                quote = kwargs.get('quotechar', '"')
                print(f"Trying to read with encoding: {encoding}, delimiter: '{delim}', quote: '{quote}'")
                df = pd.read_csv(csv_file, encoding=encoding, sep=delim, quotechar=quote)
                print(f"Successfully read with encoding: {encoding}")
                return df
            except UnicodeDecodeError as e:
                print(f"Failed with encoding {encoding}: {e}")
                continue
            except Exception as e:
                print(f"Error with encoding {encoding}: {e}")
                continue

        raise ValueError(f"Could not read file with any of the attempted encodings: {encodings_to_try}")

    def try_read_excel(self, excel_file, **kwargs):
        encodings_to_try = [kwargs.get('encoding')] if kwargs.get('encoding') else self.ENCODINGS_TO_TRY
        detected_encoding = self.detect_encoding(excel_file)
        if detected_encoding and detected_encoding not in encodings_to_try:
            encodings_to_try.insert(0, detected_encoding)
        for encoding in encodings_to_try:
            try:
                print(f"Trying to read Excel with encoding: {encoding}")
                df = pd.read_excel(excel_file, encoding=encoding)
                print(f"Successfully read Excel with encoding: {encoding}")
                return df
            except Exception as e:
                print(f"Error reading Excel with encoding {encoding}: {e}")
                continue
        raise ValueError(f"Could not read Excel file with any of the attempted encodings: {encodings_to_try}")

    def find_smiles_column(self, columns):
        for candidate in self.SMILES_CANDIDATES:
            for col in columns:
                if col.strip().lower() == candidate:
                    print(f"Auto-detected SMILES column: '{col}'")
                    return col
        print(f"Could not auto-detect SMILES column. Columns found: {list(columns)}")
        return None

    def read_csv(self, csv_file, smiles_column=None, id_column=None, name_column=None, separator=',', quotechar='"', encoding=None):
        """Read molecules from CSV file containing SMILES with robust encoding handling."""
        try:
            print(f"Reading CSV/TSV file: {csv_file}")
            df = self.try_read_csv(csv_file, sep=separator, quotechar=quotechar, encoding=encoding)
            print(f"Found {len(df)} rows in file")
            print(f"Columns: {list(df.columns)}")
            if not smiles_column:
                smiles_column = self.find_smiles_column(df.columns)
            if not smiles_column or smiles_column not in df.columns:
                print(f"Error: Could not find SMILES column. Please specify with --smiles-col. Columns: {list(df.columns)}")
                return False
            df = df.dropna(subset=[smiles_column])
            df = df[df[smiles_column].astype(str).str.strip() != '']
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
                    try:
                        AllChem.Compute2DCoords(mol)
                    except Exception as e:
                        print(f"Warning: Could not compute 2D coordinates for row {idx + 1}: {e}")
                    props = {}
                    for col in df.columns:
                        if pd.notna(row[col]) and str(row[col]).strip() != '':
                            props[col] = str(row[col])
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
            print(f"Error reading CSV/TSV file: {e}")
            return False

    def read_excel(self, excel_file, smiles_column=None, id_column=None, name_column=None, encoding=None):
        try:
            print(f"Reading Excel file: {excel_file}")
            df = self.try_read_excel(excel_file, encoding=encoding)
            print(f"Found {len(df)} rows in file")
            print(f"Columns: {list(df.columns)}")
            if not smiles_column:
                smiles_column = self.find_smiles_column(df.columns)
            if not smiles_column or smiles_column not in df.columns:
                print(f"Error: Could not find SMILES column. Please specify with --smiles-col. Columns: {list(df.columns)}")
                return False
            df = df.dropna(subset=[smiles_column])
            df = df[df[smiles_column].astype(str).str.strip() != '']
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
                    try:
                        AllChem.Compute2DCoords(mol)
                    except Exception as e:
                        print(f"Warning: Could not compute 2D coordinates for row {idx + 1}: {e}")
                    props = {}
                    for col in df.columns:
                        if pd.notna(row[col]) and str(row[col]).strip() != '':
                            props[col] = str(row[col])
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
            print(f"Error reading Excel file: {e}")
            return False

    def read_sdf(self, sdf_file):
        """Read molecules from SDF file with error handling."""
        try:
            print(f"Reading SDF file: {sdf_file}")
            suppl = Chem.SDMolSupplier(sdf_file, removeHs=False, sanitize=True)

            successful = 0
            failed = 0

            for i, mol in enumerate(suppl):
                if mol is not None:
                    self.molecules.append(mol)
                    # Extract properties
                    props = {}
                    for prop_name in mol.GetPropNames():
                        props[prop_name] = mol.GetProp(prop_name)
                    self.properties.append(props)
                    successful += 1
                else:
                    print(f"Warning: Could not read molecule {i + 1} from SDF")
                    failed += 1

            print(f"Successfully loaded: {successful} molecules")
            print(f"Failed conversions: {failed}")
            return True

        except Exception as e:
            print(f"Error reading SDF file: {e}")
            return False

    def read_mol(self, mol_file):
        """Read single molecule from MOL file with error handling."""
        try:
            print(f"Reading MOL file: {mol_file}")
            mol = Chem.MolFromMolFile(mol_file, removeHs=False, sanitize=True)

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
        """Read molecules from MOL2 file with error handling."""
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
        """Write molecules to SDF file with error handling."""
        try:
            writer = Chem.SDWriter(output_file)
            written = 0

            for i, mol in enumerate(self.molecules):
                try:
                    # Add properties back
                    if i < len(self.properties):
                        for key, value in self.properties[i].items():
                            mol.SetProp(key, value)
                    writer.write(mol)
                    written += 1
                except Exception as e:
                    print(f"Warning: Could not write molecule {i + 1}: {e}")
                    continue

            writer.close()
            print(f"Successfully wrote {written} molecules to {output_file}")
            return True
        except Exception as e:
            print(f"Error writing SDF file: {e}")
            return False

    def write_mol(self, output_file):
        """Write molecules to individual MOL files with error handling."""
        try:
            base_name = Path(output_file).stem
            base_dir = Path(output_file).parent
            written = 0

            for i, mol in enumerate(self.molecules):
                try:
                    mol_file = base_dir / f"{base_name}_{i+1}.mol"
                    Chem.MolToMolFile(mol, str(mol_file))
                    written += 1
                except Exception as e:
                    print(f"Warning: Could not write molecule {i + 1}: {e}")
                    continue

            print(f"Successfully wrote {written} molecules to individual MOL files")
            return True
        except Exception as e:
            print(f"Error writing MOL files: {e}")
            return False

    def write_mol2(self, output_file):
        """Write molecules to MOL2 file with error handling."""
        try:
            written = 0
            with open(output_file, 'w') as f:
                for i, mol in enumerate(self.molecules):
                    try:
                        mol2_block = Chem.MolToMol2Block(mol)
                        f.write(mol2_block)
                        f.write('\n')
                        written += 1
                    except Exception as e:
                        print(f"Warning: Could not write molecule {i + 1}: {e}")
                        continue

            print(f"Successfully wrote {written} molecules to {output_file}")
            return True
        except Exception as e:
            print(f"Error writing MOL2 file: {e}")
            return False

    def write_pdb(self, output_file):
        """Write molecules to PDB file with error handling."""
        try:
            written = 0
            with open(output_file, 'w') as f:
                for i, mol in enumerate(self.molecules):
                    try:
                        # Add 3D coordinates if not present
                        if mol.GetNumConformers() == 0:
                            try:
                                AllChem.EmbedMolecule(mol, randomSeed=42)
                                AllChem.MMFFOptimizeMolecule(mol)
                            except Exception as e:
                                print(f"Warning: Could not generate 3D coordinates for molecule {i + 1}: {e}")

                        pdb_block = Chem.MolToPDBBlock(mol)
                        f.write(pdb_block)
                        f.write('\n')
                        written += 1
                    except Exception as e:
                        print(f"Warning: Could not write molecule {i + 1}: {e}")
                        continue

            print(f"Successfully wrote {written} molecules to {output_file}")
            return True
        except Exception as e:
            print(f"Error writing PDB file: {e}")
            return False

    def write_xyz(self, output_file):
        """Write molecules to XYZ file with error handling."""
        try:
            written = 0
            with open(output_file, 'w') as f:
                for i, mol in enumerate(self.molecules):
                    try:
                        # Add 3D coordinates if not present
                        if mol.GetNumConformers() == 0:
                            try:
                                AllChem.EmbedMolecule(mol, randomSeed=42)
                                AllChem.MMFFOptimizeMolecule(mol)
                            except Exception as e:
                                print(f"Warning: Could not generate 3D coordinates for molecule {i + 1}: {e}")

                        xyz_block = Chem.MolToXYZBlock(mol)
                        f.write(xyz_block)
                        f.write('\n')
                        written += 1
                    except Exception as e:
                        print(f"Warning: Could not write molecule {i + 1}: {e}")
                        continue

            print(f"Successfully wrote {written} molecules to {output_file}")
            return True
        except Exception as e:
            print(f"Error writing XYZ file: {e}")
            return False

    def write_csv(self, output_file, smiles_column="SMILES"):
        """Write molecules to CSV file with SMILES and error handling."""
        try:
            data = []
            for i, mol in enumerate(self.molecules):
                try:
                    row = {smiles_column: Chem.MolToSmiles(mol)}

                    # Add properties
                    if i < len(self.properties):
                        row.update(self.properties[i])

                    data.append(row)
                except Exception as e:
                    print(f"Warning: Could not process molecule {i + 1}: {e}")
                    continue

            df = pd.DataFrame(data)
            df.to_csv(output_file, index=False, encoding='utf-8')
            print(f"Successfully wrote {len(data)} molecules to {output_file}")
            return True
        except Exception as e:
            print(f"Error writing CSV file: {e}")
            return False

    def convert(self, input_file, output_file, input_format=None, output_format=None, csv_options=None):
        """Main conversion method with enhanced error handling."""
        try:
            if input_format is None:
                ext = Path(input_file).suffix.lower().lstrip('.')
                if ext in ['csv', 'tsv', 'xls', 'xlsx']:
                    input_format = ext
                else:
                    input_format = Path(input_file).suffix.lower().lstrip('.')
            if output_format is None:
                output_format = Path(output_file).suffix.lower().lstrip('.')
            if input_format not in self.SUPPORTED_INPUT_FORMATS:
                print(f"Error: Unsupported input format '{input_format}'")
                print(f"Supported formats: {self.SUPPORTED_INPUT_FORMATS}")
                return False
            if output_format not in self.SUPPORTED_OUTPUT_FORMATS:
                print(f"Error: Unsupported output format '{output_format}'")
                print(f"Supported formats: {self.SUPPORTED_OUTPUT_FORMATS}")
                return False
            if not os.path.exists(input_file):
                print(f"Error: Input file '{input_file}' does not exist")
                return False
            if input_format in ['csv', 'tsv']:
                if csv_options is None:
                    csv_options = {}
                # For TSV, set default separator to tab
                if input_format == 'tsv' and 'separator' not in csv_options:
                    csv_options['separator'] = '\t'
                success = self.read_csv(input_file, **csv_options)
            elif input_format in ['xls', 'xlsx']:
                if csv_options is None:
                    csv_options = {}
                success = self.read_excel(input_file, **csv_options)
            elif input_format == 'sdf':
                success = self.read_sdf(input_file)
            elif input_format == 'mol':
                success = self.read_mol(input_file)
            elif input_format == 'mol2':
                success = self.read_mol2(input_file)
            if not success:
                return False
            output_dir = Path(output_file).parent
            if output_dir and not output_dir.exists():
                output_dir.mkdir(parents=True, exist_ok=True)
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
        except Exception as e:
            print(f"Unexpected error during conversion: {e}")
            return False

def main():
    parser = argparse.ArgumentParser(
        description='Convert between molecular file formats (CSV, TSV, Excel, SDF, MOL, MOL2, etc.)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Convert CSV to SDF (user must specify --separator and --quotechar if not standard)
  python molecular_converter.py input.csv output.sdf --separator ';' --quotechar '"'

  # Convert Excel to SDF
  python molecular_converter.py input.xlsx output.sdf

  # Convert SDF to CSV
  python molecular_converter.py input.sdf output.csv

  # Let the converter auto-detect the SMILES column
  python molecular_converter.py input.csv output.sdf

  # Specify SMILES column explicitly
  python molecular_converter.py input.csv output.sdf --smiles-col SMILES
        """
    )
    parser.add_argument('input_file', help='Input file path')
    parser.add_argument('output_file', help='Output file path')
    parser.add_argument('--input-format', choices=['csv', 'tsv', 'xls', 'xlsx', 'sdf', 'mol', 'mol2'], help='Input file format (auto-detected from extension if not specified)')
    parser.add_argument('--output-format', choices=['sdf', 'mol', 'mol2', 'pdb', 'xyz', 'csv'], help='Output file format (auto-detected from extension if not specified)')
    parser.add_argument('--smiles-col', help='SMILES column name (if not specified, will try to auto-detect)')
    parser.add_argument('--id-col', help='ID column name')
    parser.add_argument('--name-col', help='Name column name')
    parser.add_argument('--separator', default=',', help='CSV/TSV separator (default: ,)')
    parser.add_argument('--quotechar', default='"', help='CSV/TSV quote character (default: ")')
    parser.add_argument('--encoding', help='File encoding (auto-detected if not specified)')
    args = parser.parse_args()
    csv_options = {
        'smiles_column': args.smiles_col,
        'id_column': args.id_col,
        'name_column': args.name_col,
        'separator': args.separator,
        'quotechar': args.quotechar,
        'encoding': args.encoding
    }
    csv_options = {k: v for k, v in csv_options.items() if v is not None}
    converter = MolecularConverter()
    success = converter.convert(
        args.input_file,
        args.output_file,
        args.input_format,
        args.output_format,
        csv_options=csv_options
    )
    if success:
        print(f"\n✅ Successfully converted {args.input_file} to {args.output_file}")
    else:
        print(f"\n❌ Failed to convert {args.input_file}")
        sys.exit(1)

if __name__ == "__main__":
    main()