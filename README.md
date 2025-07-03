# Molecular Format Converter

A general-purpose tool for converting between various molecular file formats including CSV with SMILES, SDF, MOL, MOL2, PDB, and XYZ formats.

## Features

- **Multiple Input Formats**: CSV (with SMILES), SDF, MOL, MOL2
- **Multiple Output Formats**: SDF, MOL, MOL2, PDB, XYZ, CSV
- **Auto-format Detection**: Automatically detects file formats from extensions
- **Property Preservation**: Maintains molecular properties across conversions
- **2D/3D Coordinates**: Generates 2D coordinates for 2D formats, 3D for 3D formats
- **Flexible CSV Handling**: Customizable column names and separators
- **Error Handling**: Robust error handling with detailed reporting
- **Progress Reporting**: Shows conversion progress and statistics

## Installation

1. Install the required dependencies:
```bash
pip install -r requirements.txt
```

## Usage

### Basic Usage
```bash
# Convert CSV to SDF
python molecular_converter.py input.csv output.sdf

# Convert SDF to MOL2
python molecular_converter.py input.sdf output.mol2

# Convert CSV to individual MOL files
python molecular_converter.py input.csv output.mol

# Convert SDF to CSV with SMILES
python molecular_converter.py input.sdf output.csv
```

### Advanced Usage
```bash
# Specify formats explicitly
python molecular_converter.py input.csv output.mol2 --input-format csv --output-format mol2

# Custom CSV column names
python molecular_converter.py input.csv output.sdf --smiles-col "SMILES" --id-col "ID" --name-col "Name"

# Custom CSV separator (for semicolon-separated files)
python molecular_converter.py input.csv output.sdf --separator ";"
```

## Command Line Options

### Required Arguments
- `input_file`: Input file path
- `output_file`: Output file path

### Optional Arguments
- `--input-format`: Input file format (csv, sdf, mol, mol2) - auto-detected from extension
- `--output-format`: Output file format (sdf, mol, mol2, pdb, xyz, csv) - auto-detected from extension

### CSV-Specific Options
- `--smiles-col`: SMILES column name (default: "SMILES")
- `--id-col`: ID column name
- `--name-col`: Name column name
- `--separator`: CSV separator (default: ",")
- `--quotechar`: CSV quote character (default: '"')

## Supported Formats

### Input Formats
- **CSV**: Files containing SMILES strings with molecular properties
- **SDF**: Structure Data Files (multi-molecule)
- **MOL**: MDL Molfile format (single molecule)
- **MOL2**: Tripos MOL2 format

### Output Formats
- **SDF**: Structure Data Files with properties
- **MOL**: Individual MDL Molfiles (creates multiple files)
- **MOL2**: Tripos MOL2 format
- **PDB**: Protein Data Bank format (with 3D coordinates)
- **XYZ**: XYZ coordinate format (with 3D coordinates)
- **CSV**: CSV with SMILES and properties

## Examples

### Converting ChemBL CSV to SDF
```bash
# For jiang.csv with semicolon separator
python molecular_converter.py jiang.csv jiang.sdf --separator ";"
```

### Converting between different formats
```bash
# SDF to MOL2
python molecular_converter.py molecules.sdf molecules.mol2

# CSV to individual MOL files
python molecular_converter.py compounds.csv compounds.mol

# SDF to PDB (with 3D coordinates)
python molecular_converter.py molecules.sdf molecules.pdb

# Extract SMILES from SDF
python molecular_converter.py molecules.sdf molecules.csv
```

## Requirements

- Python 3.6+
- pandas
- RDKit (rdkit-pypi)

## Troubleshooting

- Make sure the CSV file uses semicolon (;) as separator
- Check that the SMILES column contains valid SMILES strings
- Ensure RDKit is properly installed
- For large files, the conversion may take some time