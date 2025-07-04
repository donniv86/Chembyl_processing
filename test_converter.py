#!/usr/bin/env python3
"""
Test script for the enhanced molecular converter
Demonstrates error handling and encoding detection features
"""

import os
import sys
from molecular_converter import MolecularConverter

def test_encoding_detection():
    """Test encoding detection functionality"""
    print("=== Testing Encoding Detection ===")

    converter = MolecularConverter()

    # Test with jiang.csv
    if os.path.exists('jiang.csv'):
        print("\nTesting encoding detection on jiang.csv:")
        encoding = converter.detect_encoding('jiang.csv')
        print(f"Detected encoding: {encoding}")
    else:
        print("jiang.csv not found, skipping encoding detection test")

def test_csv_conversion():
    """Test CSV to SDF conversion with error handling"""
    print("\n=== Testing CSV to SDF Conversion ===")

    converter = MolecularConverter()

    if os.path.exists('jiang.csv'):
        print("Converting jiang.csv to jiang_enhanced.sdf...")
        success = converter.convert(
            'jiang.csv',
            'jiang_enhanced.sdf',
            input_format='csv',
            output_format='sdf',
            csv_options={
                'smiles_column': 'SMILES',
                'name_column': 'Name'
            }
        )

        if success:
            print("✅ CSV to SDF conversion successful!")
        else:
            print("❌ CSV to SDF conversion failed!")
    else:
        print("jiang.csv not found, skipping conversion test")

def test_error_handling():
    """Test error handling with invalid files"""
    print("\n=== Testing Error Handling ===")

    converter = MolecularConverter()

    # Test with non-existent file
    print("Testing with non-existent file:")
    success = converter.convert('nonexistent.csv', 'output.sdf')
    if not success:
        print("✅ Correctly handled non-existent file")

    # Test with invalid SMILES
    print("\nTesting with invalid SMILES:")
    import pandas as pd

    # Create a test CSV with invalid SMILES
    test_data = {
        'SMILES': ['CCO', 'invalid_smiles', 'CC(C)C', 'another_invalid'],
        'Name': ['Ethanol', 'Invalid1', 'Isobutane', 'Invalid2']
    }
    test_df = pd.DataFrame(test_data)
    test_df.to_csv('test_invalid.csv', index=False)

    success = converter.convert('test_invalid.csv', 'test_output.sdf')
    if success:
        print("✅ Successfully handled invalid SMILES (converted valid ones)")

    # Clean up test file
    if os.path.exists('test_invalid.csv'):
        os.remove('test_invalid.csv')

def test_multiple_formats():
    """Test conversion to multiple output formats"""
    print("\n=== Testing Multiple Output Formats ===")

    converter = MolecularConverter()

    if os.path.exists('jiang.csv'):
        formats = ['sdf', 'mol2', 'csv']

        for fmt in formats:
            output_file = f'jiang_test.{fmt}'
            print(f"Converting to {fmt.upper()}...")

            success = converter.convert(
                'jiang.csv',
                output_file,
                input_format='csv',
                output_format=fmt
            )

            if success:
                print(f"✅ Successfully converted to {fmt.upper()}")
            else:
                print(f"❌ Failed to convert to {fmt.upper()}")

def main():
    """Run all tests"""
    print("🧪 Testing Enhanced Molecular Converter")
    print("=" * 50)

    try:
        test_encoding_detection()
        test_csv_conversion()
        test_error_handling()
        test_multiple_formats()

        print("\n" + "=" * 50)
        print("✅ All tests completed!")

    except Exception as e:
        print(f"\n❌ Test failed with error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()