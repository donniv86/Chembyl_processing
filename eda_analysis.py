#!/usr/bin/env python3
"""
Exploratory Data Analysis (EDA) for ChemBL Molecular Data
Analyzes molecular properties across different targets and provides insights.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
import warnings
warnings.filterwarnings('ignore')

class ChemBLEDA:
    """Exploratory Data Analysis for ChemBL molecular data."""

    def __init__(self, csv_file, separator=";", smiles_column="Smiles"):
        """Initialize EDA with CSV file."""
        self.csv_file = csv_file
        self.separator = separator
        self.smiles_column = smiles_column
        self.df = None
        self.molecular_properties = None

    def load_data(self):
        """Load and preprocess the CSV data."""
        print("Loading data...")
        self.df = pd.read_csv(self.csv_file, sep=self.separator, quotechar='"')
        print(f"Loaded {len(self.df)} rows and {len(self.df.columns)} columns")

        # Remove rows with empty SMILES
        self.df = self.df.dropna(subset=[self.smiles_column])
        self.df = self.df[self.df[self.smiles_column].str.strip() != '']
        print(f"After filtering empty SMILES: {len(self.df)} rows")

        return self.df

    def calculate_molecular_properties(self):
        """Calculate molecular properties from SMILES."""
        print("Calculating molecular properties...")

        properties = []
        for idx, row in self.df.iterrows():
            try:
                smiles = str(row[self.smiles_column]).strip()
                mol = Chem.MolFromSmiles(smiles)

                if mol is not None:
                    prop = {
                        'MolWt': Descriptors.MolWt(mol),
                        'LogP': Descriptors.MolLogP(mol),
                        'NumHDonors': Descriptors.NumHDonors(mol),
                        'NumHAcceptors': Descriptors.NumHAcceptors(mol),
                        'NumRotatableBonds': Descriptors.NumRotatableBonds(mol),
                        'TPSA': Descriptors.TPSA(mol),
                        'NumAtoms': mol.GetNumAtoms(),
                        'NumHeavyAtoms': mol.GetNumHeavyAtoms(),
                        'NumRings': Descriptors.RingCount(mol),
                        'AromaticRings': Descriptors.NumAromaticRings(mol),
                        'SaturatedRings': Descriptors.NumSaturatedRings(mol),
                        'SMILES': smiles
                    }

                    # Lipinski Rule of 5 violations
                    violations = 0
                    if prop['MolWt'] > 500: violations += 1
                    if prop['LogP'] > 5: violations += 1
                    if prop['NumHDonors'] > 5: violations += 1
                    if prop['NumHAcceptors'] > 10: violations += 1
                    prop['LipinskiViolations'] = violations

                    properties.append(prop)
                else:
                    print(f"Warning: Could not parse SMILES at row {idx + 1}")

            except Exception as e:
                print(f"Error processing row {idx + 1}: {e}")
                continue

        self.molecular_properties = pd.DataFrame(properties)
        print(f"Calculated properties for {len(self.molecular_properties)} molecules")

        return self.molecular_properties

    def merge_data(self):
        """Merge molecular properties with original data."""
        if self.molecular_properties is not None and len(self.molecular_properties) > 0:
            # Merge on SMILES
            self.df = self.df.merge(self.molecular_properties,
                                  left_on=self.smiles_column,
                                  right_on='SMILES',
                                  how='left')
            print(f"Merged data shape: {self.df.shape}")
        else:
            print("No molecular properties to merge")
        return self.df

    def analyze_targets(self):
        """Analyze data by different targets."""
        print("\n=== TARGET ANALYSIS ===")

        # Target distribution
        if 'Target Name' in self.df.columns:
            target_counts = self.df['Target Name'].value_counts()
            print(f"\nNumber of unique targets: {len(target_counts)}")
            print("\nTop 10 targets by number of compounds:")
            print(target_counts.head(10))

            # Target organisms
            if 'Target Organism' in self.df.columns:
                organism_counts = self.df['Target Organism'].value_counts()
                print(f"\nNumber of unique organisms: {len(organism_counts)}")
                print("\nOrganism distribution:")
                print(organism_counts)

        return target_counts if 'Target Name' in self.df.columns else None

    def analyze_activity(self):
        """Analyze activity data."""
        print("\n=== ACTIVITY ANALYSIS ===")

        # Activity status
        if 'Comment' in self.df.columns:
            activity_counts = self.df['Comment'].value_counts()
            print("\nActivity distribution:")
            print(activity_counts)

            # Activity by target
            if 'Target Name' in self.df.columns:
                activity_by_target = self.df.groupby(['Target Name', 'Comment']).size().unstack(fill_value=0)
                print("\nActivity by target (first 10 targets):")
                print(activity_by_target.head(10))

        return activity_counts if 'Comment' in self.df.columns else None

    def analyze_molecular_properties(self):
        """Analyze molecular properties."""
        print("\n=== MOLECULAR PROPERTIES ANALYSIS ===")

        if self.molecular_properties is not None:
            # Basic statistics
            print("\nMolecular properties statistics:")
            print(self.molecular_properties.describe())

            # Lipinski Rule of 5 analysis
            lipinski_violations = self.molecular_properties['LipinskiViolations'].value_counts().sort_index()
            print(f"\nLipinski Rule of 5 violations:")
            print(lipinski_violations)

            # Property ranges
            print(f"\nMolecular weight range: {self.molecular_properties['MolWt'].min():.2f} - {self.molecular_properties['MolWt'].max():.2f}")
            print(f"LogP range: {self.molecular_properties['LogP'].min():.2f} - {self.molecular_properties['LogP'].max():.2f}")
            print(f"TPSA range: {self.molecular_properties['TPSA'].min():.2f} - {self.molecular_properties['TPSA'].max():.2f}")

        return self.molecular_properties

    def analyze_by_target(self, target_name=None):
        """Analyze properties for a specific target or all targets."""
        print(f"\n=== TARGET-SPECIFIC ANALYSIS ===")

        if target_name:
            target_data = self.df[self.df['Target Name'] == target_name]
            print(f"Analyzing target: {target_name}")
            print(f"Number of compounds: {len(target_data)}")
        else:
            target_data = self.df
            print("Analyzing all targets")

        if len(target_data) > 0 and self.molecular_properties is not None:
            # Get molecular properties for target compounds
            target_smiles = target_data[self.smiles_column].tolist()
            target_props = self.molecular_properties[self.molecular_properties['SMILES'].isin(target_smiles)]

            if len(target_props) > 0:
                print(f"\nMolecular properties for target compounds:")
                print(target_props.describe())

                # Activity analysis for target
                if 'Comment' in target_data.columns:
                    target_activity = target_data['Comment'].value_counts()
                    print(f"\nActivity distribution for target:")
                    print(target_activity)

        return target_data

    def create_visualizations(self, output_dir="eda_plots"):
        """Create visualization plots."""
        import os
        os.makedirs(output_dir, exist_ok=True)

        print(f"\n=== CREATING VISUALIZATIONS ===")
        print(f"Plots will be saved to: {output_dir}")

        # Set style
        plt.style.use('default')
        sns.set_palette("husl")

        if self.molecular_properties is not None:
            # 1. Molecular weight distribution
            plt.figure(figsize=(10, 6))
            plt.hist(self.molecular_properties['MolWt'], bins=30, alpha=0.7, edgecolor='black')
            plt.xlabel('Molecular Weight')
            plt.ylabel('Frequency')
            plt.title('Distribution of Molecular Weights')
            plt.axvline(x=500, color='red', linestyle='--', label='Lipinski MW limit')
            plt.legend()
            plt.tight_layout()
            plt.savefig(f'{output_dir}/molecular_weight_distribution.png', dpi=300, bbox_inches='tight')
            plt.close()

            # 2. LogP distribution
            plt.figure(figsize=(10, 6))
            plt.hist(self.molecular_properties['LogP'], bins=30, alpha=0.7, edgecolor='black')
            plt.xlabel('LogP')
            plt.ylabel('Frequency')
            plt.title('Distribution of LogP Values')
            plt.axvline(x=5, color='red', linestyle='--', label='Lipinski LogP limit')
            plt.legend()
            plt.tight_layout()
            plt.savefig(f'{output_dir}/logp_distribution.png', dpi=300, bbox_inches='tight')
            plt.close()

            # 3. MW vs LogP scatter plot
            plt.figure(figsize=(10, 8))
            plt.scatter(self.molecular_properties['MolWt'], self.molecular_properties['LogP'],
                       alpha=0.6, s=50)
            plt.xlabel('Molecular Weight')
            plt.ylabel('LogP')
            plt.title('Molecular Weight vs LogP')
            plt.axhline(y=5, color='red', linestyle='--', alpha=0.7)
            plt.axvline(x=500, color='red', linestyle='--', alpha=0.7)
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(f'{output_dir}/mw_vs_logp_scatter.png', dpi=300, bbox_inches='tight')
            plt.close()

            # 4. Lipinski violations
            plt.figure(figsize=(10, 6))
            violations = self.molecular_properties['LipinskiViolations'].value_counts().sort_index()
            plt.bar(violations.index, violations.values, alpha=0.7, edgecolor='black')
            plt.xlabel('Number of Lipinski Violations')
            plt.ylabel('Number of Compounds')
            plt.title('Distribution of Lipinski Rule of 5 Violations')
            plt.xticks(range(len(violations)))
            plt.tight_layout()
            plt.savefig(f'{output_dir}/lipinski_violations.png', dpi=300, bbox_inches='tight')
            plt.close()

            # 5. TPSA distribution
            plt.figure(figsize=(10, 6))
            plt.hist(self.molecular_properties['TPSA'], bins=30, alpha=0.7, edgecolor='black')
            plt.xlabel('Topological Polar Surface Area (TPSA)')
            plt.ylabel('Frequency')
            plt.title('Distribution of TPSA Values')
            plt.tight_layout()
            plt.savefig(f'{output_dir}/tpsa_distribution.png', dpi=300, bbox_inches='tight')
            plt.close()

        # 6. Target distribution (if available)
        if 'Target Name' in self.df.columns:
            plt.figure(figsize=(12, 8))
            target_counts = self.df['Target Name'].value_counts().head(15)
            plt.barh(range(len(target_counts)), target_counts.values, alpha=0.7, edgecolor='black')
            plt.yticks(range(len(target_counts)), target_counts.index)
            plt.xlabel('Number of Compounds')
            plt.title('Top 15 Targets by Number of Compounds')
            plt.gca().invert_yaxis()
            plt.tight_layout()
            plt.savefig(f'{output_dir}/target_distribution.png', dpi=300, bbox_inches='tight')
            plt.close()

            # 7. Activity by target heatmap
            if 'Comment' in self.df.columns:
                activity_by_target = self.df.groupby(['Target Name', 'Comment']).size().unstack(fill_value=0)
                top_targets = activity_by_target.sum(axis=1).sort_values(ascending=False).head(10).index
                activity_subset = activity_by_target.loc[top_targets]

                plt.figure(figsize=(12, 8))
                sns.heatmap(activity_subset, annot=True, fmt='d', cmap='YlOrRd', cbar_kws={'label': 'Number of Compounds'})
                plt.title('Activity Distribution by Target (Top 10 Targets)')
                plt.xlabel('Activity Status')
                plt.ylabel('Target Name')
                plt.tight_layout()
                plt.savefig(f'{output_dir}/activity_by_target_heatmap.png', dpi=300, bbox_inches='tight')
                plt.close()

        print(f"Visualizations saved to {output_dir}/")

    def generate_report(self, output_file="eda_report.txt"):
        """Generate a comprehensive EDA report."""
        print(f"\n=== GENERATING EDA REPORT ===")

        with open(output_file, 'w') as f:
            f.write("CHEMBL EXPLORATORY DATA ANALYSIS REPORT\n")
            f.write("=" * 50 + "\n\n")

            # Basic data info
            f.write("1. DATA OVERVIEW\n")
            f.write("-" * 20 + "\n")
            f.write(f"Total compounds: {len(self.df)}\n")
            f.write(f"Total columns: {len(self.df.columns)}\n")
            f.write(f"Columns: {', '.join(self.df.columns)}\n\n")

            # Target analysis
            if 'Target Name' in self.df.columns:
                f.write("2. TARGET ANALYSIS\n")
                f.write("-" * 20 + "\n")
                target_counts = self.df['Target Name'].value_counts()
                f.write(f"Number of unique targets: {len(target_counts)}\n")
                f.write("Top 10 targets:\n")
                for target, count in target_counts.head(10).items():
                    f.write(f"  {target}: {count} compounds\n")
                f.write("\n")

                if 'Target Organism' in self.df.columns:
                    organism_counts = self.df['Target Organism'].value_counts()
                    f.write(f"Number of unique organisms: {len(organism_counts)}\n")
                    f.write("Organism distribution:\n")
                    for org, count in organism_counts.items():
                        f.write(f"  {org}: {count} compounds\n")
                    f.write("\n")

            # Activity analysis
            if 'Comment' in self.df.columns:
                f.write("3. ACTIVITY ANALYSIS\n")
                f.write("-" * 20 + "\n")
                activity_counts = self.df['Comment'].value_counts()
                f.write("Activity distribution:\n")
                for activity, count in activity_counts.items():
                    f.write(f"  {activity}: {count} compounds\n")
                f.write("\n")

            # Molecular properties
            if self.molecular_properties is not None:
                f.write("4. MOLECULAR PROPERTIES\n")
                f.write("-" * 20 + "\n")
                f.write("Property statistics:\n")
                f.write(self.molecular_properties.describe().to_string())
                f.write("\n\n")

                f.write("Lipinski Rule of 5 violations:\n")
                violations = self.molecular_properties['LipinskiViolations'].value_counts().sort_index()
                for v, count in violations.items():
                    f.write(f"  {v} violations: {count} compounds\n")
                f.write("\n")

                # Drug-likeness analysis
                drug_like = len(self.molecular_properties[self.molecular_properties['LipinskiViolations'] <= 1])
                total = len(self.molecular_properties)
                f.write(f"Drug-likeness (≤1 Lipinski violation): {drug_like}/{total} ({drug_like/total*100:.1f}%)\n\n")

            # Target-specific analysis
            if 'Target Name' in self.df.columns and self.molecular_properties is not None:
                f.write("5. TARGET-SPECIFIC ANALYSIS\n")
                f.write("-" * 20 + "\n")

                for target in self.df['Target Name'].value_counts().head(5).index:
                    target_data = self.df[self.df['Target Name'] == target]
                    target_smiles = target_data[self.smiles_column].tolist()
                    target_props = self.molecular_properties[self.molecular_properties['SMILES'].isin(target_smiles)]

                    f.write(f"Target: {target}\n")
                    f.write(f"  Number of compounds: {len(target_data)}\n")

                    if len(target_props) > 0:
                        f.write(f"  Average MW: {target_props['MolWt'].mean():.2f}\n")
                        f.write(f"  Average LogP: {target_props['LogP'].mean():.2f}\n")
                        f.write(f"  Average TPSA: {target_props['TPSA'].mean():.2f}\n")

                        if 'Comment' in target_data.columns:
                            activity = target_data['Comment'].value_counts()
                            f.write(f"  Activity: {dict(activity)}\n")
                    f.write("\n")

        print(f"Report saved to {output_file}")

    def run_full_analysis(self, output_dir="eda_results"):
        """Run complete EDA analysis."""
        print("=== CHEMBL EXPLORATORY DATA ANALYSIS ===")

        # Load data
        self.load_data()

        # Calculate molecular properties
        self.calculate_molecular_properties()

        # Merge data
        self.merge_data()

        # Run analyses
        self.analyze_targets()
        self.analyze_activity()
        self.analyze_molecular_properties()

        # Create visualizations
        self.create_visualizations(output_dir)

        # Generate report
        self.generate_report(f"{output_dir}/eda_report.txt")

        print(f"\n=== ANALYSIS COMPLETE ===")
        print(f"Results saved to: {output_dir}/")

        return self.df

def main():
    """Main function to run EDA analysis."""
    import argparse

    parser = argparse.ArgumentParser(description='Exploratory Data Analysis for ChemBL data')
    parser.add_argument('input_file', help='Input CSV file path')
    parser.add_argument('--separator', default=';', help='CSV separator (default: ;)')
    parser.add_argument('--smiles-col', default='Smiles', help='SMILES column name (default: Smiles)')
    parser.add_argument('--output-dir', default='eda_results', help='Output directory (default: eda_results)')

    args = parser.parse_args()

    # Create EDA object and run analysis
    eda = ChemBLEDA(args.input_file, args.separator, args.smiles_col)
    eda.run_full_analysis(args.output_dir)

if __name__ == "__main__":
    main()