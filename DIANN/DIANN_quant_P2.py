#!/usr/bin/env python3
"""
Quick DIA-NN MaxLFQ Extraction
Simple script to extract MaxLFQ values from filtered DIA-NN output
"""

import pandas as pd
import numpy as np
import sys
from pathlib import Path

def extract_maxlfq(input_file):
    """Extract MaxLFQ values from DIA-NN parquet file"""
    
    print("Quick DIA-NN MaxLFQ Extraction")
    print("-" * 50)
    
    # Load data
    print(f"Loading: {input_file}")
    df = pd.read_parquet(input_file)
    
    print(f"Loaded {len(df):,} rows")
    print(f"Proteins: {df['Protein.Group'].nunique():,}")
    print(f"Samples: {df['Run'].nunique()}")
    
    # Get sample names
    sample_mapping = {}
    for run in df['Run'].unique():
        sample_name = Path(run).stem
        # Clean up common suffixes
        for suffix in ['_DIA', '_dia', '.raw', '.mzML', '.d']:
            sample_name = sample_name.replace(suffix, '')
        sample_mapping[run] = sample_name
    
    print("\nSample names:")
    for orig, clean in list(sample_mapping.items())[:5]:
        print(f"  {clean}")
    if len(sample_mapping) > 5:
        print(f"  ... and {len(sample_mapping) - 5} more")
    
    # Extract MaxLFQ values
    print("\nExtracting PG.MaxLFQ values...")
    
    # Check if PG.MaxLFQ exists
    if 'PG.MaxLFQ' not in df.columns:
        print("ERROR: PG.MaxLFQ column not found!")
        print("Available columns:", df.columns.tolist())
        return
    
    # Pivot to get protein x sample matrix
    result = df.pivot_table(
        index='Protein.Group',
        columns='Run', 
        values='PG.MaxLFQ',
        aggfunc='first'  # Take first value if duplicates exist
    )
    
    # Rename columns to clean sample names
    result.columns = [sample_mapping[col] for col in result.columns]
    
    # Replace zeros with NaN
    result = result.replace(0, np.nan)
    
    # Get protein annotations
    print("Adding protein annotations...")
    annotations = df[['Protein.Group', 'Protein.Names', 'Genes']].drop_duplicates('Protein.Group')
    
    # Merge annotations with quantification
    final = annotations.merge(result, left_on='Protein.Group', right_index=True, how='right')
    final = final.set_index('Protein.Group')
    
    # Move annotation columns to the front
    info_cols = ['Protein.Names', 'Genes']
    sample_cols = [col for col in final.columns if col not in info_cols]
    final = final[info_cols + sample_cols]
    
    # Calculate statistics
    numeric_cols = final.select_dtypes(include=[np.number]).columns
    final['N_Samples'] = final[numeric_cols].notna().sum(axis=1)
    final['Mean_Intensity'] = final[numeric_cols].mean(axis=1)
    final['CV_%'] = (final[numeric_cols].std(axis=1) / final[numeric_cols].mean(axis=1) * 100)
    
    # Save output
    output_file = input_file.replace('.parquet', '_MaxLFQ.csv')
    print(f"\nSaving to: {output_file}")
    final.to_csv(output_file)
    
    # Summary statistics
    print("\n" + "=" * 50)
    print("SUMMARY")
    print("=" * 50)
    print(f"Quantified proteins: {len(final):,}")
    print(f"Samples: {len(sample_cols)}")
    print(f"Missing values: {final[numeric_cols].isna().sum().sum():,} ({final[numeric_cols].isna().sum().sum() / (len(final) * len(numeric_cols)) * 100:.1f}%)")
    print(f"Median CV: {final['CV_%'].median():.1f}%")
    print(f"Proteins in all samples: {(final['N_Samples'] == len(sample_cols)).sum():,}")
    
    print("\nDone!")
    
def main():
    if len(sys.argv) < 2:
        print("Usage: python quick_maxlfq.py <filtered_parquet_file>")
        print("\nExample:")
        print("  python quick_maxlfq.py DR12_merged_experiment_report_filtered.parquet")
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    if not input_file.endswith('.parquet'):
        print(f"Warning: Input file '{input_file}' doesn't have .parquet extension")
    
    try:
        extract_maxlfq(input_file)
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()