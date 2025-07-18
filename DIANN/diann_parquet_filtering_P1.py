#!/usr/bin/env python3
"""
DIA-NN Parquet Filtering Script with Peak Width and Intensity Filters
Filters DIA-NN output based on quality criteria including peak width and intensity
"""

import pandas as pd
import numpy as np
import sys
import os
from datetime import datetime
from pathlib import Path

def filter_diann_parquet(input_file, output_file=None, report_file=None, 
                        min_peak_width=None, min_precursor_quantity=None,
                        min_pg_maxlfq=None, min_pg_quantity=None,
                        metadata_file=None, require_all_conditions=False):
    """
    Filter DIA-NN parquet file based on standard quality criteria
    
    Parameters:
    -----------
    input_file : str
        Path to input parquet file
    output_file : str, optional
        Path to output filtered parquet file
    report_file : str, optional
        Path to filtering report text file
    min_peak_width : float, optional
        Minimum peak width (RT.Stop - RT.Start) in minutes. Default is 0.05 (3 seconds)
    min_precursor_quantity : float, optional
        Minimum Precursor.Quantity intensity threshold (precursor-level)
    min_pg_maxlfq : float, optional
        Minimum PG.MaxLFQ intensity threshold (protein-level, normalized)
    min_pg_quantity : float, optional
        Minimum PG.Quantity intensity threshold (protein-level, raw)
    metadata_file : str, optional
        Path to metadata CSV file for condition-based filtering
    require_all_conditions : bool, optional
        If True, keep only proteins present in at least 1 replicate of ALL conditions
    """
    
    # Set default min peak width if not specified
    if min_peak_width is None:
        min_peak_width = 0.05  # 3 seconds default
    
    # Generate default output filenames if not provided
    input_path = Path(input_file)
    input_dir = input_path.parent  # Get the directory of the input file
    
    if output_file is None:
        # Save in the same directory as input file
        output_file = input_dir / (input_path.stem + "_filtered.parquet")
    else:
        # If output file is just a filename (no path), save in input directory
        output_path = Path(output_file)
        if not output_path.is_absolute() and output_path.parent == Path('.'):
            output_file = input_dir / output_file
            
    if report_file is None:
        # Save in the same directory as input file
        report_file = input_dir / (input_path.stem + "_filter_report.txt")
    else:
        # If report file is just a filename (no path), save in input directory
        report_path = Path(report_file)
        if not report_path.is_absolute() and report_path.parent == Path('.'):
            report_file = input_dir / report_file
    
    # Print header
    print("=" * 60)
    print("DIA-NN FILTERING (DIA-SPECIFIC)")
    print("=" * 60)
    print(f"Input: {input_file}")
    print(f"Output: {output_file}")
    print(f"Report: {report_file}")
    print(f"Min Peak Width: {min_peak_width} minutes ({min_peak_width * 60} seconds)")
    if min_precursor_quantity:
        print(f"Min Precursor Quantity: {min_precursor_quantity}")
    if min_pg_maxlfq:
        print(f"Min PG MaxLFQ: {min_pg_maxlfq}")
    if min_pg_quantity:
        print(f"Min PG Quantity: {min_pg_quantity}")
    print()
    
    # Load data
    print("Loading parquet file...")
    df = pd.read_parquet(input_file)
    initial_rows = len(df)
    print(f"Loaded {initial_rows:,} rows")
    
    # Check for contaminants
    if 'Genes' in df.columns:
        contaminant_count = df['Genes'].fillna('').str.startswith('##').sum()
        print(f"Found {contaminant_count:,} contaminant entries (genes starting with ##)")
    
    # Calculate peak width
    if 'RT.Start' in df.columns and 'RT.Stop' in df.columns:
        df['Peak.Width'] = df['RT.Stop'] - df['RT.Start']
        print(f"\nPeak width calculated (RT.Stop - RT.Start)")
        print(f"  Mean: {df['Peak.Width'].mean():.3f} min")
        print(f"  Median: {df['Peak.Width'].median():.3f} min")
        print(f"  Min: {df['Peak.Width'].min():.3f} min")
        print(f"  Max: {df['Peak.Width'].max():.3f} min")
    else:
        print("WARNING: RT.Start or RT.Stop columns not found - skipping peak width filter")
    
    # Print DIA-relevant intensity statistics
    intensity_cols = ['Precursor.Quantity', 'PG.MaxLFQ', 'PG.Quantity', 'Fragment.Quant.Corrected']
    print("\nDIA Intensity Statistics:")
    for col in intensity_cols:
        if col in df.columns:
            print(f"\n{col}:")
            print(f"  Mean: {df[col].mean():.2e}")
            print(f"  Median: {df[col].median():.2e}")
            print(f"  Min: {df[col].min():.2e}")
            print(f"  Max: {df[col].max():.2e}")
            print(f"  10th percentile: {df[col].quantile(0.10):.2e}")
            print(f"  90th percentile: {df[col].quantile(0.90):.2e}")
    print()
    
    # Initialize report
    report = []
    report.append("DIA-NN Filtering Report (DIA-Specific)")
    report.append("=" * 60)
    report.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report.append(f"Input file: {input_file}")
    report.append(f"Initial rows: {initial_rows:,}")
    if 'Genes' in df.columns:
        contaminant_count = df['Genes'].fillna('').str.startswith('##').sum()
        report.append(f"Contaminant entries found: {contaminant_count:,}")
    report.append(f"Min peak width threshold: {min_peak_width} minutes ({min_peak_width * 60} seconds)")
    if min_precursor_quantity:
        report.append(f"Min Precursor Quantity: {min_precursor_quantity}")
    if min_pg_maxlfq:
        report.append(f"Min PG MaxLFQ: {min_pg_maxlfq}")
    if min_pg_quantity:
        report.append(f"Min PG Quantity: {min_pg_quantity}")
    report.append("")
    
    # Add peak width statistics to report
    if 'Peak.Width' in df.columns:
        report.append("Initial Peak Width Statistics:")
        report.append(f"  Mean: {df['Peak.Width'].mean():.3f} min")
        report.append(f"  Median: {df['Peak.Width'].median():.3f} min")
        report.append(f"  Std Dev: {df['Peak.Width'].std():.3f} min")
        report.append(f"  Min: {df['Peak.Width'].min():.3f} min")
        report.append(f"  Max: {df['Peak.Width'].max():.3f} min")
        report.append(f"  95th percentile: {df['Peak.Width'].quantile(0.95):.3f} min")
        report.append("")
    
    # Add DIA intensity statistics to report
    report.append("Initial DIA Intensity Statistics:")
    for col in intensity_cols:
        if col in df.columns:
            report.append(f"\n{col}:")
            report.append(f"  Mean: {df[col].mean():.2e}")
            report.append(f"  Median: {df[col].median():.2e}")
            report.append(f"  10th percentile: {df[col].quantile(0.10):.2e}")
            report.append(f"  90th percentile: {df[col].quantile(0.90):.2e}")
    report.append("")
    
    # ============================================================
    # FILTER DEFINITIONS - DIA-SPECIFIC
    # ============================================================
    filters = [
        # CONTAMINANT FILTER (remove proteins with genes starting with ##)
        ("Remove contaminants (##)", "Genes", lambda x: ~x.fillna('').str.startswith('##')),
        
        # Q-VALUE FILTERS (FDR control)
        ("Q.Value <= 0.01", "Q.Value", lambda x: x <= 0.01),
        ("Lib.Q.Value < 0.01", "Lib.Q.Value", lambda x: x < 0.01),
        ("Lib.PG.Q.Value < 0.01", "Lib.PG.Q.Value", lambda x: x < 0.01),
        ("Lib.Peptidoform.Q.Value < 0.01", "Lib.Peptidoform.Q.Value", lambda x: x < 0.01),
        ("PG.Q.Value < 0.05", "PG.Q.Value", lambda x: x < 0.05),
        
        # QUALITY SCORE FILTERS (DIA-NN confidence metrics)
        ("Quantity.Quality > 0.5", "Quantity.Quality", lambda x: x > 0.5),
        ("PG.MaxLFQ.Quality > 0.7", "PG.MaxLFQ.Quality", lambda x: x > 0.7),
    ]
    
    # PEAK WIDTH FILTER (chromatographic quality)
    if 'Peak.Width' in df.columns:
        filters.append((f"Peak.Width >= {min_peak_width}", "Peak.Width", lambda x: x >= min_peak_width))
    
    # DIA-SPECIFIC INTENSITY FILTERS (only add if thresholds are specified)
    # Precursor.Quantity - main DIA quantification at precursor level
    if min_precursor_quantity is not None and 'Precursor.Quantity' in df.columns:
        filters.append((f"Precursor.Quantity > {min_precursor_quantity}", "Precursor.Quantity", 
                       lambda x: x > min_precursor_quantity))
    
    # PG.MaxLFQ - normalized protein-level quantification (recommended for comparisons)
    if min_pg_maxlfq is not None and 'PG.MaxLFQ' in df.columns:
        filters.append((f"PG.MaxLFQ > {min_pg_maxlfq}", "PG.MaxLFQ", 
                       lambda x: x > min_pg_maxlfq))
    
    # PG.Quantity - raw protein-level quantification
    if min_pg_quantity is not None and 'PG.Quantity' in df.columns:
        filters.append((f"PG.Quantity > {min_pg_quantity}", "PG.Quantity", 
                       lambda x: x > min_pg_quantity))
    
    # Apply filters
    print("Applying filters:")
    report.append("Filters Applied:")
    report.append("-" * 60)
    
    df_filtered = df.copy()
    cumulative_removed = 0
    
    for filter_name, column, condition in filters:
        if column not in df_filtered.columns:
            print(f"  ⚠ Skipping {filter_name} - column not found")
            report.append(f"⚠ {filter_name}: SKIPPED - column '{column}' not found")
            continue
        
        # Get current data for this column
        col_data = df_filtered[column]
        
        # Create mask on current filtered dataframe
        mask = condition(col_data)
        rows_before = len(df_filtered)
        
        # Apply filter
        df_filtered = df_filtered[mask]
        rows_after = len(df_filtered)
        removed = rows_before - rows_after
        pct_removed = (removed / rows_before * 100) if rows_before > 0 else 0
        
        print(f"  {filter_name}: {rows_before:,} → {rows_after:,} (removed {removed:,}, {pct_removed:.1f}%)")
        
        # Add to report with statistics
        report.append(f"\n{filter_name}:")
        report.append(f"  Applied to: {rows_before:,} rows")
        report.append(f"  Passed: {rows_after:,} rows")
        report.append(f"  Removed: {removed:,} rows ({pct_removed:.1f}%)")
        
        # Calculate statistics on the removed rows (before filtering)
        if removed > 0:
            removed_data = col_data[~mask]
            if len(removed_data) > 0:
                report.append(f"  Statistics for removed rows:")
                
                # Check if column contains numeric data
                if column == 'Genes' or removed_data.dtype == 'object':
                    # For string columns, just report count
                    report.append(f"    Count: {len(removed_data)}")
                else:
                    # For numeric columns, report statistics
                    report.append(f"    Min: {removed_data.min():.4f}")
                    report.append(f"    Max: {removed_data.max():.4f}")
                    report.append(f"    Mean: {removed_data.mean():.4f}")
                    report.append(f"    Median: {removed_data.median():.4f}")
                    if column == 'Peak.Width':
                        report.append(f"    90th percentile: {removed_data.quantile(0.90):.4f}")
                    if column in intensity_cols:
                        report.append(f"    Values in scientific notation:")
                        report.append(f"      Mean: {removed_data.mean():.2e}")
                        report.append(f"      Median: {removed_data.median():.2e}")
    
    # Summary
    total_removed = initial_rows - len(df_filtered)
    total_pct = (total_removed / initial_rows * 100) if initial_rows > 0 else 0
    
    print(f"\nAfter quality filters: {initial_rows:,} → {len(df_filtered):,} (removed {total_removed:,}, {total_pct:.1f}%)")
    
    # ============================================================
    # CONDITION-BASED FILTERING (requires metadata)
    # ============================================================
    if metadata_file and require_all_conditions and len(df_filtered) > 0:
        print("\nApplying condition-based filtering...")
        report.append("\n" + "=" * 60)
        report.append("CONDITION-BASED FILTERING")
        report.append("=" * 60)
        
        try:
            # Load metadata
            metadata = pd.read_csv(metadata_file)
            print(f"Loaded metadata with {len(metadata)} samples")
            
            # Create run to group mapping
            run_to_group = {}
            for _, row in metadata.iterrows():
                run_name = row['Run']
                group = row['Group']
                run_to_group[run_name] = group
            
            # Get unique groups
            groups = sorted(metadata['Group'].unique())
            print(f"Found {len(groups)} conditions: {', '.join(groups)}")
            report.append(f"Conditions found: {', '.join(groups)}")
            
            # Group proteins by Protein.Group
            rows_before_condition = len(df_filtered)
            
            # For each protein group, check if it's present in all conditions
            protein_groups_to_keep = set()
            
            # Get all unique protein groups
            all_protein_groups = df_filtered['Protein.Group'].unique()
            print(f"Checking {len(all_protein_groups)} protein groups...")
            
            # Check each protein group
            proteins_in_all_conditions = 0
            for pg in all_protein_groups:
                # Get all runs for this protein group
                pg_data = df_filtered[df_filtered['Protein.Group'] == pg]
                runs_with_pg = pg_data['Run'].unique()
                
                # Check which groups have this protein
                groups_with_pg = set()
                for run in runs_with_pg:
                    if run in run_to_group:
                        groups_with_pg.add(run_to_group[run])
                
                # Keep if present in all groups
                if len(groups_with_pg) == len(groups):
                    protein_groups_to_keep.add(pg)
                    proteins_in_all_conditions += 1
            
            # Filter to keep only proteins in all conditions
            df_filtered = df_filtered[df_filtered['Protein.Group'].isin(protein_groups_to_keep)]
            
            rows_after_condition = len(df_filtered)
            removed_by_condition = rows_before_condition - rows_after_condition
            pct_removed_condition = (removed_by_condition / rows_before_condition * 100) if rows_before_condition > 0 else 0
            
            print(f"Protein groups in all {len(groups)} conditions: {proteins_in_all_conditions:,} / {len(all_protein_groups):,}")
            print(f"Condition filter: {rows_before_condition:,} → {rows_after_condition:,} (removed {removed_by_condition:,}, {pct_removed_condition:.1f}%)")
            
            report.append(f"Protein groups checked: {len(all_protein_groups):,}")
            report.append(f"Protein groups in all conditions: {proteins_in_all_conditions:,}")
            report.append(f"Rows before condition filter: {rows_before_condition:,}")
            report.append(f"Rows after condition filter: {rows_after_condition:,}")
            report.append(f"Removed by condition filter: {removed_by_condition:,} ({pct_removed_condition:.1f}%)")
            
        except Exception as e:
            print(f"WARNING: Could not apply condition filter: {e}")
            report.append(f"WARNING: Condition filter failed: {e}")
    
    # Final summary
    final_removed = initial_rows - len(df_filtered)
    final_pct = (final_removed / initial_rows * 100) if initial_rows > 0 else 0
    
    print(f"\nTotal: {initial_rows:,} → {len(df_filtered):,} (removed {final_removed:,}, {final_pct:.1f}%)")
    
    report.append("\n" + "=" * 60)
    report.append("SUMMARY")
    report.append("=" * 60)
    report.append(f"Initial rows: {initial_rows:,}")
    report.append(f"Final rows: {len(df_filtered):,}")
    report.append(f"Total removed: {final_removed:,} ({final_pct:.1f}%)")
    report.append(f"Retention rate: {100 - final_pct:.1f}%")
    
    # Final statistics
    if 'Peak.Width' in df_filtered.columns and len(df_filtered) > 0:
        report.append("\nFinal Peak Width Statistics:")
        report.append(f"  Mean: {df_filtered['Peak.Width'].mean():.3f} min")
        report.append(f"  Median: {df_filtered['Peak.Width'].median():.3f} min")
        report.append(f"  Std Dev: {df_filtered['Peak.Width'].std():.3f} min")
        report.append(f"  Min: {df_filtered['Peak.Width'].min():.3f} min")
        report.append(f"  Max: {df_filtered['Peak.Width'].max():.3f} min")
        report.append(f"  95th percentile: {df_filtered['Peak.Width'].quantile(0.95):.3f} min")
    
    if len(df_filtered) > 0:
        report.append("\nFinal DIA Intensity Statistics:")
        for col in intensity_cols:
            if col in df_filtered.columns:
                report.append(f"\n{col}:")
                report.append(f"  Mean: {df_filtered[col].mean():.2e}")
                report.append(f"  Median: {df_filtered[col].median():.2e}")
                report.append(f"  10th percentile: {df_filtered[col].quantile(0.10):.2e}")
                report.append(f"  90th percentile: {df_filtered[col].quantile(0.90):.2e}")
    
    # Column-wise summary
    report.append("\nColumn Summary:")
    report.append("-" * 60)
    
    # Check which important columns exist
    important_cols = ['Protein.Group', 'Protein.Names', 'Genes', 'Modified.Sequence', 
                     'Precursor.Id', 'Q.Value', 'PG.Q.Value', 'Quantity.Quality']
    
    for col in important_cols:
        if col in df_filtered.columns:
            unique_count = df_filtered[col].nunique()
            report.append(f"{col}: {unique_count:,} unique values")
    
    # Remove the Peak.Width column before saving (it was calculated, not original)
    if 'Peak.Width' in df_filtered.columns:
        df_filtered = df_filtered.drop('Peak.Width', axis=1)
    
    # Save filtered data
    print(f"\nSaving filtered data to {output_file}")
    df_filtered.to_parquet(output_file, index=False)
    
    # Save report
    print(f"Saving report to {report_file}")
    with open(report_file, 'w') as f:
        f.write('\n'.join(report))
    
    # Print summary statistics
    print("\nFiltering complete!")
    print(f"Retained {len(df_filtered):,} out of {initial_rows:,} rows ({100 - total_pct:.1f}%)")
    
    return df_filtered

def main():
    """Main function to handle command line arguments"""
    if len(sys.argv) < 2:
        print("Usage: python diann_parquet_filtering.py <input_file> [options]")
        print("\nOptions:")
        print("  --output OUTPUT_FILE           Filtered parquet file (default: input_filtered.parquet)")
        print("  --report REPORT_FILE           Filter report file (default: input_filter_report.txt)")
        print("  --min-peak-width WIDTH         Minimum peak width in minutes (default: 0.05 = 3 seconds)")
        print("\nDIA-specific intensity filters:")
        print("  --min-precursor-quantity VAL   Minimum Precursor.Quantity (precursor-level)")
        print("  --min-pg-maxlfq VAL           Minimum PG.MaxLFQ (normalized protein-level)")  
        print("  --min-pg-quantity VAL         Minimum PG.Quantity (raw protein-level)")
        print("\nCondition-based filtering:")
        print("  --metadata METADATA_FILE      Metadata CSV file with Run, Sample_Name, Group columns")
        print("  --require-all-conditions      Keep only proteins in at least 1 replicate of ALL conditions")
        print("\nExamples:")
        print("  # Basic filtering (Q-values and peak width only)")
        print("  python diann_parquet_filtering.py data.parquet")
        print("\n  # Add intensity filter at precursor level")
        print("  python diann_parquet_filtering.py data.parquet --min-precursor-quantity 1000")
        print("\n  # Add intensity filter at protein level (recommended)")
        print("  python diann_parquet_filtering.py data.parquet --min-pg-maxlfq 1000")
        print("\n  # Filter for proteins in all conditions")
        print("  python diann_parquet_filtering.py data.parquet --metadata metadata.csv --require-all-conditions")
        print("\n  # Multiple filters")
        print("  python diann_parquet_filtering.py data.parquet --min-peak-width 0.0833 --min-pg-maxlfq 1000 --metadata metadata.csv --require-all-conditions")
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    # Parse command line arguments
    output_file = None
    report_file = None
    min_peak_width = None
    min_precursor_quantity = None
    min_pg_maxlfq = None
    min_pg_quantity = None
    metadata_file = None
    require_all_conditions = False
    
    i = 2
    while i < len(sys.argv):
        if sys.argv[i] == '--output' and i + 1 < len(sys.argv):
            output_file = sys.argv[i + 1]
            i += 2
        elif sys.argv[i] == '--report' and i + 1 < len(sys.argv):
            report_file = sys.argv[i + 1]
            i += 2
        elif sys.argv[i] == '--min-peak-width' and i + 1 < len(sys.argv):
            min_peak_width = float(sys.argv[i + 1])
            i += 2
        elif sys.argv[i] == '--min-precursor-quantity' and i + 1 < len(sys.argv):
            min_precursor_quantity = float(sys.argv[i + 1])
            i += 2
        elif sys.argv[i] == '--min-pg-maxlfq' and i + 1 < len(sys.argv):
            min_pg_maxlfq = float(sys.argv[i + 1])
            i += 2
        elif sys.argv[i] == '--min-pg-quantity' and i + 1 < len(sys.argv):
            min_pg_quantity = float(sys.argv[i + 1])
            i += 2
        elif sys.argv[i] == '--metadata' and i + 1 < len(sys.argv):
            metadata_file = sys.argv[i + 1]
            i += 2
        elif sys.argv[i] == '--require-all-conditions':
            require_all_conditions = True
            i += 1
        else:
            # For backward compatibility - if just filenames are provided
            if i == 2 and not output_file:
                output_file = sys.argv[i]
            elif i == 3 and not report_file:
                report_file = sys.argv[i]
            i += 1
    
    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found")
        sys.exit(1)
    
    # Ensure output file has .parquet extension if specified
    if output_file and not output_file.endswith('.parquet'):
        print(f"Warning: Output file '{output_file}' should have .parquet extension")
        response = input("Add .parquet extension? (y/n): ")
        if response.lower() == 'y':
            output_file = output_file + '.parquet'
    
    # Check if using condition filter without metadata
    if require_all_conditions and not metadata_file:
        print("Error: --require-all-conditions requires --metadata to be specified")
        sys.exit(1)
    
    # Check if metadata file exists when specified
    if metadata_file and not os.path.exists(metadata_file):
        print(f"Error: Metadata file '{metadata_file}' not found")
        sys.exit(1)
    
    # Run filtering
    try:
        filter_diann_parquet(input_file, output_file, report_file, min_peak_width,
                           min_precursor_quantity, min_pg_maxlfq, min_pg_quantity,
                           metadata_file, require_all_conditions)
    except Exception as e:
        print(f"\nError during filtering: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()