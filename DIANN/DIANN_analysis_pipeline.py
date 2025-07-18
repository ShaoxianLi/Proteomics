#!/usr/bin/env python3
"""
DIA-NN Master Pipeline
Orchestrates the complete DIA-NN analysis workflow from raw parquet to final results

Pipeline steps:
1. Filter parquet file (diann_parquet_filtering_P1.py)
2. Extract MaxLFQ values (DIANN_quant_P2.py)
3. Perform statistical analysis and visualization (DIANN_view_P3.py)
"""

import subprocess
import sys
import os
from pathlib import Path
import argparse
import pandas as pd
from datetime import datetime
import shutil

class DIANNPipeline:
    def __init__(self, input_parquet, metadata_file, output_dir=None, min_peak_width=0.05,
                 min_precursor_quantity=None, min_pg_maxlfq=None, min_pg_quantity=None):
        """
        Initialize the DIA-NN analysis pipeline
        
        Parameters:
        -----------
        input_parquet : str
            Path to raw DIA-NN parquet file
        metadata_file : str
            Path to metadata CSV file with columns: Run, Sample_Name, Group, Replicate
        output_dir : str, optional
            Output directory for all results (default: creates timestamped folder)
        min_peak_width : float
            Minimum peak width in minutes for filtering (default: 0.05 = 3 seconds)
        min_precursor_quantity : float, optional
            Minimum Precursor.Quantity intensity threshold
        min_pg_maxlfq : float, optional
            Minimum PG.MaxLFQ intensity threshold (normalized)
        min_pg_quantity : float, optional
            Minimum PG.Quantity intensity threshold (raw)
        """
        self.input_parquet = Path(input_parquet)
        self.metadata_file = Path(metadata_file)
        self.min_peak_width = min_peak_width
        self.min_precursor_quantity = min_precursor_quantity
        self.min_pg_maxlfq = min_pg_maxlfq
        self.min_pg_quantity = min_pg_quantity
        
        # Set up output directory
        if output_dir is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            self.output_dir = self.input_parquet.parent / f"DIANN_Analysis_{timestamp}"
        else:
            self.output_dir = Path(output_dir)
        
        self.output_dir.mkdir(exist_ok=True, parents=True)
        
        # Define paths for intermediate files
        self.filtered_parquet = self.output_dir / f"{self.input_parquet.stem}_filtered.parquet"
        self.maxlfq_csv = self.output_dir / f"{self.input_parquet.stem}_filtered_MaxLFQ.csv"
        self.filter_report = self.output_dir / "filter_report.txt"
        
        # Copy metadata to output directory
        self.metadata_copy = self.output_dir / "metadata.csv"
        shutil.copy2(self.metadata_file, self.metadata_copy)
        
        # Script paths (assuming they're in the same directory as this script)
        script_dir = Path(__file__).parent
        self.filter_script = script_dir / "diann_parquet_filtering_P1.py"
        self.quant_script = script_dir / "DIANN_quant_P2.py"
        self.analysis_script = script_dir / "DIANN_view_P3.py"
        
        # Log file
        self.log_file = self.output_dir / "pipeline_log.txt"
        
    def log(self, message):
        """Log messages to both console and file"""
        print(message)
        with open(self.log_file, 'a') as f:
            f.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - {message}\n")
            
    def check_scripts(self):
        """Check if all required scripts exist"""
        scripts = [self.filter_script, self.quant_script, self.analysis_script]
        missing = []
        
        for script in scripts:
            if not script.exists():
                missing.append(script.name)
                
        if missing:
            self.log(f"ERROR: Missing scripts: {', '.join(missing)}")
            self.log("Please ensure all scripts are in the same directory as this master script")
            return False
            
        return True
        
    def validate_metadata(self):
        """Validate metadata file format"""
        try:
            df = pd.read_csv(self.metadata_file)
            required_cols = ['Run', 'Sample_Name', 'Group', 'Replicate']
            missing_cols = [col for col in required_cols if col not in df.columns]
            
            if missing_cols:
                self.log(f"ERROR: Metadata file missing required columns: {', '.join(missing_cols)}")
                self.log(f"Required columns: {', '.join(required_cols)}")
                return False
                
            # Check for at least 2 groups
            if df['Group'].nunique() < 2:
                self.log("WARNING: Less than 2 groups found in metadata. Statistical comparisons may be limited.")
                
            self.log(f"Metadata validated: {len(df)} samples, {df['Group'].nunique()} groups")
            return True
            
        except Exception as e:
            self.log(f"ERROR: Failed to read metadata file: {e}")
            return False
            
    def run_command(self, cmd, step_name):
        """Run a command and capture output"""
        self.log(f"\n{'='*60}")
        self.log(f"Running: {step_name}")
        self.log(f"Command: {' '.join(cmd)}")
        self.log(f"{'='*60}")
        
        try:
            # Run command and capture output
            result = subprocess.run(cmd, 
                                  capture_output=True, 
                                  text=True, 
                                  check=True)
            
            # Log output
            if result.stdout:
                self.log("Output:")
                self.log(result.stdout)
                
            if result.stderr:
                self.log("Warnings/Info:")
                self.log(result.stderr)
                
            self.log(f"✓ {step_name} completed successfully")
            return True
            
        except subprocess.CalledProcessError as e:
            self.log(f"✗ {step_name} failed with error code {e.returncode}")
            if e.stdout:
                self.log("Output:")
                self.log(e.stdout)
            if e.stderr:
                self.log("Error:")
                self.log(e.stderr)
            return False
            
    def step1_filter_parquet(self):
        """Step 1: Filter the parquet file"""
        cmd = [
            sys.executable,
            str(self.filter_script),
            str(self.input_parquet),
            "--output", str(self.filtered_parquet),
            "--report", str(self.filter_report),
            "--min-peak-width", str(self.min_peak_width)
        ]
        
        # Add intensity filters if specified
        if self.min_precursor_quantity is not None:
            cmd.extend(["--min-precursor-quantity", str(self.min_precursor_quantity)])
        if self.min_pg_maxlfq is not None:
            cmd.extend(["--min-pg-maxlfq", str(self.min_pg_maxlfq)])
        if self.min_pg_quantity is not None:
            cmd.extend(["--min-pg-quantity", str(self.min_pg_quantity)])
        
        return self.run_command(cmd, "Step 1: Parquet Filtering")
        
    def step2_extract_maxlfq(self):
        """Step 2: Extract MaxLFQ values"""
        cmd = [
            sys.executable,
            str(self.quant_script),
            str(self.filtered_parquet)
        ]
        
        return self.run_command(cmd, "Step 2: MaxLFQ Extraction")
        
    def step3_statistical_analysis(self):
        """Step 3: Run statistical analysis and visualization"""
        cmd = [
            sys.executable,
            str(self.analysis_script),
            str(self.maxlfq_csv),
            str(self.metadata_copy)
        ]
        
        return self.run_command(cmd, "Step 3: Statistical Analysis")
        
    def create_summary_report(self):
        """Create a final summary report"""
        self.log("\nCreating pipeline summary report...")
        
        summary = []
        summary.append("DIA-NN PIPELINE SUMMARY REPORT")
        summary.append("=" * 60)
        summary.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        summary.append(f"\nInput Files:")
        summary.append(f"  - Parquet: {self.input_parquet}")
        summary.append(f"  - Metadata: {self.metadata_file}")
        summary.append(f"\nOutput Directory: {self.output_dir}")
        summary.append(f"\nParameters:")
        summary.append(f"  - Min Peak Width: {self.min_peak_width} minutes")
        if self.min_precursor_quantity:
            summary.append(f"  - Min Precursor Quantity: {self.min_precursor_quantity}")
        if self.min_pg_maxlfq:
            summary.append(f"  - Min PG MaxLFQ: {self.min_pg_maxlfq}")
        if self.min_pg_quantity:
            summary.append(f"  - Min PG Quantity: {self.min_pg_quantity}")
        
        # Add filtering summary if available
        if self.filter_report.exists():
            summary.append(f"\nFiltering Summary:")
            with open(self.filter_report, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    if "Initial rows:" in line or "Final rows:" in line or "Retention rate:" in line:
                        summary.append(f"  {line.strip()}")
        
        # Add quantification summary if available
        if self.maxlfq_csv.exists():
            df = pd.read_csv(self.maxlfq_csv, index_col=0)
            summary.append(f"\nQuantification Summary:")
            summary.append(f"  - Total proteins: {len(df)}")
            summary.append(f"  - Samples: {len([col for col in df.columns if col not in ['Protein.Names', 'Genes', 'N_Samples', 'Mean_Intensity', 'CV_%']])}")
            
        # List key output files
        summary.append(f"\nKey Output Files:")
        key_files = [
            ("Filtered data", self.filtered_parquet),
            ("MaxLFQ matrix", self.maxlfq_csv),
            ("Filter report", self.filter_report),
            ("Analysis results", self.output_dir / "analysis_results"),
            ("Pipeline log", self.log_file)
        ]
        
        for desc, path in key_files:
            if Path(path).exists():
                summary.append(f"  - {desc}: {path.name}")
                
        # Add analysis results summary
        analysis_dir = self.output_dir / "analysis_results"
        if analysis_dir.exists():
            summary.append(f"\nAnalysis Results Generated:")
            for file in sorted(analysis_dir.iterdir()):
                if file.is_file():
                    summary.append(f"  - {file.name}")
                    
        # Save summary
        summary_file = self.output_dir / "pipeline_summary.txt"
        with open(summary_file, 'w') as f:
            f.write('\n'.join(summary))
            
        self.log(f"Summary report saved to: {summary_file}")
        
    def run_pipeline(self):
        """Run the complete pipeline"""
        self.log("=" * 80)
        self.log("DIA-NN MASTER PIPELINE")
        self.log("=" * 80)
        self.log(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        self.log(f"Output directory: {self.output_dir}")
        
        # Check prerequisites
        if not self.check_scripts():
            self.log("Pipeline aborted: Missing required scripts")
            return False
            
        if not self.validate_metadata():
            self.log("Pipeline aborted: Invalid metadata file")
            return False
            
        # Run pipeline steps
        steps = [
            ("Filtering parquet file", self.step1_filter_parquet),
            ("Extracting MaxLFQ values", self.step2_extract_maxlfq),
            ("Statistical analysis", self.step3_statistical_analysis)
        ]
        
        for step_name, step_func in steps:
            if not step_func():
                self.log(f"\nPipeline failed at: {step_name}")
                self.log("Check the log file for details")
                return False
                
        # Create summary report
        self.create_summary_report()
        
        self.log("\n" + "=" * 80)
        self.log("PIPELINE COMPLETED SUCCESSFULLY!")
        self.log(f"All results saved to: {self.output_dir}")
        self.log(f"End time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        self.log("=" * 80)
        
        return True

def create_example_metadata():
    """Create an example metadata file"""
    example = """Run,Sample_Name,Group,Replicate
20240301_Sample1_DIA.raw,Sample1,Control,1
20240301_Sample2_DIA.raw,Sample2,Control,2
20240301_Sample3_DIA.raw,Sample3,Control,3
20240301_Sample4_DIA.raw,Sample4,Treatment,1
20240301_Sample5_DIA.raw,Sample5,Treatment,2
20240301_Sample6_DIA.raw,Sample6,Treatment,3"""
    
    with open("metadata_example.csv", 'w') as f:
        f.write(example)
    print("Created example metadata file: metadata_example.csv")
    print("\nMetadata file should contain the following columns:")
    print("- Run: The run/file name (should match what's in the parquet file)")
    print("- Sample_Name: A clean sample name for display")
    print("- Group: The experimental group/condition")
    print("- Replicate: The replicate number within each group")

def main():
    parser = argparse.ArgumentParser(
        description="DIA-NN Master Pipeline - Orchestrates filtering, quantification, and analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage (Q-values and peak width only)
  python diann_master_pipeline.py data.parquet metadata.csv
  
  # Specify output directory
  python diann_master_pipeline.py data.parquet metadata.csv --output my_analysis
  
  # Add intensity filter at protein level (recommended for DIA)
  python diann_master_pipeline.py data.parquet metadata.csv --min-pg-maxlfq 1000
  
  # Multiple filters
  python diann_master_pipeline.py data.parquet metadata.csv --min-peak-width 0.0833 --min-pg-maxlfq 1000
  
  # Create example metadata file
  python diann_master_pipeline.py --create-metadata-example
        """
    )
    
    parser.add_argument('input_parquet', nargs='?', help='Input DIA-NN parquet file')
    parser.add_argument('metadata_file', nargs='?', help='Metadata CSV file')
    parser.add_argument('--output', '-o', dest='output_dir', 
                       help='Output directory (default: timestamped folder)')
    parser.add_argument('--min-peak-width', type=float, default=0.05,
                       help='Minimum peak width in minutes (default: 0.05 = 3 seconds)')
    parser.add_argument('--min-precursor-quantity', type=float, default=None,
                       help='Minimum Precursor.Quantity (precursor-level intensity)')
    parser.add_argument('--min-pg-maxlfq', type=float, default=None,
                       help='Minimum PG.MaxLFQ (normalized protein intensity, recommended)')
    parser.add_argument('--min-pg-quantity', type=float, default=None,
                       help='Minimum PG.Quantity (raw protein intensity)')
    parser.add_argument('--create-metadata-example', action='store_true',
                       help='Create an example metadata file and exit')
    
    args = parser.parse_args()
    
    # Handle metadata example creation
    if args.create_metadata_example:
        create_example_metadata()
        return
        
    # Check required arguments
    if not args.input_parquet or not args.metadata_file:
        parser.print_help()
        sys.exit(1)
        
    # Check input files exist
    if not Path(args.input_parquet).exists():
        print(f"Error: Input parquet file not found: {args.input_parquet}")
        sys.exit(1)
        
    if not Path(args.metadata_file).exists():
        print(f"Error: Metadata file not found: {args.metadata_file}")
        sys.exit(1)
        
    # Run pipeline
    pipeline = DIANNPipeline(
        args.input_parquet,
        args.metadata_file,
        args.output_dir,
        args.min_peak_width,
        args.min_precursor_quantity,
        args.min_pg_maxlfq,
        args.min_pg_quantity
    )
    
    success = pipeline.run_pipeline()
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()