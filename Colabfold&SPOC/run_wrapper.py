# run_wrapper.py
import os
import sys
import glob
import subprocess
import tempfile
import re
def run_prediction(input_folder):
    input_folder = os.path.abspath(input_folder)
    folder_name = os.path.basename(input_folder)
    # Match only top-level files
    a3m_files = glob.glob(os.path.join(input_folder, '*.a3m'))
    json_files = sorted(glob.glob(os.path.join(input_folder, '*_scores_rank_*_model_*_seed_000.json')))
    pdb_files = sorted(glob.glob(os.path.join(input_folder, '*_unrelaxed_rank_*_model_*_seed_000.pdb')))
    all_files = a3m_files + json_files + pdb_files

    if not all_files:
        print(f"[{input_folder}] No valid input files found.")
        return

    # Temporary clean directory with only relevant files
    with tempfile.TemporaryDirectory() as tmpdir:
        for file_path in all_files:
            link_path = os.path.join(tmpdir, os.path.basename(file_path))
            os.symlink(os.path.abspath(file_path), link_path)

        print(f"[{input_folder}] Running prediction...")
        output_file = f"{folder_name}_SPOC_output.csv"
        subprocess.run(["python3", "run.py", tmpdir,"--output",output_file], cwd=os.path.dirname(__file__), check=True)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python run_wrapper.py <folder1> [<folder2> ...]")
        sys.exit(1)

    for folder in sys.argv[1:]:
        run_prediction(folder)
