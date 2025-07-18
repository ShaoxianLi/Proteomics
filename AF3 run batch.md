### STEP 1: Generate JSON files in batch and file List
Make sure this script is in the same folder with the CSV file download from TMT viewer. The UniprotID should be in column 2 and GeneSymbol should be in column 3 
```
import os
import csv
import json
from Bio import ExPASy, SwissProt
# Define protein A as MAPLC3B
protein_A_uniprot = "Q9GZQ8"  # MAPLC3B的UniProt ID
protein_A_id = "A"

# Output folder
folder_name = "MAPLC3B_protein_pairs"
os.makedirs(folder_name, exist_ok=True)

# Input CSV with protein UniProt IDs and Gene Symbols
csv_file = "SL18_Live_pq_153 modified.csv"  # 你的实际CSV文件名

# CSV列索引配置
# 根据你的CSV格式: EPPK1__P58107, P58107, EPPK1, Epiplakin, ...
UNIPROT_COLUMN = 1  # UniProt ID所在列（第2列，0-based索引为1）
GENE_SYMBOL_COLUMN = 2  # Gene Symbol所在列（第3列，0-based索引为2）

# First, get MAPLC3B sequence
try:
    handle_A = ExPASy.get_sprot_raw(protein_A_uniprot)
    record_A = SwissProt.read(handle_A)
    seq_A = record_A.sequence
    entry_name_A = record_A.entry_name
    print(f"Successfully retrieved MAPLC3B sequence: {entry_name_A}")
except Exception as e:
    print(f"Failed to retrieve MAPLC3B sequence: {e}")
    exit(1)

# Process each protein target (chainB) 
processed_count = 0
failed_count = 0

with open(csv_file, "r", encoding="utf-8-sig") as file:
    reader = csv.reader(file)
    header = next(reader)  # Skip header line
    print(f"CSV header: {header}")
    
    # 验证列索引是否有效
    if len(header) <= max(UNIPROT_COLUMN, GENE_SYMBOL_COLUMN):
        print(f"Error: CSV doesn't have enough columns. Expected at least {max(UNIPROT_COLUMN, GENE_SYMBOL_COLUMN)+1} columns.")
        exit(1)
    
    print(f"Using column {UNIPROT_COLUMN} for UniProt ID: {header[UNIPROT_COLUMN]}")
    print(f"Using column {GENE_SYMBOL_COLUMN} for Gene Symbol: {header[GENE_SYMBOL_COLUMN]}")
    
    for row_num, row in enumerate(reader, start=2):
        if len(row) <= max(UNIPROT_COLUMN, GENE_SYMBOL_COLUMN):
            print(f"Row {row_num}: Insufficient columns")
            continue
            
        chainB_id = row[UNIPROT_COLUMN].strip()  # UniProt ID
        gene_symbol = row[GENE_SYMBOL_COLUMN].strip()  # Gene Symbol
        
        if not chainB_id:
            print(f"Row {row_num}: Empty UniProt ID")
            continue
            
        if not gene_symbol:
            print(f"Row {row_num}: Empty Gene Symbol for {chainB_id}")
            continue
            
        try:
            # Fetch protein sequence for chainB
            print(f"Processing {chainB_id} ({gene_symbol})...")
            handle_B = ExPASy.get_sprot_raw(chainB_id)
            record_B = SwissProt.read(handle_B)
            seq_B = record_B.sequence
            
            # Build JSON structure for protein-protein interaction
            json_data = {
                "name": f"MAPLC3B_{gene_symbol}",
                "sequences": [
                    {
                        "protein": {
                            "id": protein_A_id,
                            "sequence": seq_A
                        }
                    },
                    {
                        "protein": {
                            "id": "B", 
                            "sequence": seq_B
                        }
                    }
                ],
                "modelSeeds": [1],
                "dialect": "alphafold3",
                "version": 1
            }
            
            # Save JSON file using gene symbol
            output_filename = f"pair_MAPLC3B_{gene_symbol}.json"
            output_path = os.path.join(folder_name, output_filename)
            
            # 检查文件是否已存在（避免重复基因名的冲突）
            counter = 1
            while os.path.exists(output_path):
                output_filename = f"pair_MAPLC3B_{gene_symbol}_{counter}.json"
                output_path = os.path.join(folder_name, output_filename)
                counter += 1
            
            with open(output_path, "w") as json_file:
                json.dump(json_data, json_file, indent=2)
                
            print(f"✓ Wrote: {output_path}")
            processed_count += 1
            
        except Exception as e:
            print(f"✗ Failed to process {chainB_id} ({gene_symbol}): {e}")
            failed_count += 1

# Create file list for batch submission
list_file_path = f"{folder_name}_list.txt"
json_files = [f for f in os.listdir(folder_name) if f.endswith('.json')]
json_files.sort()  # Sort for consistent ordering

with open(list_file_path, 'w') as list_file:
    for json_file in json_files:
        list_file.write(f"{json_file}\n")

print(f"\n=== Processing Complete ===")
print(f"Successfully processed: {processed_count}")
print(f"Failed: {failed_count}")
print(f"Files saved in '{folder_name}' folder.")
print(f"Created file list: {list_file_path}")
print(f"Total JSON files: {len(json_files)}")
print(f"Remember to update the job array size in your batch script to [1-{len(json_files)}]")
```
### STEP 3: submit the job
- Before submitting this script, run the following command once:
```
ls $HOME/af3_input_files/citrate_hits/*.json > $HOME/af3_input_files/citrate_hits_list.txt
```

- Then confirm the number of JSONs (for updating the #BSUB -J "af3_json[1-1293]" line):
```
wc -l $HOME/af3_input_files/citrate_hits_list.txt
```

- Save this as something like ```submit_af3_json_jobs.bash```

```
#!/bin/bash
#BSUB -q gpu
#BSUB -J "af3_json[1-1293]"   # ⬅️ Update to match number of JSON files
#BSUB -gpu "num=1:mode=exclusive_process:j_exclusive=yes"
#BSUB -R "rusage[mem=2G] span[hosts=1]"
#BSUB -n 8
#BSUB -W 8:00

# Path to AF3 container
export ALPHAIMG=/share/pkg/containers/alphafold/alphafold_unified-mem-3.0.0.sif

# JSON input list and folders
JSON_LIST="$HOME/af3_input_files/citrate_hits_demo_list.txt" # ⬅️ change
INPUT_DIR="$HOME/af3_input_files/citrate_hits_demo" # ⬅️ change
OUTPUT_DIR="$HOME/af3_output_files/citrate_hits_demo" # ⬅️ change

# Make sure output folder exists
mkdir -p "$OUTPUT_DIR"

# Get JSON file for this job index
JSON_PATH=$(sed -n "${LSB_JOBINDEX}p" "$JSON_LIST")
BASENAME=$(basename "$JSON_PATH" .json)

# Redirect logs to output folder
exec > "$OUTPUT_DIR/${BASENAME}.out" 2> "$OUTPUT_DIR/${BASENAME}.err"

# Run AlphaFold3 prediction
singularity exec --nv --bind \
/home/qing.yu-umw/public_databases:/databases,\
/home/qing.yu-umw/models:/models,\
${INPUT_DIR}:/input,\
${OUTPUT_DIR}:/output \
${ALPHAIMG} \
env XLA_FLAGS="--xla_gpu_cuda_data_dir=/usr/local/cuda-12.2 --xla_disable_hlo_passes=custom-kernel-fusion-rewriter" \
python /app/alphafold/run_alphafold.py \
--db_dir=/databases \
--model_dir=/models \
--json_path="/input/$(basename "$JSON_PATH")" \
--output_dir=/output \
--flash_attention_implementation xla

```
