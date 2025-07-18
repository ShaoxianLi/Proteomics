## Batch AF3 Script
### STEP 1: Generate JSON files in batch
- Create 2 folders in the home directory named:
  ```
  af3_input_files
  ```
  af3_iutput_files is where all JSON files are stored
  
  ```
  af3_output_files
  ```
### STEP 2: Generate Your JSON File List
- Before submitting this script, run the following command once:
```
ls $HOME/af3_input_files/citrate_hits/*.json > $HOME/af3_input_files/citrate_hits_list.txt
```

- Then confirm the number of JSONs (for updating the #BSUB -J "af3_json[1-1293]" line):
```
wc -l $HOME/af3_input_files/citrate_hits_list.txt
```

- Save this as something like ```submit_af3_json_jobs.bash```

### STEP 3: submit the job
```
#!/bin/bash
#BSUB -q gpu
#BSUB -J "af3_json[1-1293]"   # ⬅️ Update to match number of JSON files
#BSUB -gpu "num=1:mode=exclusive_process:j_exclusive=yes"
#BSUB -R "rusage[mem=4G] span[hosts=1]"
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
