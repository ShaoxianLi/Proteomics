#!/bin/bash
#BSUB -q gpu
#BSUB -J "af3_json[1-1341]"   # ⬅️ Update to match number of JSON files
#BSUB -gpu "num=1:mode=exclusive_process:j_exclusive=yes"
#BSUB -R "rusage[mem=2G] span[hosts=1]"
#BSUB -n 8
#BSUB -W 8:00

# Path to AF3 container
export ALPHAIMG=/share/pkg/containers/alphafold/alphafold_unified-mem-3.0.0.sif

# JSON input list and folders
JSON_LIST="$HOME/af3_input_files/MAPLC3B_protein_pairs_list.txt" # ⬅️ change
INPUT_DIR="$HOME/af3_input_files/MAPLC3B_protein_pairs" # ⬅️ change
OUTPUT_DIR="$HOME/af3_output_files/MAPLC3B_protein" # ⬅️ change

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

