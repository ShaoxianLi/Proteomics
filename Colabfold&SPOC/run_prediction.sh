#!/bin/bash
#BSUB -q large
#BSUB -n 40
#BSUB -J spoc_predict
#BSUB -W 96:00
#BSUB -oo spoc_predict.out
#BSUB -eo spoc_predict.err

# Setup micromamba
eval "$(micromamba shell hook --shell=bash)"
micromamba activate spoc_venv

# Define the path to run_wrapper.py
RUN_WRAPPER_PATH="/home/shaoxian.li11-umw/SPOC/run_wrapper.py"

# Go to the directory that contains all the predictions
cd /home/shaoxian.li11-umw/colabfold/All_multimer_MAPLC3B

for folder in */; do
    folder=${folder%/}
    echo "Running on $folder"
    python3 "$RUN_WRAPPER_PATH" "$folder"
done
