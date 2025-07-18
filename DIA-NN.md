# How to use DIA-NN in HPC?
DIA-NN is an automated software suite for data-independent acquisition (DIA) proteomics data processing.
You can found more information here: https://github.com/vdemichev/DiaNN

Now you are ready to work with diann!

## Load diann module
### diann-2.1.0 is available under module control as module 'diann/2.1.0'
- (Optional) Check all the modules available: ``` module avail```
Use ```q``` to exit out the module viewing
- Load the module: ```module load diann/2.1.0```
- Run diann with the following command
For an interactive shell:
```bash
bsub -q interactive -n16 -W 4:00 -R "rusage[mem=4G] span[hosts=1]" -Is singularity shell $DIANNIMG /bin/bash
```
For a non-interactive job:
``` bash
   bsub -q short singularity exec $DIANNIMG /diann-2.1.0/diann-linux -h
```

Now you are ready to run diann!

## Step 1: Generate the spectral library
```bash
diann-linux \
  --verbose 4 \
  --threads 16 \
  --predictor \
  --gen-spec-lib \
  --fasta-search \
  --fasta uniprotkb_proteome_UP000005640_2025_04_26.fasta \
  --cut K*,R* \
  --var-mods 1 \
  --var-mod UniMod:35,15.9949146221,M \
  --fixed-mod UniMod:4,57.021464,C \
  --missed-cleavages 2 \
  --max-pep-len 30 \
  --min-pep-len 7 \
  --max-pr-charge 5 \
  --min-pr-charge 2 \
  --max-pr-charge 980 \
  --min-pr-mz 380 \
  --pg-level 1 \
  --peptidoforms \
  --out-lib uniprot_predicted.speclib
```

For generating PTM spectral library (e.g. phosphorylation)
```bash
diann-linux \
  --verbose 4 \
  --threads 16 \
  --predictor \
  --gen-spec-lib \
  --fasta-search \
  --fasta uniprotkb_proteome_UP000005640_2025_04_26.fasta \
  --cut K*,R* \
  --var-mods 3 \
  --var-mod UniMod:35,15.9949146221,M \
  --var-mod UniMod:21,79.966331,STY \
  --fixed-mod UniMod:4,57.021464,C \
  --missed-cleavages 1 \
  --max-pep-len 30 \
  --min-pep-len 7 \
  --max-pr-charge 5 \
  --min-pr-charge 2 \
  --max-pr-charge 980 \
  --min-pr-mz 380 \
  --pg-level 1 \
  --peptidoforms \
  --out-lib uniprot_phos.speclib
```

## Step 2: Analyze the raw data using the generated spectral library
```bash
diann-linux \
  --verbose 4 \
  --threads 16 \
  --f qymp1_00254_qy_dia_8min_rpe1_3.raw \
  --fasta uniprotkb_proteome_UP000005640_2025_04_26.fasta \
  --lib uniprot_predicted.predicted.speclib \
  --peptidoforms \
  --export-quant \
  --reannotate \
  --out qymp1_00254_qy_dia_8min_rpe1_3_report.parquet
```
## Write all the command into a script for simplier execution (for generating phospho-library)
- Create a .sh file to submit the job and run the command: ```nano diann_job.sh```
- It will open the nano text editor. You can paste the following example command into the text. This is the example for generating a spectral_library.
``` bash
   #!/bin/sh
   #BSUB -q short
   #BSUB -n 32
   #BSUB -R "rusage[mem=8G] span[hosts=1]"
   #BSUB -W 8:00
   #BSUB -u fangyi.zhai@umassmed.edu
   #BSUB -N
   #BSUB -o "%phospho.out"
   #BSUB -e "%phospho.err"

   # Go to the diann directory
   cd diann

   # Load DIANN module
   module load diann/2.1.0
   
   # Load into Singularity container and open shell
   singularity exec $DIANNIMG /diann-2.1.0/diann-linux \
     --verbose 4 \
     --threads 32 \
     --predictor \
     --gen-spec-lib \
     --fasta-search \
     --fasta uniprotkb_proteome_UP000005640_2025_04_26.fasta \
     --cut K*,R* \
     --var-mods 3 \
     --var-mod UniMod:35,15.9949146221,M \
     --var-mod UniMod:21,79.966331,STY \
     --fixed-mod UniMod:4,57.021464,C \
     --missed-cleavages 1 \
     --max-pep-len 30 \
     --min-pep-len 7 \
     --max-pr-charge 5 \
     --min-pr-charge 2 \
     --max-pr-charge 980 \
     --min-pr-mz 380 \
     --pg-level 1 \
     --peptidoforms \
     --out-lib uniprot_phos.speclib
```
## If you want to run multiple files in parellel
```bash
#!/bin/bash

# List of input .raw files
FILES=(*.raw)

# Build the --f flags for all files
FILE_FLAGS=""
for FILE in "${FILES[@]}"; do
  FILE_FLAGS+="--f ${FILE} "
done

# Submit a single job for all files
bsub -q short \
     -J "diann" \
     -n 32 \
     -R "rusage[mem=8G] span[hosts=1]" \
     -W 4:00 \
     -o "quant.out" \
     -e "quant.err" \
     "module load diann/2.1.0 && \
      singularity exec \$DIANNIMG /diann-2.1.0/diann-linux \
        --verbose 4 \
        --threads 32 \
        $FILE_FLAGS \
        --fasta ../human_yeast_proteome.fasta \
        --lib ../uniprot_human_yeast.predicted.speclib \
        --matrices \
        --matrix-qvalue 0.01 \
        --matrix-spec-q 0.05 \
        --peptidoforms \
        --export-quant \
        --reannotate \
        --reanalyse \
        --out report.parquet"

```
# Troubleshooting
## How to find files in the Singularity named diann and find the executive line
- (Optional) Look for the all files in the Singularity named diann:
```bash
find / -name "*diann*" 2>/dev/null
```
You should see something like
```
/home/fangyi.zhai-umw/DIA-NN-2/DIA-NN-2/diann-2.1.0
/home/fangyi.zhai-umw/DIA-NN-2/DIA-NN-2/diann-2.1.0/diann-linux
/home/fangyi.zhai-umw/DIA-NN-2/DIA-NN-2/diann-2.1.0/diann-stats.py
/home/fangyi.zhai-umw/micromamba/envs/diann
```
The executive line for diann is ```diann-linux```
