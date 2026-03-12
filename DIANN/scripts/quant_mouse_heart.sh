#!/bin/bash
#BSUB -J DIANN_mouse_heart                 # Job name
#BSUB -q long                              # Queue: up to 8 h
#BSUB -n 32                                # 32 cores
#BSUB -R "rusage[mem=8192] span[hosts=1]"  # 8 GB per core, single node
#BSUB -W 8:00                              # Wall-time limit
#BSUB -u shaoxian.li11@umassmed.edu        # Email notifications
#BSUB -N                                   # Send email when done
#BSUB -o mouse_heart_%J.out                # stdout
#BSUB -e mouse_heart_%J.err                # stderr

set -euo pipefail

# ============================================================
# 用户配置区域
# ============================================================
RAW_DATA_DIR="/home/shaoxian.li11-umw/DIANN/raw/JZ"
RAW_FILE="SL_Heart_test_DIA_260311_aa01679.raw"
FASTA_DIR="/home/shaoxian.li11-umw/DIANN/fasta"
FASTA_FILE="mouse_uniprotkb_proteome_UP000000589_2026_03_12.fasta"
SPECLIB_DIR="/pi/qing.yu-umw/home/shaoxian.li11-umw/DIANN/lib"
SPECLIB_FILE="mouse_predicted.predicted.speclib"
RESULTS_DEST="/pi/qing.yu-umw/home/shaoxian.li11-umw/DIANN/results/JZ_heart"
# ============================================================

export TMPDIR=${TMPDIR:-/tmp}
export DIANNIMG=/share/pkg/containers/diann/diann_2.1.0.sif
echo "Scratch TMPDIR: $TMPDIR"
echo "Using container: $DIANNIMG"

module load diann/2.1.0

# 准备临时工作目录
mkdir -p "$TMPDIR"/{raw,lib,fasta,results}
cp "${RAW_DATA_DIR}/${RAW_FILE}" "$TMPDIR/raw/"
cp "${SPECLIB_DIR}/${SPECLIB_FILE}" "$TMPDIR/lib/"
cp "${FASTA_DIR}/${FASTA_FILE}" "$TMPDIR/fasta/"

cd "$TMPDIR"

echo ">>> Starting DIA-NN analysis..."
echo ">>> RAW file: ${RAW_FILE}"
echo ">>> Spectral library: ${SPECLIB_FILE}"
echo ">>> FASTA: ${FASTA_FILE}"

singularity exec "$DIANNIMG" /diann-2.1.0/diann-linux \
  --verbose 2 \
  --threads 32 \
  --f "raw/${RAW_FILE}" \
  --lib "lib/${SPECLIB_FILE}" \
  --fasta "fasta/${FASTA_FILE}" \
  --qvalue 0.01 \
  --peptidoforms \
  --export-quant \
  --matrices \
  --matrix-qvalue 0.01 \
  --matrix-spec-q 0.05 \
  --protein-inference \
  --pg-level 1 \
  --var-mods 1 \
  --var-mod "UniMod:35,15.994915,M" \
  --fixed-mod "UniMod:4,57.021464,C" \
  --met-excision \
  --reannotate \
  --out results/SL_Heart_test_DIA_report.parquet

echo ">>> Contents of results/:"
ls -lh results/

# 拷贝结果到永久存储
mkdir -p "$RESULTS_DEST"
rsync -av results/SL_Heart_test_DIA_report.parquet "$RESULTS_DEST/"
rsync -av results/*_matrix.tsv "$RESULTS_DEST/"

echo "Done! Results saved to ${RESULTS_DEST}/"
