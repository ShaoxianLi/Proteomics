#!/bin/bash
#BSUB -J DIANN_mouse_speclib               # Job name
#BSUB -q long                              # Queue: up to 8 h
#BSUB -n 32                                # 32 cores
#BSUB -R "rusage[mem=8192] span[hosts=1]"  # 8 GB per core, single node
#BSUB -W 8:00                              # Wall-time limit
#BSUB -u shaoxian.li11@umassmed.edu        # Email notifications
#BSUB -N                                   # Send email when done
#BSUB -o mouse_speclib_%J.out              # stdout
#BSUB -e mouse_speclib_%J.err              # stderr

set -euo pipefail

# ============================================================
# 用户配置区域
# ============================================================
FASTA_DIR="/home/shaoxian.li11-umw/DIANN/fasta"
FASTA_FILE="mouse_uniprotkb_proteome_UP000000589_2026_03_12.fasta"
OUTPUT_LIB="mouse_predicted.speclib"
LIB_DEST="/pi/qing.yu-umw/home/shaoxian.li11-umw/DIANN/lib"
# ============================================================

export TMPDIR=${TMPDIR:-/tmp}
export DIANNIMG=/share/pkg/containers/diann/diann_2.1.0.sif
echo "Scratch TMPDIR: $TMPDIR"
echo "Using container: $DIANNIMG"

module load diann/2.1.0

# 准备临时工作目录
mkdir -p "$TMPDIR"/{fasta,lib}
cp "${FASTA_DIR}/${FASTA_FILE}" "$TMPDIR/fasta/"

cd "$TMPDIR"

echo ">>> Starting spectral library generation..."
echo ">>> FASTA: ${FASTA_FILE}"

singularity exec "$DIANNIMG" /diann-2.1.0/diann-linux \
  --verbose 4 \
  --threads 32 \
  --predictor \
  --gen-spec-lib \
  --fasta-search \
  --fasta "fasta/${FASTA_FILE}" \
  --cut K*,R* \
  --var-mods 1 \
  --var-mod "UniMod:35,15.994915,M" \
  --fixed-mod "UniMod:4,57.021464,C" \
  --missed-cleavages 1 \
  --max-pep-len 30 \
  --min-pep-len 7 \
  --max-pr-charge 5 \
  --min-pr-charge 2 \
  --max-pr-mz 980 \
  --min-pr-mz 380 \
  --pg-level 1 \
  --peptidoforms \
  --met-excision \
  --out-lib "lib/${OUTPUT_LIB}"

echo ">>> Spectral library generation complete!"
echo ">>> Contents of lib/:"
ls -lh lib/

# 拷贝结果到永久存储
mkdir -p "$LIB_DEST"
rsync -av lib/*.speclib "$LIB_DEST/"

echo "Done! Spectral library saved to ${LIB_DEST}/"
