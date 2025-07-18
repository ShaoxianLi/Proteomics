DIANN HPC 提交指南
🔧 1. 快速配置
修改脚本中的路径：
# 必须修改的部分（在脚本的"用户配置区域"）：

RAW_DATA_DIR="/path/to/your/raw/files"                    # 你的.raw文件目录
SPECLIB_FILE="/path/to/your/library.speclib"              # 你的光谱库文件路径
RESULTS_DIR="/path/to/your/results"                       # 结果保存目录
OUTPUT_PREFIX="your_experiment"                           # 实验名称前缀

# 还要修改邮箱：
#BSUB -u your.email@umassmed.edu
🚀 2. 提交作业
# 1. 保存脚本
nano diann_annotation.sh
# 复制脚本内容，修改路径，保存退出

# 2. 设置执行权限
chmod +x diann_annotation.sh

# 3. 提交作业
bsub < diann_annotation.sh

# 4. 检查作业状态
bjobs
bjobs -l JOBID  # 查看详细信息
📊 3. 监控作业
# 查看输出（实时）
tail -f diann_JOBID.out

# 查看错误日志
tail -f diann_JOBID.err

# 取消作业
bkill JOBID
🔍 4. 关键改进点
✅ Speclib优化搜索：
•	--lib - 使用光谱库搜索
•	--export-quant - 导出定量信息
•	--protein-inference - 蛋白质推断
•	--pg-level 1 - 蛋白质组级别分析
•	--peptidoforms - 肽形式分析
✅ 更好的输出检查：
•	自动检查 Protein.Names 和 Genes 列
•	统计注释覆盖率
•	详细的文件验证
✅ 增强的错误处理：
•	文件存在性检查
•	运行状态验证
•	更详细的日志记录
🎯 5. 仅使用Speclib搜索
这个脚本已经简化为只使用speclib文件进行搜索：
1.	纯Library搜索：
2.	--lib lib/"$SPECLIB_BASENAME"
3.	优化的输出格式：
4.	--matrices \
5.	--export-quant \
6.	--protein-inference \
7.	质量控制：
8.	--qvalue 0.01 \
9.	--matrix-qvalue 0.01 \
10.	--matrix-spec-q 0.05 \
📋 6. 检查清单
提交前确认：
•	[ ] 修改了所有路径
•	[ ] 修改了邮箱地址
•	[ ] RAW文件目录正确
•	[ ] Speclib文件路径正确
•	[ ] 有足够的磁盘空间（建议100GB+）
🆘 7. 常见问题
如果作业失败：
# 检查错误日志
cat diann_JOBID.err

# 检查输出日志
cat diann_JOBID.out

# 检查临时目录空间
df -h $TMPDIR
如果注释仍然缺失：
1.	检查生成的 pg_matrix.tsv 文件
2.	确认speclib文件包含完整的蛋白质注释
3.	考虑使用包含注释的FASTA重新生成speclib
4.	检查DIANN版本兼容性
💡 8. 优化建议
对于大数据集：
•	增加内存：#BSUB -R "rusage[mem=12288]"
•	增加时间：#BSUB -W 12:00
•	使用更多核心：#BSUB -n 64
对于小数据集：
•	减少资源：#BSUB -n 16 和 #BSUB -R "rusage[mem=4096]"
•	使用短队列：#BSUB -q short

有desthiobiotin mod：
#!/bin/bash
#BSUB -J DIANN_quant                      # Job name
#BSUB -q long                            # Queue: up to 8 h
#BSUB -n 32                              # 32 cores
#BSUB -R "rusage[mem=8192] span[hosts=1]" # 8 GB per core (256 GB total), single node
#BSUB -W 8:00                            # Wall-time limit: 8 h
#BSUB -u shaoxian.li11@umassmed.edu      # Email notifications
#BSUB -o diann_%J.out                    # stdout
#BSUB -e diann_%J.err                    # stderr

set -euo pipefail

export TMPDIR=${TMPDIR:-/tmp}
export DIANNIMG=/share/pkg/containers/diann/diann_2.1.0.sif
echo "Scratch TMPDIR: $TMPDIR"
echo "Using container: $DIANNIMG"

module load diann/2.1.0

mkdir -p "$TMPDIR"/{raw,lib,results}
cp /pi/qing.yu-umw/home/shaoxian.li11-umw/DIANN/raw/SL14N15_250620_RAPL2_timestarv_proteinlevel_nobiotinvsbiotin/*.raw "$TMPDIR/raw/"
cp /pi/qing.yu-umw/home/shaoxian.li11-umw/DIANN/lib/Apex_libmiss1var3_withapex_new2.predicted.speclib "$TMPDIR/lib/"

RAW_ARGS=""
for f in "$TMPDIR"/raw/*.raw; do
  RAW_ARGS+=" --f $f"
done
echo "Will process RAW files:$RAW_ARGS"

cd "$TMPDIR"
singularity exec "$DIANNIMG" /diann-2.1.0/diann-linux \
  --verbose 2 \
  --threads 32 \
  $RAW_ARGS \
  --lib   lib/Apex_libmiss1var3_withapex_new2.predicted.speclib \
  --qvalue 0.01 \
  --peptidoforms \
  --export-quant \
  --matrices \
  --matrix-qvalue 0.01 \
  --matrix-spec-q 0.05 \
  --protein-inference \
  --pg-level 1 \
  --var-mods 2 \
  --var-mod "UniMod:35,15.994915,M" \
  --var-mod "UniMod:2127,331.1896,Y" \
  --fixed-mod "UniMod:4,57.021464,C" \
  --met-excision \
  --out results/SL14N15_experiment_report.parquet

echo ">>> Contents of results/:"
ls -lh results/

DEST=/pi/qing.yu-umw/home/shaoxian.li11-umw/DIANN/results/S14N15BvsNB
mkdir -p "$DEST"
rsync -av results/SL14N15_experiment_report.parquet "$DEST/"
rsync -av results/*_matrix.tsv                "$DEST/"

echo "Done! Results now in $DEST/"


无DBPmod：
#!/bin/bash
#BSUB -J DIANN_quant                      # Job name
#BSUB -q long                            # Queue: up to 8 h
#BSUB -n 32                              # 32 cores
#BSUB -R "rusage[mem=8192] span[hosts=1]" # 8 GB per core (256 GB total), single node
#BSUB -W 8:00                            # Wall-time limit: 8 h
#BSUB -u shaoxian.li11@umassmed.edu      # Email notifications
#BSUB -o diann_%J.out                    # stdout
#BSUB -e diann_%J.err                    # stderr

set -euo pipefail

export TMPDIR=${TMPDIR:-/tmp}
export DIANNIMG=/share/pkg/containers/diann/diann_2.1.0.sif
echo "Scratch TMPDIR: $TMPDIR"
echo "Using container: $DIANNIMG"

module load diann/2.1.0

mkdir -p "$TMPDIR"/{raw,lib,results}
cp /pi/qing.yu-umw/home/shaoxian.li11-umw/DIANN/raw/QC/20250708_8cells_aa000253_20250708135045.raw "$TMPDIR/raw/"
cp /pi/qing.yu-umw/home/shaoxian.li11-umw/DIANN/lib/lib1missc1varmod.speclib "$TMPDIR/lib/"

RAW_ARGS=""
for f in "$TMPDIR"/raw/*.raw; do
  RAW_ARGS+=" --f $f"
done
echo "Will process RAW files:$RAW_ARGS"

cd "$TMPDIR"
singularity exec "$DIANNIMG" /diann-2.1.0/diann-linux \
  --verbose 2 \
  --threads 32 \
  $RAW_ARGS \
  --lib   lib/lib1missc1varmod.speclib \
  --qvalue 0.01 \
  --peptidoforms \
  --export-quant \
  --matrices \
  --matrix-qvalue 0.01 \
  --matrix-spec-q 0.05 \
  --protein-inference \
  --pg-level 1 \
  --var-mods 2 \
  --var-mod "UniMod:35,15.994915,M" \
  --fixed-mod "UniMod:4,57.021464,C" \
  --met-excision \
  --out results/QC_8cells_DIA_report.parquet

echo ">>> Contents of results/:"
ls -lh results/

DEST=/pi/qing.yu-umw/home/shaoxian.li11-umw/DIANN/results/QC/250709_8cells
mkdir -p "$DEST"
rsync -av results/250709_8cells.parquet "$DEST/"
rsync -av results/*_matrix.tsv                "$DEST/"

echo "Done! Results now in $DEST/"

