#!/usr/bin/env bash
# 旧版 DMPK 基线作业：保留 1 kb 锚点和全 read NM 比例过滤，便于结果对照。
#SBATCH -J dmpk_extract
#SBATCH -p CPUqueue_fcp1
#SBATCH -q stzhu_1_qos
#SBATCH --account=stzhu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH -o logs/%x.%j.out
#SBATCH -e logs/%x.%j.err

set -euo pipefail

# 固定软件与数据路径，确保计算节点使用同一输入和环境。
PROJECT='/home/data/user/stzhu/projects/重复区域分析/experiment-trs-map'
BAM='/home/data/user/stzhu/测序data/ChineseQuartet.LCL5.GRCh38.HiFi.minimap2.bam'
SAMTOOLS='/home/data/user/stzhu/miniconda3/envs/samtools_env/bin/samtools'
PYTHON='/home/data/user/stzhu/miniconda3/bin/python3'
OUTPUT_DIR="${PROJECT}/results/DMPK_test"

mkdir -p "${OUTPUT_DIR}"

# 将作业身份、节点和时间写入标准输出，便于排查调度与运行问题。
echo "job_id=${SLURM_JOB_ID}"
echo "node=$(hostname)"
echo "started=$(date --iso-8601=seconds)"

# 旧版参数用于复现实验基线，不代表当前推荐配置。
srun "${PYTHON}" "${PROJECT}/scripts/extract_candidate_reads.py" \
  --bam "${BAM}" \
  --chrom chr19 \
  --start 45770205 \
  --end 45770264 \
  --flank 1000 \
  --min-mapq 20 \
  --max-nm-fraction 0.01 \
  --samtools "${SAMTOOLS}" \
  --output-dir "${OUTPUT_DIR}" \
  --prefix DMPK

echo "finished=$(date --iso-8601=seconds)"
