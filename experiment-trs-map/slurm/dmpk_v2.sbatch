#!/usr/bin/env bash
# 当前 DMPK v2 作业：从目录读取位点，采用可配置锚点与 motif 证据筛选。
#SBATCH -J dmpk_select_v2
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
OUTPUT_DIR="${PROJECT}/results/DMPK_v2"

mkdir -p "${OUTPUT_DIR}"

# 将作业身份、节点和时间写入标准输出，便于排查调度与运行问题。
echo "job_id=${SLURM_JOB_ID}"
echo "node=$(hostname)"
echo "started=$(date --iso-8601=seconds)"

# 通过 srun 继承 Slurm 分配的 CPU 与内存资源，线程数与 cpus-per-task 一致。
srun "${PYTHON}" "${PROJECT}/scripts/extract_candidate_reads.py" \
  --bam "${BAM}" \
  --repeats "${PROJECT}/catalogs/example.repeats.bed" \
  --repeat-id DMPK \
  --anchor-length 100 \
  --fetch-padding 1000 \
  --min-mapq 20 \
  --min-motif-fraction 0.70 \
  --min-repeat-bases 12 \
  --samtools "${SAMTOOLS}" \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --output-dir "${OUTPUT_DIR}" \
  --prefix DMPK_v2

echo "finished=$(date --iso-8601=seconds)"
