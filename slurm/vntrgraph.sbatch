#!/usr/bin/env bash
#SBATCH -J vntrgraph
#SBATCH -p CPUqueue_fcp1
#SBATCH -q stzhu_1_qos
#SBATCH --account=stzhu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH -o vntrgraph.%j.out
#SBATCH -e vntrgraph.%j.err

set -euo pipefail
PROJECT='/home/data/user/stzhu/重复区域分析/vntr-point-graph'
cd "$PROJECT"
ARGS=()
if [[ -n "${LOCUS:-}" ]]; then
  ARGS+=(--locus "$LOCUS")
fi
/home/data/user/stzhu/miniconda3/bin/python3 scripts/vntrgraph.py run \
  --bam-dir bamdata \
  --reference reference/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  --catalog catalogs/VNTR_candidates.bed \
  --out results/vntrgraph \
  --samtools /home/data/user/stzhu/miniconda3/envs/samtools_env/bin/samtools \
  "${ARGS[@]}"
