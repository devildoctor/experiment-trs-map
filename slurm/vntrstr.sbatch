#!/usr/bin/env bash
#SBATCH -J str-ATXN2
#SBATCH -p CPUqueue_fcp1
#SBATCH -q stzhu_1_qos
#SBATCH --account=stzhu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH -o str-ATXN2.%j.out
#SBATCH -e str-ATXN2.%j.err
set -euo pipefail
PROJECT='/home/data/user/stzhu/重复区域分析/vntr-point-graph'
RUN_DIR="${RUN_DIR:?Set RUN_DIR to a new versioned source snapshot}"
cd "$RUN_DIR"
sha256sum --quiet -c SHA256SUMS
export PYTHONUNBUFFERED=1
export PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="$PROJECT/.deps/vntr-anchors${PYTHONPATH:+:$PYTHONPATH}"
PY=/home/data/user/stzhu/miniconda3/bin/python3
ARGS=(--bam-dir "$PROJECT/bamdata_by_locus"
      --reference "$PROJECT/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa"
      --catalog "$RUN_DIR/catalogs/grch38_str_expression_atxn2_v1.bed"
      --locus ATXN2
      --samtools /home/data/user/stzhu/miniconda3/envs/samtools_env/bin/samtools)
"$PY" scripts/vntrgraph.py run "${ARGS[@]}" --out "$RUN_DIR/results/pointgraph"
"$PY" scripts/vntrreads.py "${ARGS[@]}" --out "$RUN_DIR/results/readview" --workers 2
"$PY" scripts/vntrstructure.py "${ARGS[@]}" --out "$RUN_DIR/results/structure" --workers 2 --gene-strand '-'
