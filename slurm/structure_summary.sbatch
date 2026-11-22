#!/usr/bin/env bash
#SBATCH -J vntr-summary
#SBATCH -p CPUqueue_fcp1
#SBATCH -q stzhu_1_qos
#SBATCH --account=stzhu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
set -euo pipefail
PROJECT='/home/data/user/stzhu/重复区域分析/vntr-point-graph'
RUN_DIR="${RUN_DIR:?}"
LOCUS="${LOCUS:?}"
case "$LOCUS" in ATXN2) STRAND=-;; ABCA7|D1S80) STRAND=unknown;; *) exit 2;; esac
cd "$RUN_DIR"
export PYTHONUNBUFFERED=1 OMP_NUM_THREADS=1
PYTHON=/home/data/user/stzhu/miniconda3/bin/python3
SAMTOOLS=/home/data/user/stzhu/miniconda3/envs/samtools_env/bin/samtools
"$PYTHON" scripts/vntrstructure.py --bam-dir "$PROJECT/bamdata_by_locus" \
  --reference "$PROJECT/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa" \
  --catalog "$RUN_DIR/catalogs/grch38_structure_summary_atxn2_d1s80_abca7_v1.bed" \
  --out "$RUN_DIR/results" --locus "$LOCUS" --samtools "$SAMTOOLS" \
  --workers "${SLURM_CPUS_PER_TASK:-4}" --gene-strand "$STRAND" \
  --variant-support 3 --max-variants 256 --min-sample-fraction 0.10 --discovery-rounds 3
"$PYTHON" tests/verify_structure_results.py "$RUN_DIR/results/$LOCUS"
"$PYTHON" tests/verify_summary_results.py "$RUN_DIR/results/$LOCUS"
