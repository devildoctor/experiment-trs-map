#!/usr/bin/env bash
#SBATCH -J vntrstructure
#SBATCH -p CPUqueue_fcp1
#SBATCH -q stzhu_1_qos
#SBATCH --account=stzhu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=06:00:00
set -euo pipefail
PROJECT='/home/data/user/stzhu/重复区域分析/vntr-point-graph'
RUN_DIR="${RUN_DIR:?Set RUN_DIR}"
LOCUS="${LOCUS:?Set LOCUS}"
case "$LOCUS" in IL1RN|PER3|ABCA7) ;; *) exit 2;; esac
cd "$RUN_DIR"
sha256sum --quiet -c SHA256SUMS
export PYTHONUNBUFFERED=1 OMP_NUM_THREADS=1
PYTHON=/home/data/user/stzhu/miniconda3/bin/python3
SAMTOOLS=/home/data/user/stzhu/miniconda3/envs/samtools_env/bin/samtools
"$PYTHON" -m unittest discover -s tests -v
"$PYTHON" tests/integration_structure.py "$SAMTOOLS"
"$PYTHON" scripts/vntrgraph.py --structure-view \
  --bam-dir "$PROJECT/bamdata_by_locus" \
  --reference "$PROJECT/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa" \
  --catalog "$RUN_DIR/catalogs/grch38_structure_seeds_il1rn_per3_abca7_v1.bed" \
  --out "$RUN_DIR/results" --locus "$LOCUS" --samtools "$SAMTOOLS" \
  --workers "${SLURM_CPUS_PER_TASK:-2}" --variant-support 3 --max-variants 12
"$PYTHON" tests/verify_structure_results.py "$RUN_DIR/results/$LOCUS"
