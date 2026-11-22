#!/usr/bin/env bash
#SBATCH -J vntr-partition
#SBATCH -p CPUqueue_fcp1
#SBATCH -q stzhu_1_qos
#SBATCH --account=stzhu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
set -euo pipefail
PROJECT='/home/data/user/stzhu/重复区域分析/vntr-point-graph'
RUN_DIR="${RUN_DIR:?Set RUN_DIR to a frozen source snapshot}"
cd "$RUN_DIR"
sha256sum --quiet -c SHA256SUMS
export PYTHONUNBUFFERED=1 OMP_NUM_THREADS=1
export PYTHONPATH="$PROJECT/.deps/vntr-anchors${PYTHONPATH:+:$PYTHONPATH}"
PYTHON=/home/data/user/stzhu/miniconda3/bin/python3
SAMTOOLS=/home/data/user/stzhu/miniconda3/envs/samtools_env/bin/samtools
"$PYTHON" -m unittest discover -s tests -v
"$PYTHON" tests/integration_structure.py "$SAMTOOLS"
for LOCUS in ATXN2 ABCA7; do
  "$PYTHON" scripts/vntrstructure.py \
    --bam-dir "$PROJECT/bamdata_by_locus" \
    --reference "$PROJECT/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa" \
    --catalog catalogs/grch38_partition_pilot_v1.bed \
    --locus-model "catalogs/$LOCUS.partition.v1.json" --core-output \
    --out results --locus "$LOCUS" --samtools "$SAMTOOLS" \
    --workers "${SLURM_CPUS_PER_TASK:-2}"
  "$PYTHON" tests/verify_partition_results.py "results/$LOCUS"
done
sha256sum --quiet -c SHA256SUMS
