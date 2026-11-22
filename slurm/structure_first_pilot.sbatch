#!/usr/bin/env bash
#SBATCH -J vntr-structure-v2
#SBATCH -p CPUqueue_fcp1
#SBATCH -q stzhu_1_qos
#SBATCH --account=stzhu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:40:00
set -euo pipefail
PROJECT='/home/data/user/stzhu/重复区域分析/vntr-point-graph'
RUN_DIR="${RUN_DIR:?Set RUN_DIR to the frozen source snapshot}"
cd "$RUN_DIR"
sha256sum --quiet -c SHA256SUMS
export PYTHONUNBUFFERED=1 OMP_NUM_THREADS=1
export PYTHONPATH="$PROJECT/.deps/vntr-anchors${PYTHONPATH:+:$PYTHONPATH}"
PYTHON=/home/data/user/stzhu/miniconda3/bin/python3
SAMTOOLS=/home/data/user/stzhu/miniconda3/envs/samtools_env/bin/samtools
"$PYTHON" -m unittest discover -s tests -p test_pathgraph.py -v
"$PYTHON" -m unittest discover -s tests -p test_partition.py -v
"$PYTHON" -m unittest discover -s tests -p test_pathstructure.py -v
for LOCUS in ATXN2 D1S80; do
  "$PYTHON" scripts/vntrpathgraph_cli.py \
    --bam-dir "$PROJECT/bamdata_by_locus" \
    --reference "$PROJECT/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa" \
    --catalog catalogs/grch38_structure_summary_atxn2_d1s80_abca7_v1.bed \
    --out results --locus "$LOCUS" --samtools "$SAMTOOLS" --preview-paths 10 \
    --graph-mode structure --locus-model "catalogs/$LOCUS.structure.v1.json" \
    --structure-policy "catalogs/$LOCUS.structure-policy.v1.json" --max-loss-bp 2
  "$PYTHON" tests/verify_structure_first.py "results/$LOCUS" \
    "$PROJECT/runs-result/20261105_path_structure_061916/results/$LOCUS"
done
sha256sum --quiet -c SHA256SUMS
printf '%s\n' 'STRUCTURE_FIRST_PILOT_PASS'
