#!/usr/bin/env bash
#SBATCH -J str-gfa-ATXN2
#SBATCH -p CPUqueue_fcp1
#SBATCH -q stzhu_1_qos
#SBATCH --account=stzhu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH -o str-gfa-ATXN2.%j.out
#SBATCH -e str-gfa-ATXN2.%j.err
set -euo pipefail
PROJECT='/home/data/user/stzhu/重复区域分析/vntr-point-graph'
RUN_DIR="${RUN_DIR:?Set RUN_DIR to an immutable source snapshot}"
BASELINE_DIR="${BASELINE_DIR:?Set BASELINE_DIR to the pre-change source snapshot}"
BEFORE_RESULT="${BEFORE_RESULT:?Set BEFORE_RESULT to the previous ATXN2 pointgraph directory}"
OUT_DIR="${OUT_DIR:-$RUN_DIR/results}"
cd "$RUN_DIR"
sha256sum --quiet -c SHA256SUMS
export PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="$PROJECT/.deps/vntr-anchors${PYTHONPATH:+:$PYTHONPATH}"
PY=/home/data/user/stzhu/miniconda3/bin/python3
mkdir -p "$OUT_DIR"
"$PY" -m unittest discover -s tests -v > "$OUT_DIR/unit-tests.log" 2>&1
"$PY" tests/integration_synthetic.py > "$OUT_DIR/integration.log" 2>&1
"$PY" tests/regression_vntr_graph.py --baseline "$BASELINE_DIR" > "$OUT_DIR/vntr-regression.log" 2>&1
"$PY" scripts/vntrgraph.py run \
    --bam-dir "$PROJECT/bamdata_by_locus" \
    --reference "$PROJECT/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa" \
    --catalog "$RUN_DIR/catalogs/grch38_str_expression_atxn2_v1.bed" \
    --locus ATXN2 \
    --samtools /home/data/user/stzhu/miniconda3/envs/samtools_env/bin/samtools \
    --out "$OUT_DIR/pointgraph"
"$PY" tests/validate_str_graph.py --result "$OUT_DIR/pointgraph/ATXN2" \
    --before "$BEFORE_RESULT" > "$OUT_DIR/graph-validation.log" 2>&1
echo STR_GFA_ACCEPTANCE_PASS
