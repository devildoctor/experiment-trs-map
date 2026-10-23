#!/usr/bin/env bash
#SBATCH -J d1s80-terminal-units
#SBATCH -p CPUqueue_fcp1
#SBATCH -q stzhu_1_qos
#SBATCH --account=stzhu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=01:00:00
set -euo pipefail
RUN_DIR="${RUN_DIR:?Set a new versioned run directory}"
BASE_RUN="${BASE_RUN:?Set the frozen D1S80 graph run}"
OLD_RUN="${OLD_RUN:?Set the original joint-unit analysis run}"
REFERENCE="${REFERENCE:?Set the unchanged original reference}"
PY=/home/data/user/stzhu/miniconda3/bin/python3
export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
export PYTHONPATH="$RUN_DIR/source/scripts:$RUN_DIR/source/.deps/vntr-anchors"
cd "$RUN_DIR/source"
sha256sum --quiet -c SHA256SUMS
mkdir -p "$RUN_DIR/acceptance" "$RUN_DIR/results"
find "$BASE_RUN/results" -type f -print0 | sort -z | xargs -0 sha256sum > "$RUN_DIR/acceptance/frozen-before.sha256"
"$PY" -m unittest discover -s tests -p 'test_*.py' -v > "$RUN_DIR/acceptance/unit-tests.log" 2>&1
"$PY" scripts/vntrdbg_analysis.py --dbg-dir "$BASE_RUN/results/dbg/D1S80" \
  --pointgraph-dir "$BASE_RUN/results/pointgraph/D1S80" \
  --out "$RUN_DIR/results/baseline_recheck" --k 21 33 > "$RUN_DIR/acceptance/baseline.log" 2>&1
for table in joint_queries events joint_support joint_molecules molecule_events primitive_units unit_summary motif_counts; do
  cmp "$RUN_DIR/results/baseline_recheck/$table.tsv" "$OLD_RUN/results/D1S80/$table.tsv"
done
"$PY" scripts/vntrdbg_analysis.py --dbg-dir "$BASE_RUN/results/dbg/D1S80" \
  --pointgraph-dir "$BASE_RUN/results/pointgraph/D1S80" \
  --out "$RUN_DIR/results/D1S80" --k 21 33 --terminal-units --terminal-min-support 3 \
  --catalog "$RUN_DIR/source/catalogs/grch38_pointgraph_10_loci_v1.bed" \
  --reference "$REFERENCE" > "$RUN_DIR/acceptance/analysis.log" 2>&1
"$PY" tests/validate_vntrdbg_analysis.py "$RUN_DIR/results/D1S80" \
  "$BASE_RUN/results/dbg/D1S80/input_paths.json" > "$RUN_DIR/acceptance/base-validation.log" 2>&1
"$PY" tests/validate_vntrterminal_run.py "$RUN_DIR/results/D1S80" \
  "$OLD_RUN/results/D1S80" "$BASE_RUN/results/dbg/D1S80/input_paths.json" \
  "$BASE_RUN/results/pointgraph/D1S80" > "$RUN_DIR/acceptance/terminal-validation.log" 2>&1
sha256sum --quiet -c "$RUN_DIR/acceptance/frozen-before.sha256"
sha256sum --quiet -c SHA256SUMS
find "$RUN_DIR/results" -type f -print0 | sort -z | xargs -0 sha256sum > "$RUN_DIR/acceptance/results.sha256"
echo D1S80_TERMINAL_UNITS_PASS
