#!/usr/bin/env bash
#SBATCH -J d1s80-joint-units
#SBATCH -p CPUqueue_fcp1
#SBATCH -q stzhu_1_qos
#SBATCH --account=stzhu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=01:00:00
set -euo pipefail

RUN_DIR="${RUN_DIR:?Set RUN_DIR to a new, versioned run directory}"
BASE_RUN="${BASE_RUN:?Set BASE_RUN to the frozen D1S80 graph run}"
PY=/home/data/user/stzhu/miniconda3/bin/python3
export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
export PYTHONPATH="$RUN_DIR/source/scripts:$RUN_DIR/source/.deps/vntr-anchors"

cd "$RUN_DIR/source"
sha256sum --quiet -c SHA256SUMS
mkdir -p "$RUN_DIR/acceptance" "$RUN_DIR/results"
"$PY" -m unittest discover -s tests -p 'test_vntrdbg_analysis.py' -v \
  > "$RUN_DIR/acceptance/unit-tests.log" 2>&1
"$PY" scripts/vntrdbg_analysis.py \
  --dbg-dir "$BASE_RUN/results/dbg/D1S80" \
  --pointgraph-dir "$BASE_RUN/results/pointgraph/D1S80" \
  --out "$RUN_DIR/results/D1S80" --k 21 33 \
  > "$RUN_DIR/acceptance/analysis.log" 2>&1
"$PY" tests/validate_vntrdbg_analysis.py "$RUN_DIR/results/D1S80" \
  "$BASE_RUN/results/dbg/D1S80/input_paths.json" \
  > "$RUN_DIR/acceptance/validation.log" 2>&1
echo D1S80_JOINT_UNITS_PASS
