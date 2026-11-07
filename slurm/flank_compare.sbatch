#!/usr/bin/env bash
#SBATCH -J vntr-flank-compare
#SBATCH -p CPUqueue_fcp1
#SBATCH -q stzhu_1_qos
#SBATCH --account=stzhu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=03:00:00
set -euo pipefail
PROJECT='/home/data/user/stzhu/重复区域分析/vntr-point-graph'
RUN="$PROJECT/runs-result/20261108_flank_v2"
SOURCE="${SOURCE:?}"
LOCUS="${LOCUS:?}"
cd "$SOURCE"
export PYTHONUNBUFFERED=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONPATH="$RUN/deps:$PROJECT/.deps/vntr-anchors"
export NUMBA_CACHE_DIR="$RUN/numba-cache-$LOCUS"
PYTHON=/home/data/user/stzhu/miniconda3/bin/python3
"$PYTHON" -m unittest discover -s tests -v
"$PYTHON" scripts/run_flank_comparison.py --project "$PROJECT" --run "$RUN" --locus "$LOCUS"
