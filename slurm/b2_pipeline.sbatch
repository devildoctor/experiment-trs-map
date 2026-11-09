#!/usr/bin/env bash
#SBATCH -J vntr-b2
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
RUN="${RUN:?}"
cd "${SOURCE:?}"
export PYTHONUNBUFFERED=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONPATH="$PROJECT/runs-result/20261108_flank_v2/deps:$PROJECT/.deps/vntr-anchors"
export NUMBA_CACHE_DIR="$RUN/numba-cache"
SAMTOOLS=/home/data/user/stzhu/miniconda3/envs/samtools_env/bin/samtools
PYTHON=/home/data/user/stzhu/miniconda3/bin/python3
"$PYTHON" -m unittest discover -s tests -v
for check in integration_flanks integration_all_reads integration_pathgraph integration_structure integration_review_namespaces integration_synthetic; do
    "$PYTHON" "tests/$check.py" "$SAMTOOLS"
done
echo B2_CODE_TESTS_PASS
"$PYTHON" scripts/run_b2_pipeline.py --project "$PROJECT" --out "$RUN/results" \
    --previous-b2 "$PROJECT/runs-result/20261108_flank_v2/final_replay/comparison/B2" --samtools "$SAMTOOLS"
echo B2_ALL_DONE
