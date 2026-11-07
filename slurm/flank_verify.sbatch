#!/usr/bin/env bash
#SBATCH -J vntr-flank-verify
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
cd "${SOURCE:?}"
export PYTHONUNBUFFERED=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONPATH="$RUN/deps:$PROJECT/.deps/vntr-anchors"
export NUMBA_CACHE_DIR="$RUN/numba-cache-verify"
export SAMTOOLS=/home/data/user/stzhu/miniconda3/envs/samtools_env/bin/samtools
PYTHON=/home/data/user/stzhu/miniconda3/bin/python3
"$PYTHON" -m unittest discover -s tests -v
for check in integration_flanks integration_all_reads integration_pathgraph integration_structure integration_review_namespaces integration_synthetic; do
    "$PYTHON" "tests/$check.py" "$SAMTOOLS"
done
echo FINAL_CODE_TESTS_PASS
for LOCUS in ATXN2 D1S80; do
    "$PYTHON" scripts/run_flank_comparison.py --project "$PROJECT" --run "$RUN/final_replay" --locus "$LOCUS"
    "$PYTHON" scripts/audit_flank_results.py --project "$PROJECT" --run "$RUN/final_replay" --locus "$LOCUS"
done
echo FINAL_REAL_REPLAY_PASS
