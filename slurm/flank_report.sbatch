#!/usr/bin/env bash
#SBATCH -J vntr-flank-report
#SBATCH -p CPUqueue_fcp1
#SBATCH -q stzhu_1_qos
#SBATCH --account=stzhu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
set -euo pipefail
PROJECT='/home/data/user/stzhu/重复区域分析/vntr-point-graph'
RUN="$PROJECT/runs-result/20261108_flank_v2"
export PYTHONPATH="$RUN/source_trial03/scripts:$RUN/deps:$PROJECT/.deps/vntr-anchors"
export NUMBA_CACHE_DIR="$RUN/numba-cache-report"
/home/data/user/stzhu/miniconda3/bin/python3 "$RUN/report_flank_evidence.py" --project "$PROJECT" --run "$RUN"
