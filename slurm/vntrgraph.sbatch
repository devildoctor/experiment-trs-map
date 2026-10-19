#!/usr/bin/env bash
#SBATCH -J vntrgraph
#SBATCH -p CPUqueue_fcp1
#SBATCH -q stzhu_1_qos
#SBATCH --account=stzhu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH -o vntrgraph.%j.out
#SBATCH -e vntrgraph.%j.err

set -euo pipefail
PROJECT='/home/data/user/stzhu/重复区域分析/vntr-point-graph'
RUN_DIR="${RUN_DIR:?Set RUN_DIR to a versioned source snapshot}"
cd "$RUN_DIR"
sha256sum --quiet -c SHA256SUMS
ARGS=()
if [[ -n "${LOCUS:-}" ]]; then
  ARGS+=(--locus "$LOCUS")
fi
if [[ -f "$RUN_DIR/sources.tsv" ]]; then
  ARGS+=(--sources "$RUN_DIR/sources.tsv")
fi
if [[ -n "${REFERENCE_EQUIVALENCE:-}" ]]; then
  ARGS+=(--reference-equivalence "$REFERENCE_EQUIVALENCE")
fi
export PYTHONUNBUFFERED=1
export PYTHONPATH="$PROJECT/.deps/vntr-anchors${PYTHONPATH:+:$PYTHONPATH}"
/home/data/user/stzhu/miniconda3/bin/python3 "$RUN_DIR/scripts/vntrgraph.py" run \
  --bam-dir "$PROJECT/bamdata_by_locus" \
  --reference "$PROJECT/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa" \
  --catalog "$RUN_DIR/catalogs/grch38_pointgraph_10_loci_v1.bed" \
  --out "$RUN_DIR/results" \
  --samtools /home/data/user/stzhu/miniconda3/envs/samtools_env/bin/samtools \
  "${ARGS[@]}"
