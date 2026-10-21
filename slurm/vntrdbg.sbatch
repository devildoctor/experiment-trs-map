#!/usr/bin/env bash
#SBATCH -J vntr-dbg-D1S80
#SBATCH -p CPUqueue_fcp1
#SBATCH -q stzhu_1_qos
#SBATCH --account=stzhu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
set -euo pipefail
RUN_DIR="${RUN_DIR:?Set RUN_DIR to the versioned run directory}"
PROJECT='/home/data/user/stzhu/重复区域分析/vntr-point-graph'
PY=/home/data/user/stzhu/miniconda3/bin/python3
SAMTOOLS=/home/data/user/stzhu/miniconda3/envs/samtools_env/bin/samtools
export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
export PYTHONPATH="$PROJECT/.deps/vntr-anchors"
cd "$RUN_DIR/source"
sha256sum --quiet -c SHA256SUMS
sha256sum --quiet -c "$RUN_DIR/results/preparation/bam-before.sha256"
mkdir -p "$RUN_DIR/results/acceptance"
"$PY" -m unittest discover -s tests -v > "$RUN_DIR/results/acceptance/unit-tests.log" 2>&1
"$PY" tests/integration_synthetic.py > "$RUN_DIR/results/acceptance/integration.log" 2>&1
"$PY" scripts/vntrgraph.py --dbg-view \
 --evidence-dir "$RUN_DIR/results/pointgraph/D1S80" \
 --bam-dir "$PROJECT/bamdata_by_locus" \
 --reference "$PROJECT/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa" \
 --catalog catalogs/grch38_pointgraph_10_loci_v1.bed \
 --models catalogs/grch38_graph_models_v1.tsv --locus D1S80 \
 --k 17 21 31 33 49 65 --samtools "$SAMTOOLS" \
 --out "$RUN_DIR/results/dbg/D1S80" > "$RUN_DIR/results/acceptance/build.log" 2>&1
"$PY" tests/validate_dbg.py "$RUN_DIR/results/dbg/D1S80" > "$RUN_DIR/results/acceptance/dbg-validation.log" 2>&1
"$PY" scripts/vntrgraph.py run --bam-dir "$PROJECT/bamdata_by_locus" \
 --reference "$PROJECT/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa" \
 --catalog catalogs/grch38_str_expression_atxn2_v1.bed --locus ATXN2 \
 --samtools "$SAMTOOLS" --out "$RUN_DIR/results/str-regression" > "$RUN_DIR/results/acceptance/str-regression.log" 2>&1
"$PY" - "$RUN_DIR" "$PROJECT" <<'PY'
from pathlib import Path
import hashlib,json,sys
r,p=map(Path,sys.argv[1:])
old=p/'runs-result/20261020T082940Z_str-gfa_ATXN2/results/pointgraph/ATXN2'
new=r/'results/str-regression/ATXN2'
names=['candidates.tsv','assignments.tsv','genotypes.tsv','str_paths.json','str_paths.tsv','tokens.tsv','motifs.tsv','status.json','molecules.tsv','ATXN2.gfa','ATXN2.json','ATXN2.svg','ATXN2.paths.svg']
for name in names:
    assert (old/name).read_bytes()==(new/name).read_bytes(),name
for file in (old/'paths-svg').glob('*.svg'):
    assert file.read_bytes()==(new/'paths-svg'/file.name).read_bytes(),file.name
result=dict(status='PASS',byte_identical=names,individual_path_svgs=len(list((new/'paths-svg').glob('*.svg'))))
(r/'results/acceptance/str-regression.json').write_text(json.dumps(result,indent=2)+'\n')
PY
sha256sum --quiet -c "$RUN_DIR/results/preparation/bam-before.sha256"
echo VNTR_DBG_ACCEPTANCE_PASS
