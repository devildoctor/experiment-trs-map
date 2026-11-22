#!/usr/bin/env bash
# Submit ATXN2 pilot, then nine loci from the same immutable snapshot.
set -euo pipefail
run_dir="$(realpath "${1:?Usage: submit.sh RUN_DIR}")"
cd "$run_dir"
sha256sum --quiet -c SHA256SUMS
test ! -e jobs.tsv || { echo 'jobs.tsv exists; use a fresh run directory' >&2; exit 1; }
mkdir -p logs results
submitted=()
rollback() {
  if ((${#submitted[@]})); then scancel "${submitted[@]}"; fi
}
trap rollback ERR
printf 'locus\tjob_id\tdependency\n' > jobs.tsv
pilot=''
for locus in ATXN2 D1S80 PER3 IL1RN ACAN PLIN4 WDR7 CEL ABCA7 MUC1; do
  options=()
  dependency='none'
  if [[ -n "$pilot" ]]; then
    dependency="afterok:$pilot"
    options+=(--dependency="$dependency" --kill-on-invalid-dep=yes)
  fi
  id="$(sbatch --parsable --job-name="vntr2-$locus" \
    --output="$run_dir/logs/$locus.%j.out" --error="$run_dir/logs/$locus.%j.err" \
    --export="ALL,RUN_DIR=$run_dir,LOCUS=$locus" \
    "${options[@]}" "$run_dir/slurm/vntrgraph.sbatch")"
  id="${id%%;*}"
  [[ "$id" =~ ^[0-9]+$ ]]
  submitted+=("$id")
  printf '%s\t%s\t%s\n' "$locus" "$id" "$dependency" >> jobs.tsv
  if [[ -z "$pilot" ]]; then pilot="$id"; fi
done
trap - ERR
cat jobs.tsv
