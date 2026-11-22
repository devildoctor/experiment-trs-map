#!/usr/bin/env bash
# Submit exactly the three user-requested loci from an immutable source snapshot.
set -euo pipefail
run_dir="$(realpath "${1:?Usage: submit_reads.sh RUN_DIR}")"
cd "$run_dir"
sha256sum --quiet -c SHA256SUMS
test ! -e reads_jobs.tsv || { echo 'Already submitted; use a new snapshot' >&2; exit 1; }
mkdir -p logs results
printf 'locus\tjob_id\n' > reads_jobs.tsv
for locus in IL1RN PER3 ABCA7; do
  id="$(sbatch --parsable --job-name="vntr3-$locus" \
    --output="$run_dir/logs/$locus.%j.out" --error="$run_dir/logs/$locus.%j.err" \
    --export="ALL,RUN_DIR=$run_dir,LOCUS=$locus" "$run_dir/slurm/vntrreads.sbatch")"
  id="${id%%;*}"
  [[ "$id" =~ ^[0-9]+$ ]]
  printf '%s\t%s\n' "$locus" "$id" >> reads_jobs.tsv
done
cat reads_jobs.tsv
