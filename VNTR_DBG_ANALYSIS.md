# D1S80 ordered branch witnesses and 16-bp unit decomposition

This is an additive analysis of a frozen `vntrdbg` run. It does not modify the
existing GFA, candidate sequences, read assignments, or genotype status.

`scripts/vntrdbg_analysis.py` reads `input_paths.json` and the matching
point-graph `motifs.tsv`/`tokens.tsv`. For k=21 and k=33, it finds branch-edge
blocks from exact observed k+1-mers. An observed edge must have three distinct
molecules in one sample to define a supported branch; reference and candidate
edges are retained as explicit hypotheses. Adjacent, nonoverlapping branch
blocks define an ordered interval in the candidate's actual base sequence.

`joint_queries.tsv` gives the exact interval and candidate coordinates. For an
interval, `joint_molecules.tsv` records every contiguous exact occurrence in
one observed path. Multiple hits or paths from one `molecule_key` still count
once in `joint_support.tsv`. Local-only paths can witness a contained interval,
but never a whole candidate. A `FULL` query requires a two-flank-anchored OBS
whose entire sequence equals the candidate. Candidate paths never vote for
themselves, and a pair of separate local fragments is never joined.

This exact test is a **lower bound**, sensitive to sequencing errors. A missing
exact match is `uninformative`, not a contradiction. This version does not
estimate a support fraction or make genotype calls. `candidate_occurrences`
exposes repeated query sequences and therefore possible placement ambiguity.
The source sequence, per-molecule occurrence coordinates, and upstream source
hashes are retained for audit.

The unit analysis runs the existing fixed-phase dynamic-programming
decomposer with only usable, exactly 16-bp D1S80 base templates. It creates a
new layer beside the legacy annotations. `primitive_units.tsv` retains the
actual bases, template, half-open sequence coordinates, edit operations,
primary/alternative scheme, and overlapping legacy-token indexes.
`unit_summary.tsv` reports full 16-bp units, partial units, unknown bases,
length/16 and a conservative count status separately. A 48-bp legacy match
can become three 16-bp units, but only if the sequence-level decomposition
finds them; `legacy_full_tokens` remains unchanged. `motif_counts.tsv` groups
the primary unit layer without discarding its exact sequence rows.

Run on a compute node with `slurm/vntrdbg_analysis.sbatch`, setting `RUN_DIR`
to a new run directory containing a source snapshot and `BASE_RUN` to the
frozen D1S80 graph run. The batch job verifies source checksums, executes the
new unit tests, computes outputs, and independently validates molecule counts,
source slices, full-length scope and sequence reconstruction.
