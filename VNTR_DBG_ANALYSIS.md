# D1S80 ordered branch witnesses and 16-bp unit decomposition

This is an additive analysis of a frozen `vntrdbg` run. It does not modify the
existing GFA, candidate sequences, read assignments, or genotype status.

`scripts/vntrdbg_analysis.py` reads `input_paths.json` and the matching
point-graph `motifs.tsv`/`tokens.tsv`. It verifies D1S80 identity, each sequence
hash and the point-graph hashes against frozen `input_summary.json` before
writing results. For k=21 and k=33, it finds branch-edge
blocks from exact observed k+1-mers. An observed edge must have `--branch-min-support` distinct
molecules in one sample (default 3) to define a supported branch; reference and candidate
edges are retained as explicit hypotheses. Adjacent, nonoverlapping branch
blocks define an ordered interval in the candidate's actual base sequence.

The branch threshold belongs to this analysis stage, is recorded in `manifest.json`, and is not implicitly inherited from graph filtering. Set it explicitly when comparing analyses with different thresholds. Candidates without same-sample OBS produce zero tested/support/conflict counts, `support_fraction=NA`, and header-only molecule evidence tables rather than an export error.

`joint_queries.tsv` gives the exact interval and candidate coordinates. For an
interval, `joint_molecules.tsv` records every contiguous exact occurrence in
one observed path. Multiple hits or paths from one `molecule_key` still count
once in `joint_support.tsv`. Local-only paths can witness a contained interval,
but never a whole candidate. A `FULL` query requires a two-flank-anchored OBS
whose entire sequence equals the candidate. Candidate paths never vote for
themselves, and a pair of separate local fragments is never joined.

The raw `exact_support_molecules` count is a **lower bound** on contiguous
sequence witnesses, sensitive to sequencing errors and repeat placement.
`branch_positions_json` records actual branching edge positions; merged
nonbranch gaps do not inflate `branch_edge_count`.

`events.tsv` defines query intervals and explicit sequence states from
same-sample candidates. These are graph intervals, not normalized biological
variants. An adjacent-pair query references its two ordered block IDs. Its
complete intervening sequence is tested, not just two endpoint matches.
`molecule_events.tsv` records every individual OBS assessment and coordinates.
`joint_molecules.tsv` additionally aggregates assessments per molecule:

- FULL requires two-flank anchoring and an exact whole-sequence match to the
  target (support) or a different same-sample candidate state (conflict).
  Identical candidate sequences are equivalent states, not competitors.
- A local interval requires unique, nonoverlapping, correctly ordered k-base
  left and right contexts in both the target candidate and OBS. The bounded
  OBS sequence must equal a defined target or alternative state. An OBS
  identical to the entire anchored target also supports its known intervals.
- Repeated contexts, missing coverage, unmodeled mismatches and sequencing
  errors remain uninformative. Exact matches without unique placement remain
  in the raw witness columns, separate from resolved `support_molecules`.
- A molecule with both supporting and conflicting resolved records is
  uninformative (`mixed_states_in_one_molecule`), never two votes.

`informative_molecules = support_molecules + conflict_molecules`.
`support_fraction = support_molecules / informative_molecules`, or `NA` for a
zero denominator. This is conditional on exact modeled states; it is **not**
an allele frequency, genotype probability, or complete error-aware estimate.
Candidate assignment is retained as metadata, never used to select witnesses.

The unit analysis runs the existing fixed-phase dynamic-programming
decomposer with only usable D1S80 templates whose exact primitive period is
16 bp. It processes every REF, OBS and CAND path, including OBS shorter than k.
It creates a
new layer beside the legacy annotations. `primitive_units.tsv` retains the
actual bases, template, half-open sequence coordinates, edit operations,
primary/alternative scheme, and overlapping legacy-token indexes.
Every unit table carries `path_id`, kind, coverage, molecule identity and
sequence hash; candidate IDs and legacy annotations apply only to CAND.
`unit_summary.tsv` reports full 16-bp template units, partial units, unknown bases,
length/16 and a conservative count status separately. A 48-bp legacy match
can become three 16-bp units, but only if the sequence-level decomposition
finds them; `legacy_full_tokens` remains unchanged. `motif_counts.tsv` groups
the primary unit layer without discarding its exact sequence rows. A full
template unit can consume 15/17 observed bases when the allowed edit is an
indel; `observed_unit_bp` and the edit CIGAR preserve this difference.

Local OBS always have `count_scope=observed_fragment` and `primitive_cn=NA`.
The primary and one second-best annotation are retained with scores, score
gap, fixed phase policy and `alternative_is_tied`. A second-best annotation
is not necessarily tied. Any tied annotation keeps CN unresolved even if
the two displayed counts agree, because the algorithm does not enumerate
all tied solutions. Unknown/partial sequence also prevents a resolved CN.
These counts describe the supplied sequence and do not establish a genotype.

Run on a compute node with `slurm/vntrdbg_analysis.sbatch`, setting `RUN_DIR`
to a new run directory containing a source snapshot and `BASE_RUN` to the
frozen D1S80 graph run. The batch job verifies source checksums, executes the
new unit tests, computes outputs, and independently validates molecule counts,
source slices, full-length scope, graph branch positions, edit operations,
all-path coverage and sequence reconstruction. The original graph and STR
modules remain unchanged. Existing results are not overwritten.
