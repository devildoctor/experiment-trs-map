# D1S80 local A-E display, 2026-11-10

This opt-in entry point implements the user-approved rules without changing flank
admission, B2 cuts, original P IDs, molecule provenance, or the frozen motif catalog.
Run `scripts/run_local_display.py --input SOURCE_MEMBERS --model catalogs/D1S80.structure.v1.json --out OUTPUT`.
The old v1/v2 entry points and result directories remain available.
The regular `vntrpathgraph_cli.py --graph-mode structure` also accepts the v3
policy with `--loss-budget` in [0, 0.10] and `--preview-paths 0`.

## Changed scope

- `vntrstructure.segment` optionally exports every primary-optimal candidate edge
  and mandatory boundary. Its default output and scoring are unchanged.
- `vntrpathstructure.annotate_paths` forwards that opt-in export flag.
- `vntrlocal.py` implements local comparison, per-visit variants, count parameters,
  shared unresolved fragments, regional interpretation selection, and local edits.
- `vntrlocal_export.py` produces fixed-string original/represented GFA, all-P TXT,
  a variable-node display graph, exact visit vectors, and separate support tables.
- The new runner, v3 policy, and tests are scoped to this display implementation.
- The existing CLI/export dispatcher explicitly routes schema 3 to these functions.

## Approved behavior and implementation limits

A: Same stable type and comparison position can share a node. Each occurrence
retains its ordered actual sequence choice, even if unit lengths differ.

B: Consecutive compatible units become one run with exact observed counts. No
unobserved intermediate count or arbitrary cross-path variant combination is asserted.

C: Stable internal subtype changes are eligible, including U01/U02. Each original
unit is compared directly with a frozen actual observation at the same comparison
slot. Maximum two edit operations per unit; total edits for each P must be at most
floor(original_length * fraction). The reference is never edited. UNKNOWN,
partial, terminal-template, weak-approximate and ambiguous occurrences retain
their bases. Unit counts and order of occurrences remain fixed; approved type
changes are explicitly recorded. A represented whole string may be synthetic,
so no other observed P is falsely assigned as its full representative.

Candidates use a deterministic observed representative per type and local slot.
Same-source-type/same-target batches are evaluated against the parameterized,
chain-compressed graph. Only positive N+E reductions are accepted. This is an
audited greedy search, not a claim of globally optimal lossy compression. Stable
same-type microvariants already have zero-cost parameterization; their bases are
not rewritten merely to reduce the variant dictionary.

D: Short UNKNOWN fragments at the same selected comparison slot and with the
same neighboring position instances/types/statuses share a display family,
regardless of length. Every P retains the actual sequence, length and ambiguity.
A block at least one catalog unit long remains a separate literal block to avoid
hiding an entire additional structure. Merely matching neighbor motif names is
insufficient. Alignment options remain available; selected slots are layout
coordinates, not claims of unique homology. UNKNOWN is never automatically filled,
deleted or relabeled as a known motif. Different uncertain candidate sets remain
in the per-P annotation rather than being declared equivalent.

E: Enumerate all alternatives between consecutive mandatory boundaries in the
primary-optimal DAG. Holding the other paths fixed, choose the alternative with
the fewest nodes in the full compressed display, then fewest edges, then a stable
lexical rule. Iterate over regions and P IDs until no choice changes. All regional
alternatives and their dependencies are retained through the original DAG. This
is a regional fixed point, not a proof of a joint global optimum. P0001's annotation
is frozen as the common coordinate origin. Alternative counts are model-derived
annotation counts, not biological copy-number truth.

## Outputs and comparison

Run a zero-loss control plus independent 1%, 3%, 5%, 10% budgets. Select the lowest
actual edit burden among requested budgets attaining 90% of the largest N+E
reduction over the v3 zero-loss control. If several are identical, choose the
smallest budget. Report the old v2 result separately so zero-loss parameterization
and actual base rewriting do not get conflated.

Every output preserves all P. Display nodes contain count/variant/length options;
their per-P visit vectors expand exactly into the literal represented GFA and
string. A reversible edit ledger restores the unchanged input. Parent-node route
support is not the exact support of any single representative template.

The graph is an observation display. Its free walks are not additional alleles,
and interpretation choices do not add P records or molecular support.
