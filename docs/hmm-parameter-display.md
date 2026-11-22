# Frozen-loss HMM display parameterization

Scope: move weak-match status and short-skip neighbor context from node identity
into per-visit parameters; mark `LOSS*` only for changes inside that node.

The runner consumes a completed single-Viterbi HMM result. It preserves its
annotations, comparison slots, represented strings and exact 1% loss ledger.
It does not run the HMM, graph-based annotation selection or loss optimizer.

```sh
python scripts/run_hmm_parameters.py --source PREVIOUS_RESULT_DIR --out NEW_RESULT_DIR
```

Same-position, same-type units share a family even when their weak status differs;
consecutive compatible units can form one run. Each occurrence retains its actual
string, `edit_distance` and `weak_approximate`. The main label uses `~` when any
instance differs from its template; detailed weak status is in the visit parameters.

Short skips retain their own comparison position and occurrence order. Their
left/right original neighbor context is saved in `skip_context` for each visit,
rather than splitting the node family. Skips at least one template unit long
retain their existing literal-block grouping.

Local loss membership compares each occurrence's original and represented
strings, so substitutions and boundary insertions/deletions are counted at the
affected node. A path edited elsewhere does not color this node. Only node loss
markers change; edge support calculations are unchanged.

The runner checks byte-identical original/represented FASTA, GFA, sequence TXT
and loss ledger against its input result. The existing per-P loss maximum remains
1%; no additional edits or new edit targets are selected.
