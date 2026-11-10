# D1S80 HMM single-path display

`scripts/run_hmm_display.py` annotates every frozen P sequence with an independent
profile HMM and one Viterbi trace. The motif/skip concept follows TRGT; this Python
implementation has its own topology and is not expected to reproduce TRGT scores.
M/I/D states describe motif variants; a background skip state emits literal
UNKNOWN blocks. The nine forward 16-bp templates and existing exact terminal
templates are unchanged. All original strings survive tokenization.

Equal-score paths are resolved by fixed motif/state order. There is no alternate
path enumeration and no E graph-based annotation selection. An absent ambiguity
flag means only that one trace was selected, not that annotation is unique.

The existing A/B/D parameterization and C local sequence replacement are reused.
The requested output has a 1% maximum edit budget per original P and at most two
edits per eligible unit. UNKNOWN, weak matches and terminal templates retain
their strings. HMM alignment differences themselves do not consume this budget.
Only graph-simplifying replacements consume it. The reference is never rewritten.

Run from the repository:

```sh
python scripts/run_hmm_display.py --input SOURCE_MEMBERS.json \
  --model catalogs/D1S80.structure.v1.json --out NEW_RESULT_DIR
```

The runner writes a 0% control and the requested 1% graph, with SVG, literal GFA,
full original/represented sequences and per-path parameters. The full BAM CLI
also accepts `--structure-policy catalogs/D1S80.structure-policy.hmm.json`
and `--loss-budget .01` in structure mode.
