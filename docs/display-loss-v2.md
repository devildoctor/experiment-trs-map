# D1S80 display loss v2

This opt-in policy preserves all original P IDs and molecule memberships.
The existing v1 policy and default partition entry points remain available.

Run the four-budget comparison on an already verified B2 input:

```sh
python scripts/run_display_loss.py --input INPUT/source_members.json \
  --model catalogs/D1S80.structure.v1.json \
  --policy catalogs/D1S80.structure-policy.v2.json --out NEW_OUTPUT
```

The original B2 cut points are not extracted again. The input file hash, frozen
model hash, policy hash, and code commit are recorded. Every original string,
reference P0001 and all molecule provenance are preserved. A single-budget
full BAM workflow is also available through vntrpathgraph_cli.py with
`--graph-mode structure --structure-policy ...v2.json --loss-budget 0.01`.

Partition selection uses the whole primary-optimal DAG: maximize exact full
template coverage, minimize genuinely unknown bases, then fixed orientation,
motif ID, boundaries and operation order. Its original optimal-DAG ambiguity
intervals remain. Unit correspondence minimizes primary type edit cost, then
maximizes exact paired bases, minimizes gap runs, and chooses left reference
gap positions before operation order. Literal bases at the same comparison
position can share independently of RT/internal annotation labels.

Protected signatures retain ordered unit identities and orientations. Unknown
content is retained literally. A terminal partial within 2 bp of its template
can occupy the same unit slot as its complete version. Each corresponding
unit permits at most 2 base edits; whole-path loss is floored from ORIGINAL
length. The original and represented coordinates are both in losses.jsonl.
Terminal completion is billed once as I or D, not as an additional event.

Observed representatives are prioritized by complete unit explanation, HiFi
support, total support, template deviation and hash. The reference is immutable
and excluded as a lossy target. Once used as a target, a representative is
locked; original members bind directly to it, never through a chain.

Candidate batches come from representative groups and current node/edge owner
sets. Each is evaluated by actual active display nodes plus edges AFTER the
same run parameterization and unbranched chain compression. Only positive
complexity gain is accepted. Greedy order is gain/edit bp, larger gain, fewer
edits, target priority, stable pair ordering. This is not a global optimizer.

GFA S records always contain literal fixed bases; links add no bases. Display
repeat counts are per path visit and only observed count sets are recorded.
Modules compress unbranched chains without hiding variable repeat parameters.
Both parameter expansion and recorded GFA expansion must reconstruct the
represented string. Local exact molecule support is projected from each edit
script, not copied from whole-representative support.

0%, 1%, 3%, 5% are independently evaluated. Dominated results are excluded;
select the lowest actual total edit loss reaching 90% of maximum tested
complexity reduction, then the lower budget on a tie. This is an explicit
display convention, not a biological neutrality threshold or genotype model.

Outputs include original/represented GFA and FASTA, full two-sequence TXT,
parameterized display_graph.json, annotation and correspondence audits, direct
merge batches, per-P reversible edit scripts, molecule memberships, four-tier
comparison, and a zoomable single SVG canvas. Only registered full walks are
observations; arbitrary graph recombinations are not additional sequences.
