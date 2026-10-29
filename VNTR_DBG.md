# VNTR compact de Bruijn graph pilot

The explicit `--dbg-view` entry point is independent of the existing point,
read, structure and STR branches. ATXN2 remains an STR unit graph. D1S80
explicitly selects `compact_dbg` in `catalogs/grch38_graph_models_v1.tsv`;
other loci retain their existing paths pending a structural assessment.

## Semantics

* Raw vertices are exact, reference-oriented k-mers. An edge must be observed
  as a contiguous k+1-mer in an input sequence, never inferred from overlap.
* Nonbranching node-disjoint chains are compacted. Every path start/end is a
  compaction boundary. Each S contains the actual chain sequence; each L has
  an exact k-1 overlap. Node visits do not directly count repeat units.
* P is a standard ordered list with explicit k-1 CIGAR overlaps. Reference,
  observed and candidate paths are distinguished by `pt` and `cv` user tags.
  Reconstruct by taking the first S then removing k-1 bases from each next S.
* Complete observations reuse current two-flank extraction decisions. Other
  fragments require a primary or supplementary target alignment, MAPQ>=20,
  consistent molecule SEQ, and no excluded anchor conflict. They remain
  `alignment_local_only`, not flank-unique or full-allele evidence. N CIGAR
  operations and ambiguous bases split paths; supplementary SEQs are never
  stitched. Each stored slice and source record are retained.
* REF is the catalog repeat interval only, with no synthetic flanks appended
  to observations. CAND paths retain the previous inference and uncertainty;
  reference and candidates add zero observed-molecule support.
* Molecule support is a set, not traversal count. Platform counts may overlap.
  Every node/edge carries sample/platform support in the JSON/TSV sidecar.
* `.supported` retains raw nodes/edges supported by at least `--min-support` (default 3) independent
  molecules within one sample. REF is retained as labeled context. Observed
  and candidate paths are split at removed nodes/edges; no bridge or complete
  path is fabricated. All evidence remains in the unsuffixed graph.

## Output naming

Runs: `runs-result/<UTC>_vntr-dbg_D1S80/`, with baseline/source snapshots,
SHA256SUMS, logs, jobs and acceptance evidence. Data remain server-side.

`results/dbg/D1S80/k33/D1S80.k33.gfa` is the complete evidence graph.
The paired `.json`, `.svg` and `.dot` represent exactly the same topology.
`.paths.svg` includes every exported P; `.summary.paths.svg` shows only REF
and CAND for inspection. `paths-svg/` contains one image per exported P.
The `.supported.*` files and `paths-svg.supported/` follow the same contract.
`*.paths.tsv` records path provenance, and `*.support.tsv` maps nodes/edges
to molecule sets. `motif_annotations.json` preserves existing candidate token
coordinates, exact bases and ambiguity. Annotation colors identify exact
slices; a motif can cross a node boundary. OBS paths are not newly motif-called.

Each graph records its actual `min_support` and per-sample support rule. Edge color compares the maximum within-sample molecule support with this threshold; line width and molecule labels still describe total distinct molecules. REF context can remain below threshold.

All SVG node sizes are schematic. Labels show sequence length and molecule
support; full sequence is in GFA/JSON and SVG tooltips. Raw graphs above 450
nodes use Graphviz neato with VPSC overlap removal, straight edges and 100
layout iterations; smaller graphs use dot. This works with the installed
Graphviz without optional triangulation support. The layout choice does not
change the graph. Per-path boxes are in exact traversal order, wrap after 8
visits, and show newly spelled bases after overlap subtraction.

## Execution and acceptance

Use `slurm/vntrdbg.sbatch` after fresh pointgraph extraction from the same BAM,
catalog, reference and explicit ALT-equivalence table. Run all computation,
tests and rendering through Slurm. The six k values are 17,21,31,33,49,65.
Existing extraction files are checked against catalog, reference metadata and
BAM inventory; the run driver verifies BAM checksums before and after work.

Acceptance independently parses GFA and verifies every sequence, overlap,
source k-mer edge, unique-molecule support, topology SVG and per-path SVG.
The normal STR branch is rerun on actual ATXN2 BAMs and compared byte-for-byte
with the prior STR results. Input BAM/BAI checksums must remain unchanged.

Compare k using supported branch counts, merge-and-split ambiguity,
candidate edge retention, repeated candidate k-mers, and platform-stratified
full observation retention. These are descriptive readability and evidence
metrics, not validated genotype accuracy. Smallest/cleanest graphs are not
automatically best. ALT equivalence remains an explicit location assumption;
changing graph k does not replace flank validation.

Format: https://gfa-spec.github.io/GFA-spec/GFA1.html
