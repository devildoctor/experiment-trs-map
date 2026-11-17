> 2026-11-17：当前仅保留单元注释和终端证据库；下方旧实验结果为历史记录。

# D1S80 right-terminal templates (opt-in)

The shared D1S80 last 16 bases were absent from the ordinary dictionary. Its
decomposition could consume one base as an insertion into the preceding unit,
leaving 15 unknown bases. This feature nominates exact oriented terminal 16-mers
from REF and full OBS and validates them using independent molecules. It does
not modify repeat boundaries, sequences, the ordinary dictionary, the GFA, or
joint-support rules. Candidate paths never nominate or vote; REF only nominates.

## Current library interface (2026-11-17)

The former analysis CLI has been removed. Terminal evidence validation remains
in `scripts/vntrterminal.py:validate_terminals`; sequence-only annotation remains
in `scripts/vntrunits.py:primitive_units`. These are Python library functions,
not standalone command-line entry points. The caller must supply the required
records, molecule evidence, anchor decisions and unchanged reference context.
The former CLI's input/provenance checks are not automatically supplied by these
libraries. See their signatures and the terminal evidence/integration tests.

For each sample and exact oriented terminal template, a supporting molecule must:

1. Have a frozen two-flank-anchored OBS and a `complete_support` anchor decision.
2. Match the stored source records, original/stored query intervals, sequence hash
   and actual raw query slice. Conflicting observations cannot vote.
3. Match the terminal 16 bases exactly; its immediately preceding 16 bases must
   fit an original primitive template with at most one edit. New terminal templates
   cannot validate this inner context.
4. Match the next 16 bases of the unchanged reference in its actual raw query.
5. Meet at least three independent molecule keys in that sample. Duplicate paths
   are not independent; platform subtotals may overlap. Local paths cannot vote.

The fixed-phase primitive period must be 16. Existing exact oriented templates
need no new transition. Other nominations, failed contexts and rejection reasons
remain in `terminal_unit_candidates.tsv` and `terminal_molecules.tsv`. These files
contain genomic bases and molecule identities and remain on the analysis server.

Accepted templates are supplied to `decompose(right_terminal_motifs=...)`, which
only allows their full/partial matches to end at the original sequence end. It
does not freely rotate or reverse the supplied terminal template. Existing edit
thresholds, scores, exact slices and genuine indels remain. The default empty
argument leaves all existing callers unchanged. Only the sample's accepted set
is used for its candidates/full OBS; REF uses the accepted union. Local OBS keep
the old dictionary. A terminal template is never made globally available internally.

The new run preserves primary/alternative annotations and identifies terminal
templates in unit rows. `annotated_full_units` is a descriptive count. All
`primitive_cn` values are explicitly `NA` under `withheld_terminal_validation`,
including otherwise resolved counts. No genotype or main graph is updated.

The former Slurm wrapper and its graph-dependent acceptance validator were
removed with the analysis CLI. Their original workflow is available in the
[historical document](https://github.com/devildoctor/experiment-trs-map/blob/9c5739f18adf31add56b9fb36752a5eeb527ad89/D1S80_TERMINAL_UNITS.md).
The results below describe that historical run, not a current runnable command.

## Validated frozen D1S80 run (2026-10-23)

- Implementation commit: `cc26593d277bb5b0d6013711ed9da03698a57f96`.
- Run: `runs-result/20261023T141043Z_d1s80-terminal-units`.
- Baseline: `runs-result/20261023T124501Z_d1s80-joint-units-complete`.
- Input: `runs-result/20261021T090055Z_vntr-dbg_D1S80` (199 paths).
- Slurm 163620: cn06, COMPLETED 0:0, 67 seconds. All 186 unit tests passed.
- Default mode's eight TSVs are byte-identical to baseline. Opt-in mode's five
  joint-evidence TSVs are byte-identical. Frozen graph artifacts passed all
  before/after hashes. Source snapshot hashes passed before/after.
- 47 distinct terminal nominations; only the shared terminal template passed
  the same-sample threshold. Independent molecule support: LCL5 14 (HiFi9/ONT5),
  LCL6 16 (HiFi13/ONT3), LCL7 17 (HiFi14/ONT3), LCL8 12 (HiFi9/ONT3), total59.
- All eight candidates and REF: terminal unknown15 -> 0, partial0, full units+1.
  Primary annotation before the final32 bases is unchanged; sequence and BED
  are unchanged. All 199 primitive_cn values are NA.
- C10 stays432bp, annotated full units26 -> 27, score40.0 -> 9.9. Its final
  intervals are [400,416) and [416,432), both16=. Other annotation ties remain;
  this is not a released CN27 or genotype call.
- Candidate full-unit changes: LCL5_C10 26->27; LCL5_C4 22->23;
  LCL5_C5 22->23; LCL6_C2 22->23; LCL6_C6 26->27; LCL7_C6 22->23;
  LCL8_C2 22->23; LCL8_C5 26->27. REF28->29.
- Review `results/D1S80/comparison.tsv` and `verification.json` under the run.
  Raw terminal-template/molecule evidence stays server-side. No GitHub upload.
