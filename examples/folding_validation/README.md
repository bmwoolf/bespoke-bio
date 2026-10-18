# Folding validation scientific handoff

## Start here

- **Find a structure:** [clickable structure index](STRUCTURE_INDEX.md) or [spreadsheet-friendly CSV](STRUCTURE_INDEX.csv). These map every CIF to its target, binder, condition, seed, chain roles and native filter status.
- **Understand the results:** [scientific findings](FINAL_FINDINGS.md) and [summary figures](summary_figures/README.md).
- **Get exact scores and provenance:** [science.json](science.json). Each index entry has a lookup key into this file; paths are relative to this directory.

## What is in each directory?

| Directory | Contents |
|---|---|
| `structures/controls/` | 36 target-only folds testing MSA and template effects. |
| `structures/candidates/` | 36 standard reevaluations of earlier binders and controls, including complex and binder-alone cases. This is not the new-design shortlist. |
| `structures/positive_control/` | 3 predictions of the known natural IGF1–IGFBP4 interaction using paired/unpaired MSAs. |
| `structures/variants/` | 9 SERPINA1 reference/S/Z folds across three seeds. Variant backgrounds are described in the findings. |
| `structures/focused_validation/` | 18 standard predictions for the 3 new native-filter-passing binders: complex and binder-alone, each at seeds 11, 22 and 33. |
| `native/FGF2/`, `native/IL1RN/` | Original BoltzGen binder/complex outputs for all 64 new candidates, including native filter failures. These differ from the standard refolded structures above. |
| `evidence/` | Associated numerical confidence arrays; the index CSV links PAE and pLDDT files where exported. |
| `summary_figures/` | PNG/PDF plots explaining the comparisons. |

## Worked example: `structures/candidates/0000.cif`

This is a **predicted HMGB1–binder complex**, not a standalone binder or the best-ranked result.

| Field | Value |
|---|---|
| Target | HMGB1 experimental two-box monomer construct |
| Binder ID | `binder_08f2d30ba3c804fb` |
| Target chain | A, 166 residues |
| Binder chain | Z, 75 residues |
| Prediction seed | 11 |
| Evaluation | Fixed-sequence reevaluation of an earlier designed candidate |
| Original native filter | Pass; does not establish experimental binding |
| Exact lookup | `science.json` → `reports.candidates.rows[0]` |
| Structure | [0000.cif](structures/candidates/0000.cif) |
| Confidence arrays | [PAE](evidence/candidates/0000_pae.npz), [pLDDT](evidence/candidates/0000_plddt.npz) |

`0000` is only an export index within that report. Use the **full relative path**, never the basename alone, as the file key. A seed labels a repeated model prediction, not a different sequence or an experimental replicate. Open CIF files in a molecular viewer such as ChimeraX or PyMOL; use the chain mapping in the index to distinguish target and binder. Chain Z is the binder in this example, but do not assume this for every native output.

## Reading the scores

- **pLDDT (0–100):** local folding confidence; higher is more confident.
- **Directional chain-pair ipTM (0–1):** interface confidence in each reported direction; higher is more confident, not a binding probability or affinity.
- **PAE (Å):** predicted relative-position uncertainty; lower is more confident.
- **RMSD (Å):** geometric difference; read the comparison label. Agreement with a design is different from agreement with an experimental structure.
- **Native pass/fail:** the original BoltzGen filter result, separate from subsequent validation and experimental binding.

## Interpretation limits

Completed computational evidence; scientific interpretation remains required.

See science.json for scores and structure links. Native filters are preserved separately from standard fixed-sequence validation. All native bundle candidates are retained, including failures; counts are not restricted to successful examples.

- Computational predictions only; no experimental binding, selectivity, safety, or benefit demonstrated.
- Native design filters and independent fixed-sequence validation are separate evidence.
- When provided, posthoc target reference geometry is evaluated after prediction; those reference coordinates are not conditioning inputs to the predictor.
- Template-assisted recovery is not independent fold recovery; historical structures may be training data.
- Three seeds measure limited prediction variability, not equilibrium ensembles or statistical significance.
- SERPINA1 primary comparison uses shared experimental core positions, excluding variant sites; all-residue comparisons can be dominated by the unresolved amino-terminal tail.
- SERPINA1 variants are single substitutions on the canonical background, not full natural haplotypes.
- The experimental positive binder and shuffled decoy calibrate this workflow only; the decoy is not a validated nonbinder and may differ in folding.
- Oligomer context and competing biological partners require separate interpretation; a confidence score does not establish binding.
