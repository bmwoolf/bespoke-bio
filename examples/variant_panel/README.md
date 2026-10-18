# Public reference / variant panel

Seven sequence states across four human proteins: GDF15 reference/H202D, factor VII reference/R413Q, APOE3/APOE4, and VEGFA165 reference. These are **mature research constructs**, not full precursor sequences or individual patient records.

- `targets.jsonl`: one upstream-compatible JSON object per line, using `target_id`, `gene_symbol`, `phase_scenario`, `amino_acid_sequence`, and its SHA256. Preparing these inputs does not require GPU inference.
- `panel.json`: complete sequences, construct coordinates, variant mappings and source metadata.
- `SOURCE_NOTES.md`: primary/database sources and biological limitations.

VEGF has no fabricated coding variant: commonly studied regulatory SNPs do not change the amino-acid sequence. Every coding pair differs at one recorded residue. These targets cover stress signaling, coagulation, vascular signaling and lipid biology; their inclusion does not establish a therapeutic benefit from binding or inhibition.

See [comparison workflow](../../docs/variant-comparison.md). The local campaign's complete native outputs and review report are stored in `outputs/variant_campaign`, which is excluded from Git.

## Completed campaign — 17 October 2026

Generated 28 distinct candidate records across seven target states, then ran 75 fixed-sequence comparison predictions: 21 target-only predictions and 54 complexes (nine selected binders × two target states × three seeds). All 75 completed and passed sequence/coordinate/output validation.

| Target | Reference candidates passing native filters | Variant candidates passing native filters |
|---|---:|---:|
| GDF15 / H202D | 2/4 | 1/4 |
| Factor VII / R413Q | 0/4 | 0/4 |
| APOE3 / APOE4 | 0/4 | 0/4 |
| VEGFA165 | 0/4 | Not tested |

**No convincing variant-selective binder emerged.** Native structural filter passes are not experimental binding validation. GDF15's three passing candidates had low interface confidence in the subsequent no-template comparison. The most persistent APOE candidate (`binder_fd9bee82d286d9c3`) scored mean interface ipTM 0.418 on reference versus 0.430 on the variant, with overlapping replicate variation; it also failed the original structural filters. These are model confidence scores, not affinities.

Target-only predictions were often low confidence, and between-variant structural differences were not clearly separated from variation among repeated predictions of the same sequence. They do not establish biological conformational effects. Native generation evaluation uses target templates while the fixed-sequence comparisons do not; absolute scores across those protocols are not directly comparable.

- `binders.fasta`: all 28 candidate sequences with native pass/fail labels.
- `candidates.csv`: candidate scores and native filter outcomes.
- [Fixed-binder comparison](fixed_binder_confidence.png).
- [Target structural variation](target_structure_variation.png).

The complete local archive includes structures, native outputs, logs, inputs, provenance and verification receipts. All 835 downloaded files passed checksum verification. The campaign GPU was deleted and deletion verified; estimated GPU usage was $5.51 excluding storage, within the authorized $100 ceiling. The team API host and persistent model cache were retained.

Next scientific priority: reliable structures in the correct biological assembly, followed by larger candidate batches and independent complex assessment. This pilot does not support recommending any candidate for treatment.
