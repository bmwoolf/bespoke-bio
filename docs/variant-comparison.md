# Reference and variant comparison

The comparison tools separate three questions:

1. Do independently generated candidates pass BoltzGen's native structural filters?
2. How much do repeated target-only predictions vary within a sequence and between reference and variant sequences?
3. Does the **same fixed binder sequence** receive different predicted interface scores against the two target sequences?

Different generated sequences or a native filter pass do not establish binding, affinity, allele specificity, agonism, or antagonism. Target-only repeat predictions are model samples, not equilibrium molecular ensembles.

## Inputs and execution

`scripts/refold_fixed_sequences.py` runs in the pinned BoltzGen environment. Supply a manifest with `tasks`, each containing a unique `id` and `target_sequence`. Add `binder_sequence` for a complex; omit it for a target-only fold. For analysis, retain `target_id` and `scored_on`; complexes additionally require `binder_id` and `designed_on`.

```json
{"tasks": [{
  "id": "target_reference_alone",
  "target_id": "target",
  "scored_on": "reference",
  "target_sequence": "ACDEFGHIKLMNPQRSTVWY"
}]}
```

The sequence above is a format example, not a biological benchmark.

```bash
python scripts/refold_fixed_sequences.py tasks.json \
  --checkpoint /cache/boltz2_conf_final.ckpt \
  --moldir /cache/mols.zip \
  --output crossfold --seeds 11 22 33 --validate-only
```

Remove `--validate-only` to run GPU predictions. The helper preserves fixed sequences, uses empty MSAs and no target templates, and retains one prediction per task per seed. It fails if any expected prediction archive is missing. Lightning/Torch seeds are recorded, but native feature construction initializes a separate NumPy RNG, so bitwise reproducibility and identical random draws across sequence pairs are not claimed.

## Reports

A campaign directory contains `generation_manifest.json` with `cohorts`; each entry identifies `target_id`, `state`, `sequence`, `directory`, and execution `status`. Completed candidate bundles live at `<directory>/live/bundle/bundle.json`. Refold outputs retain `receipt.json` and `seed_N/refold_cif/*.cif`.

```bash
python scripts/build_variant_campaign_report.py outputs/variant_campaign \
  --refold outputs/variant_campaign/crossfold \
  --output outputs/variant_campaign/report
python scripts/render_variant_campaign_report.py outputs/variant_campaign
```

Native design evaluation uses target templates; the fixed-sequence comparison deliberately uses none. Absolute scores across those protocols are not interchangeable. Low-confidence target-only folds cannot establish biological conformational changes.

The report includes native pass/fail and failure reasons in linked bundles, Boltz and ESMFold2 score distributions, exact sequence diversity, target C-alpha RMSD, fixed-binder cross-scores, and mutation-to-binder heavy-atom distances. Chains are located by exact sequence; ambiguous matches are reported unavailable. The 5 Å contact cutoff is a descriptive geometric convention, not an energetic criterion.

Candidate selection is recorded before cross-scoring: prefer native-filter passes, then higher Boltz interface confidence; also include the best mutation-contact candidate when distinct. If none pass, any selected candidate remains explicitly labeled failed and exploratory. No epitope constraint is used during generation.

The [public input panel](../examples/variant_panel/README.md) comprises mature GDF15 reference/H202D, factor VII reference/R413Q, APOE3/APOE4, and VEGFA165 reference. Common regulatory VEGF SNPs are not converted into fabricated protein substitutions. Native dimers, lipidation, clotting modifications and cofactors are omitted in these single-chain inputs; biological interpretation requires appropriate assembly context and experimental validation.
