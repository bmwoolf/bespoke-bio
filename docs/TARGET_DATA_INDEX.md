# Target and binder data index

The broad campaigns and focused validation are separate datasets. Start with [the target-to-folder CSV](TARGET_DATA_INDEX.csv): it maps all 30 completed broad/focused generation cohorts to readable local folders and recorded Runpod S3 job prefixes. S3 paths come from execution receipts, not a fresh S3 availability audit. The persistent volume is retained after GPU deletion.

- Original panel: GDF15 (reference/H202D), F7 (reference/R413Q), APOE (E3/E4; precursor C130R), VEGFA165. 7 states, 28 candidates. GitHub sequences/scores: `examples/variant_panel/`. Full local results: `outputs/variant_campaign/`; archive `variant_campaign_results.tar.gz`.
- Expanded panel: IGF1, PF4, SERPINA1 (reference/E288V/E366K), FGB (reference/R478K), FGF2, IL1RN, CXCL12, TFPI, IL6, MSTN, THBD, FGF21, CCL2, CXCL8, IL18, LEP, HMGB1, C5. 21 states, 104 candidates. GitHub sequences/scores: `examples/expanded_variant_panel/`. Full local results: `outputs/variant_expansion/`; archive `variant_expansion_results.tar.gz`.
- Focused redesign: FGF2 and IL1RN, 32 candidates each. GitHub structures and evidence: `examples/folding_validation/`. Full local results: `outputs/folding_validation/`; archive `scientific_campaign.tar.gz`.
- Earlier PD-L1 smoke tests: `examples/boltzgen_results/`, `examples/boltzgen_api_results/`, `examples/shared_service_smoke/` include bundles/structures.

Within local generation `live/bundle/`, read `bundle.json` for target/candidate paths; candidate folders contain `binder.cif` and `complex.cif`, and `target/` contains prepared/source structures when supplied. `live/downloaded/` retains raw native artifacts. Remote jobs contain equivalent native outputs; list the mapped prefix to inspect its exact layout. Fixed-sequence repeat predictions are in each campaign's `crossfold/`, or the latest campaign's validation outputs.

Use Runpod S3 credentials with region US-MO-1 and endpoint https://s3api-us-mo-1.runpod.io. The API-state volume hj1wqrrktq is not the scientific-results volume. Cache files are model weights, not binder outputs.

These are 22 distinct broad-panel target proteins plus PD-L1. Counts are generated candidate records, not experimentally validated binders. Many targets use domains/truncated constructs; see each input manifest. Native filter failures are intentionally preserved. IL10/EPO were prepared but not generated; CCL26 was an intake test, not a completed binder run.
