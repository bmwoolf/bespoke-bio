# Bespoke Bio

This repository implements protein translation/packaging and a target-to-binder design
pipeline. The active design engine is **BoltzGen**, with no required epitope or hotspot.
Targets arrive as amino-acid sequences, PDB/mmCIF structures, or an upstream target bundle.
Results include designed sequences, predicted complexes, filter outcomes, confidence scores
and downloadable provenance.

The `biocad` package and CLI name are retained in the code. A CAD website, interdigitated
assembly design and folded multi-module constructs are outside the active delivery.
Earlier experiments are preserved under [research history](docs/LEGACY.md).

## Run the combined API

```sh
uv sync --extra design --extra dev --locked
uv run uvicorn biocad.service:app --host 127.0.0.1 --port 8000
```

Open `http://127.0.0.1:8000/docs`. The combined app includes translation and design routes.
Preparation and saved-result replay work locally. In shared managed mode, explicit job
submission starts a disposable GPU, collects results and deletes it under a persistent
$100 phase budget. The authenticated gateway keeps job state on separate persistent storage.
Local static-worker mode remains available. See [shared service operation](docs/SHARED_SERVICE.md).

Try the real saved PD-L1 demonstration:

```sh
curl -X POST http://127.0.0.1:8000/v1/designs/demo \
  -H 'Content-Type: application/json' --data '{}'
```

Use the returned run ID with `GET /v1/designs/{run_id}/results` for candidate JSON or
`GET /v1/designs/{run_id}/download` for the full ZIP. The demonstration is explicitly
labeled as cached replay; it does not generate new candidates for the request.

## Connect the upstream feed

The proposed team contract uses `targets.jsonl`, `provenance.json`, and declared structure
files. The client embeds structure text and preserves target, variant and RNA metadata.
This is a proposal for review, not an agreed external standard.

```sh
# Prepare independent jobs from an upstream package; no inference starts.
uv run biocad submit-bundle path/to/upstream-package --output outputs/team-intake

# Explicitly execute on a configured worker, wait and download results.
uv run biocad submit-bundle path/to/upstream-package --output outputs/team-live --execute --wait
```

Start with the [feed contract and Python client](docs/FEED_CONTRACT.md),
[API walkthrough](docs/HACKATHON_HANDOFF.md), or [OpenAPI schema](examples/design_api_requests/openapi.json).
For producing the upstream package, see [translation and packaging](docs/TRANSLATION.md).

## What is verified

- Full seven-stage BoltzGen runs completed for public PD-L1 structure and sequence inputs.
- The fresh sequence-only HTTP job returned a 75-residue candidate that passes native
  filters, with complex recovery RMSD 0.83 Å and ipTM 0.41.
- A fresh GPU restored the persistent runtime/model cache in 76 seconds. Both test GPU pods
  were deleted; the persistent cache remains billable.
- Automated tests cover upstream package/client integration, authenticated shared access,
  persistent queue/budget rules, GPU cleanup, native result export and failure handling.

A native filter pass is computational evidence, not measured binding. Interface confidence
is modest in the retained examples. See [fresh API results](examples/boltzgen_api_results/summary.json)
and [the earlier two-candidate report](examples/boltzgen_results/review.html).

The [expanded target and variant study](examples/expanded_variant_panel/README.md) includes
public sequences, candidate scores, repeated comparisons and assembly checks. Across both
research panels, 132 designs yielded 15 native-filter passes; no variant-selective binding
was established. See the [findings and limitations](examples/expanded_variant_panel/FINDINGS.md).

## Documentation

- [Current implementation and limits](docs/STATUS.md)
- [Active pipeline and API architecture](docs/INTEGRATION.md)
- [GPU operation and persistent cache](docs/COMPUTE.md)
- [Optional bounded Foldseek template retrieval](docs/RETRIEVAL_DEMO.md)
- [Remaining work](docs/PUNCHLIST.md)
- [Upstream system boundaries](docs/architecture.md)

## Verify

```sh
uv run pytest
uv run ruff check src tests
```

This software produces research predictions, not clinical results.
