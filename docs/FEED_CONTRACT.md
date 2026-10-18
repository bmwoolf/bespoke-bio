# Proposed upstream → design contract for team review

This is a working proposal based on the repository's existing upstream output, not an agreed external standard. The service and client already execute this format. The team can review the fields here and adapt the boundary without changing the BoltzGen scientific pipeline.

## Simplest starting point

For one target, POST `/v1/designs` with `sequence`, `binder_length` (default 75) and `num_designs` (default 2). Or send `structure_text`, `format`, and the author `chain`. No epitope is required. A request prepares inputs; a separate explicit `/run` request starts a managed GPU job on the shared service, or uses a configured static worker in local mode.

For the existing upstream feed, use a directory containing:

- `targets.jsonl`: one translated `Protein` record per line, as emitted by `bespoke_bio.api.build_bundle`.
- `provenance.json`: `case_id`, `reference` (`genome_build`, `annotation`), plus source provenance.
- Declared relative structure files, if available; unresolved proteins retain their sequence.

The supplied client embeds declared structure text and computes record counts. It checks sequence hashes and rejects paths escaping the package. It preserves variant and expression metadata. RNA expression is not used as a binding objective or interpreted as protein abundance.

## Submit without hand-assembling JSON

```python
from biocad.client import DesignClient

client = DesignClient("http://127.0.0.1:8000")
intake = client.submit_bundle("path/to/upstream-package", num_designs=2, binder_length=75)
for record in intake["records"]:
    child_id = record.get("child_run_id")
    if not child_id:
        print("Intake failed:", record)
        continue
    print(child_id, client.status(child_id))
    # Explicit paid inference on the managed service (or configured local worker):
    # client.start(child_id)
    # status = client.wait(child_id, timeout=3600)
    # if status["status"] == "completed":
    #     results = client.results(child_id)
    #     client.download(child_id, f"results/{child_id}.zip")
```

Equivalent CLI preparation:

```sh
uv run biocad submit-bundle path/to/upstream-package --output outputs/team-intake
# Explicitly execute, wait, and download completed child jobs:
uv run biocad submit-bundle path/to/upstream-package --output outputs/team-live --execute --wait
```

Each prepared upstream record has an independent child run. Parent status/results indexes refer to those children. One failure does not overwrite another target's results. An unavailable managed queue or missing static-worker configuration returns 503; pending candidate results return 409. A client wait timeout does not cancel remote inference.

## Downstream return values

`GET /v1/designs/{child_run_id}/results` returns typed JSON with:

- Run ID, status, mode, and target-state IDs.
- Candidate ID and exact binder sequence.
- Native filter outcome, failed-filter names, flags and metrics.
- Interface confidence metrics, units, and evidence classification.
- Binder-only and complex structure download URLs, the portable JSON bundle URL, and ZIP URL.

Before completion, inspect status/progress instead. Cached replay uses `actual_cached_replay`; fresh execution uses `live_model_outputs`. Native filter pass and measured binding remain separate claims. Stable component IDs and target sequence hashes are preserved in the full bundle.

## Review decisions

1. Confirm whether the upstream team will emit the current `Protein` record shape or needs a thinner adapter.
2. Confirm sequence-only versus structure-first handoff, including chain selection for multichain inputs.
3. Confirm the shared API address/token with the team. Managed mode owns GPU startup/shutdown, conservative spend reservations and persistent status; local static-worker mode remains available.
4. Confirm which fields the downstream team needs in addition to candidate sequence, structures, scores and pass/fail.

The proposed wire schema and examples are in `examples/design_api_requests/openapi.json` and the adjacent request fixtures. No new GPU pod is needed to review preparation or saved-result integration.
