# Shared hackathon design API

The deployment uses a small, continuously reachable CPU gateway and starts GPU workers only for explicitly submitted jobs. It is one API, not another website. The gateway retains the API job queue and artifacts under `/workspace/shared-service`; GPU model weights remain on their existing **separate** network volume so subsequent workers reuse them. The CPU gateway uses its own state volume; GPU workers never mount that CPU volume.

## Operator deployment

The operator first provisions the gateway and registers the SSH public key. CPU gateway selection must be based on the live Runpod CPU catalog and a compatible data center. The current CLI supports CPU creation but does not expose CPU flavor/vCPU sizing flags; the REST v2 CPU catalog and create schema supply those fields. The REST schema supports a network volume for CPU pods, while a persistent pod volume is not supported for CPU pods. A dedicated CPU-only network volume is used for gateway state. It must be distinct from the GPU model-cache volume; never mount the CPU state volume on GPU workers.

Create a **private, mode-0600 JSON file outside tracked artifacts** containing environment configuration. Generate the service bearer token locally using a cryptographic random generator; never put the token or Runpod key in Git, URLs or public command arguments. Required entries:

- `BIOCAD_MANAGED_WORKERS`: `"1"`
- `BIOCAD_API_TOKEN`: private shared service token
- `RUNPOD_API_KEY`: private Runpod API key. A dedicated, restricted key may be supplied here instead of the existing account key; the helper does not create keys, replace the account key or broaden permissions. Its permissions must support the exact required pod/catalog reads and lifecycle actions.
- `BIOCAD_GATEWAY_POD_ID`: the exact gateway pod ID for its 24-hour lifetime guard
- `BIOCAD_SERVICE_EXPIRES_AT`: absolute Unix expiry timestamp, at most 24 hours ahead, shared with the service admission guard

The deployment helper sets `BIOCAD_REPO_DIR` (for installed-wheel fixture lookup), `BIOCAD_RUNS_DIR`, `BIOCAD_CLOUD_STATE_DIR` and `RUNPODCTL`. It sets `BIOCAD_SSH_KEY` and `BIOCAD_WORKER_KEY` to **CPU-host-local** `/root/.config/biocad-shared/worker_key`; credentials are never placed on the network volume shared with GPU workers. The SSH private key supplied as `--worker-key` must correspond to the public key registered with Runpod for GPU workers. Any additional cloud-manager limits should be set to the managed worker module's documented configuration; the overall phase cap remains $100.

```sh
python3 scripts/deploy_shared.py \
  --host GATEWAY_IP --port GATEWAY_SSH_PORT \
  --ssh-key /absolute/path/to/gateway-key \
  --env-file /absolute/path/to/private-environment.json \
  --worker-key /absolute/path/to/registered-gpu-key \
  --binary-path /absolute/path/to/runpodctl-linux-amd64
```

`--binary-path` is optional: absent it, the gateway downloads the pinned public Linux AMD64 CLI itself. Version 3.0.1 is checksum-verified either way. The helper creates a source/artifact allowlist archive, separately transfers the private environment and key, applies restrictive permissions and starts a **single** Uvicorn process under a local supervisor. The supervisor restarts an exited API process until the persisted deadline; SQLite startup recovery reconciles interrupted jobs. Redeploy replaces only its own supervisor/child and retains the independent lifetime guard. Managed dispatch uses a durable SQLite queue; do not launch multiple independent Uvicorn workers.

No cloud resources are created by this deployment helper. Invoking it without `--host` only builds `outputs/shared_deploy/service.tar.gz`. That archive excludes secrets, vendor installations, Python bytecode and old compressed bundles.

Gateway paths:

- `runs/`: SQLite queue, prepared requests, statuses and downloaded results
- `cloud_leases/`: managed GPU lifecycle receipts
- `state/gateway_deadline`: non-secret lifetime deadline
- `logs/`: API and lifetime guard output

Authentication material is **outside `/workspace`**: `/root/.config/biocad-shared/environment.json` and `worker_key`, directory mode 0700 and files mode 0600. The Runpod key is loaded only into CPU gateway process environments, never passed to GPU worker environment configuration or bundled with jobs. GPU workers receive only prepared design inputs and retain model caches.

Executable code, the Python environment and the pinned CLI are CPU-local under `/opt/biocad-service/{app,runtime,bin}`. They are not loaded from the GPU cache volume.

## Teammate access

Use the operator-provided HTTPS base URL (typically Runpod's port-8000 HTTP proxy) and bearer token. Swagger's **Authorize** control accepts the bearer credential. For shell use, load the token into an environment variable without placing its value in shell history:

```sh
curl -sS "$BIOCAD_BASE_URL/openapi.json" \
  -H "Authorization: Bearer $BIOCAD_API_TOKEN"
curl -sS -X POST "$BIOCAD_BASE_URL/v1/designs/demo" \
  -H "Authorization: Bearer $BIOCAD_API_TOKEN" \
  -H 'Content-Type: application/json' --data '{}'
```

For fresh designs:

1. POST `/v1/designs` with the checked-in sequence fixture or structure text. No hotspot/epitope selection is required.
2. POST `{}` to the returned run's `/run` URL. In managed mode, this queues a real GPU job; it can create a worker and incur charges.
3. Poll the status URL until `completed` or `failed`. Download all artifacts using `/v1/designs/{run_id}/download`.
4. For upstream intake, submit each record's `child_run_url` independently; use its child status and download URLs.
5. Cancel a job through its `/cancel` endpoint when needed. Check the resulting cleanup receipts; do not assume losing the HTTP connection cancels remote execution.

The demo endpoint replays genuine saved PD-L1 GPU outputs immediately. It does not produce a new design. Fresh jobs are labeled separately. A native model filter pass is a computational screening result, not measured affinity or experimental binding.

## Cost, lifetime and cleanup

The service lifetime is capped at 24 hours from its first start; its persisted deadline survives process restarts. The gateway starts an independent local guard that stops the service and reconciles outstanding GPU leases, including provisioning records whose pod IDs were not returned, confirms their deletion, then deletes **only its configured gateway pod ID** at that deadline. This is a best-effort guard: the process cannot enforce a deadline if its host fails. The operator must maintain an external deadline and verify cloud cleanup.

GPU managers should serialize work, enforce the phase's $100 cap, download results and terminate task-created GPU pods after use. Gateway CPU, disk and network-volume storage continue to cost money outside GPU execution. Stopping a pod is not deleting its retained storage. Do not delete the shared model-cache volume merely to stop a GPU worker.

A network-mounted gateway data directory survives gateway pod replacement. If gateway deployment instead falls back to non-persistent container storage, back up `runs/` and `cloud_leases/` before deletion. Keep secret backups separate and private. The gateway refuses new work near its configured expiry so a two-hour GPU job cannot be intentionally admitted immediately before teardown. Before final teardown, download the result ZIPs, inspect managed GPU cleanup receipts, stop/delete the gateway and account for any deliberately retained storage.

Replacing the CPU pod requires reinstalling its local runtime, code and credentials through an explicitly authorized deployment; a network volume restores job data, not CPU-local secrets. Do not migrate the local secret directory onto the shared volume for convenience.

A network volume is a shared trust boundary: root on any pod mounting it can read or modify its files. Isolation therefore uses **separate CPU-state and GPU-cache volumes**, plus CPU-local credentials and executable code. Root permissions alone on a shared volume would not isolate secrets. This is a scoped hackathon service, not a hardened multi-tenant deployment.

The deployment is not verified merely because the process starts. The operator must make an authenticated request from outside the pod, test the cached demo, then complete one new GPU job and verify its downloaded bundle and worker cleanup.

## Python client

Set `BIOCAD_API_TOKEN` in your environment, then use the same client against the shared URL:

```python
from biocad.client import DesignClient
client = DesignClient("https://YOUR-GATEWAY-8000.proxy.runpod.net")
print(client.service_status())  # readiness, phase budget and lifetime
# Submit a bundle, start child jobs, poll status, retrieve results and download ZIPs.
# client.cancel(run_id) requests cancellation and GPU cleanup.
```

All coworkers sharing the service token can access the same project jobs. The token does not expose the Runpod API key or grant direct cloud-account access. The application cap is a conservative spend policy, not a provider-enforced account cap or protection against a stolen cloud credential.

## Simple target JSON

For a JSON file containing `sequence` (or `structure_text`) and optional `binder_length` / `num_designs`, use:

```sh
# BIOCAD_API_TOKEN must contain the separate service token.
uv run biocad submit-json target.json --url "$BIOCAD_BASE_URL" --output outputs/first-target --execute --wait
```

The same command routes a proposed upstream JSON with a `targets` list to upstream intake. Declared structure files must be embedded as `artifact_texts`; for a package directory, `submit-bundle` embeds them automatically. Omitting `--execute` performs preparation only. Completed results and ZIPs are saved to the chosen local receipt directory.

## Verified deployment behavior

A real shared-HTTPS sequence-only PD-L1 job completed all seven native stages, returned a candidate, downloaded a verified ZIP and automatically deleted its GPU (independent read-back returned 404). Estimated GPU cost was $0.56, not finalized billing. The candidate passes native filters but has low interface confidence (ipTM 0.20); this is execution evidence, not a validated binder. Queued cancellation released its reservation without another GPU. Completed/cancelled statuses and the $10 retained reservation survived a service process restart. See [sanitized proof](../examples/shared_service_smoke/summary.json).

The current deployment retains a small CPU gateway until its 24-hour expiry; new GPU jobs are rejected within two hours of that expiry. The gateway's 10 GB state volume and the separate 100 GB model cache survive host deletion and remain billable (approximately $0.70/month and $7/month respectively).
