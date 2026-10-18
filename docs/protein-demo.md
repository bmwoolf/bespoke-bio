# Protein catalog and design demo

This demo connects a verified Stage 1 protein package to Jacob's design API. It preserves
every biological target and sends the selected record unchanged, plus `binder_length` and
`num_designs`. Expression stays null when unmeasured. The complete sequence and both hashes
remain in the request. The standalone CDS-based Step 2 is bypassed. An optional, independently published structure
library supplies exact-sequence previews; the design service predicts candidate complexes.
The preview is not sent as a structural constraint. This route does not demonstrate the
original three modules running in sequence.

The integration targets [design service PR #1](https://github.com/bmwoolf/bespoke-bio/pull/1),
commit `ee1db463e1cbc18450bbcbb0491dae040fbd024f`. Its executable API/schema is the upstream
contract. Run that revision separately; this branch does not import its Step 2 changes.

## Run on the data VM

Keep the package, SQLite index, job receipts, downloads, credentials and model files outside
the checkout. Only source, documentation and small synthetic tests belong in Git.

```bash
uv sync --extra stage1 --extra dev
# Node 22.12+; web/.nvmrc pins the tested Node release.
npm --prefix web ci
npm --prefix web run build
uv run python -m bespoke_bio.demo /data/proteins /data/demo/catalog.sqlite3

export BESPOKE_DEMO_INDEX=/data/demo/catalog.sqlite3
export BESPOKE_DEMO_RUNS=/data/demo/jobs
export BIOCAD_BASE_URL=https://your-design-service.example
export BIOCAD_API_TOKEN_FILE=/private/service-token
uv run uvicorn bespoke_bio.demo.service:app --host 127.0.0.1 --port 8787 --workers 1
```

Use a private token file with mode `600`. The token remains on the server. For an isolated
local design API, use `http://127.0.0.1:8788` and omit the token setting if that API is configured
without authentication. It still needs Jacob's source, dependencies and cached example files.

Access the VM through an SSH tunnel:

```bash
ssh -N -L 8787:127.0.0.1:8787 your-vm
```

Open `http://127.0.0.1:8787`. Keep one gateway worker and bind to loopback; this is a single-user
demo. For protected HTTPS access, use the operator proxy below. Existing raw-read jobs can continue. Limit
indexing/testing to spare CPU capacity; browsing and preparation do not need a GPU.
Use a new private runs directory when switching to another design service; receipts belong
to the service that created them.

## Public website and data API

The demo is available at [bespoke-bio-demo.westus2.cloudapp.azure.com](https://bespoke-bio-demo.westus2.cloudapp.azure.com/).
It serves the full public HG002 catalog and explicitly published completed results from the
VM. Browsing needs no token or tunnel. **Generate keys** opens `/operator/` with the selected
protein. The browser asks for a separate operator username/password; Jacob's service token
never enters the browser. After signing in, **Generate keys** asks for confirmation, prepares
the selected record, and starts its design job. Progress and returned key/complex models
appear in the same workspace. The public process remains read-only and is not configured
with service credentials or private receipts. Operator jobs are not automatically published to visitors.

To publish saved results, export selected completed job IDs from the private gateway:

```bash
uv run python -m bespoke_bio.demo.published \
  --gateway http://127.0.0.1:8787 --job COMPLETED_LOCAL_JOB_ID \
  --output /data/demo/published-release-1

# In a separate process with BIOCAD_* and BESPOKE_DEMO_RUNS unset:
export BESPOKE_DEMO_INDEX=/data/demo/catalog.sqlite3
export BESPOKE_DEMO_PUBLISHED=/data/demo/published-release-1
export BESPOKE_DEMO_HOSTS=demo.example.org
uv run uvicorn bespoke_bio.demo.service:app --host 127.0.0.1 --port 8789 --workers 2
```

The snapshot contains normalized results and checksummed copies of their structures, bundle
JSON and ZIP. Startup verifies the catalog identity and every artifact. It never copies
private job receipts; export fails for incomplete jobs or an existing destination. The
public process rejects live-service configuration and every write request. Its OpenAPI
schema contains only reads. GET requests support cross-origin access without credentials.

| Endpoint | Contents |
| --- | --- |
| `/api/catalog` | Full counts, source hashes and availability flags |
| `/api/targets?q=TP53&changed=true&limit=25&offset=0` | Paged summaries; continue through `total` |
| `/api/targets/{target_id}` | Original complete protein record, sequence and provenance |
| `/api/targets/{target_id}/structure` | Exact-sequence preview metadata, or null |
| `/api/structures/{file}` | Published coordinates, full confidence data or provenance |
| `/api/jobs` | Explicitly published jobs only |
| `/api/jobs/{job_id}/results` | All exported candidates and artifact/download links |
| `/docs`, `/openapi.json` | Executable API documentation and schema |

Run the public process as a persistent systemd service, with the data mount available,
automatic restart and bounded CPU/memory. Keep both Python services on loopback. Put an
HTTPS reverse proxy in front of port 8789, permit public TCP 80/443, and retain private SSH
access. The deployed VM uses an [Azure public DNS label](https://learn.microsoft.com/en-us/azure/virtual-network/ip-services/public-ip-addresses)
and [Caddy automatic HTTPS](https://caddyserver.com/docs/automatic-https); certificates renew
automatically while the domain resolves and the service is reachable. A minimal Caddyfile is:

```caddyfile
demo.example.org {
    encode zstd gzip
    redir /operator /operator/ 308
    handle_path /operator/* {
        basic_auth {
            operator REPLACE_WITH_CADDY_PASSWORD_HASH
        }
        header {
            Cache-Control no-store
            defer
        }
        reverse_proxy 127.0.0.1:8787 {
            header_up -Authorization
        }
    }
    handle {
        reverse_proxy 127.0.0.1:8789
    }
}
```

Use [Caddy's password hashing](https://caddyserver.com/docs/caddyfile/directives/basic_auth)
and HTTPS. Protect the entire `/operator/*` path, including APIs, artifacts and documentation;
never expose the gateway port directly. The operator client prefixes all requests and
download links. Enable `BESPOKE_DEMO_OPERATOR_ACCESS=1` on the public instance only after
the protected route exists. Configure the public hostname on both processes. Browser Basic
authentication is suitable for this shared demo; per-person access and audit require an
identity provider for broader use. Browsers may remember the operator login until closed.

### Key Vault credentials

Install the optional `vault` extra in the VM runtime. Store one JSON secret containing
`base_url`, `api_token`, `expires_at` (with timezone), and optionally `docs_url`. Grant the VM's
managed identity **Key Vault Secrets User** on this secret; the public process needs no
secret file. Follow Microsoft's [VM identity guidance](https://learn.microsoft.com/en-us/azure/key-vault/general/tutorial-python-virtual-machine).
In the private gateway's systemd unit, use `RuntimeDirectory=bespoke-demo-operator`,
`RuntimeDirectoryMode=0700`, and this startup command:

```bash
python -m bespoke_bio.demo.vault --vault-url https://YOUR-VAULT.vault.azure.net \
  --secret biocad-service --output /run/bespoke-demo-operator/service.json
```

Set `BIOCAD_SERVICE_CONFIG_FILE=/run/bespoke-demo-operator/service.json`. The helper writes
an atomic mode-600 runtime file, and the gateway enforces expiry before contacting the
service. Set the Vault secret's expiry metadata as well. After rotating the secret, restart
the operator service to retrieve its latest version. This intentionally avoids a Vault call
for every UI request. Never put secret values in commands, Git, frontend variables or logs.
The separate operator password may also be stored in Vault; Caddy receives only its hash.

Public reads use a 60-second HTTP cache and bounded per-worker catalog caches. Saved files
are served directly from disk, without contacting Jacob's API. Hashed assets are immutable;
retain previous asset files during deployments so open tabs can finish loading. When switching
to an installed release, clear any old `PYTHONPATH` override and verify the served catalog
metadata before completing the switch. The index
and snapshot are immutable releases: create new paths and restart workers when updating.
The VM must remain running for the URL to work; no development tunnel is involved.

## Publish target structures

A partial preview library is separate from the complete Step 2 bundle contract. It never
marks unresolved sequences as folded or removes targets from the catalog. Identical full
sequences share coordinates; every biological record keeps its identity and original JSON.

The publication input is a version-1 JSON manifest with `catalog_sha256` and a `structures`
array. Each entry supplies `amino_acid_sequence`, `sequence_sha256`, `file`, `file_sha256`,
`confidence: {file, sha256}`, `source`, `mean_plddt` (0–100 or null), and a `provenance` object.
Paths are relative to that manifest. Coordinates must be mmCIF; confidence is JSON or NPZ.
Retain model revision, inference settings, source URLs and validation evidence in provenance.

```bash
uv run python -m bespoke_bio.demo.structures \
  --index /data/demo/catalog.sqlite3 \
  --source /data/reviewed-structures/manifest.json \
  --output /data/demo/structures-release-1
export BESPOKE_DEMO_STRUCTURES=/data/demo/structures-release-1
# Restart the public and/or private service with this setting.
```

Publication checks full sequence hashes, ordered residues, backbone coverage, finite
coordinates and artifact hashes, then builds an immutable target index. Startup rechecks
catalog identity, index cardinality and every artifact. This verifies data integrity; it is
not an independent measurement of prediction accuracy. The optional provenance field
`peptide_geometry_screen_flags` exposes a nonnegative count of flagged bond distances in the
UI; null means no screen result was supplied. Review confidence semantics and structural
quality before publication. Source artifacts and generated releases stay outside Git.

The **Proteins** tab defaults to **Structure ready** when a library exists. Uncheck it to browse the
full catalog, including records without coordinates. `/api/targets?ready=true` exposes the
same filter. `/api/catalog` reports both complete and structure-ready counts. These counts
measure coverage, not confidence or biological significance. List order is alphabetical,
not a scientific ranking. Downloads retain full confidence arrays and provenance; the browser
loads only the selected model. Structures use immutable, content-hashed URLs.

## React application

`web/` contains Vite, React and strict TypeScript source. A single workspace places the
protein catalog beside the molecular viewer, with candidates below it. On mobile the catalog
opens in a drawer. History, sequences, provenance, metrics and design settings stay available
without occupying the main view. Components separate catalog, target preparation, job status,
candidates and the Mol* viewer. Zod validates API responses;
TanStack Query handles cancellation, caching and polling. Mutations never retry automatically.
`npm --prefix web run dev` serves the app on port 5173 and proxies `/api` to the gateway.

The production build writes ignored assets into `src/bespoke_bio/demo/static/`. Build before
packaging a wheel or deploying. The gateway compresses and caches hashed assets. Private API
responses remain uncached; public snapshots cache reads. Mol* loads on demand from the same server;
no CDN is needed for the application.

Existing mmCIF/PDB coordinates supply surface, ribbon and atom representations. Gaussian
surfaces request GPU computation and fall back to CPU when Mol* determines the device or
structure is unsuitable. Appearance reports the selected path and coordinate counts. WebGL renders the resulting
geometry. Metallic finish controls material lighting only; it does not need another dataset
and is not a measured density map. Sequence-only records need coordinates before 3D viewing.

The public workspace opens on **Saved keys**, with evaluated benchmark constructs kept
separate from donor records. List order is not a scientific ranking. Selecting a donor
record replaces the benchmark with that donor's own context. The sidebar reopens that exact target's completed,
active or prepared job when available. Otherwise it displays that exact sequence's saved
structure, or explicitly says that no structure has been published.
Selection is preserved in the URL across refresh and back/forward navigation. Cross-target
job links cannot display a model under the wrong protein's heading. A delayed
preparation response cannot switch the active selection. The Target/Key/Complex buttons choose
the artifact; Target is the separate preview or refold when available, while Complex shows
the design service's predicted bound geometry; Surface/Ribbon/Atoms changes its representation. Reset and fullscreen controls
stay with the viewer. Chain labels/colors come from Mol*'s actual structure color legend;
they do not infer which chain is the target or binder.

## What is retained and what runs

- Indexing first verifies the complete portable Stage 1 package, then checks its hashes again.
  Each biological target retains its original JSON; identical sequences remain distinct
  targets. An existing index is never overwritten. Build a new index after a package changes.
- The browser pages through all records and fetches a full sequence only when selected.
  Prefix search accepts gene or transcript; exact target IDs also work. The non-reference
  filter compares sequence hashes. Pagination does not cap the catalog.
- **Prepare target** calls `POST /v1/designs`. It validates and prepares one protein;
  it does not run inference. Local receipts bind the job to its submitted record and package
  hash. Repeating the same selection/settings reuses the local receipt. After an ambiguous
  timeout, inspect the upstream job before retrying; do not delete receipts blindly.
- **Example** opens the existing cached receipt, or calls `POST /v1/designs/demo` once when
  none exists on the private gateway. Public snapshots never submit this request.
  This is the only automatic POST; it never starts inference or retries after a
  failure. It displays previously generated results, explicitly separate from donor records.
  All exported candidates remain available in the selector,
  including failed computational filters. Structures and ZIP downloads come from that job.
- GPU inference is disabled by default. Only set `BESPOKE_DEMO_ALLOW_RUN=1` after configuring
  and authorizing live execution with the service operator. Starting then requires an explicit
  **Generate keys** action and confirmation. For a new target this prepares and dispatches
  once; for an existing prepared job it dispatches that job. Ambiguous submissions are never
  automatically retried. The protein stays visible while progress is polled, and candidates
  appear only when returned. Viewing cached models never invokes `/run`.

This UI prepares one chosen target at a time; it does not dispatch the entire catalog to
GPUs. The service's 2-design, 75-residue defaults are small demo settings, not a validated
design campaign. Confidence scores and filter passes are model evidence, not measured binding
or affinity. An empty result is shown as empty, never replaced with a cached donor result.

## Verification and limits

```bash
uv run pytest
uv run ruff check .
npm --prefix web run test
npm --prefix web run build
npm --prefix web run format:check
# Against the running HG002 demo, without starting inference:
npm --prefix web run test:browser
# Read-only checks against the public deployment:
BESPOKE_DEMO_PUBLIC=1 BESPOKE_DEMO_TEST_URL=https://your-public-host \
npm --prefix web run test:browser
```

The tests cover source corruption, biological cardinality, exact handoff, unknown submissions,
authentication failures, disabled execution, cross-origin requests and job-bound artifacts.
The gateway's own schema is generated at `/openapi.json`. Test the deployed service separately
by preparing a real non-reference record, comparing its saved target sequence/hash, opening
the cached candidates and downloading their structures. Do not call `/run` in integration tests.

Authenticated preparation and two explicitly authorized A1BG design jobs completed on
2026-10-18. Their four candidates retain exact input identity and native filter outcomes.
This is computational integration evidence, not experimental binding validation. Browser
tests intercept paid actions and exercise the UI without starting additional GPU jobs.

HTTP requests have a 60-second timeout; JSON responses are limited to 8 MiB and individual
artifact/ZIP responses to 128 MiB. Larger jobs need streaming downloads. Unknown submissions
and restarts during dispatch require operator reconciliation; the gateway never retries POSTs
automatically. The [Mol* viewer](https://molstar.org/docs/v5/plugin/instance/) is pinned in
`web/package-lock.json`. Without WebGL, original structure downloads remain available.
Browser tests cover the real cached mmCIF, representation and candidate changes, viewer
teardown, fullscreen, mobile drawers, selection identity, late preparation responses, empty
results, refresh/back navigation and absence of inference calls. Public snapshot tests cover
corruption, offline upstream operation, exact record retention and rejection of writes.
Headless tests can exercise the GPU algorithm through
software WebGL; they do not measure physical GPU performance.

For the demo, a verified subset is a coverage decision, not evidence that every protein was
modeled or that the models are equally reliable. HG002 is public benchmark data, not a soldier
or a sensor stream. Candidate isoforms/phase scenarios are not proof of protein expression.
The current flow demonstrates sequence identity → predicted target → candidate key/complex;
variant-specific binding, off-target selectivity and therapeutic benefit remain unestablished.
ESMFold2's [model card](https://huggingface.co/biohub/ESMFold2) explicitly distinguishes reduced
inference compute from maximal accuracy. BoltzGen's [pipeline outputs](https://github.com/HannesStark/boltzgen#pipeline-output)
provide computational design evidence, not measured binding.

## Saved key batches

The public site opens on **Saved keys**. Each sidebar row represents one published job;
identical gene names never merge different constructs or jobs. The selected complex opens in
3D with both target and key visible; Target/Key/Complex buttons switch structures. **Proteins** retains the
complete donor catalog and links completed jobs by exact target ID. A fold alone is labeled
“Protein structure · no saved keys”. Public browsing never starts inference.

`GET /api/designs?q=...&offset=...&limit=25` pages completed published jobs.
`GET /api/designs/{job_id}` returns display metadata and provenance. Native failures remain
visible. Native candidate `target_structure_url` is a chain extracted from its original complex.
Evaluated `predictions` instead link independent target and binder refolds plus a separately
predicted complex, indexed by exact target, candidate, context and seed.

Import Jacob's allowlisted public scientific export after copying it outside Git:

```bash
uv run python -m bespoke_bio.demo.import_science \
  --source /data/jacob/examples/folding_validation \
  --snapshot /data/current-published --output /data/new-published \
  --revision FULL_SOURCE_COMMIT
```

The importer validates every indexed CIF and confidence artifact checksum, coordinate
finiteness, backbone completeness, exact binder sequences, and consistent target constructs.
It imports every native candidate, including failures, and preserves the complete public
study in the ZIP download. Standard refolds, controls and seed repeats remain separate
study evidence; they are not counted as new candidate keys. Native pass status is not proof
of binding, and target processing, oligomers and competing partners still matter.

The website serves this immutable snapshot from Azure; no browser token or live GPU is
needed to view imported results. Generated files stay outside Git. Replace the public
snapshot setting only after validation, keeping the previous directory for rollback.

For the three reevaluated candidates, **Standard refold** is selected by default with seed
11; seeds 22/33 and **Original design** remain available in **Prediction**. The Target view uses the matching
MSA/no-template control; Key uses the independent binder refold. Complex loads exactly one
CIF. The single Target/Key/Complex selector switches complete artifacts; it never assembles
separate coordinate frames or changes their positions. The native filter label stays separate from
view-specific computational confidence; the complete evidence JSON retains all metrics.

Study IDs such as FGF2 and IL1RN are benchmark constructs, not HG002 target IDs. The importer
checks exact sequences across native and evaluated structures and rejects ambiguous or
unmatched prediction joins. It does not attach a benchmark to a donor based on gene name.
