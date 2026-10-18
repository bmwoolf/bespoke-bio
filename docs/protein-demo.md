# Sequence-only protein design demo

This demo connects a verified Stage 1 protein package to Jacob's design API. It preserves
every biological target and sends the selected record unchanged, plus `binder_length` and
`num_designs`. Expression stays null when unmeasured. The complete sequence and both hashes
remain in the request. The standalone CDS-based Step 2 is bypassed: structural prediction
belongs to the design service in this route. This does not demonstrate the original three
modules running in sequence.

The integration targets [design service PR #1](https://github.com/bmwoolf/bespoke-bio/pull/1),
commit `95507e664491bb6f1986704e2da1771cf817e525`. Its executable API/schema is the upstream
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
demo, not an authenticated public application. Existing raw-read jobs can continue. Limit
indexing/testing to spare CPU capacity; browsing and preparation do not need a GPU.
Use a new private runs directory when switching to another design service; receipts belong
to the service that created them.

## React application

`web/` contains Vite, React and strict TypeScript source. A single workspace places the
protein catalog beside the molecular viewer, with candidates below it. On mobile the catalog
opens in a drawer. History, sequences, provenance, metrics and design settings stay available
without occupying the main view. Components separate catalog, target preparation, job status,
candidates and the Mol* viewer. Zod validates API responses;
TanStack Query handles cancellation, caching and polling. Mutations never retry automatically.
`npm --prefix web run dev` serves the app on port 5173 and proxies `/api` to the gateway.

The production build writes ignored assets into `src/bespoke_bio/demo/static/`. Build before
packaging a wheel or deploying. The gateway compresses and caches hashed assets, while API
responses remain uncached. Mol* loads on demand from the same server; no CDN is needed.

Existing mmCIF/PDB coordinates supply surface, ribbon and atom representations. Gaussian
surfaces request GPU computation and fall back to CPU when Mol* determines the device or
structure is unsuitable. Appearance reports the selected path and coordinate counts. WebGL renders the resulting
geometry. Metallic finish controls material lighting only; it does not need another dataset
and is not a measured density map. Sequence-only records need coordinates before 3D viewing.

The workspace opens on the explicitly labeled saved PD-L1 example and its first exported
candidate; list order is not a scientific ranking. Selecting a donor record replaces that
example with the donor's own context and an empty viewer until results exist. A delayed
preparation response cannot switch the active selection. The Complex/Binder selector chooses
the artifact; Surface/Ribbon/Atoms changes its representation. Reset and fullscreen controls
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
  none exists. This is the only automatic POST; it never starts inference or retries after a
  failure. It displays previously generated results, explicitly separate from donor records.
  All exported candidates remain available in the selector,
  including failed computational filters. Structures and ZIP downloads come from that job.
- GPU inference is disabled by default. Only set `BESPOKE_DEMO_ALLOW_RUN=1` after configuring
  and authorizing live execution with the service operator. Starting then requires an explicit
  **Generate binders** action and confirmation. Neither preparing nor viewing a cached
  example invokes `/run`.

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
```

The tests cover source corruption, biological cardinality, exact handoff, unknown submissions,
authentication failures, disabled execution, cross-origin requests and job-bound artifacts.
The gateway's own schema is generated at `/openapi.json`. Test the deployed service separately
by preparing a real non-reference record, comparing its saved target sequence/hash, opening
the cached candidates and downloading their structures. Do not call `/run` in integration tests.

Without the shared service token, authenticated shared-service behavior remains unverified.
Without a fresh GPU run, deployment, inference and cleanup remain unverified. Cached replay
demonstrates the result path only. No timing here predicts the cost of processing all proteins.

HTTP requests have a 60-second timeout; JSON responses are limited to 8 MiB and individual
artifact/ZIP responses to 128 MiB. Larger jobs need streaming downloads. Unknown submissions
and restarts during dispatch require operator reconciliation; the gateway never retries POSTs
automatically. The [Mol* viewer](https://molstar.org/docs/v5/plugin/instance/) is pinned in
`web/package-lock.json`. Without WebGL, original structure downloads remain available.
Browser tests cover the real cached mmCIF, representation and candidate changes, viewer
teardown, fullscreen, mobile drawers, selection identity, late preparation responses, empty
results and absence of inference calls. Headless tests can exercise the GPU algorithm through
software WebGL; they do not measure physical GPU performance.
