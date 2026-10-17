# Bespoke Bio agent guide

Read this file first, then read `docs/architecture.md`, `README.md`, and
`examples/case.json`. Treat the code and generated `/openapi.json` schema as the
authoritative executable contract. Update the documentation and tests whenever that contract
changes.

## System boundaries

Bespoke Bio contains three steps:

1. Step 1 processes DNA and RNA. It owns alignment, variant calling and filtering,
   transcript/CDS resolution, ORF selection, expression quantification, and phase scenarios.
2. Step 2 is implemented in this repository. It validates Step 1 output, applies variants to
   resolved CDS sequences, translates proteins, resolves or predicts structures, and packages
   the results.
3. Step 3 consumes folded targets to select binding surfaces and design candidate keys.

Do not move Step 1 cleanup or Step 3 protein-design behavior into Step 2.

## Step 2 invariants

- Step 1 supplies coding-orientation CDS sequences and one-based CDS-relative variants.
- Variant `ref` must exactly match the supplied reference CDS.
- Do not infer phase. A variant-bearing record with an unknown haplotype needs an explicit
  `phase_scenario`.
- Translate every supplied transcript record regardless of expression level. Expression is
  provenance and ranking evidence, not a direct abundance measurement.
- Deduplicate folding work by the SHA-256 hash of the complete amino-acid sequence while
  retaining every biological target record.
- Reuse RCSB PDB or AlphaFold DB structures only after an exact full-sequence match.
- Every unique amino-acid sequence must end with a structure. Fail the run if any remains
  unresolved.
- Folding is provider-agnostic. Supported backends are `boltz2`, `esmfold`, and `openfold`.
- Do not hardcode cloud providers, GPU types, worker counts, or private-data assumptions into
  the biological pipeline.
- This software produces research predictions, not clinical results.

## Shared execution path

Both transports call `build_bundle()` in `src/bespoke_bio/api.py`:

- HTTP: `POST /v1/build?folding_backend=...`
- CLI: `python -m bespoke_bio.api INPUT OUTPUT --folding-backend ...`

The HTTP response includes a `bundle_id`. Package files are available from
`GET /v1/bundles/{bundle_id}/{file}`. FastAPI publishes the live schema at
`GET /openapi.json`.

The handoff package contains:

```text
proteins.fasta
targets.jsonl
provenance.json
targets/<sequence-hash>.cif or .pdb
targets/<sequence-hash>.*.json
```

Step 3 should consume `targets.jsonl` and the referenced structure files. It must not rely on
temporary model-specific working files.

## Folding backends

Each backend accepts unique `Protein` records and returns one normalized `Structure` result per
sequence hash. Backend-specific model packages, weights, databases, MSA services, hardware, and
parallel execution belong to the runtime environment.

- Boltz-2 uses `boltz` and requires an explicitly configured MSA service.
- ESMFold uses `esm-fold` and accepts the generated multi-FASTA batch directly.
- OpenFold uses the official inference script configured by
  `BESPOKE_BIO_OPENFOLD_SCRIPT` and `BESPOKE_BIO_OPENFOLD_TEMPLATE_MMCIF_DIR`.

Add a new model by implementing the existing folding-runner contract and registering it in
`FOLDING_BACKENDS`. Do not fork the translation or packaging pipeline for a model.

## Development workflow

Keep the implementation small, readable, and dependency-light. Avoid speculative abstractions,
fixed biological cutoffs, and infrastructure that is not required by the current contract.

Before handing off a change, run:

```bash
uv run pytest
uv run ruff check .
```

Do not commit model weights, sequence databases, generated bundles, credentials, or local
runtime configuration.
