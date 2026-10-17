# Bespoke Bio

Bespoke Bio's middle layer converts expressed, variant-bearing coding transcripts into
amino-acid sequences. Those target proteins become the inputs to structure prediction
and key generation.

The first version is intentionally small:

1. Accept a versioned coding sequence, expression measurement, and normalized variants.
2. Validate and apply the variants using 1-based CDS coordinates.
3. Translate the personalized CDS and deduplicate identical protein sequences.
4. Return stable sequence hashes and explicit translation warnings.

Folding, binding-surface selection, and key generation belong to the downstream Step 3
system and are intentionally absent from this service's response contract.

## Run locally

```bash
uv sync --extra dev
uv run uvicorn bespoke_bio.api:app --reload
```

Open `http://127.0.0.1:8000/docs`, or submit the included example:

```bash
curl -X POST http://127.0.0.1:8000/v1/translate \
  -H 'content-type: application/json' \
  --data @examples/case.json
```

## Input boundary

The upstream data pipeline owns read alignment, variant calling, rarity filtering, and
transcript quantification. It sends this service:

- normalized rare variants in genomic coordinates;
- transcript-level expression findings;
- a genome build and annotation release;
- haplotype/phase information when available; and
- assembled transcript sequences for novel isoforms, when available.

This service owns transcript/CDS resolution, application of variants to the correct CDS,
and translation into protein sequences. The current `/v1/translate` endpoint is the tested
translation kernel and accepts a resolved CDS; the Step 1 bundle adapter will sit in front
of it.

Inside the translation kernel, variants use VCF-style replacement semantics: at
`cds_position`, replace `ref` with `alt`. Coordinates refer to the original reference CDS,
even when earlier indels change the personalized sequence length.

## Test

```bash
uv run pytest
uv run ruff check .
```

This software produces research predictions, not clinical results.

See [docs/architecture.md](docs/architecture.md) for the end-to-end system,
software choices, databases, and ownership boundaries.

