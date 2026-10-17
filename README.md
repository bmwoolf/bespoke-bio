# Bespoke Bio

Bespoke Bio's middle layer converts expressed, variant-bearing coding transcripts into
amino-acid sequences. Those target proteins become the inputs to structure prediction
and key generation.

The first version is intentionally small:

1. Accept a versioned coding sequence, expression measurement, and normalized variants.
2. Validate and apply the variants using 1-based CDS coordinates.
3. Translate the personalized CDS and deduplicate identical protein sequences.
4. Return stable sequence hashes and explicit translation warnings.

Folding and key-generation results are represented in the response contract, but are
not fabricated when a model backend is unavailable.

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

The upstream data pipeline owns read alignment, variant calling, transcript quantification,
and conversion from genomic coordinates to transcript/CDS coordinates. It sends this service:

- a genome build and annotation release;
- a versioned transcript ID and spliced coding sequence;
- transcript-level expression;
- variants expressed against that coding sequence; and
- haplotype/phase information when available.

Variants use VCF-style replacement semantics: at `cds_position`, replace `ref` with `alt`.
Coordinates refer to the original reference CDS, even when earlier indels change the
personalized sequence length.

## Test

```bash
uv run pytest
uv run ruff check .
```

This software produces research predictions, not clinical results.

