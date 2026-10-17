# Bespoke Bio

Bespoke Bio is the complete DNA/RNA-to-designed-protein system. The current implementation
focuses on Step 2, the translation and folding module: it converts expressed, variant-bearing
transcripts into amino-acid sequences and folded target structures for Step 3.

The first version is intentionally small:

1. Accept a versioned coding sequence, expression measurement, and normalized variants.
2. Validate and apply the variants using 1-based CDS coordinates.
3. Translate the personalized CDS and deduplicate identical protein sequences.
4. Return stable sequence hashes and explicit translation warnings.
5. Fold each unique personalized protein and preserve structure confidence and provenance.

Binding-surface selection and key generation belong to Bespoke Bio's Step 3 module. Folding
belongs to Step 2 and its structures form part of the Step 2 response contract.

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

Bespoke Bio's Step 1 module performs read alignment, variant calling, rarity filtering,
and transcript quantification. It sends the Step 2 module:

- normalized rare variants in genomic coordinates;
- transcript-level expression findings;
- a genome build and annotation release;
- haplotype/phase information when available; and
- assembled transcript sequences for novel isoforms, when available.

The Step 2 module performs transcript/CDS resolution, application of variants to the correct
CDS, translation into protein sequences, and folding those proteins. The current `/v1/translate`
endpoint is the tested translation kernel and accepts a resolved CDS; the Step 1 bundle
adapter and folding worker will complete the Step 2 pipeline around it.

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
software choices, databases, and step boundaries.

