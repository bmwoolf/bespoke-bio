# Bespoke Bio

Bespoke Bio is the complete DNA/RNA-to-designed-protein system. This repository currently
implements Step 2: translate the resolved DNA/RNA findings from Step 1, resolve or predict
one structure per unique protein, and package the results for Step 3.

The implementation is intentionally small. It:

1. validates resolved coding sequences, expression measurements, and CDS-relative variants;
2. generates possible proteins when heterozygous variants are unphased;
3. translates, hashes, and deduplicates personalized amino-acid sequences;
4. optionally retrieves exact structures from RCSB PDB or AlphaFold DB;
5. optionally folds unresolved sequences with a local `colabfold_batch`; and
6. writes `proteins.fasta`, `targets.jsonl`, structure files, and `provenance.json` for Step 3.

Binding-surface selection and key generation belong to Step 3.

## Run the translation API

```bash
uv sync --extra dev
uv run uvicorn bespoke_bio.api:app --reload
```

Open `http://127.0.0.1:8000/docs`, or submit the example:

```bash
curl -X POST http://127.0.0.1:8000/v1/translate \
  -H 'content-type: application/json' \
  --data @examples/case.json
```

## Build the Step 3 handoff

The default command stays offline and marks structures as unresolved:

```bash
uv run python -m bespoke_bio.api examples/case.json outputs/demo
```

Add exact-sequence lookup or fold unresolved proteins with an installed ColabFold worker:

```bash
uv run python -m bespoke_bio.api examples/case.json outputs/demo --lookup
uv run python -m bespoke_bio.api examples/case.json outputs/demo --lookup --fold-missing
```

`--lookup` sends amino-acid sequences to RCSB PDB. A default ColabFold installation may
send sequences to its public MSA service. Do not enable either option for sensitive data
unless the relevant data-handling policy permits it; use local databases instead.

## Input boundary

Step 1 is responsible for alignment, variant calling and filtering, transcript/CDS
resolution, ORF selection, expression quantification, and phasing when available. It sends
Step 2 one JSON document containing:

- the pinned genome build and annotation release;
- the resolved coding sequence for each expressed transcript;
- variants expressed in 1-based CDS coordinates;
- transcript-level expression;
- haplotype/phase information when available; and
- an optional UniProt accession for exact AlphaFold DB lookup.

Variants use replacement semantics: at `cds_position`, replace `ref` with `alt`. Coordinates
refer to the original CDS even when an earlier indel changes the personalized sequence length.

## Verify

```bash
uv run pytest
uv run ruff check .
```

This software produces research predictions, not clinical results.

See [docs/architecture.md](docs/architecture.md) for the system contracts and step boundaries.
