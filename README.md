# Bespoke Bio

Bespoke Bio is the complete DNA/RNA-to-designed-protein system. This repository currently
implements Step 2: translate the resolved DNA/RNA findings from Step 1, resolve or predict
one structure per unique protein, and package the results for Step 3.

The implementation is intentionally small. It:

1. validates resolved coding sequences, expression measurements, and CDS-relative variants;
2. accepts only explicit haplotypes or phase scenarios supplied by Step 1;
3. translates, hashes, and deduplicates personalized amino-acid sequences;
4. reuses exact structures from RCSB PDB or AlphaFold DB when available;
5. folds every remaining sequence with `colabfold_batch`; and
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

```bash
uv run python -m bespoke_bio.api examples/case.json outputs/demo
```

The command looks for exact existing structures first, then sends every unresolved unique
sequence to `colabfold_batch`. The command fails if any target still lacks a structure.
ColabFold controls its own CPU, GPU, or remote-MSA configuration; this program does not select
hardware. A compatible `colabfold_batch` executable must be available on `PATH` whenever an
exact existing structure cannot be reused.

## Input boundary

Step 1 is responsible for alignment, variant calling and filtering, transcript/CDS
resolution, ORF selection, expression quantification, and haplotype/scenario grouping. It
sends Step 2 one JSON document containing:

- the pinned genome build and annotation release;
- the resolved coding sequence for each expressed transcript;
- variants expressed in 1-based CDS coordinates;
- transcript-level expression;
- an explicit haplotype or phase scenario for every variant set; and
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
