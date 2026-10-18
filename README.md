# Bespoke Bio

Bespoke Bio is the complete DNA/RNA-to-designed-protein system. This is the shared repository
for Steps 1, 2, and 3. Step 2 is the first implemented contribution: it translates the resolved
DNA/RNA findings from Step 1, resolves or predicts one structure per unique protein, and
packages the results for Step 3.

The Step 2 implementation:

1. validates resolved coding sequences, expression measurements, and CDS-relative variants;
2. accepts only explicit haplotypes or phase scenarios supplied by Step 1;
3. translates, hashes, and deduplicates personalized amino-acid sequences;
4. reuses exact structures from RCSB PDB or AlphaFold DB when available;
5. folds every remaining sequence with a selected Boltz-2, ESMFold, or OpenFold backend; and
6. writes `proteins.fasta`, `targets.jsonl`, structure files, and `provenance.json` for Step 3.

Binding-surface selection and key generation belong to Step 3.

## Run Step 1

Step 1 imports precomputed human variant evidence and exports complete protein sequences in
`protein-targets.jsonl` and a deduplicated `proteins.fasta`. It retains every supported target
and records why other transcripts were deferred. Expression remains null until measured.

```bash
uv sync --extra stage1 --extra dev
uv run bespoke-stage1 process-precomputed --settings /work/run.json --output /work/attempt-001
uv run bespoke-stage1 validate-proteins --root /work/attempt-001/proteins
```

See [the Stage 1 guide](docs/stage-one.md) for inputs, settings, container requirements and the
protein JSON contract. Keep runtime settings and generated data outside Git. The HG002
precomputed-evidence route has passed a full export and independent replay. The raw DNA/RNA
workflow adapters still require full-case qualification. Step 2's current API accepts CDS and
variants; direct protein ingestion requires a separate interface update.

## Run the sequence-only design demo

The Vite/React/TypeScript [protein demo](docs/protein-demo.md), with a Mol* viewer, browses
a complete verified Step 1 export and sends one unchanged protein record to the separate
design service. It shows preparation status,
cached or newly generated candidates, and structures. GPU inference is disabled by default.
Keep the full dataset and runtime files on the VM, outside Git. This demo bypasses standalone
Step 2; the original CDS-to-structure API below remains available.

## Run the API

```bash
uv sync --extra dev
BESPOKE_BIO_OUTPUT_DIR=outputs \
uv run uvicorn bespoke_bio.api:app --reload
```

Open `http://127.0.0.1:8000/docs`, or run the complete Step 2 pipeline with ESMFold:

```bash
curl -X POST 'http://127.0.0.1:8000/v1/build?folding_backend=esmfold' \
  -H 'content-type: application/json' \
  --data @examples/case.json
```

The response includes translated proteins, completed structures, and a `bundle_id`. Any package
file can then be downloaded from
`GET /v1/bundles/{bundle_id}/{file}`. The API and CLI both call the same `build_bundle()` core
function.

For API deployments, set model-specific arguments once in `BESPOKE_BIO_BACKEND_ARGS`. For
example, a Boltz-2 server can use
`--use_msa_server --msa_server_url https://msa.example.org`. The CLI accepts the same arguments
explicitly with `--backend-arg`.

## Build the Step 3 handoff

```bash
uv run python -m bespoke_bio.api examples/case.json outputs/demo \
  --folding-backend esmfold
```

The command looks for exact existing structures first, then sends every unresolved unique
sequence to the selected model. It fails if any target still lacks a structure. Model-specific
arguments can be repeated with `--backend-arg=VALUE`; hardware configuration remains outside
the biological pipeline.

Supported backends:

| Backend | Required runtime | Input prepared by Bespoke Bio | Output collected |
| --- | --- | --- | --- |
| `boltz2` | `boltz` on `PATH` | one YAML file per sequence | mmCIF and confidence JSON |
| `esmfold` | `esm-fold` on `PATH` | one deduplicated FASTA batch | PDB with mean pLDDT |
| `openfold` | an OpenFold checkout and databases | one FASTA file per sequence | PDB with mean pLDDT |

The Boltz-2 adapter requires an MSA service explicitly. The public service is suitable for a
small demo; a 1,000-protein run should point at a separately operated service:

```bash
uv run python -m bespoke_bio.api examples/case.json outputs/demo \
  --folding-backend boltz2 \
  --backend-arg=--use_msa_server \
  --backend-arg=--msa_server_url \
  --backend-arg=https://msa.example.org
```

OpenFold's official inference script also requires its template and sequence databases. The
adapter locates the script and template directory from the environment and forwards the
remaining official OpenFold arguments:

```bash
export BESPOKE_BIO_OPENFOLD_SCRIPT=/opt/openfold/run_pretrained_openfold.py
export BESPOKE_BIO_OPENFOLD_TEMPLATE_MMCIF_DIR=/data/pdb_mmcif/mmcif_files

uv run python -m bespoke_bio.api examples/case.json outputs/demo \
  --folding-backend openfold \
  --backend-arg=--use_precomputed_alignments \
  --backend-arg=/data/alignments \
  --backend-arg=--model_device \
  --backend-arg=cuda:0
```

The model packages and their large weights or databases are intentionally not Python
dependencies of this small service. A local environment or GPU container supplies the selected
runtime.

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
Coding assistants should start with [AGENTS.md](AGENTS.md).
