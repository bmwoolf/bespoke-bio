# System architecture

## Goal

Convert every expressed rare-variant transcript received from the upstream data pipeline
into a personalized target protein, predict its structure, and generate ranked candidate
protein binders ("keys").

The system has three ownership boundaries:

```text
Step 1: sequence processing
    -> expressed, variant-bearing CDS records

Step 2: Bespoke Bio middle layer
    -> personalized proteins, structures, and candidate keys

Step 3: catalog UI and experimental workflow
    -> searchable results and lab feedback
```

## Step 1: upstream sequence processing

The data team owns alignment, variant calling, rarity filtering, transcript selection,
expression quantification, phasing when available, and preparation of spliced coding
sequences.

Suggested reference data:

| Purpose | Database |
| --- | --- |
| Human reference assembly | GRCh38 |
| Gene, transcript, and CDS annotation | GENCODE, with an exact release pinned |
| Default representative transcript | MANE Select |
| Population allele frequency | gnomAD v4.1 GRCh38 |

Suggested software:

| Task | Software |
| --- | --- |
| DNA alignment | BWA-MEM2 |
| Germline variant calling | DeepVariant |
| Variant normalization | bcftools |
| Transcript/protein consequence annotation | Ensembl VEP |
| RNA alignment | STAR |
| Transcript-level expression | Salmon |

The handoff is one validated JSON record per case. Each transcript includes a versioned
transcript ID, spliced CDS, TPM, variants in one-based CDS coordinates, and haplotype
information when available. Raw reads and BAM files do not cross this boundary.

```text
POST /v1/translate
Content-Type: application/json

examples/case.json
```

## Step 2A: variant to personalized protein

Software already implemented:

- Python 3.12
- FastAPI for the service boundary
- Pydantic for strict input and output validation
- A small tested standard-code translator in `bespoke_bio.translation`
- pytest and Ruff for verification

Processing:

1. Validate transcript, expression, and variant fields.
2. Apply non-overlapping variants using original CDS coordinates.
3. Translate the personalized CDS to amino acids.
4. Report missing starts, partial codons, ambiguous codons, internal stops, and unknown
   phasing explicitly.
5. Deduplicate identical proteins by SHA-256 sequence hash.

Expression is retained as evidence and ranking metadata. It is not interpreted as protein
abundance.

## Step 2B: target structure

Resolve structures in this order:

1. **RCSB PDB** for a suitable experimental structure with a verified sequence mapping.
2. **AlphaFold Protein Structure Database** for an unchanged reference protein mapped by
   UniProt accession.
3. **ColabFold 1.6.3** for a novel or variant amino-acid sequence.

All lookups and model runs are cached by amino-acid sequence hash. Store coordinates as
mmCIF, plus provenance, model version, pLDDT, PAE, and the exact input sequence. Structure
predictions are research outputs, not experimentally validated structures.

## Step 2C: lock selection

For the first implementation, the lock is a variant-centered surface patch:

1. Map the changed amino-acid position onto the structure.
2. Confirm that the residue or nearby residues are surface-accessible.
3. Select nearby exposed residues as candidate hotspot residues.
4. Preserve domain and active-site annotations from UniProt/InterPro when available.

This keeps lock selection deterministic and inspectable. More sophisticated pocket and
functional-site models can be added later.

## Step 2D: key generation

Use **BindCraft** for the first complete key-generation worker. It accepts a target
structure and target residues and provides an integrated binder-design pipeline based on
AlphaFold2 optimization, ProteinMPNN, and PyRosetta.

Each design run must record:

- target sequence and structure hashes;
- selected hotspot residues;
- software commit and model weights;
- random seed and generation settings;
- generated key sequence;
- target-key complex structure; and
- confidence and interface metrics.

Do not represent an in-silico score as measured binding affinity. The API returns ranked
candidate keys for experimental testing.

If BindCraft proves too restrictive, replace that worker with the more modular sequence:

```text
RFdiffusion -> ProteinMPNN -> ColabFold-Multimer -> interface filtering
```

The public API contract does not change when the key-generation backend changes.

## Step 2 output

```text
result.json
targets.fasta
keys.fasta
targets/<target-id>.cif
keys/<key-id>.cif
complexes/<target-id>_<key-id>.cif
provenance.json
```

The JSON response carries expression, source variants, sequence hashes, model confidence,
interface scores, warnings, and ranked key IDs.

## Step 3: catalog and lab feedback

Suggested software:

- Next.js for the catalog UI
- Mol* for interactive structures
- The FastAPI service for result retrieval

For the hackathon, keep metadata in JSON and structure files on disk. Add PostgreSQL and
S3-compatible object storage only when concurrent users or durable hosted jobs require it.

Experimental measurements should be appended as new evidence; they must not overwrite the
original predictions.

## Runtime shape

Translation remains synchronous and fast. Folding and key design run on a Linux NVIDIA GPU
worker and expose job status:

```text
POST /v1/jobs
GET  /v1/jobs/{job_id}
GET  /v1/jobs/{job_id}/result
```

Start with a single worker and a file-backed job directory. Add Redis or a distributed queue
only after one worker becomes a demonstrated limitation.

## Later pathogen mode

Pathogen support changes the input adapter, not the middle-layer contract. Add a taxon ID,
reference accession, organism-specific genetic-code table, and pathogen protein databases
from NCBI RefSeq/GenBank or UniProt Proteomes. Translation, structure resolution, lock
selection, key generation, and ranking remain the same.

Add host-selectivity checks before advancing pathogen keys: compare candidate keys against
the human proteome and explicitly track predicted human off-targets.

## Primary resources

- GENCODE: https://www.gencodegenes.org/
- MANE: https://www.ncbi.nlm.nih.gov/refseq/MANE/
- gnomAD: https://gnomad.broadinstitute.org/
- Ensembl VEP: https://www.ensembl.org/info/docs/tools/vep/
- UniProt: https://www.uniprot.org/
- RCSB PDB APIs: https://www.rcsb.org/docs/programmatic-access/web-apis-overview
- AlphaFold DB: https://alphafold.ebi.ac.uk/
- ColabFold: https://github.com/sokrypton/ColabFold
- BindCraft: https://github.com/martinpacesa/BindCraft
- RFdiffusion: https://github.com/RosettaCommons/RFdiffusion
- ProteinMPNN: https://github.com/dauparas/ProteinMPNN
- Mol*: https://molstar.org/
