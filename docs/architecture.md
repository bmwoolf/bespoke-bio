# System architecture

## Step boundaries

```text
Step 1: DNA + RNA processing
    -> resolved transcript CDS, CDS-relative variants, expression, and phase

Step 2: translation + folding
    -> personalized proteins, structures, confidence, and provenance

Step 3: protein design
    -> lock selection and candidate keys
```

Bespoke Bio comprises all three steps. Step 1 owns genomic and transcript cleanup. Step 2
starts from resolved coding sequences and ends with folded target structures. Step 3 consumes
those structures and generates candidate keys.

## Step 1 -> Step 2

Step 1 performs alignment, variant calling and rarity filtering, transcript/CDS resolution,
novel-transcript ORF selection, expression quantification, and haplotype/scenario grouping.

Suggested Step 1 references and software:

| Responsibility | Software/database |
| --- | --- |
| Human reference | GRCh38 |
| Gene, transcript, and CDS annotation | GENCODE, pinned release |
| Representative transcripts | MANE Select |
| Population frequency | gnomAD, pinned release |
| DNA alignment and variant calling | BWA-MEM2, DeepVariant, bcftools |
| RNA alignment and expression | STAR, Salmon |
| Variant consequences and CDS coordinates | Ensembl VEP, pinned offline cache |
| Novel-transcript ORFs | TransDecoder |

Step 1 emits one JSON document. `examples/case.json` is executable against the current API.

```json
{
  "case_id": "patient-001",
  "reference": {
    "genome_build": "GRCh38",
    "annotation": "GENCODE-50"
  },
  "transcripts": [
    {
      "gene_id": "ENSG00000141510.19",
      "gene_symbol": "TP53",
      "transcript_id": "ENST00000269305.9",
      "uniprot_id": "P04637",
      "haplotype": "unknown",
      "phase_scenario": "candidate-1",
      "coding_sequence": "ATG...TAA",
      "expression": {
        "tpm": 42.7,
        "raw_count": 1834
      },
      "variants": [
        {
          "variant_id": "17:7675088:C:T",
          "cds_position": 524,
          "ref": "G",
          "alt": "A",
          "genotype": "0/1",
          "phase_set": null,
          "quality": "PASS"
        }
      ]
    }
  ]
}
```

Contract requirements:

- `coding_sequence` is the resolved CDS in coding orientation, including its start codon.
- `cds_position` is one-based against the supplied reference CDS.
- `ref` must match the reference CDS at that position.
- One record is supplied per known haplotype or explicit candidate phase scenario.
- A variant-bearing record with `haplotype: unknown` must include `phase_scenario`.
- Novel transcripts arrive with an already selected CDS; Step 2 does not predict their ORFs.
- `uniprot_id` is optional and is used only for exact AlphaFold DB lookup.

## Step 2: translation and folding

### Software

| Responsibility | Software/database |
| --- | --- |
| API and validation | Python, FastAPI, Pydantic |
| Translation and packaging | `bespoke_bio.api` |
| Experimental structure lookup | RCSB PDB Search and Data APIs |
| Reference prediction lookup | AlphaFold DB API |
| Missing-structure prediction | `colabfold_batch` |
| Verification | pytest, Ruff |

### Processing flow

```text
Step 1 JSON
    -> validate CDS, variants, expression, and reference versions
    -> apply each explicit haplotype/scenario to the reference CDS
    -> translate to amino-acid sequences
    -> deduplicate by SHA-256 sequence hash
    -> reuse an exact existing structure when available
    -> fold every unresolved unique sequence with ColabFold
    -> write the Step 3 handoff package
```

Each transcript is retained regardless of expression level. Expression is evidence and
ranking metadata, not a direct protein-abundance measurement.

### Translation and phasing

The service rejects reference mismatches, overlapping variants, variants beyond the CDS, and
variant sets without an explicit haplotype or phase scenario. It reports missing starts,
partial or ambiguous codons, and internal stops. Step 2 does not infer phase or generate
combinations; it translates each variant set supplied by Step 1 exactly once.

### Structures

Structure resolution runs once per unique amino-acid sequence:

1. Search RCSB PDB at 100% identity and verify the complete polymer sequence exactly.
2. If a UniProt accession was supplied, accept AlphaFold DB only when its complete sequence
   exactly matches the personalized protein.
3. Run `colabfold_batch` for every unresolved sequence.

RCSB and AlphaFold structures are stored as mmCIF. ColabFold structures are stored as PDB,
with its score JSON retained when available. Every target records the source, source ID,
file, confidence information, and sequence hash. Step 2 fails if any sequence remains without
a structure. ColabFold owns the CPU, GPU, or remote execution configuration; Step 2 only passes
it a deduplicated FASTA batch and collects its results.

## Step 2 -> Step 3

The package contains only the files Step 3 needs:

```text
proteins.fasta
targets.jsonl
provenance.json
targets/<sequence-hash>.cif          # exact RCSB or AlphaFold match, when found
targets/<sequence-hash>.pdb          # ColabFold result, when run
targets/<sequence-hash>.*.json       # confidence data, when available
```

Each `targets.jsonl` record contains the target identity, transcript and variant provenance,
expression, reference and personalized amino-acid sequences, phase scenario, sequence hash,
translation warnings, and the completed structure result:

```json
{
  "target_id": "target_f82a91",
  "gene_symbol": "TP53",
  "transcript_id": "ENST00000269305.9",
  "haplotype": "unknown",
  "phase_scenario": "possible_2_of_2",
  "source_variant_ids": ["17:7675088:C:T"],
  "expression": {"tpm": 42.7, "raw_count": 1834},
  "amino_acid_sequence": "MEEPQ...",
  "sequence_hash": "sha256:...",
  "warnings": [],
  "structure": {
    "status": "complete",
    "source": "ColabFold",
    "file": "targets/<sequence-hash>.pdb",
    "mean_plddt": 86.2,
    "confidence_file": "targets/<sequence-hash>.scores.json"
  }
}
```

Binding surfaces, key sequences, target-key complexes, interface scores, and rankings are
Step 3 outputs and are intentionally absent from the Step 2 package.

## Current implementation

The repository implements the complete Step 2 boundary described above:

- `POST /v1/translate` validates and translates explicit Step 1 scenarios.
- `python -m bespoke_bio.api INPUT OUTPUT` resolves or folds every unique sequence and writes
  the completed handoff package.

## Later pathogen mode

Pathogen support changes the Step 1 reference adapter, not the Step 2 contract. Step 1 must
also supply a taxon ID and use the organism's reference accession, transcript/CDS annotation,
and genetic code. Step 2 currently uses the standard nuclear genetic code.

## Primary resources

- GENCODE: https://www.gencodegenes.org/
- MANE: https://www.ncbi.nlm.nih.gov/refseq/MANE/
- gnomAD: https://gnomad.broadinstitute.org/
- Ensembl VEP: https://www.ensembl.org/info/docs/tools/vep/
- TransDecoder: https://github.com/TransDecoder/TransDecoder
- RCSB PDB: https://www.rcsb.org/
- AlphaFold DB: https://alphafold.ebi.ac.uk/
- ColabFold: https://github.com/sokrypton/ColabFold
