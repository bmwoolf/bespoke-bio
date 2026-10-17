# System architecture

## Step boundaries

```text
Step 1: DNA + RNA processing module
    -> rare genomic variants and transcript-expression findings

Step 2: Translation module + folding
    -> personalized amino-acid sequences with provenance + folded proteins

Step 3: Structure and protein-design module
    -> lock selection, and key generation
```

Bespoke Bio comprises all three modules. Step 1 produces DNA- and RNA-derived findings.
Step 2 converts those findings into validated protein sequences and folded target
structures. Step 3 consumes those structures to identify binding surfaces and generate
candidate keys.

## Step 1 -> Step 2

Step 1 performs read alignment, variant calling, rarity filtering, RNA alignment,
transcript/isoform quantification, and phasing when available.

Suggested upstream references:

| Purpose | Database |
| --- | --- |
| Human reference assembly | GRCh38 |
| Gene and transcript annotation | GENCODE, with an exact release pinned |
| Representative transcript labels | MANE Select |
| Population allele frequency | gnomAD v4.1 GRCh38 |

Suggested upstream software:

| Task | Software |
| --- | --- |
| DNA alignment | BWA-MEM2 |
| Germline variant calling | DeepVariant |
| Variant normalization | bcftools |
| RNA alignment | STAR |
| Transcript-level expression | Salmon |

Step 1 emits a processed-data bundle rather than raw reads:

```text
manifest.json
variants.vcf.gz
variants.vcf.gz.tbi
expression.tsv
novel-transcripts.fasta     # optional
novel-transcripts.gtf       # optional
```

### `manifest.json`

```json
{
  "case_id": "sample-001",
  "genome_build": "GRCh38",
  "annotation": "GENCODE-<pinned-release>",
  "variants": "variants.vcf.gz",
  "expression": "expression.tsv",
  "novel_transcripts_fasta": null,
  "novel_transcripts_gtf": null
}
```

### `variants.vcf.gz`

Required information:

- chromosome, position, reference allele, and alternate allele;
- sample genotype;
- phase set or haplotype when available;
- quality/filter status; and
- stable variant identifier when available.

All variants should be normalized against the exact GRCh38 reference named in the
manifest. A bare rsID or protein-change label is not enough to reconstruct a protein.

### `expression.tsv`

```text
transcript_id        gene_id             gene_symbol  tpm   raw_count
ENST00000269305.9    ENSG00000141510.19  TP53         42.7  1834
```

Transcript-level values are required whenever available. Gene-only expression cannot
determine which splice isoform, and therefore which protein, was expressed.

Novel RNA isoforms should include both their assembled sequence and transcript model.
The manifest must state whether an assembled sequence is already sample-specific so that
the same DNA variants are not applied twice.

### `1_2_api.json`
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
      "strand": "-",
      "expression": {
        "tpm": 42.7,
        "raw_count": 1834
      },
      "variants": [
        {
          "variant_id": "17:7675088:C:T",
          "chrom": "17",
          "position": 7675088,
          "ref": "C",
          "alt": "T",
          "genotype": "0/1",
          "haplotype": null,
          "phase_set": null,
          "population_af": 0.00002,
          "quality": "PASS"
        }
      ]
    }
  ]
}
```


## Step 2: Translation and folding module

### Software

| Responsibility | Software/database |
| --- | --- |
| API and validation | Python 3.12, FastAPI, Pydantic |
| Variant/transcript consequences | Ensembl VEP, run against a pinned offline cache |
| Known transcript and CDS sequences | GENCODE FASTA/GTF |
| Representative-transcript metadata | MANE Select |
| Novel-transcript ORF prediction | TransDecoder |
| Translation | Tested `/v1/translate` API kernel |
| Protein identity and metadata | UniProt |
| Experimental structure lookup | RCSB PDB |
| Reference structure predictions | AlphaFold DB |
| Variant and novel protein folding | ColabFold |
| Structure interchange format | mmCIF |
| Verification | pytest and Ruff |

VEP and GENCODE are part of Step 2 because Step 1 sends genomic variants rather than
pre-resolved CDS variants. The exact VEP cache, genome FASTA, GENCODE release, and MANE
release must be recorded for every run. Every folded structure must also record its source,
software version, model weights, exact amino-acid input, and confidence outputs.

### Processing flow

```text
variants.vcf.gz ------------------------┐
                                       v
                           VEP transcript mapping
                                       |
expression.tsv -> expressed transcripts+
                                       |
GENCODE CDS FASTA ----------------------+
                                       v
                       group variants by transcript/haplotype
                                       |
                                       v
                         apply variants to reference CDS
                                       |
                                       v
                           translate CDS into protein
                                       |
                                       v
                     deduplicate by amino-acid sequence hash
                                       |
                                       v
                    exact-sequence structure lookup/cache
                                       |
                                       v
                        fold unresolved sequences with ColabFold
                                       |
                                       v
                  targets.jsonl + proteins.fasta + targets/*.cif
```

For a known transcript, use its annotated GENCODE CDS. For a novel assembled transcript,
use TransDecoder to identify candidate coding regions and clearly mark the result as a
predicted ORF.

Each record is processed even when expression is low. Expression remains attached as
evidence and Step 3 ranking metadata; it is not interpreted as protein abundance.

### Translation checks

The service reports rather than hides:

- reference-allele mismatches;
- overlapping or out-of-range variants;
- unknown phasing;
- missing start codons;
- ambiguous codons;
- frameshifts and partial codons; and
- premature/internal stop codons.

Variants on a known haplotype are applied together before translation. Unphased variants
must be represented as uncertain rather than silently asserted to occur on one protein.

### Folding

Resolve a structure for each unique amino-acid sequence in this order:

1. Use RCSB PDB when an experimental structure has an exact, verified sequence match.
2. Use AlphaFold DB when the personalized sequence exactly matches a reference protein.
3. Run ColabFold when the sequence contains a variant or is novel.

Cache structures by amino-acid sequence hash so identical proteins are folded once. Store
the structure as mmCIF and retain pLDDT, PAE, model/version information, and the exact input
sequence. A predicted structure remains a research prediction rather than an experimental
measurement.

### `2_3_api.json`
```json
{
  "case_id": "patient-001",
  "targets": [
    {
      "target_id": "target_f82a91",
      "gene_symbol": "TP53",
      "transcript_id": "ENST00000269305.9",
      "haplotype": "unknown",
      "source_variants": ["17:7675088:C:T"],
      "reference_sequence_hash": "sha256:...",
      "variant_sequence": "MEEPQ...",
      "protein_change": "p.R175H",
      "expression_tpm": 42.7,
      "translation_status": "complete",
      "structure": {
        "status": "complete",
        "file": "targets/target_f82a91.cif",
        "source": "ColabFold",
        "mean_plddt": 86.2,
        "pae_file": "targets/target_f82a91.pae.json"
      }
    }
  ],
  "failures": []
}
```

## Step 2 -> Step 3

Step 2 returns personalized sequences, folded target structures, and provenance:

```text
proteins.fasta
reference-proteins.fasta
targets.jsonl
targets/<target-id>.cif
targets/<target-id>.pae.json
rejected-records.jsonl
provenance.json
```

Example protein record:

```json
{
  "target_id": "target_f82a91",
  "case_id": "sample-001",
  "gene_id": "ENSG00000141510.19",
  "gene_symbol": "TP53",
  "transcript_id": "ENST00000269305.9",
  "haplotype": "unknown",
  "source_variants": ["17:7675088:C:T"],
  "expression_tpm": 42.7,
  "reference_amino_acid_sequence": "...",
  "amino_acid_sequence": "...",
  "changed_residues": [175],
  "sequence_hash": "sha256:...",
  "orf_source": "GENCODE",
  "warnings": ["unphased_variants"],
  "structure": {
    "status": "complete",
    "file": "targets/target_f82a91.cif",
    "source": "ColabFold",
    "mean_plddt": 86.2,
    "pae_file": "targets/target_f82a91.pae.json"
  }
}
```

The Bespoke Bio Step 3 module ingests `targets.jsonl` and the target mmCIF files for lock
selection and key generation. Key sequences, target-key complexes, interface scores, and
key rankings are Step 3 outputs and are intentionally absent from `2_3_api.json`.

## Current implementation status

The current Step 2 implementation contains the inner translation kernel:

```text
resolved CDS + CDS-relative variants -> personalized amino-acid sequence
```

The next implementation slices are:

1. Accept the manifest, VCF, and expression table.
2. Run or consume VEP transcript mappings.
3. Resolve reference CDS sequences from a pinned GENCODE dataset.
4. Convert genomic variants into the existing translation-kernel input.
5. Add the exact-sequence structure resolver and ColabFold worker.
6. Export the Step 3 FASTA, JSONL, and mmCIF bundle.

The sequence adapter does not require a GPU; the ColabFold worker does. Both must use pinned
reference/model versions and emit provenance.

## Later pathogen mode

Pathogen support changes the reference adapter rather than the translation contract. Add
a taxon ID, reference accession, organism-specific genetic-code table, and pathogen
transcript/CDS annotations from NCBI RefSeq/GenBank or UniProt Proteomes.

The Step 2 output remains amino-acid FASTA, folded target structures, and provenance for
Step 3.

## Primary resources

- GENCODE: https://www.gencodegenes.org/
- MANE: https://www.ncbi.nlm.nih.gov/refseq/MANE/
- gnomAD: https://gnomad.broadinstitute.org/
- Ensembl VEP: https://www.ensembl.org/info/docs/tools/vep/
- TransDecoder: https://github.com/TransDecoder/TransDecoder
- NCBI RefSeq: https://www.ncbi.nlm.nih.gov/refseq/
- UniProt: https://www.uniprot.org/
- RCSB PDB: https://www.rcsb.org/
- AlphaFold DB: https://alphafold.ebi.ac.uk/
- ColabFold: https://github.com/sokrypton/ColabFold
