from hashlib import sha256
from typing import Annotated, Literal

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, ConfigDict, Field, field_validator

# API schema

Identifier = Annotated[str, Field(min_length=1, max_length=200)]


def normalize_dna(value: str) -> str:
    value = value.upper().replace("U", "T")
    invalid = sorted(set(value) - set("ACGTN"))
    if invalid:
        raise ValueError(f"unsupported nucleotide symbols: {', '.join(invalid)}")
    return value


class Model(BaseModel):
    model_config = ConfigDict(extra="forbid")


class Reference(Model):
    genome_build: Identifier
    annotation: Identifier


class Expression(Model):
    tpm: float = Field(ge=0)
    raw_count: float | None = Field(default=None, ge=0)


class Variant(Model):
    variant_id: Identifier
    cds_position: int = Field(ge=1, description="One-based position in the reference CDS")
    ref: str = Field(min_length=1)
    alt: str = ""
    genotype: str | None = None
    phase_set: str | None = None
    quality: str | None = None

    @field_validator("ref", "alt")
    @classmethod
    def validate_allele(cls, value: str) -> str:
        return normalize_dna(value)


class Transcript(Model):
    gene_id: Identifier
    gene_symbol: Identifier
    transcript_id: Identifier
    haplotype: Identifier = "unknown"
    coding_sequence: str = Field(min_length=3)
    expression: Expression
    variants: list[Variant] = Field(default_factory=list)

    @field_validator("coding_sequence")
    @classmethod
    def validate_sequence(cls, value: str) -> str:
        return normalize_dna(value)


class TranslateRequest(Model):
    case_id: Identifier
    reference: Reference
    transcripts: list[Transcript] = Field(min_length=1)


class Protein(Model):
    target_id: Identifier
    gene_id: Identifier
    gene_symbol: Identifier
    transcript_id: Identifier
    haplotype: Identifier
    source_variant_ids: list[str]
    expression: Expression
    reference_amino_acid_sequence: str
    amino_acid_sequence: str
    amino_acid_length: int
    sequence_hash: str
    first_changed_residue: int | None
    translation_status: Literal["complete", "complete_with_warnings"]
    warnings: list[str]


class TranslateResponse(Model):
    case_id: Identifier
    reference: Reference
    targets: list[Protein]
    target_count: int
    unique_sequence_count: int


# Translation

# Standard genetic code in T, C, A, G base order.
_CODONS = (first + second + third for first in "TCAG" for second in "TCAG" for third in "TCAG")
_AMINO_ACIDS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = dict(zip(_CODONS, _AMINO_ACIDS, strict=True))


class VariantError(ValueError):
    pass


def apply_variants(reference: str, variants: list[Variant]) -> str:
    """Apply non-overlapping variants using positions from the reference CDS."""
    pieces: list[str] = []
    cursor = 0

    for variant in sorted(variants, key=lambda item: item.cds_position):
        start = variant.cds_position - 1
        end = start + len(variant.ref)
        if start < cursor:
            raise VariantError(f"{variant.variant_id} overlaps an earlier variant")
        if end > len(reference):
            raise VariantError(f"{variant.variant_id} extends beyond the reference CDS")
        observed = reference[start:end]
        if observed != variant.ref:
            raise VariantError(
                f"{variant.variant_id} expected {variant.ref} at CDS position "
                f"{variant.cds_position}, found {observed}"
            )
        pieces.extend((reference[cursor:start], variant.alt))
        cursor = end

    pieces.append(reference[cursor:])
    return "".join(pieces)


def translate_cds(cds: str) -> tuple[str, list[str]]:
    warnings: list[str] = []
    if not cds.startswith("ATG"):
        warnings.append("missing_start")
    if len(cds) % 3:
        warnings.append("partial_codon")

    protein: list[str] = []
    codon_count = len(cds) // 3
    for index in range(codon_count):
        amino_acid = CODON_TABLE.get(cds[index * 3 : index * 3 + 3], "X")
        if amino_acid == "X":
            warnings.append("ambiguous_codon")
        if amino_acid == "*":
            if index < codon_count - 1:
                warnings.append("internal_stop")
            break
        protein.append(amino_acid)

    return "".join(protein), warnings


def sequence_hash(sequence: str) -> str:
    return f"sha256:{sha256(sequence.encode()).hexdigest()}"


def first_change(reference: str, alternate: str) -> int | None:
    for index, (left, right) in enumerate(zip(reference, alternate), start=1):
        if left != right:
            return index
    return min(len(reference), len(alternate)) + 1 if len(reference) != len(alternate) else None


def translate_case(request: TranslateRequest) -> TranslateResponse:
    targets: list[Protein] = []

    for transcript in request.transcripts:
        reference_protein, _ = translate_cds(transcript.coding_sequence)
        personalized_cds = apply_variants(transcript.coding_sequence, transcript.variants)
        protein, warnings = translate_cds(personalized_cds)
        if transcript.haplotype == "unknown" and transcript.variants:
            warnings.append("unphased_variants")

        protein_hash = sequence_hash(protein)
        identity = (
            f"{request.case_id}|{transcript.transcript_id}|{transcript.haplotype}|{protein_hash}"
        )
        targets.append(
            Protein(
                target_id=f"target_{sha256(identity.encode()).hexdigest()[:12]}",
                gene_id=transcript.gene_id,
                gene_symbol=transcript.gene_symbol,
                transcript_id=transcript.transcript_id,
                haplotype=transcript.haplotype,
                source_variant_ids=[variant.variant_id for variant in transcript.variants],
                expression=transcript.expression,
                reference_amino_acid_sequence=reference_protein,
                amino_acid_sequence=protein,
                amino_acid_length=len(protein),
                sequence_hash=protein_hash,
                first_changed_residue=first_change(reference_protein, protein),
                translation_status="complete_with_warnings" if warnings else "complete",
                warnings=warnings,
            )
        )

    return TranslateResponse(
        case_id=request.case_id,
        reference=request.reference,
        targets=targets,
        target_count=len(targets),
        unique_sequence_count=len({target.sequence_hash for target in targets}),
    )


# HTTP API

app = FastAPI(
    title="Bespoke Bio",
    summary="Translate expressed rare variants into target proteins",
)


@app.get("/healthz")
def healthcheck() -> dict[str, str]:
    return {"status": "ok"}


@app.post("/v1/translate", response_model=TranslateResponse)
def translate(request: TranslateRequest) -> TranslateResponse:
    try:
        return translate_case(request)
    except VariantError as error:
        raise HTTPException(status_code=422, detail=str(error)) from error
