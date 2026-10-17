from __future__ import annotations

from typing import Annotated, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator

Identifier = Annotated[str, Field(min_length=1, max_length=200)]


def normalize_nucleotides(value: str) -> str:
    normalized = value.upper().replace("U", "T")
    invalid = sorted(set(normalized) - set("ACGTN"))
    if invalid:
        raise ValueError(f"unsupported nucleotide symbols: {', '.join(invalid)}")
    return normalized


class StrictModel(BaseModel):
    model_config = ConfigDict(extra="forbid")


class ReferenceBuild(StrictModel):
    genome_build: Identifier
    annotation: Identifier


class Expression(StrictModel):
    tpm: float = Field(ge=0)
    raw_count: float | None = Field(default=None, ge=0)


class CdsVariant(StrictModel):
    variant_id: Identifier
    cds_position: int = Field(ge=1, description="One-based position in the reference CDS")
    ref: str = Field(min_length=1)
    alt: str = Field(min_length=0)
    genotype: str | None = None
    phase_set: str | None = None
    quality: str | None = None

    @field_validator("ref", "alt")
    @classmethod
    def validate_allele(cls, value: str) -> str:
        return normalize_nucleotides(value)


class TranscriptInput(StrictModel):
    gene_id: Identifier
    gene_symbol: Identifier
    transcript_id: Identifier
    haplotype: Identifier = "unknown"
    coding_sequence: str = Field(min_length=3)
    expression: Expression
    variants: list[CdsVariant] = Field(default_factory=list)

    @field_validator("coding_sequence")
    @classmethod
    def validate_coding_sequence(cls, value: str) -> str:
        return normalize_nucleotides(value)


class TranslateRequest(StrictModel):
    case_id: Identifier
    reference: ReferenceBuild
    transcripts: list[TranscriptInput] = Field(min_length=1)


class TranslationWarning(StrictModel):
    code: Literal[
        "ambiguous_codon",
        "internal_stop",
        "missing_start",
        "partial_codon",
        "unphased_variants",
    ]
    message: str


class TargetProtein(StrictModel):
    target_id: Identifier
    gene_id: Identifier
    gene_symbol: Identifier
    transcript_id: Identifier
    haplotype: Identifier
    source_variant_ids: list[str]
    expression: Expression
    reference_amino_acid_sequence: str
    amino_acid_sequence: str
    amino_acid_length: int = Field(ge=0)
    sequence_hash: str
    first_changed_residue: int | None = Field(default=None, ge=1)
    translation_status: Literal["complete", "complete_with_warnings"]
    warnings: list[TranslationWarning]


class TranslateResponse(StrictModel):
    case_id: Identifier
    reference: ReferenceBuild
    targets: list[TargetProtein]
    target_count: int = Field(ge=0)
    unique_sequence_count: int = Field(ge=0)
