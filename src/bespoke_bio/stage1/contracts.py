"""Portable records. Runtime paths belong in manifests, never in public source locks."""

from typing import Annotated, Literal
from urllib.parse import urlsplit

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

Digest = Annotated[str, Field(pattern=r"^[0-9a-f]{64}$")]
Identifier = Annotated[str, Field(min_length=1, max_length=200)]


class Record(BaseModel):
    """Reject unknown fields and non-finite numbers at every external boundary."""

    model_config = ConfigDict(extra="forbid", allow_inf_nan=False)


class Artifact(Record):
    """One immutable file; bytes and SHA-256 describe its exact stored representation."""

    role: Identifier
    location: str
    bytes: int = Field(ge=0)
    sha256: Digest
    format: str
    index_role: str | None = None
    tool: str | None = None

    @field_validator("location")
    @classmethod
    def unsigned_location(cls, value: str) -> str:
        parsed = urlsplit(value)
        if parsed.query or parsed.fragment or parsed.username or parsed.password:
            raise ValueError("artifact locations must not contain credentials or query strings")
        return value


class ReadInput(Record):
    """A delivered FASTQ mate, with specimen identity independent of donor identity."""

    donor: Identifier
    specimen: Identifier
    library: Identifier
    assay: Literal["germline_wgs", "bulk_rna"]
    lane: Identifier
    mate: Literal[1, 2]
    protocol: str
    passage: str | None = None
    aliquot: str | None = None
    artifact: Artifact


class Intake(Record):
    """Complete pilot delivery. Each library/lane has exactly one ordered mate pair."""

    schema_version: Literal["1.0.0"] = "1.0.0"
    case_id: Identifier
    sex: Literal["XX", "XY"]
    reads: list[ReadInput] = Field(min_length=4)

    @model_validator(mode="after")
    def pairs(self):
        keys = [(r.library, r.lane, r.mate) for r in self.reads]
        if len(keys) != len(set(keys)):
            raise ValueError("duplicate library/lane/mate")
        if len({r.donor for r in self.reads}) != 1:
            raise ValueError("pilot requires one donor")
        for assay in ("germline_wgs", "bulk_rna"):
            group = [r for r in self.reads if r.assay == assay]
            if len({(r.library, r.specimen, r.protocol) for r in group}) != 1:
                raise ValueError(f"pilot requires one {assay} library and protocol")
            for lane in {r.lane for r in group}:
                pair = [r for r in group if r.lane == lane]
                if sorted(r.mate for r in pair) != [1, 2]:
                    raise ValueError(f"incomplete pair: {assay}/{lane}")
        if len({r.library for r in self.reads}) != 2:
            raise ValueError("DNA and RNA require distinct library identities")
        return self


class Policy(Record):
    """Versioned pilot gates; changes alter the run identity and require qualification."""

    schema_version: Literal["1.0.0"] = "1.0.0"
    name: str = "human-germline-bulk-rna-v1"
    min_depth: int = Field(default=10, ge=1)
    min_genotype_quality: int = Field(default=20, ge=0)
    min_base_quality: int = Field(default=20, ge=0)
    min_mapping_quality: int = Field(default=20, ge=0)
    identity_lod: float = 5


class ExpressionRow(Record):
    """Salmon transcript-total abundance; NumReads is an estimated fragment count."""

    transcript_id: Identifier
    length: int = Field(gt=0)
    effective_length: float = Field(ge=0)
    tpm: float = Field(ge=0)
    num_reads: float = Field(ge=0)


class Disposition(Record):
    """One result for each transcript in the primary GTF, including unselected records."""

    transcript_id: Identifier
    status: Literal["out_of_scope", "deferred", "partially_resolved", "selected"]
    reasons: list[str]
    primary_reason: str | None = None
    not_evaluated: list[str] = Field(default_factory=list)


class Check(Record):
    """A named acceptance result linked to the command and evidence that established it."""

    name: Identifier
    status: Literal["passed", "failed", "not_run"]
    command: list[str]
    evidence_roles: list[str]
    detail: str | None = None


RELEASE_CHECKS = frozenset(
    {
        "paired_reads",
        "reference_compatibility",
        "tool_locks",
        "dna_integrity",
        "rna_integrity",
        "sample_identity",
        "strandedness",
        "callability",
        "variant_reconciliation",
        "reference_translation",
        "handoff_translation",
        "disposition_accounting",
        "artifact_integrity",
    }
)


class Manifest(Record):
    """Publication record. Only ready manifests may be submitted to Step 2."""

    schema_version: Literal["1.0.0"] = "1.0.0"
    run_id: Digest
    attempt_id: Identifier
    case_id: Identifier
    status: Literal["ready", "no_targets", "failed"]
    case_sha256: Digest | None
    adapter_sha256: Digest
    contract_sha256: Digest
    reference_sha256: Digest
    config_sha256: Digest
    intake_sha256: Digest
    artifacts: list[Artifact]
    checks: list[Check]
    transcript_counts: dict[str, Annotated[int, Field(ge=0)]]
    target_count: int = Field(ge=0)
    unique_sequence_count: int = Field(ge=0)

    @model_validator(mode="after")
    def ready_requires_evidence(self):
        roles = [a.role for a in self.artifacts]
        if len(roles) != len(set(roles)):
            raise ValueError("artifact roles must be unique")
        names = [c.name for c in self.checks]
        if len(names) != len(set(names)):
            raise ValueError("acceptance check names must be unique")
        if self.unique_sequence_count > self.target_count:
            raise ValueError("unique sequences cannot exceed biological targets")
        for check in self.checks:
            if not set(check.evidence_roles).issubset(roles):
                raise ValueError(f"check references missing evidence: {check.name}")
            if check.status == "passed" and (not check.command or not check.evidence_roles):
                raise ValueError("passed checks require their command and evidence")
        if self.status == "ready":
            if not self.case_sha256 or self.target_count == 0 or "case" not in roles:
                raise ValueError("ready requires a nonempty hashed case")
            if not RELEASE_CHECKS.issubset(names) or any(c.status != "passed" for c in self.checks):
                raise ValueError("ready requires all required and declared release checks to pass")
        if self.status == "no_targets" and self.target_count:
            raise ValueError("no_targets cannot contain targets")
        return self
