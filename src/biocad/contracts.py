from typing import Literal

import numpy as np
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator


class Record(BaseModel):
    model_config = ConfigDict(extra="forbid")


class ResidueRef(Record):
    chain: str
    number: int
    insertion: str = ""


class Transform(Record):
    """Column-vector convention: destination = rotation @ source + translation, Å."""

    rotation: list[list[float]] = Field(default_factory=lambda: np.eye(3).tolist())
    translation: list[float] = Field(default_factory=lambda: [0.0, 0.0, 0.0])

    @model_validator(mode="after")
    def proper(self):
        r, t = np.array(self.rotation), np.array(self.translation)
        if (
            r.shape != (3, 3)
            or t.shape != (3,)
            or not np.isfinite(r).all()
            or not np.isfinite(t).all()
        ):
            raise ValueError("Transform must contain finite 3x3 rotation and 3-vector translation")
        if not np.allclose(r.T @ r, np.eye(3), atol=1e-6) or not np.isclose(
            np.linalg.det(r), 1.0, atol=1e-6
        ):
            raise ValueError("Rotation must be proper and orthonormal; reflections rejected")
        return self


class TargetState(Record):
    id: str
    source: str
    prepared: str
    source_sha256: str
    chains: dict[str, str]
    residue_mapping: list[dict]
    unresolved: list[dict]
    context: dict
    evidence: Literal["experimental", "predicted", "synthetic"]
    variants: list[str] = []
    measurements: list[dict] = []


class BindingSite(Record):
    id: str
    target_state_id: str
    residues: list[ResidueRef]
    origin: Literal["observed_interface", "mapped_analogy", "operator_selection", "proposed_patch"]
    frame: Transform | None = None
    evidence: dict


class Port(Record):
    terminus: Literal["N", "C"]
    residue: ResidueRef
    anchors: list[str]
    frame: Transform | None
    connection_type: Literal["peptide_fusion"] = "peptide_fusion"


class Component(Record):
    id: str
    structure: str
    sequence: str
    sequence_sha256: str
    role: str
    ports: list[Port]
    fixed_residues: list[ResidueRef]
    evidence: dict


class Instance(Record):
    id: str
    component_id: str
    output_chain: str
    transform: Transform


class Assembly(Record):
    id: str
    instances: list[Instance]
    connections: list[dict]
    required_interfaces: list[tuple[str, str]]
    state: str
    structure: str
    hypothesis: str
    target_placements: list[dict] = []


class Evaluation(Record):
    id: str
    subject_id: str
    target_state_id: str | None = None
    method: str
    evidence: str
    status: Literal["pass", "fail", "missing", "assessed", "unsupported"]
    metrics: dict[str, float | int | None]
    units: dict[str, str]
    details: dict = {}


class DesignRequest(Record):
    target_state_ids: list[str] = Field(min_length=1)
    recognition_mode: Literal["cover", "selective", "evaluate_existing"] = "evaluate_existing"
    exclusion_state_ids: list[str] = []
    binder_length: int = Field(default=75, ge=25, le=200)
    site: list[ResidueRef] = Field(min_length=1)
    attachment: Literal["N", "C"] = "C"
    terminal_clearance_angstrom: float = Field(default=8.0, gt=0)
    linker_sequence: str = "GGGGSGGGGSGGGGS"
    second_module_offset_angstrom: float = Field(default=45.0, gt=0)
    second_module_axis: Literal["x", "y", "z"] = "y"

    @field_validator("linker_sequence")
    @classmethod
    def valid_linker(cls, value):
        if not value or len(value) > 100 or set(value) - set("ARNDCQEGHILKMFPSTWYV"):
            raise ValueError("Linker must contain 1..100 canonical amino acids")
        return value

    seed: int = Field(default=7, ge=0)
    steps: int = Field(default=50, ge=1, le=500)


class DesignResponse(Record):
    schema_version: str = "0.1.0"
    request: DesignRequest
    target_states: list[TargetState]
    sites: list[BindingSite]
    components: list[Component]
    assemblies: list[Assembly]
    evaluations: list[Evaluation]
    stages: dict[str, dict]
    artifacts: dict[str, str]


class PlacementSearchRequest(Record):
    """Bounded rigid packing of two complete bound modules; no backbone editing."""

    offsets_angstrom: list[float] = Field(
        default_factory=lambda: [30.0, 40.0, 50.0, 65.0], min_length=1, max_length=8
    )
    z_rotations_degrees: list[float] = Field(
        default_factory=lambda: [0.0, 90.0], min_length=1, max_length=8
    )
    linker_sequence: str = "GGGGSGGGGSGGGGS"
    attachment: Literal["N", "C"] = "C"
    vdw_overlap_tolerance_angstrom: float = Field(default=0.4, ge=0.0, le=1.0)
    top_k: int = Field(default=3, ge=1, le=10)

    @field_validator("offsets_angstrom")
    @classmethod
    def valid_offsets(cls, values):
        if any(not np.isfinite(v) or v <= 0 or v > 250 for v in values):
            raise ValueError("Offsets must be finite and in (0, 250] Å")
        return values

    @field_validator("z_rotations_degrees")
    @classmethod
    def valid_angles(cls, values):
        if not np.isfinite(values).all():
            raise ValueError("Angles must be finite")
        return values

    @field_validator("linker_sequence")
    @classmethod
    def valid_search_linker(cls, value):
        return DesignRequest.valid_linker(value)

    @model_validator(mode="after")
    def bounded(self):
        if len(self.offsets_angstrom) * len(self.z_rotations_degrees) * 6 > 192:
            raise ValueError("Search capped at 192 sampled layouts")
        return self


class MotifPlanRequest(Record):
    target_scaffold_length: int | None = Field(default=None, ge=30, le=300)
    flank_residues: int = Field(default=2, ge=0, le=5)
    preserve_native_disulfides: bool = True
    bridge_min_residues: int = Field(default=4, ge=1, le=40)
    bridge_max_residues: int = Field(default=24, ge=1, le=60)
    bridge_slack_residues: int = Field(default=3, ge=0, le=10)
    terminal_residues: int = Field(default=8, ge=1, le=30)
    maximum_scaffold_length: int = Field(default=180, ge=30, le=300)

    @model_validator(mode="after")
    def ordered_budget(self):
        if self.bridge_min_residues > self.bridge_max_residues:
            raise ValueError("Bridge minimum exceeds maximum")
        return self
