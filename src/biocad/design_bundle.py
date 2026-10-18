"""Export actual BoltzGen outputs as reusable, target-linked molecular components."""

import math
import shutil
from pathlib import Path

import gemmi
import numpy as np
from pydantic import BaseModel, ConfigDict

from .boltzgen_adapter import COMMIT, summarize
from .contracts import Component, Evaluation, TargetState
from .pipeline import save
from .structures import atoms, ports, protein_residues, sequence, sha, write_chains


class DesignBundle(BaseModel):
    model_config = ConfigDict(extra="forbid")
    schema_version: str = "0.1.0"
    target_states: list[TargetState]
    components: list[Component]
    evaluations: list[Evaluation]
    provenance: dict


def _number(value):
    try:
        number = float(value)
        return number if math.isfinite(number) else None
    except (TypeError, ValueError):
        return None


def export_bundle(results, target_state_path, output, *, mode="actual_cached_model_outputs"):
    """Export native candidates; refuse mismatched target sequences or ambiguous chains."""
    if mode not in {"actual_cached_model_outputs", "live_model_outputs"}:
        raise ValueError("Unsupported model output mode")
    results, target_state_path, output = map(Path, (results, target_state_path, output))
    output.mkdir(parents=True, exist_ok=True)
    state = TargetState.model_validate_json(target_state_path.read_text())
    if len(state.chains) != 1:
        raise ValueError("A single target chain is required")
    target_seq = next(iter(state.chains.values()))
    metrics = results / "all_designs_metrics.csv"
    if metrics.exists():
        import tempfile

        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp) / "final_ranked_designs"
            folder.mkdir()
            shutil.copy2(metrics, folder / metrics.name)
            report = summarize(temp, output / "native_summary.json")
    else:
        report = summarize(results, output / "native_summary.json")
    components, evaluations, fasta = [], [], []
    target_out = output / "target"
    target_out.mkdir(exist_ok=True)
    for field in ("source", "prepared"):
        path = getattr(state, field)
        if not path:
            continue
        source = Path(path) if Path(path).is_absolute() else target_state_path.parent / path
        if not source.is_file():
            raise FileNotFoundError(source)
        dest = target_out / (field + "".join(source.suffixes))
        shutil.copy2(source, dest)
        setattr(state, field, str(dest.relative_to(output)))
    for candidate in report["candidates"]:
        filename = Path(candidate["file_name"])
        if filename.name != str(filename):
            raise ValueError("Candidate file_name must be a basename")
        # Native filtering copies ranked names into nested folders; the unranked
        # refolded complex is retained here for every candidate (including failures).
        # Never fall back to a before_refolding backbone for these sequence metrics.
        locations = [
            results / "intermediate_designs_inverse_folded" / "refold_cif" / filename,
            results / "candidates" / filename,
            results / "final_ranked_designs" / filename,
        ]
        source = next((p for p in locations if p.is_file()), None)
        if source is None:
            raise FileNotFoundError(filename)
        structure = gemmi.read_structure(str(source))
        if len(structure) != 1:
            raise ValueError("Candidate must have exactly one model")
        binder_seq = candidate["sequence"]
        binders = [c for c in structure[0] if sequence(c) == binder_seq]
        targets = [c for c in structure[0] if sequence(c) == target_seq]
        if len(binders) != 1 or len(targets) != 1 or binders[0].name == targets[0].name:
            raise ValueError("Candidate target or binder sequence mismatches supplied target state")
        binder, target = binders[0], targets[0]
        component_id = "binder_" + sha((sha(source.read_bytes()) + binder_seq).encode())[:16]
        folder = output / component_id
        folder.mkdir(exist_ok=True)
        shutil.copy2(source, folder / "complex.cif")
        write_chains([binder], folder / "binder.cif")
        component = Component(
            id=component_id,
            structure=str((folder / "binder.cif").relative_to(output)),
            sequence=binder_seq,
            sequence_sha256=sha(binder_seq.encode()),
            role="computational_binder_candidate",
            ports=ports(binder),
            fixed_residues=[],
            evidence={
                "target_state_id": state.id,
                "native_candidate_id": candidate["id"],
                "complex_structure": str((folder / "complex.cif").relative_to(output)),
                "source_sha256": sha(source.read_bytes()),
                "engine_commit": COMMIT,
                "measured_binding": False,
            },
        )
        components.append(component)
        fasta.append(f">{component_id} native={candidate['id']} target={state.id}\n{binder_seq}\n")
        evaluations.append(
            Evaluation(
                id=component_id + "_filters",
                subject_id=component_id,
                target_state_id=state.id,
                method="BoltzGen native default filters",
                evidence="computational",
                status="pass" if candidate["passed"] else "fail",
                metrics={k: _number(v) for k, v in candidate["metrics"].items()},
                units={},
                details={
                    "failed_filters": candidate["failed_filters"],
                    "flags": candidate["upstream_filter_flags"],
                    "acceptance_scope": "default filters; does not establish binding",
                },
            )
        )
        evaluations.append(
            Evaluation(
                id=component_id + "_binding",
                subject_id=component_id,
                target_state_id=state.id,
                method="predicted interface confidence",
                evidence="computational prediction; no experimental binding",
                status="assessed",
                metrics={
                    k: _number(candidate["metrics"].get(k))
                    for k in ("design_to_target_iptm", "esmfold2_ipsae_min")
                },
                units={},
                details={
                    "selection_gate": "not calibrated; report values independently of filter pass"
                },
            )
        )
        target_xyz, _ = atoms(target)
        clearance = {}
        for port in component.ports:
            residue = protein_residues(binder)[0 if port.terminus == "N" else -1]
            anchor = residue.find_atom("N" if port.terminus == "N" else "C", "*")
            if anchor is not None:
                point = np.array([anchor.pos.x, anchor.pos.y, anchor.pos.z])
                clearance[port.terminus + "_target_clearance_angstrom"] = float(
                    np.linalg.norm(target_xyz - point, axis=1).min()
                )
        evaluations.append(
            Evaluation(
                id=component_id + "_ports",
                subject_id=component_id,
                target_state_id=state.id,
                method="terminal atom distance to bound target",
                evidence="coordinate geometry",
                status="assessed",
                metrics=clearance,
                units={k: "angstrom" for k in clearance},
                details={
                    "limitation": "point clearance only; linker path, solvent exposure and fusion feasibility untested"
                },
            )
        )
    bundle = DesignBundle(
        target_states=[state],
        components=components,
        evaluations=evaluations,
        provenance={
            "engine": "BoltzGen",
            "commit": COMMIT,
            "mode": mode,
            "source_metrics_sha256": sha(
                (
                    metrics
                    if metrics.exists()
                    else results / "final_ranked_designs/all_designs_metrics.csv"
                ).read_bytes()
            ),
            "binding_evidence": "no measured binding or affinity",
        },
    )
    (output / "binders.fasta").write_text("".join(fasta))
    save(output / "bundle.json", bundle.model_dump())
    save(
        output / "manifest.json",
        {
            str(p.relative_to(output)): sha(p.read_bytes())
            for p in sorted(output.rglob("*"))
            if p.is_file() and p.name != "manifest.json"
        },
    )
    return bundle.model_dump()
