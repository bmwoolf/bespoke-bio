import importlib.metadata
import json
import platform
import shutil
import time
from pathlib import Path

import gemmi
import numpy as np
from pydantic import BaseModel

from .contracts import (
    Assembly,
    BindingSite,
    Component,
    DesignRequest,
    DesignResponse,
    Evaluation,
    Instance,
    Transform,
)
from .geometry import distances, frame, pair_metrics
from .retrieval import move, search, transfer
from .structures import (
    atoms,
    ports,
    prepare,
    protein_residues,
    ref,
    select_indices,
    sequence,
    sha,
    write_chains,
)


def save(path, value):
    Path(path).write_text(
        json.dumps(
            value.model_dump() if isinstance(value, BaseModel) else value, indent=2, allow_nan=False
        )
        + "\n"
    )


def interface(subject, state_id, a, b, evidence="geometry_only"):
    metrics = pair_metrics(atoms(a)[0], atoms(b)[0])
    return Evaluation(
        id="eval_" + sha((subject + state_id + json.dumps(metrics)).encode())[:16],
        subject_id=subject,
        target_state_id=state_id,
        method="heavy_atom_geometry_v1",
        evidence=evidence,
        status="fail" if metrics["severe_clash_pairs"] else "assessed",
        metrics=metrics,
        units={
            "min_heavy_atom_distance": "angstrom",
            "severe_clash_pairs": "count (<1.5 Å)",
            "contact_atom_pairs": "count (<5 Å)",
        },
        details={
            "confidence": None,
            "affinity": None,
            "function": None,
            "interpretation": "Contact is geometry, not binding validation",
        },
    )


def run_example(fixture_dir, output, request=None, foldseek="foldseek"):
    start = time.monotonic()
    out = Path(output)
    if out.exists() and any(out.iterdir()):
        raise ValueError("Use an empty output directory to preserve prior runs")
    out.mkdir(parents=True, exist_ok=True)
    fixtures = Path(fixture_dir)
    for name in ["5JDS.cif", "5J89.cif"]:
        shutil.copyfile(fixtures / name, out / name)
    source = gemmi.read_structure(str(out / "5JDS.cif"))
    target = source[0]["A"].clone()
    binder = source[0]["B"].clone()
    alternate = gemmi.read_structure(str(out / "5J89.cif"))[0]["A"].clone()
    s1 = prepare(
        out / "5JDS.cif",
        out / "target.cif",
        ["A"],
        context={
            "pdb": "5JDS",
            "assembly": "deposited asymmetric unit chains A/B; inspect original assembly metadata",
            "condition": "nanobody-bound PD-L1 IgV domain",
        },
    )
    s2 = prepare(
        out / "5J89.cif",
        out / "alternate.cif",
        ["A"],
        context={
            "pdb": "5J89",
            "condition": "small-molecule-bound PD-L1 crystal conformation; not an allelic variant",
            "omitted_context": "other PD-L1 chains, ligand and waters remain in original 5J89.cif",
        },
    )
    target = gemmi.read_structure(s1.prepared)[0]["A"].clone()
    alternate = gemmi.read_structure(s2.prepared)[0]["A"].clone()
    # Reparse binder via the same explicit conformer policy.
    binder_st = source.clone()
    binder_st.remove_alternative_conformations()
    binder = binder_st[0]["B"].clone()
    bind_atoms = atoms(binder)[0]
    target_res = protein_residues(target)
    selected = [
        ref("A", r)
        for r in target_res
        if distances(atoms(gemmi_chain(r))[0], bind_atoms).min() < 5.0
    ]
    req = request or DesignRequest(target_state_ids=[s1.id, s2.id], site=selected)
    if set(req.target_state_ids) != {s1.id, s2.id}:
        raise ValueError("Demo request must name its two exact fixture states")
    idx = select_indices(s1, req.site, "A")
    anchors = [target_res[i].find_atom("CA", "*") for i in idx[:3]]
    site_frame = (
        frame(*[[a.pos.x, a.pos.y, a.pos.z] for a in anchors])
        if len(anchors) == 3 and all(a is not None for a in anchors)
        else None
    )
    site = BindingSite(
        id="site_" + sha(json.dumps([r.model_dump() for r in req.site]).encode())[:16],
        target_state_id=s1.id,
        residues=req.site,
        origin="operator_selection" if request else "observed_interface",
        frame=site_frame,
        evidence={
            "source": "5JDS A/B",
            "definition": "selected by request"
            if request
            else "target residues with any heavy atom within 5 Å of nanobody",
            "frame_anchors": [r.model_dump() for r in req.site[:3]],
            "surface_accessibility": "not computed",
        },
    )
    write_chains([binder], out / "component.cif")
    seq = sequence(binder)
    comp = Component(
        id="component_" + sha((seq + sha((out / "component.cif").read_bytes())).encode())[:16],
        structure="component.cif",
        sequence=seq,
        sequence_sha256=sha(seq.encode()),
        role="recognition",
        ports=ports(binder),
        fixed_residues=[ref("B", r) for r in protein_residues(binder)],
        evidence={
            "type": "experimental_structure",
            "pdb": "5JDS",
            "chain": "B",
            "construct_scope": "observed-only computational construct; unresolved native termini omitted, not experimentally validated as this truncation",
            "de_novo_cached_output": False,
        },
    )
    library = out / "library"
    library.mkdir()
    write_chains([target], library / "5JDS_target.cif")
    retrieval = search(out / "alternate.cif", library, out / "retrieval", executable=foldseek)
    if not retrieval["hits"]:
        raise RuntimeError("No Foldseek hit; see retrieval log")
    partner, pose = transfer(retrieval["hits"][0], alternate, target, binder)
    save(out / "retrieval.json", retrieval)
    save(out / "transferred_pose.json", pose)
    write_chains([alternate, renamed(partner, "K")], out / "transferred_complex.cif")
    query_res = protein_residues(alternate)
    mapped = [
        ref(alternate.name, query_res[p["query_index"]])
        for p in pose["residue_pairs"]
        if p["template_index"] in idx
    ]
    mapped_site = BindingSite(
        id=site.id + "_mapped",
        target_state_id=s2.id,
        residues=mapped,
        origin="mapped_analogy",
        evidence={
            "source_site_id": site.id,
            "ca_rmsd_angstrom": pose["ca_rmsd_angstrom"],
            "unmapped_source_residues": len(idx) - len(mapped),
        },
    )
    evals = [
        interface(comp.id, s1.id, binder, target, "experimental_pose_geometry"),
        interface(comp.id, s2.id, partner, alternate, "transferred_analogy"),
    ]
    # One candidate, two states. Absence of prediction is explicit in each matrix cell.
    for s in [s1, s2]:
        evals.append(
            Evaluation(
                id="missing_confidence_" + s.id,
                subject_id=comp.id,
                target_state_id=s.id,
                method="independent_full_atom_reprediction",
                evidence="missing",
                status="missing",
                metrics={"chain_pair_iptm": None},
                units={"chain_pair_iptm": "0..1"},
            )
        )
    terminal = next(p for p in comp.ports if p.terminus == req.attachment)
    clearance = (
        float(distances([terminal.frame.translation], atoms(target)[0]).min())
        if terminal.frame
        else None
    )
    evals.append(
        Evaluation(
            id="attachment_" + comp.id,
            subject_id=comp.id,
            target_state_id=s1.id,
            method="terminal_clearance_v1",
            evidence="geometry_only",
            status="missing"
            if clearance is None
            else ("pass" if clearance >= req.terminal_clearance_angstrom else "fail"),
            metrics={"clearance": clearance, "threshold": req.terminal_clearance_angstrom},
            units={"clearance": "angstrom", "threshold": "angstrom"},
            details={
                "terminus": req.attachment,
                "scope": "nearest target-context atom; does not establish linker folding",
            },
        )
    )
    assemblies = []
    for label, offset in [
        ("proposed", req.second_module_offset_angstrom),
        ("deliberate_overlap", 0.0),
    ]:
        vector = [0.0, 0.0, 0.0]
        vector["xyz".index(req.second_module_axis)] = offset
        tr = Transform(translation=vector)
        second = move(binder, tr)
        second_target = move(target, tr)
        aid = (
            "assembly_"
            + sha(
                (
                    comp.id
                    + label
                    + json.dumps(tr.model_dump(), sort_keys=True)
                    + req.linker_sequence
                    + req.attachment
                ).encode()
            )[:16]
        )
        path = label + ".cif"
        write_chains(
            [
                renamed(target, "T"),
                renamed(binder, "K"),
                renamed(second_target, "U"),
                renamed(second, "L"),
            ],
            out / path,
        )
        first_res = protein_residues(binder)
        second_res = protein_residues(second)
        # C->N or second C->first N; carbon/nitrogen anchors, not CA distance.
        ends = (
            (first_res[-1].find_atom("C", "*"), second_res[0].find_atom("N", "*"))
            if req.attachment == "C"
            else (second_res[-1].find_atom("C", "*"), first_res[0].find_atom("N", "*"))
        )
        gap = float(np.linalg.norm(np.array(ends[0].pos.tolist()) - ends[1].pos.tolist()))
        # Necessary upper bound only: linker residues plus joining bond lengths.
        max_reach = 3.8 * (len(req.linker_sequence) + 1)
        reach = "fail" if gap > max_reach else "assessed"
        assembly = Assembly(
            id=aid,
            instances=[
                Instance(id="key1", component_id=comp.id, output_chain="K", transform=Transform()),
                Instance(id="key2", component_id=comp.id, output_chain="L", transform=tr),
            ],
            connections=[
                {
                    "type": "peptide_fusion",
                    "from": "key1:C" if req.attachment == "C" else "key2:C",
                    "to": "key2:N" if req.attachment == "C" else "key1:N",
                    "linker_sequence": req.linker_sequence,
                    "status": "unmodeled",
                    "reach_status": reach,
                }
            ],
            required_interfaces=[("key1", "target1"), ("key2", "target2")],
            state="rigid_placement_hypothesis",
            structure=path,
            target_placements=[
                {
                    "id": "target1",
                    "state_id": s1.id,
                    "output_chain": "T",
                    "transform": Transform().model_dump(),
                },
                {
                    "id": "target2",
                    "state_id": s1.id,
                    "output_chain": "U",
                    "transform": tr.model_dump(),
                },
            ],
            hypothesis="Two PD-L1 encounters placed together artificially; relative target pose not observed. Tandem identical nanobodies, not established recruitment or switching.",
        )
        assemblies.append(assembly)
        for name, a, b in [
            ("key1:target1", binder, target),
            ("key2:target2", second, second_target),
            ("key1:key2", binder, second),
            ("target1:target2", target, second_target),
            ("key1:target2", binder, second_target),
            ("key2:target1", second, target),
        ]:
            evals.append(interface(aid + "/" + name, s1.id, a, b))
        evals.append(
            Evaluation(
                id=aid + "_linker",
                subject_id=aid,
                method="linker_upper_bound_v1",
                evidence="geometry_only",
                status=reach,
                metrics={"anchor_distance": gap, "maximum_contour_estimate": max_reach},
                units={"anchor_distance": "angstrom", "maximum_contour_estimate": "angstrom"},
                details={
                    "connection_modeled": False,
                    "sequence": seq + req.linker_sequence + seq,
                    "scope": "necessary reach bound; no positive feasibility claim",
                },
            )
        )
    for assembly in assemblies:
        failures = [
            e.id for e in evals if e.subject_id.startswith(assembly.id) and e.status == "fail"
        ]
        evals.append(
            Evaluation(
                id=assembly.id + "_summary",
                subject_id=assembly.id,
                target_state_id=s1.id,
                method="assembly_geometry_summary_v1",
                evidence="geometry_only",
                status="fail" if failures else "assessed",
                metrics={"failed_checks": len(failures)},
                units={"failed_checks": "count"},
                details={
                    "failed_evaluation_ids": failures,
                    "joint_reprediction": "missing",
                    "connected_construct": "unmodeled",
                },
            )
        )
        evals.append(
            Evaluation(
                id=assembly.id + "_alternate_missing",
                subject_id=assembly.id,
                target_state_id=s2.id,
                method="whole_assembly_state_evaluation",
                evidence="missing",
                status="missing",
                metrics={"confidence": None},
                units={"confidence": "not computed"},
            )
        )
    (out / "chains.fasta").write_text(
        "".join(
            f">{c}\n{s}\n"
            for c, s in [
                ("target_reference_observed", s1.chains["A"]),
                ("target_alternate_observed", s2.chains["A"]),
                ("key", seq),
                ("tandem_proposed_unmodeled", seq + req.linker_sequence + seq),
            ]
        )
    )
    for s in [s1, s2]:
        s.source = Path(s.source).name
        s.prepared = Path(s.prepared).name
    response = DesignResponse(
        request=req,
        target_states=[s1, s2],
        sites=[site, mapped_site],
        components=[comp],
        assemblies=assemblies,
        evaluations=evals,
        stages={
            "preparation": {"status": "complete"},
            "retrieval": {
                "status": "complete",
                "summary": "bounded Foldseek search and mapped partner pose",
            },
            "binder_design": {
                "status": "blocked",
                "summary": "No GPU/weights configured; retained experimental nanobody fixture, no generated scores",
            },
            "assembly": {
                "status": "partial",
                "summary": "Rigid geometry evaluated; fusion linker and joint reprediction missing",
            },
        },
        artifacts={
            n: n
            for n in [
                "chains.fasta",
                "component.cif",
                "target.cif",
                "alternate.cif",
                "retrieval.json",
                "transferred_pose.json",
                "transferred_complex.cif",
                "proposed.cif",
                "deliberate_overlap.cif",
            ]
        },
    )
    save(out / "request.json", req)
    save(out / "response.json", response)
    from .report import render

    render(out, response)
    save(
        out / "scene.json",
        {
            "structures": [
                {"path": n, "purpose": n.removesuffix(".cif")}
                for n in [
                    "component.cif",
                    "transferred_complex.cif",
                    "proposed.cif",
                    "deliberate_overlap.cif",
                ]
            ],
            "site": site.model_dump(),
            "units": "angstrom",
        },
    )
    manifest = {
        "schema_version": "0.1.0",
        "python": platform.python_version(),
        "packages": {p: importlib.metadata.version(p) for p in ["gemmi", "numpy", "pydantic"]},
        "elapsed_seconds": time.monotonic() - start,
        "files": {p.name: sha(p.read_bytes()) for p in out.iterdir() if p.is_file()},
        "foldseek_version": retrieval["version"],
        "evidence": "experimental fixture replay plus live geometry/retrieval; no binder model run",
    }
    save(out / "manifest.json", manifest)
    return response


def gemmi_chain(residue):
    c = gemmi.Chain("selection")
    c.add_residue(residue.clone())
    return c


def renamed(chain, name):
    c = chain.clone()
    c.name = name
    return c
