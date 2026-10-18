"""Geometric support audit: evidence for preservation/redesign, not essentiality."""

import json
from pathlib import Path

import gemmi
import numpy as np

from .pipeline import save
from .structures import protein_residues, ref


def solvent_area(chain, samples=256, probe=1.4):
    atoms = [(r, a) for r in protein_residues(chain) for a in r if not a.element.is_hydrogen]
    xyz = np.array([a.pos.tolist() for _, a in atoms])
    radii = np.array([a.element.vdw_r + probe for _, a in atoms])
    k = np.arange(samples) + 0.5
    z = 1.0 - 2 * k / samples
    angle = k * np.pi * (3.0 - np.sqrt(5.0))
    sphere = np.stack(
        [np.sqrt(1 - z * z) * np.cos(angle), np.sqrt(1 - z * z) * np.sin(angle), z], axis=1
    )
    squared = ((xyz[:, None, :] - xyz[None, :, :]) ** 2).sum(-1)
    per_residue = {}
    for i, (residue, atom) in enumerate(atoms):
        neighbors = np.flatnonzero(
            (squared[i] < (radii[i] + radii) ** 2) & (np.arange(len(atoms)) != i)
        )
        points = xyz[i] + sphere * radii[i]
        exposed = np.ones(samples, dtype=bool)
        for j in neighbors:
            exposed &= ((points - xyz[j]) ** 2).sum(-1) >= radii[j] ** 2
            if not exposed.any():
                break
        area = float(exposed.mean() * 4 * np.pi * radii[i] ** 2)
        key = (residue.seqid.num, residue.seqid.icode.strip())
        per_residue[key] = per_residue.get(key, 0.0) + area
    return per_residue


def audit_core(assembly, plan_dir, output):
    assembly, plan_dir = Path(assembly), Path(plan_dir)
    plan = json.loads((plan_dir / "motif_plan.json").read_text())
    layout = json.loads((assembly / "placement_search.json").read_text())
    source = next(r for r in layout["layouts"] if r["id"] == plan["source_layout_id"])
    st = gemmi.read_structure(str(assembly / source["artifact"]))
    modules = []
    for chain_name in ["K", "L"]:
        chain = st[0][chain_name]
        residues = protein_residues(chain)
        areas = solvent_area(chain)
        kept = {
            (r["number"], r["insertion"])
            for b in plan["blocks"]
            if b["chain"] == chain_name
            for r in b["source_residues"]
        }
        support = []
        for a in residues:
            ka = (a.seqid.num, a.seqid.icode.strip())
            if ka in kept:
                continue
            ax = np.array(
                [
                    x.pos.tolist()
                    for x in a
                    if x.name not in ["N", "CA", "C", "O"] and not x.element.is_hydrogen
                ]
            )
            if not len(ax):
                continue
            neighbors = []
            for b in residues:
                kb = (b.seqid.num, b.seqid.icode.strip())
                if kb not in kept or abs(a.seqid.num - b.seqid.num) < 3:
                    continue
                bx = np.array(
                    [
                        x.pos.tolist()
                        for x in b
                        if x.name not in ["N", "CA", "C", "O"] and not x.element.is_hydrogen
                    ]
                )
                if len(bx) and np.linalg.norm(ax[:, None, :] - bx[None, :, :], axis=-1).min() < 4.5:
                    neighbors.append(ref(chain_name, b).model_dump())
            if neighbors:
                support.append(
                    {
                        "residue": ref(chain_name, a).model_dump(),
                        "name": a.name,
                        "approximate_isolated_module_sasa_angstrom2": areas[ka],
                        "retained_motif_neighbors": neighbors,
                        "core_like": areas[ka] < 10.0
                        and a.name
                        in {"VAL", "ILE", "LEU", "MET", "PHE", "TRP", "TYR", "ALA", "CYS"},
                        "policy": "preserve this packing or require its replacement by the redesigned scaffold; deletion is not validated",
                    }
                )
        modules.append(
            {
                "chain": chain_name,
                "omitted_packing_support": support,
                "core_like_omissions": [r["residue"] for r in support if r["core_like"]],
                "retained_buried_hydrophobic_residues": [
                    ref(chain_name, r).model_dump()
                    for r in residues
                    if (r.seqid.num, r.seqid.icode.strip()) in kept
                    and areas[(r.seqid.num, r.seqid.icode.strip())] < 10.0
                    and r.name in {"VAL", "ILE", "LEU", "MET", "PHE", "TRP", "TYR", "ALA", "CYS"}
                ],
            }
        )
    report = {
        "method": "256-point Shrake–Rupley-style geometric solvent exposure plus nonlocal side-chain contacts <4.5 Å",
        "probe_radius_angstrom": 1.4,
        "core_like_rule": "isolated-module residue SASA <10 Å², hydrophobic/aromatic/Cys, contacts a retained motif nonlocally",
        "evidence": "geometry_derived; not mutational essentiality",
        "source_layout_id": plan["source_layout_id"],
        "modules": modules,
        "decision": "Scaffold replacement experiment, not validated framework deletion. The 82-residue motif mask alone is insufficient evidence of preserved folding.",
        "must_preserve_or_rebuild": [
            "Both target-facing recognition geometries",
            "Observed disulfide endpoints and compatible pairing",
            "Hydrophobic packing against the recognition loops",
            "A coherent scaffold and peptide backbone",
        ],
        "limits": [
            "SASA is approximate and calculated for isolated modules",
            "Beta-sheet hydrogen-bond support is not quantified by this side-chain screen",
            "A structure cannot establish the smallest function-preserving core without redesign/prediction and experiments",
        ],
        "primary_support": "https://pmc.ncbi.nlm.nih.gov/articles/PMC5341541/",
        "numbering_note": "Our author-chain residue IDs differ from the publication immunoglobulin numbering; use the exported source mapping",
    }
    save(output, report)
    return report
