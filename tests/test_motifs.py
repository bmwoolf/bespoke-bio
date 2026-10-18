import json
from itertools import pairwise
from pathlib import Path

import gemmi
import numpy as np
import pytest

from biocad.contracts import MotifPlanRequest, Transform
from biocad.motifs import check_motif_geometry, prepare_plan
from biocad.retrieval import move
from biocad.structures import write_chains

ROOT = Path(__file__).resolve().parents[1]


def test_motif_plan_preserves_support_and_exact_residue_maps(tmp_path):
    plan = prepare_plan(ROOT / "examples/compact_assembly", tmp_path / "plan")
    assert len(plan["native_disulfides"]) == 4
    normalized = gemmi.read_structure(str(tmp_path / "plan/model_input.cif"))
    source = gemmi.read_structure(str(ROOT / "examples/compact_assembly/compact_1.cif"))
    mapping = json.loads((tmp_path / "plan/residue_mapping.json").read_text())
    for bond in plan["native_disulfides"]:
        for end in ["first", "second"]:
            assert any(bond[end] in block["source_residues"] for block in plan["blocks"])
    for entry in mapping:
        src = entry["source"]
        a = next(
            r
            for r in source[0][src["chain"]]
            if r.seqid.num == src["number"] and r.seqid.icode.strip() == src["insertion"]
        )
        b = normalized[0][entry["input_chain"]][entry["input_number"] - 1]
        assert a.name == b.name
        assert np.array_equal(
            [[x.pos.x, x.pos.y, x.pos.z] for x in a], [[x.pos.x, x.pos.y, x.pos.z] for x in b]
        )
    picked = [t for t in plan["topologies"] if t["id"] in plan["selected_topology_ids"]]
    assert {t["provenance_order"] for t in picked} == {"tandem_control"}
    assert len({t["scaffold_length"] for t in picked}) == 1
    # A motif selected for one module appears exactly once in each chain proposal.
    for topology in plan["topologies"]:
        assert sorted(topology["order"]) == sorted(b["id"] for b in plan["blocks"])
        spans = topology["output_motif_spans"]
        assert all(a["output_last"] < b["output_first"] for a, b in pairwise(spans))
    short = prepare_plan(
        ROOT / "examples/compact_assembly",
        tmp_path / "short",
        MotifPlanRequest(bridge_min_residues=1, bridge_max_residues=1),
    )
    assert short["locally_viable_topologies"] == 0
    assert json.loads((tmp_path / "short/rfd3_inputs.draft.json").read_text()) == {}


def synthetic_motif_fixture(plan_dir, plan):
    # Geometry-only control containing fixed motifs and targets, no invented bridges.
    source = gemmi.read_structure(str(plan_dir / "model_input.cif"))
    topology = next(t for t in plan["topologies"] if t["id"] == plan["selected_topology_ids"][0])
    scaffold = gemmi.Chain("A")
    for span in topology["output_motif_spans"]:
        for src, dst in zip(
            range(span["input_first"], span["input_last"] + 1),
            range(span["output_first"], span["output_last"] + 1),
            strict=True,
        ):
            residue = source[0][span["input_chain"]][src - 1].clone()
            residue.seqid = gemmi.SeqId(dst, " ")
            scaffold.add_residue(residue)
    targets = []
    for name, dst in [("T", "B"), ("U", "C")]:
        c = source[0][name].clone()
        c.name = dst
        targets.append(c)
    return [scaffold] + targets, topology["id"]


def test_target_relative_validator_detects_movement_and_missing_motifs(tmp_path):
    plan_dir = tmp_path / "plan"
    plan = prepare_plan(ROOT / "examples/compact_assembly", plan_dir)
    chains, tid = synthetic_motif_fixture(plan_dir, plan)
    rigid = Transform(rotation=[[0, -1, 0], [1, 0, 0], [0, 0, 1]], translation=[20.0, -4.0, 12.0])
    write_chains([move(c, rigid) for c in chains], tmp_path / "rigid.cif")
    assert check_motif_geometry(plan_dir, tmp_path / "rigid.cif", tid)["status"] == "pass"
    assert check_motif_geometry(plan_dir, tmp_path / "rigid.cif", tid)["scope"].startswith(
        "Motif identity"
    )
    first = chains[0][0]
    for atom in first:
        atom.pos = gemmi.Position(atom.pos.x + 10, atom.pos.y, atom.pos.z)
    write_chains(chains, tmp_path / "moved.cif")
    assert check_motif_geometry(plan_dir, tmp_path / "moved.cif", tid)["status"] == "fail"
    del chains[0][-1]
    write_chains(chains, tmp_path / "missing.cif")
    result = check_motif_geometry(plan_dir, tmp_path / "missing.cif", tid)
    assert any(
        b["status"] == "missing" and b["fixed_atom_rmsd_angstrom"] is None for b in result["blocks"]
    )


def test_motif_budget_contract():
    with pytest.raises(ValueError):
        MotifPlanRequest(bridge_min_residues=10, bridge_max_residues=3)
