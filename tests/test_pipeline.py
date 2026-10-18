import json
from pathlib import Path

import gemmi
import numpy as np
import pytest

from biocad.contracts import DesignResponse, ResidueRef, Transform
from biocad.geometry import apply, fit, inverse, pair_metrics
from biocad.mosaic_adapter import compile_request
from biocad.pipeline import run_example
from biocad.structures import prepare, select_indices, sequence, write_chains

ROOT = Path(__file__).resolve().parents[1]


def test_proper_transforms_roundtrip_and_pose():
    p = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 2.0]])
    t = Transform(rotation=[[0, -1, 0], [1, 0, 0], [0, 0, 1]], translation=[12, 3, -1])
    assert np.allclose(apply(apply(p, t), inverse(t)), p)
    fitted, rmsd = fit(p, apply(p, t))
    assert rmsd < 1e-10
    partner = np.array([[2, 4, 7.0]])
    assert np.allclose(apply(partner, fitted), apply(partner, t))
    with pytest.raises(ValueError):
        Transform(rotation=[[-1, 0, 0], [0, 1, 0], [0, 0, 1]])
    with pytest.raises(ValueError):
        fit([[0, 0, 0], [1, 0, 0], [2, 0, 0]], [[0, 0, 0], [1, 0, 0], [2, 0, 0]])


def test_selection_survives_crop_and_insertion(tmp_path):
    st = gemmi.read_structure(str(ROOT / "examples/5JDS.cif"))
    c = st[0]["A"].clone()
    # Build a discontinuous crop with an insertion-coded residue.
    cropped = gemmi.Chain("A")
    for idx in [3, 7, 11]:
        cropped.add_residue(c[idx].clone())
    cropped[1].seqid = gemmi.SeqId(99, "A")
    path = tmp_path / "crop.cif"
    write_chains([cropped], path)
    state = prepare(path, tmp_path / "prepared.cif", ["A"])
    refs = [
        ResidueRef(chain="A", number=99, insertion="A"),
        ResidueRef(chain="A", number=c[11].seqid.num),
    ]
    assert select_indices(state, refs, "A") == [1, 2]
    with pytest.raises(ValueError):
        select_indices(state, [ResidueRef(chain="A", number=-99)], "A")


def test_full_context_can_reject_domain_fit():
    binder = [[0, 0, 0]]
    cropped = [[10, 0, 0]]
    full = [[10, 0, 0], [0.2, 0, 0]]
    assert pair_metrics(binder, cropped)["severe_clash_pairs"] == 0
    assert pair_metrics(binder, full)["severe_clash_pairs"] > 0


@pytest.fixture(scope="module")
def bundle(tmp_path_factory):
    out = tmp_path_factory.mktemp("bundle")
    r = run_example(ROOT / "examples", out, foldseek=str(ROOT / "vendor/foldseek/bin/foldseek"))
    return out, r


def test_real_retrieval_identity_transfer_and_export(bundle):
    out, r = bundle
    assert r.target_states[0].chains["A"].startswith("AFTVTV")
    assert r.components[0].sequence.startswith("QVQLQ")
    assert (
        sequence(gemmi.read_structure(str(out / "component.cif"))[0][0]) == r.components[0].sequence
    )
    h = json.loads((out / "retrieval.json").read_text())["hits"][0]
    assert (
        float(h["qtmscore"]) > 0.8
    )  # Independent crystal conformations and different observed loops
    pose = json.loads((out / "transferred_pose.json").read_text())
    assert pose["ca_rmsd_angstrom"] < 2
    assert len(pose["residue_pairs"]) > 100
    DesignResponse.model_validate_json((out / "response.json").read_text())


def test_matrix_missing_results_and_rejected_assembly(bundle):
    _, r = bundle
    cid = r.components[0].id
    assert {
        e.target_state_id
        for e in r.evaluations
        if e.subject_id == cid and e.method == "heavy_atom_geometry_v1"
    } == set(r.request.target_state_ids)
    assert all(
        e.metrics["chain_pair_iptm"] is None and e.status == "missing"
        for e in r.evaluations
        if e.method == "independent_full_atom_reprediction"
    )
    bad = r.assemblies[1]
    assert any(
        e.status == "fail" and "key1:key2" in e.subject_id
        for e in r.evaluations
        if e.subject_id.startswith(bad.id)
    )
    assert [i.component_id for i in bad.instances] == [
        i.component_id for i in r.assemblies[0].instances
    ]
    for assembly in r.assemblies:
        for a, b in assembly.required_interfaces:
            assert any(e.subject_id == f"{assembly.id}/{a}:{b}" for e in r.evaluations)


def test_actual_config_change_and_compiler(bundle, tmp_path):
    _, r = bundle
    plan = compile_request(r.request, r.target_states[0])
    assert min(plan["epitope_indices_target_relative_zero_based"]) >= 0
    modified = r.request.model_copy(
        update={"terminal_clearance_angstrom": 1000.0, "second_module_offset_angstrom": 200.0}
    )
    other = run_example(
        ROOT / "examples", tmp_path, modified, foldseek=str(ROOT / "vendor/foldseek/bin/foldseek")
    )
    first = next(e for e in r.evaluations if e.method == "terminal_clearance_v1")
    second = next(e for e in other.evaluations if e.method == "terminal_clearance_v1")
    assert first.status == "pass" and second.status == "fail"
    assert (
        next(e for e in other.evaluations if e.method == "linker_upper_bound_v1").status == "fail"
    )
    assert r.components[0].id == other.components[0].id
