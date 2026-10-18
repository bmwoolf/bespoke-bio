from pathlib import Path

import gemmi
import numpy as np
import pytest

from biocad.contracts import PlacementSearchRequest, Transform
from biocad.geometry import apply, distances
from biocad.placement import cross_screen, run_search
from biocad.structures import atoms, sequence

ROOT = Path(__file__).resolve().parents[1]


def test_vdw_overlap_is_stricter_than_coincident_atom_check():
    # Two carbons 2.8 Å apart do not trigger the old 1.5 Å severe-clash rule,
    # but overlap their vdW spheres by 0.6 Å and must fail compact packing.
    result = cross_screen(
        np.array([[0.0, 0.0, 0.0]]),
        np.array([[2.8, 0.0, 0.0]]),
        np.array([1.7]),
        np.array([1.7]),
        0.4,
    )
    assert result["vdw_overlap_pairs"] == 1
    assert result["minimum_vdw_gap_angstrom"] == pytest.approx(-0.6)


def test_rigid_search_preserves_binding_and_rejects_impossible_connection(tmp_path):
    request = PlacementSearchRequest(offsets_angstrom=[40.0], z_rotations_degrees=[0.0, 90.0])
    with np.errstate(all="raise"):
        result = run_search(ROOT / "examples/replay", tmp_path / "normal", request)
    assert result["eligible_layouts"] > 0
    assert result["layouts"][0]["kind"] == "overlap_control"
    assert result["layouts"][0]["status"] == "rejected"
    key = gemmi.read_structure(str(tmp_path / "normal/key.cif"))[0][0].clone()
    target = gemmi.read_structure(str(tmp_path / "normal/target.cif"))[0][0].clone()
    k, t = atoms(key, ca=True)[0], atoms(target, ca=True)[0]
    for row in result["layouts"]:
        tr = Transform(**row["transform"])
        assert np.allclose(distances(k, t), distances(apply(k, tr), apply(t, tr)), atol=1e-10)
        if row["id"] not in result["top_layout_ids"]:
            continue
        assert row["status"] == "passes_geometry_screen"
        assert all(m["vdw_overlap_pairs"] == 0 for m in row["forbidden_pairs"].values())
        st = gemmi.read_structure(str(tmp_path / "normal" / row["artifact"]))
        assert sequence(st[0]["K"]) == sequence(st[0]["L"]) == sequence(key)
        assert np.allclose(atoms(st[0]["L"])[0], apply(atoms(key)[0], tr), atol=0.001)
        assert (
            row["compactness"]["key_ca_radius_of_gyration_angstrom"]
            < result["spread_out_reference"]["compactness"]["key_ca_radius_of_gyration_angstrom"]
        )
    short = request.model_copy(update={"linker_sequence": "G"})
    other = run_search(ROOT / "examples/replay", tmp_path / "short", short)
    assert other["eligible_layouts"] == 0
    assert result["component_id"] == other["component_id"]
    assert all(row["confidence"] is None for row in other["layouts"])
    assert (tmp_path / "short/review.html").exists()


def test_search_budget_and_linker_validation():
    with pytest.raises(ValueError):
        PlacementSearchRequest(offsets_angstrom=[float("nan")])
    with pytest.raises(ValueError):
        PlacementSearchRequest(linker_sequence="G*G")
    with pytest.raises(ValueError):
        PlacementSearchRequest(offsets_angstrom=[10.0] * 8, z_rotations_degrees=[0.0] * 8)
