"""Synthetic geometry exercises native-output parsing without model inference."""

import importlib.util
import json
from pathlib import Path

import gemmi
import numpy as np
import pytest

spec = importlib.util.spec_from_file_location(
    "analyze_standard", Path(__file__).parents[1] / "scripts/analyze_standard_validation.py"
)
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)


def write_cif(path, chains, transform=False):
    structure = gemmi.Structure()
    structure.add_model(gemmi.Model("1"))
    xyz = [[0, 0, 0], [2, 0, 0], [0, 3, 0], [0, 0, 4]]
    for name, sequence, positions, confidence in chains:
        chain = gemmi.Chain(name)
        for i, (letter, pos) in enumerate(zip(sequence, positions)):
            residue = gemmi.Residue()
            residue.name = (
                gemmi.find_tabulated_residue(letter).name
                if len(letter) > 1
                else {"A": "ALA", "G": "GLY", "S": "SER", "T": "THR"}[letter]
            )
            residue.seqid = gemmi.SeqId(pos, " ")
            atom = gemmi.Atom()
            atom.name = "CA"
            atom.element = gemmi.Element("C")
            point = np.array(xyz[i % 4], dtype=float)
            if transform:
                point = point[[1, 0, 2]] * [-1, 1, 1] + [20, 5, -4]
            atom.pos = gemmi.Position(*point)
            atom.b_iso = confidence
            residue.add_atom(atom)
            chain.add_residue(residue)
        structure[0].add_chain(chain)
    path.parent.mkdir(parents=True, exist_ok=True)
    structure.make_mmcif_document().write_file(str(path))


def case(sequence="AGST", binder=None):
    return {
        "id": "T",
        "sequence": sequence,
        "binder_sequence": binder,
        "target_id": "T",
        "kind": "target_alone",
    }


def fixture_prediction(root, *, sequence="AGST", binder=None, plddt=None):
    chains = [("A", sequence, list(range(1, len(sequence) + 1)), 80)]
    if binder:
        chains.append(("B", binder, list(range(1, len(binder) + 1)), 60))
    write_cif(root / "T_model_0.cif", chains, transform=True)
    (root / "confidence_T_model_0.json").write_text(
        json.dumps({"ptm": 0.8, "iptm": 0.4, "chains_ptm": {"0": 0.85, "1": 0.65}})
    )
    if plddt is None:
        plddt = [0.8] * len(sequence) + ([0.6] * len(binder) if binder else [])
    np.savez(root / "plddt_T_model_0.npz", plddt=plddt)
    return chains


def reference(root, sequence="AGST", positions=(1, 2, 3, 4), role="held_out"):
    path = root / "reference.cif"
    write_cif(path, [("A", sequence, list(positions), 80)])
    return {
        "role": role,
        "pdb_id": "TEST",
        "matched_structure": str(path),
        "matched_chain": "A",
        "mapping": [{"input_position": p} for p in positions],
    }


def test_correct_chain_units_confidence_and_rigid_fit(tmp_path):
    fixture_prediction(tmp_path / "prediction", binder="AGST")
    ref = reference(tmp_path)
    result = analysis.parse_prediction(
        tmp_path / "prediction",
        case(binder="AGST"),
        {"sequence": "AGST", "references": [ref]},
        tmp_path,
    )
    assert result["status"] == "completed"
    assert result["chains"]["A"]["plddt"]["mean"] == pytest.approx(80)
    assert result["chains"]["B"]["plddt"]["mean"] == pytest.approx(60)
    assert result["chains"]["B"]["ptm"] == 0.65
    assert result["reference_comparisons"][0]["rmsd_angstrom"] < 1e-6
    assert result["reference_comparisons"][0]["compared_input_positions"] == [1, 2, 3, 4]


def test_exact_sequence_mismatch_is_invalid(tmp_path):
    fixture_prediction(tmp_path, sequence="AGSA")
    result = analysis.parse_prediction(
        tmp_path, case(), {"sequence": "AGST", "references": []}, tmp_path
    )
    assert result["status"] == "invalid_artifact"
    assert "sequence mismatch" in result["error"]


def test_wrong_npz_scale_not_silently_rescaled(tmp_path):
    fixture_prediction(tmp_path, plddt=[80] * 4)
    result = analysis.parse_prediction(
        tmp_path, case(), {"sequence": "AGST", "references": []}, tmp_path
    )
    assert result["status"] == "invalid_artifact"
    assert "0..1" in result["error"]


def test_missing_results_and_missing_confidence_explicit(tmp_path):
    result = analysis.parse_prediction(
        tmp_path, case(), {"sequence": "AGST", "references": []}, tmp_path
    )
    assert result["status"] == "missing_results"
    write_cif(tmp_path / "T_model_0.cif", [("A", "AGST", [1, 2, 3, 4], 100)])
    result = analysis.parse_prediction(
        tmp_path, case(), {"sequence": "AGST", "references": []}, tmp_path
    )
    assert result["status"] == "incomplete_results"
    assert result["chains"]["A"]["plddt"]["mean"] is None  # Default CIF B=100 is not confidence.


def test_reference_gaps_use_positions_not_repeated_sequence_search(tmp_path):
    ref = reference(tmp_path, sequence="AST", positions=[1, 3, 4])
    prediction = {
        "ca": {
            1: np.array([0, 0, 0]),
            2: np.array([999, 999, 999]),
            3: np.array([2, 0, 0]),
            4: np.array([0, 3, 0]),
        }
    }
    result = analysis.compare_reference(prediction, ref, "AGST", tmp_path)
    assert result["compared_input_positions"] == [1, 3, 4]
    assert result["rmsd_angstrom"] < 1e-6
    assert result["input_coverage"] == 0.75


def test_template_conditioned_label_never_independent(tmp_path):
    ref = reference(tmp_path, role="template")
    pred = analysis.chain_data(tmp_path / "reference.cif")["A"]
    result = analysis.compare_reference(pred, ref, "AGST", tmp_path)
    assert result["scope"] == "template-conditioned recovery"


def test_complete_grid_includes_missing_seeds_and_never_imputes(tmp_path):
    import hashlib

    controls = tmp_path / "controls.json"
    controls.write_text(
        json.dumps({"targets": [{"gene_symbol": "T", "sequence": "AGST", "references": []}]})
    )
    prepared = tmp_path / "prepared"
    prepared.mkdir()
    conditions = [
        "empty_msa__no_template",
        "empty_msa__target_template",
        "server_msa__no_template",
        "server_msa__target_template",
    ]
    inputs = []
    commands = []
    for condition in conditions:
        inputs.append(
            {
                **case(),
                "condition": condition,
                "msa": condition.split("_")[0],
                "target_template": "__target_template" in condition,
            }
        )
        for seed in [11, 22, 33]:
            commands.append(
                {
                    "condition": condition,
                    "seed": seed,
                    "expected_prediction_directory": f"results/{condition}/seed_{seed}/predictions",
                }
            )
    manifest = {
        "package": "boltz==2.2.1",
        "controls_sha256": hashlib.sha256(controls.read_bytes()).hexdigest(),
        "inputs": inputs,
        "commands": commands,
        "limitations": [],
    }
    (prepared / "manifest.json").write_text(json.dumps(manifest))
    fixture_prediction(prepared / "results" / conditions[0] / "seed_11/predictions/T")
    result = analysis.analyze(controls, prepared, tmp_path / "analysis", base=tmp_path)
    assert result["expected_predictions"] == 12
    assert result["completed_predictions"] == 1
    assert len(result["summary"]) == 4
    assert result["summary"][0]["target_plddt"]["n"] == 1
    assert (tmp_path / "analysis/predictions.csv").exists()
    assert "not independent" in (tmp_path / "analysis/README.md").read_text()


def test_alignment_does_not_hide_mirror_reflection(tmp_path):
    ref = reference(tmp_path)
    pred = analysis.chain_data(tmp_path / "reference.cif")["A"]
    pred["ca"] = {p: xyz * [-1, 1, 1] for p, xyz in pred["ca"].items()}
    result = analysis.compare_reference(pred, ref, "AGST", tmp_path)
    assert result["rmsd_angstrom"] > 0.5


def test_pae_direction_and_chain_confidence_are_distinct(tmp_path):
    fixture_prediction(tmp_path, binder="AGST")
    pae = np.ones((8, 8))
    pae[:4, 4:] = 3
    pae[4:, :4] = 7
    np.savez(tmp_path / "pae_T_model_0.npz", pae=pae)
    result = analysis.parse_prediction(
        tmp_path, case(binder="AGST"), {"sequence": "AGST", "references": []}, tmp_path
    )
    assert result["status"] == "completed"
    assert result["pae_angstrom"]["target_to_binder_mean"] == 3
    assert result["pae_angstrom"]["binder_to_target_mean"] == 7
    assert result["chains"]["A"]["ptm"] != result["chains"]["B"]["ptm"]


@pytest.mark.parametrize("contents", ["[]", '{"ptm": NaN}'])
def test_malformed_confidence_is_invalid_without_breaking_report(tmp_path, contents):
    fixture_prediction(tmp_path)
    (tmp_path / "confidence_T_model_0.json").write_text(contents)
    result = analysis.parse_prediction(
        tmp_path, case(), {"sequence": "AGST", "references": []}, tmp_path
    )
    assert result["status"] == "invalid_artifact"


def test_prior_custom_bfactor_scale_differs_from_standard_boltz(tmp_path):
    directory = (
        tmp_path
        / "T_reference/live/downloaded/results/intermediate_designs_inverse_folded/refold_cif"
    )
    write_cif(directory / "sample.cif", [("A", "AGST", [1, 2, 3, 4], 0.8)])
    config = directory.parents[1] / "config/folding.yaml"
    config.parent.mkdir(parents=True, exist_ok=True)
    config.write_text(
        "data:\n  _target_: boltzgen.task.predict.data_from_generated.FromGeneratedDataModule\n  target_templates: true\n"
    )
    rows = analysis.generation_context(
        tmp_path, {"T": {"sequence": "AGST", "references": []}}, tmp_path
    )
    assert rows[0]["status"] == "compared"
    assert rows[0]["target_CA_bfactor"]["mean"] == pytest.approx(0.8)
    assert rows[0]["target_plddt_0_100"]["mean"] == pytest.approx(80)
    assert rows[0]["target_CA_bfactor_native_scale"].startswith("0..1")


def test_actual_custom_comparison_selected_from_receipt_not_generation(tmp_path):
    root = tmp_path / "crossfold"
    root.mkdir()
    task = {
        "id": "T_reference_alone",
        "target_id": "T",
        "target_sequence": "AGST",
        "scored_on": "reference",
    }
    receipt = {
        "mode": "fixed_sequence_refolding",
        "target_templates": False,
        "sequence_redesign": False,
        "checkpoint": "custom_checkpoint",
        "tasks": [task],
        "results": [{"id": task["id"], "seed": 11}],
    }
    (root / "receipt.json").write_text(json.dumps(receipt))
    write_cif(root / "seed_11/refold_cif/T_reference_alone.cif", [("A", "AGST", [1, 2, 3, 4], 0.3)])
    archive = root / "seed_11/fold_out_npz"
    archive.mkdir()
    np.savez(archive / "T_reference_alone.npz", ptm=np.array([0.4]))
    rows = analysis.crossfold_context(root, {"T": {"sequence": "AGST", "references": []}}, tmp_path)
    assert len(rows) == 1
    assert rows[0]["kind"] == "custom_target_alone"
    assert rows[0]["target_templates"] is False
    assert rows[0]["target_plddt_0_100"]["mean"] == pytest.approx(30)
    assert rows[0]["metrics"]["ptm"] == pytest.approx(0.4)
    assert "not standard Boltz CLI" in rows[0]["protocol"]
    receipt["target_templates"] = True
    (root / "receipt.json").write_text(json.dumps(receipt))
    with pytest.raises(ValueError, match="does not verify"):
        analysis.crossfold_context(root, {"T": {"sequence": "AGST", "references": []}}, tmp_path)
