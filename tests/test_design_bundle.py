import json
from pathlib import Path

import gemmi
import pytest

from biocad.boltzgen_adapter import prepare, prepare_target
from biocad.design_bundle import DesignBundle, export_bundle

ROOT = Path(__file__).resolve().parents[1]


def test_generic_pdb_preserves_author_mapping_and_sequence(tmp_path):
    st = gemmi.read_structure(str(ROOT / "examples/5JDS.cif"))
    del st[0]["A"]
    st[0]["B"].name = "Z"
    source = tmp_path / "input.pdb"
    st.write_pdb(str(source))
    out = tmp_path / "prepared"
    plan = prepare_target(source, out, binder_length=60, num_designs=3)
    mapping = json.loads((out / "target_residue_mapping.json").read_text())
    state = json.loads((out / "target_state.json").read_text())
    assert plan["target_residues"] == 127
    assert mapping[0]["source_chain"] == "Z"
    assert mapping[0]["source_author_number"] == st[0]["Z"][0].seqid.num
    assert mapping[-1]["input_number"] == 127
    assert (out / "source.pdb").read_bytes() == source.read_bytes()
    assert state["context"]["input_mode"] == "observed_structure"
    assert len(plan["target_sequence_sha256"]) == 64
    assert plan["num_designs"] == 3


def test_sequence_only_is_explicit_and_multichain_requires_selection(tmp_path):
    fasta = tmp_path / "target.fasta"
    fasta.write_text(">de_novo\nACDEFGHIKLMNPQRSTVWY\n")
    out = tmp_path / "sequence"
    plan = prepare_target(fasta, out)
    spec = json.loads((out / "design.yaml").read_text())
    assert spec["entities"][1]["protein"]["sequence"] == "ACDEFGHIKLMNPQRSTVWY"
    assert plan["target_input_mode"] == "sequence_only"
    assert not (out / "target.cif").exists()
    assert (
        json.loads((out / "target_state.json").read_text())["context"][
            "sequence_only_prediction_performed"
        ]
        is False
    )
    with pytest.raises(ValueError, match="Select an author chain"):
        prepare_target(ROOT / "examples/5JDS.cif", tmp_path / "ambiguous")


def test_actual_candidates_export_typed_evidence_and_refuse_target_mismatch(tmp_path):
    prep = tmp_path / "prep"
    prepare(ROOT / "examples", prep)
    bundle = export_bundle(
        ROOT / "examples/boltzgen_results", prep / "target_state.json", tmp_path / "bundle"
    )
    typed = DesignBundle.model_validate(bundle)
    assert len(typed.components) == 2
    filters = [e for e in typed.evaluations if e.id.endswith("_filters")]
    assert sorted(e.status for e in filters) == ["fail", "pass"]
    bindings = [e for e in typed.evaluations if e.id.endswith("_binding")]
    assert all(e.status == "assessed" for e in bindings)
    assert all(e.metrics["design_to_target_iptm"] < 0.31 for e in bindings)
    assert all(e.metrics["esmfold2_ipsae_min"] < 0.06 for e in bindings)
    assert all(len(c.ports) == 2 and c.sequence_sha256 for c in typed.components)
    state = json.loads((prep / "target_state.json").read_text())
    state["chains"]["A"] = "A" * 115
    (prep / "target_state.json").write_text(json.dumps(state))
    with pytest.raises(ValueError, match="mismatches"):
        export_bundle(
            ROOT / "examples/boltzgen_results", prep / "target_state.json", tmp_path / "bad"
        )
