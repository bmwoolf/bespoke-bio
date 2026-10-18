"""Selection and coordinate-free validation preparation; no inference or networking."""

import importlib.util
import json
from pathlib import Path

import pytest

spec = importlib.util.spec_from_file_location(
    "focused_validation", Path(__file__).parents[1] / "scripts/prepare_focused_validation.py"
)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def fixtures(tmp_path):
    design = tmp_path / "design"
    native = design / "T/live/bundle"
    native.mkdir(parents=True)
    sequence = "ACDE"
    components = []
    evaluations = []
    for index, (binder, passed, quality) in enumerate(
        [("AGST", "pass", 0.4), ("GGST", "fail", 0.9), ("AAST", "pass", 0.8), ("AAST", "pass", 0.7)]
    ):
        identifier = f"binder_{index}"
        folder = native / identifier
        folder.mkdir()
        (folder / "binder.cif").write_text("data_binder\n")
        (folder / "complex.cif").write_text(f"data_original_pose_{index}\n")
        components.append(
            {
                "id": identifier,
                "sequence": binder,
                "sequence_sha256": module.sequence_hash(binder),
                "structure": identifier + "/binder.cif",
                "evidence": {
                    "complex_structure": identifier + "/complex.cif",
                    "source_sha256": module.sha(folder / "complex.cif"),
                },
            }
        )
        evaluations.append(
            {
                "subject_id": identifier,
                "target_state_id": "target",
                "method": "BoltzGen native default filters",
                "status": passed,
                "metrics": {"quality_score": quality},
            }
        )
    bundle = {
        "target_states": [{"id": "target", "chains": {"A": sequence}}],
        "components": components,
        "evaluations": evaluations,
    }
    (native / "bundle.json").write_text(json.dumps(bundle))
    prepared = design / "T/prepared"
    prepared.mkdir(parents=True)
    (prepared / "boltzgen_plan.json").write_text(
        json.dumps(
            {
                "target_sequence_sha256": module.sequence_hash(sequence),
                "num_designs": 32,
                "budget": 8,
            }
        )
    )
    controls = tmp_path / "controls.json"
    controls.write_text(
        json.dumps(
            {
                "targets": [
                    {
                        "gene_symbol": "T",
                        "sequence": sequence,
                        "sequence_sha256": module.sequence_hash(sequence),
                    }
                ]
            }
        )
    )
    standard = tmp_path / "standard"
    (standard / "msas").mkdir(parents=True)
    msa = standard / "msas/T.a3m"
    msa.write_text(">query\nACDE\n>homolog\nAC-E\n")
    (standard / "manifest.json").write_text(
        json.dumps(
            {
                "controls_sha256": module.sha(controls),
                "frozen_msas": {"T": {"path": "msas/T.a3m", "sha256": module.sha(msa)}},
            }
        )
    )
    return design, controls, standard, native


def test_selects_unique_passing_sequences_and_preserves_all_provenance(tmp_path):
    design, controls, standard, native = fixtures(tmp_path)
    output = tmp_path / "validation"
    result = module.prepare(design, controls, standard, output, genes=["T"])
    chosen = [r["candidate_id"] for r in result["selection_audit"] if r["selected"]]
    assert chosen == [
        "binder_0",
        "binder_2",
    ]  # audit retains original order, command selection is quality-ranked.
    assert result["prediction_count"] == 12
    assert result["source_provenance"]["T"]["all_candidates"] == 4
    assert result["source_provenance"]["T"]["native_passed"] == 3
    assert (output / "provenance/T/native_bundle/binder_1/complex.cif").read_bytes() == (
        native / "binder_1/complex.cif"
    ).read_bytes()
    for row in result["inputs"]:
        schema = json.loads((output / row["input"]).read_text())
        assert "templates" not in schema and "constraints" not in schema
        assert schema["sequences"][-1]["protein"]["id"] == "Z"
        assert schema["sequences"][-1]["protein"]["msa"] == "empty"
        assert row["source_native_filter"] == "pass"
        assert row["target_template"] is False
        if row["kind"] == "fixed_binder_complex":
            assert schema["sequences"][0]["protein"]["msa"] == "msas/T.a3m"
            assert row["target_chains"][0]["sequence"] == "ACDE"
        else:
            assert len(schema["sequences"]) == 1
    assert all("--use_msa_server" not in c["argv"] for c in result["commands"])


def test_no_native_pass_fails_before_publishing_inputs(tmp_path):
    design, controls, standard, native = fixtures(tmp_path)
    path = native / "bundle.json"
    bundle = json.loads(path.read_text())
    for evaluation in bundle["evaluations"]:
        evaluation["status"] = "fail"
    path.write_text(json.dumps(bundle))
    output = tmp_path / "out"
    with pytest.raises(ValueError, match="No native-filter-passing"):
        module.prepare(design, controls, standard, output, genes=["T"])
    assert not output.exists()


@pytest.mark.parametrize(
    "tamper", ["msa_bytes", "msa_query", "target_hash", "binder_hash", "complex_hash"]
)
def test_source_tampering_rejected_before_output(tmp_path, tamper):
    design, controls, standard, native = fixtures(tmp_path)
    if tamper in ["msa_bytes", "msa_query"]:
        msa = standard / "msas/T.a3m"
        msa.write_text(">q\nACDF\n")
        if tamper == "msa_query":
            p = standard / "manifest.json"
            d = json.loads(p.read_text())
            d["frozen_msas"]["T"]["sha256"] = module.sha(msa)
            p.write_text(json.dumps(d))
    elif tamper == "target_hash":
        p = design / "T/prepared/boltzgen_plan.json"
        d = json.loads(p.read_text())
        d["target_sequence_sha256"] = "wrong"
        p.write_text(json.dumps(d))
    elif tamper == "binder_hash":
        p = native / "bundle.json"
        d = json.loads(p.read_text())
        d["components"][0]["sequence_sha256"] = "wrong"
        p.write_text(json.dumps(d))
    else:
        (native / "binder_0/complex.cif").write_text("changed")
    with pytest.raises(ValueError):
        module.prepare(design, controls, standard, tmp_path / "out", genes=["T"])
    assert not (tmp_path / "out").exists()


def test_binder_selection_limit_and_tie_break_are_explicit(tmp_path):
    _, _, _, native = fixtures(tmp_path)
    bundle = json.loads((native / "bundle.json").read_text())
    chosen, audit = module.select_candidates(bundle, "ACDE", 1)
    assert [r["candidate_id"] for r in chosen] == ["binder_2"]
    assert len(audit) == 4
    assert sum(r["selected"] for r in audit) == 1


def test_positive_control_record_is_evidence_not_promoted_binder(tmp_path):
    path = tmp_path / "positive.json"
    path.write_text(
        json.dumps(
            {
                "cases": [
                    {
                        "condition": "positive_msa_control",
                        "status": "assessed",
                        "seed": 11,
                        "reference_contact_recovery": {
                            "contact_recall": 0.9,
                            "reference_interface_CA_jointfit_RMSD_A": 0.4,
                        },
                    }
                ]
            }
        )
    )
    result = module.positive_context(path)
    assert result["results"][0]["reference_contact_recall"] == 0.9
    assert "not proof" in result["scope"]
    assert result["sha256"] == module.sha(path)


def test_candidate_analyzer_marks_no_template_integrity_not_applicable(tmp_path):
    spec = importlib.util.spec_from_file_location(
        "candidate_analysis_focused",
        Path(__file__).parents[1] / "scripts/analyze_candidate_validation.py",
    )
    analyzer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(analyzer)
    row = {"target_chains": [{"chain": "A", "sequence": "ACDE"}], "target_template": False}
    result = analyzer.target_assembly_integrity({}, row, tmp_path / "nonexistent")
    assert result["status"] == "not_applicable"
    assert "No target template" in result["reason"]


def pose_fixture(tmp_path):
    import gemmi
    import numpy as np

    analyzer_spec = importlib.util.spec_from_file_location(
        "focused_pose_analyzer",
        Path(__file__).parents[1] / "scripts/analyze_candidate_validation.py",
    )
    analyzer = importlib.util.module_from_spec(analyzer_spec)
    analyzer_spec.loader.exec_module(analyzer)

    def save(path, chains, rotate=False):
        structure = gemmi.Structure()
        structure.add_model(gemmi.Model("1"))
        names = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "G": "GLY", "S": "SER", "T": "THR"}
        xyz = np.array([[0, 0, 0], [2, 0, 0], [0, 3, 0], [0, 0, 4]], dtype=float)
        for name, seq, offset in chains:
            chain = gemmi.Chain(name)
            for i, code in enumerate(seq):
                residue = gemmi.Residue()
                residue.name = names[code]
                residue.seqid = gemmi.SeqId(i + 1, " ")
                atom = gemmi.Atom()
                atom.name = "CA"
                atom.element = gemmi.Element("C")
                atom.b_iso = 90
                p = xyz[i] + offset
                if rotate:
                    p = p[[1, 0, 2]] * [-1, 1, 1] + [10, 5, -2]
                atom.pos = gemmi.Position(*p)
                residue.add_atom(atom)
                chain.add_residue(residue)
            structure[0].add_chain(chain)
        path.parent.mkdir(parents=True, exist_ok=True)
        structure.make_mmcif_document().write_file(str(path))

    reference = tmp_path / "provenance/complex.cif"
    save(reference, [("R", "ACDE", [0, 0, 0]), ("B", "AGST", [5, 0, 0])])
    prediction = tmp_path / "prediction.cif"
    save(prediction, [("A", "ACDE", [0, 0, 0]), ("Z", "AGST", [5, 0, 0])], True)
    row = {
        "kind": "fixed_binder_complex",
        "target_chains": [{"chain": "A", "sequence": "ACDE"}],
        "binder_chain": "Z",
        "binder_sequence": "AGST",
        "candidate_id": "new_binder",
        "source_design_complex_pose_EVALUATION_ONLY": "provenance/complex.cif",
        "source_design_complex_pose_sha256": module.sha(reference),
    }
    return analyzer, reference, prediction, row, save


def test_explicit_evaluation_pose_has_priority_and_recovers_new_candidate(tmp_path):
    analyzer, reference, prediction, row, _ = pose_fixture(tmp_path)
    chains = analyzer.read_chains(prediction, {"A": "ACDE", "Z": "AGST"})
    result = analyzer.reference_pose(chains, row, {"fixed_candidates": []}, tmp_path)
    assert result["status"] == "assessed"
    assert result["binder_CA_RMSD_after_target_alignment_A"] < 1e-6
    assert result["target_anchor_CA_RMSD_A"] < 1e-6
    assert result["reference"] == str(reference)
    assert "NOT accuracy" in result["interpretation"]


def test_explicit_evaluation_reference_cannot_escape_root_or_change_hash(tmp_path):
    analyzer, reference, prediction, row, _ = pose_fixture(tmp_path)
    chains = analyzer.read_chains(prediction, {"A": "ACDE", "Z": "AGST"})
    for unsafe in ["../outside.cif", str(reference)]:
        result = analyzer.reference_pose(
            chains, {**row, "source_design_complex_pose_EVALUATION_ONLY": unsafe}, None, tmp_path
        )
        assert result["status"] == "unavailable"
    result = analyzer.reference_pose(
        chains, {**row, "source_design_complex_pose_sha256": "bad"}, None, tmp_path
    )
    assert result["status"] == "unavailable"
    assert "SHA256" in result["reason"]


def test_binder_alone_gets_fold_agreement_without_pose_claim(tmp_path):
    analyzer, _, _, row, save = pose_fixture(tmp_path)
    path = tmp_path / "alone.cif"
    save(path, [("Z", "AGST", [15, 12, 7])], True)
    result = analyzer.reference_pose(
        analyzer.read_chains(path, {"Z": "AGST"}),
        {**row, "kind": "binder_alone", "target_chains": []},
        None,
        tmp_path,
    )
    assert result["status"] == "assessed"
    assert result["binder_own_fold_CA_RMSD_A"] < 1e-6
    assert "not_applicable" in result["pose_status"]
    assert "NOT experimental accuracy" in result["interpretation"]


def test_posthoc_target_reference_is_not_prediction_conditioning(tmp_path):
    analyzer, _, prediction, row, save = pose_fixture(tmp_path)
    ref = tmp_path / "evaluation/target.cif"
    save(ref, [("A", "ACDE", [0, 0, 0])])
    row["target_analysis_references_EVALUATION_ONLY"] = [
        {
            "role": "held_out",
            "pdb_id": "TEST",
            "path": "evaluation/target.cif",
            "sha256": module.sha(ref),
            "input_positions": [1, 2, 3, 4],
            "prediction_chain": "A",
            "reference_chain": "A",
        }
    ]
    result = analyzer.posthoc_target_references(
        analyzer.read_chains(prediction, {"A": "ACDE", "Z": "AGST"}), row, tmp_path
    )
    assert result["references"][0]["status"] == "assessed"
    assert result["references"][0]["global_CA_RMSD_A"] < 1e-6
    assert "not supplied" in result["references"][0]["conditioning"]


def test_experimental_references_copied_only_as_evaluation_metadata(tmp_path):
    design, controls, standard, _ = fixtures(tmp_path)
    _, reference, _, _, _ = pose_fixture(tmp_path)
    content = json.loads(controls.read_text())
    content['targets'][0]['references'] = [{
        'pdb_id': 'TEST', 'role': 'held_out', 'source': str(reference),
        'source_sha256': module.sha(reference), 'matched_structure': str(reference),
        'matched_chain': 'R', 'mapping': [
            {'input_position': i, 'observed_residue': aa, 'reference_residue': aa}
            for i, aa in enumerate('ACDE', 1)]}]
    controls.write_text(json.dumps(content))
    frozen = json.loads((standard / 'manifest.json').read_text())
    frozen['controls_sha256'] = module.sha(controls)
    (standard / 'manifest.json').write_text(json.dumps(frozen))
    output = tmp_path / 'validation'
    manifest = module.prepare(design, controls, standard, output, genes=['T'])
    for row in manifest['inputs']:
        refs = row['target_analysis_references_EVALUATION_ONLY']
        assert refs[0]['input_positions'] == [1, 2, 3, 4]
        assert module.sha(output / refs[0]['path']) == module.sha(reference)
        schema = json.loads((output / row['input']).read_text())
        assert 'templates' not in schema
        assert 'experimental_refs' not in json.dumps(schema)
    content['targets'][0]['references'][0]['mapping'][0]['observed_residue'] = 'G'
    with pytest.raises(ValueError, match='exact sequence'):
        module.experimental_references(content['targets'][0])
