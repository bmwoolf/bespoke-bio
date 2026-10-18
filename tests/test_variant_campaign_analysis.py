import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

SPEC = importlib.util.spec_from_file_location(
    "variant_analysis", Path(__file__).parents[1] / "scripts/analyze_variant_campaign.py"
)
analysis = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(analysis)


def test_rigid_rotation_and_translation_not_reflection():
    points = np.array([[0., 0., 0.], [2., 0., 0.], [0., 3., 0.], [0., 0., 4.]])
    rotation = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
    assert analysis.proper_rmsd(points @ rotation + 9, points) < 1e-12
    reflected = points.copy()
    reflected[:, 0] *= -1
    assert analysis.proper_rmsd(reflected, points) > .5


def test_cohorts_deduplicate_sequence_and_preserve_mutations(tmp_path):
    bundle = {"components": [{"id": f"b{i}", "sequence": "AAA",
                              "role": "computational_binder_candidate"} for i in range(2)],
              "evaluations": [{"subject_id": f"b{i}", "method": "native filters",
                               "status": "pass", "metrics": {"iptm": .2 + .2 * i}}
                              for i in range(2)]}
    (tmp_path / "bundle.json").write_text(json.dumps(bundle))
    result = analysis.analyze({"cohorts": [
        {"target_id": "x", "state": "reference", "sequence": "AAA"},
        {"target_id": "x", "state": "variant", "sequence": "AVA", "bundle": "bundle.json"}
    ]}, tmp_path)
    variant = result["cohorts"][1]
    assert variant["n_candidates"] == 2
    assert variant["n_unique_sequences"] == 1
    assert variant["metrics"]["iptm"]["mean"] == pytest.approx(.3)
    assert variant["mutations"] == [{"position": 2, "reference": "A", "variant": "V"}]
    assert result["fixed_binder_cross_scoring"] == []


def test_cross_scoring_never_pairs_different_sequences():
    rows = [{"target_id": "x", "binder_id": "same-label", "binder_sequence": seq,
             "designed_on": "reference", "scored_on": state, "metrics": {"iptm": score}}
            for seq, state, score in [("AAA", "reference", .4), ("AAA", "variant", .2),
                                      ("GGG", "variant", .9)]]
    result = analysis.analyze({"cross_scores": rows})["fixed_binder_cross_scoring"]
    assert len(result) == 2
    assert result[0]["mean_metric_delta_vs_reference"]["variant"]["iptm"] == pytest.approx(-.2)
    assert result[1]["mean_metric_delta_vs_reference"]["variant"] == {}


def test_residue_mapping_validates_identity_and_missing_ca(tmp_path):
    import gemmi

    st = gemmi.Structure()
    model = gemmi.Model("1")
    chain = gemmi.Chain("B")
    for number, name in [(10, "ALA"), (11, "VAL"), (12, "ALA")]:
        residue = gemmi.Residue()
        residue.name = name
        residue.seqid = gemmi.SeqId(number, " ")
        atom = gemmi.Atom()
        atom.name = "N" if number == 11 else "CA"
        atom.element = gemmi.Element("C")
        atom.pos = gemmi.Position(number, 0, 0)
        residue.add_atom(atom)
        chain.add_residue(residue)
    model.add_chain(chain)
    st.add_model(model)
    st.write_pdb(str(tmp_path / "target.pdb"))
    spec = {"path": "target.pdb", "chain": "B", "residue_positions": [1, 2, 3]}
    assert set(analysis.read_ca(spec, "AVA", tmp_path)) == {1, 3}
    with pytest.raises(ValueError, match="identity mismatch"):
        analysis.read_ca(spec, "AAA", tmp_path)


def test_seed_deltas_match_seed_and_fixed_sequence(monkeypatch):
    import importlib

    monkeypatch.syspath_prepend(str(Path(__file__).parents[1] / "scripts"))
    builder = importlib.import_module("build_variant_campaign_report")
    rows = [{"target_id": "x", "binder_sequence": "AAA", "designed_on": "reference",
             "scored_on": state, "seed": seed, "metrics": {"iptm": score}}
            for state, seed, score in [("reference", 11, .2), ("reference", 22, .4),
                                       ("variant", 22, .3), ("variant", 33, .9)]]
    result = builder.matched_seed_deltas(rows)[0]
    assert result["n_matched_seeds"] == 1
    assert result["paired_deltas"][0]["seed"] == 22
    assert result["delta_distributions"]["iptm"]["mean"] == pytest.approx(-.1)
    assert result["unmatched_reference_seeds"] == [11]
    assert result["unmatched_variant_seeds"] == [33]


def test_builder_keeps_incomplete_cohorts_explicit(tmp_path, monkeypatch):
    import importlib

    monkeypatch.syspath_prepend(str(Path(__file__).parents[1] / "scripts"))
    builder = importlib.import_module("build_variant_campaign_report")
    (tmp_path / "generation_manifest.json").write_text(json.dumps({"cohorts": [
        {"target_id": "x", "state": "reference", "sequence": "AAA",
         "directory": "generation/x", "status": "running"}]}))
    manifest, report = builder.build(tmp_path)
    assert len(manifest["cohorts"]) == 1
    assert not report["complete_generation"]
    assert report["availability"][0]["generation_status"] == "running"
    assert report["cohorts"][0]["n_candidates"] == 0
    assert not report["refold_evidence_supplied"]


def test_variant_contact_maps_by_sequence_and_ignores_hydrogen(tmp_path):
    import gemmi

    st = gemmi.Structure()
    model = gemmi.Model("1")
    for chain_name, residues in [("Q", [("GLY", 5), ("GLY", 30)]),
                                 ("Z", [("ALA", 0), ("VAL", 20)])]:
        chain = gemmi.Chain(chain_name)
        for i, (name, x) in enumerate(residues):
            residue = gemmi.Residue()
            residue.name = name
            residue.seqid = gemmi.SeqId(111 + i, " ")
            atom = gemmi.Atom()
            atom.name, atom.element = "CA", gemmi.Element("C")
            atom.pos = gemmi.Position(x, 0, 0)
            residue.add_atom(atom)
            if chain_name == "Q":
                hydrogen = gemmi.Atom()
                hydrogen.name, hydrogen.element = "H", gemmi.Element("H")
                hydrogen.pos = gemmi.Position(20, 0, 0)
                residue.add_atom(hydrogen)
            chain.add_residue(residue)
        model.add_chain(chain)
    st.add_model(model)
    path = tmp_path / "complex.pdb"
    st.write_pdb(str(path))
    result = analysis.variant_contacts(path, "AV", "GG", [1, 2])
    assert result["target_chain"] == "Z"
    assert result["binder_chain"] == "Q"
    assert result["residues"][0]["min_heavy_atom_distance_angstrom"] == 5
    assert result["residues"][0]["within_cutoff"]
    assert result["residues"][1]["min_heavy_atom_distance_angstrom"] == 10
    assert not result["residues"][1]["within_cutoff"]
    ambiguous = st[0]["Z"].clone()
    ambiguous.name = "R"
    st[0].add_chain(ambiguous)
    st.write_pdb(str(path))
    assert analysis.variant_contacts(path, "AV", "GG", [1])["status"] == "unavailable"


def test_native_confidence_requires_normalized_scale(tmp_path):
    import gemmi

    st = gemmi.Structure()
    model = gemmi.Model('1')
    chain = gemmi.Chain('X')
    for i, value in enumerate([.8, .4], 1):
        residue = gemmi.Residue()
        residue.name, residue.seqid = 'ALA', gemmi.SeqId(i, ' ')
        atom = gemmi.Atom()
        atom.name, atom.element, atom.b_iso = 'CA', gemmi.Element('C'), value
        residue.add_atom(atom)
        chain.add_residue(residue)
    model.add_chain(chain)
    st.add_model(model)
    path = tmp_path / 'confidence.pdb'
    st.write_pdb(str(path))
    result = analysis.native_fold_confidence(path, 'AA', positions=[2])
    assert result['chains']['target']['mean'] == pytest.approx(.6)
    assert result['variant_positions'][0]['normalized_plddt'] == pytest.approx(.4)
    st[0][0][0][0].b_iso = 80
    st.write_pdb(str(path))
    assert analysis.native_fold_confidence(path, 'AA')['status'] == 'unavailable'


def test_experimental_calibration_matches_sequences_not_chain_names(tmp_path, monkeypatch):
    import importlib

    import gemmi

    monkeypatch.syspath_prepend(str(Path(__file__).parents[1] / 'scripts'))
    builder = importlib.import_module('build_variant_campaign_report')
    prepared = tmp_path / 'generation/x/prepared'
    prepared.mkdir(parents=True)
    predicted = tmp_path / 'refold/seed_11/refold_cif'
    predicted.mkdir(parents=True)
    for path, chain_name, shift in [(prepared / 'target.cif', 'Z', 0),
                                    (predicted / 'x_alone.cif', 'A', 10)]:
        st = gemmi.Structure()
        model = gemmi.Model('1')
        chain = gemmi.Chain(chain_name)
        for i, (x, y, z) in enumerate([(0, 0, 0), (2, 0, 0), (0, 3, 0)], 1):
            residue = gemmi.Residue()
            residue.name, residue.seqid = 'ALA', gemmi.SeqId(i, ' ')
            residue.het_flag = 'A'
            atom = gemmi.Atom()
            atom.name, atom.element, atom.b_iso = 'CA', gemmi.Element('C'), .8
            atom.pos = gemmi.Position(x + shift, y, z)
            residue.add_atom(atom)
            chain.add_residue(residue)
        model.add_chain(chain)
        st.add_model(model)
        st.setup_entities()
        st.make_mmcif_document().write_file(str(path))
    (tmp_path / 'generation_manifest.json').write_text(json.dumps({'cohorts': [
        {'target_id': 'x', 'state': 'reference', 'sequence': 'AAA',
         'directory': 'generation/x', 'input_mode': 'observed_structure'}]}))
    (tmp_path / 'refold/receipt.json').write_text(json.dumps({
        'tasks': [{'id': 'x_alone', 'target_id': 'x', 'scored_on': 'reference',
                   'target_sequence': 'AAA'}], 'seeds': [11],
        'results': [{'id': 'x_alone', 'seed': 11, 'metrics': {'ptm': .8}}]}))
    _, report = builder.build(tmp_path, [tmp_path / 'refold'])
    result = report['experimental_target_calibration'][0]
    assert result['status'] == 'assessed'
    assert result['n_matched'] == 3
    assert result['ca_rmsd_angstrom'] < 1e-12
