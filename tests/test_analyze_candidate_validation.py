import importlib.util
import json
from pathlib import Path

import gemmi
import numpy as np

SPEC = importlib.util.spec_from_file_location(
    'analyze_candidates', Path(__file__).parents[1] / 'scripts/analyze_candidate_validation.py')
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_oligomer_interface_scores_cannot_mix_target_target_with_binder():
    confidence = {'pair_chains_iptm': {
        '0': {'1': .99, '2': .10}, '1': {'0': .98, '2': .20},
        '2': {'0': .30, '1': .40}}}
    binder, targets = module.pair_scores(confidence, ['A', 'B', 'Z'], 'Z')
    assert binder == {'A|Z': .10, 'B|Z': .20, 'Z|A': .30, 'Z|B': .40}
    assert targets == {'A|B': .99, 'B|A': .98}


def test_native_scale_directional_pae_and_sequence_integrity(tmp_path):
    base = tmp_path / 'predictions/case'
    base.mkdir(parents=True)
    structure, model = gemmi.Structure(), gemmi.Model('1')
    for name, code, x in [('A', 'ALA', 0), ('Z', 'CYS', 4)]:
        chain, residue = gemmi.Chain(name), gemmi.Residue()
        residue.name, residue.seqid = code, gemmi.SeqId(1, ' ')
        atom = gemmi.Atom()
        atom.name, atom.element = 'CA', gemmi.Element('C')
        atom.pos, atom.b_iso = gemmi.Position(x, 0, 0), 80
        residue.add_atom(atom)
        chain.add_residue(residue)
        model.add_chain(chain)
    structure.add_model(model)
    structure.make_mmcif_document().write_file(str(base / 'case_model_0.cif'))
    (base / 'confidence_case_model_0.json').write_text(json.dumps({
        'iptm': .9, 'pair_chains_iptm': {'0': {'1': .2}, '1': {'0': .4}}}))
    np.savez(base / 'plddt_case_model_0.npz', plddt=np.array([.8, .8]))
    np.savez(base / 'pae_case_model_0.npz', pae=np.array([[1, 2], [7, 1]]))
    row = {'id': 'case', 'target_chains': [{'chain': 'A', 'sequence': 'A'}],
           'binder_chain': 'Z', 'binder_sequence': 'C'}
    command = {'seed': 11, 'expected_prediction_directory': 'predictions'}
    result = module.analyze_case(tmp_path, row, command)
    assert result['status'] == 'analyzed'
    assert result['binder_plddt_0_100']['mean'] == 80
    interface = result['binder_target_interfaces']['A']
    assert interface['residue_contact_pairs'] == 1
    assert interface['contact_pae_binder_rows_target_columns_A']['mean'] == 7
    assert interface['contact_pae_target_rows_binder_columns_A']['mean'] == 2
    row['binder_sequence'] = 'D'
    assert module.analyze_case(tmp_path, row, command)['status'] == 'unavailable'


def test_one_global_assembly_fit_detects_subunit_rearrangement():
    a = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    reference = {'A': a, 'B': a + [10, 0, 0]}
    rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    rigid = {name: xyz @ rotation + [7, -3, 2] for name, xyz in reference.items()}
    fit = module.assembly_alignment(rigid, reference)
    assert fit['global_CA_RMSD_A'] < 1e-10
    moved = {'A': a, 'B': reference['B'] + [0, 5, 0]}
    changed = module.assembly_alignment(moved, reference)
    assert changed['global_CA_RMSD_A'] > .5
    assert changed['per_chain']['A']['CA_RMSD_after_global_fit_A'] > .5
    assert changed['per_chain']['B']['CA_RMSD_after_global_fit_A'] > .5
    # Each isolated subunit is still unchanged; separate fits would hide assembly change.
    assert module.assembly_alignment({'B': moved['B']}, {'B': reference['B']})['global_CA_RMSD_A'] < 1e-10


def test_heavy_overlap_is_strict_and_excludes_hydrogen():
    def residue(x, hydrogen=False):
        result = gemmi.Residue()
        atom = gemmi.Atom()
        atom.name = 'H' if hydrogen else 'CA'
        atom.element = gemmi.Element('H' if hydrogen else 'C')
        atom.pos = gemmi.Position(x, 0, 0)
        atom.occ = 1
        result.add_atom(atom)
        return result
    overlap = module.heavy_atom_overlap([residue(0), residue(10, True)],
                                        [residue(1.9), residue(2), residue(10)])
    assert overlap['atom_pairs_below_cutoff'] == 1
    assert overlap['binder_residues_below_cutoff'] == 1
    assert overlap['minimum_heavy_atom_distance_A'] == 1.9


def test_target_anchored_pose_detects_relocation_despite_identical_binder_fold():
    target = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    binder = target + [5, 0, 0]
    relocated = binder + [0, 12, 0]
    result = module.aligned_pose_metrics(target, target, relocated, binder)
    assert result['target_CA_RMSD_A'] < 1e-10
    assert result['binder_own_fold_CA_RMSD_A'] < 1e-10
    assert abs(result['binder_CA_RMSD_after_target_alignment_A'] - 12) < 1e-10
    rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    exact = module.aligned_pose_metrics(target @ rotation + 3, target,
                                        binder @ rotation + 3, binder)
    assert exact['binder_CA_RMSD_after_target_alignment_A'] < 1e-10


def test_reference_contact_recovery_and_strict_five_angstrom_boundary():
    def residues(code, positions):
        output = []
        for index, x in enumerate(positions, 1):
            residue = gemmi.Residue()
            residue.name, residue.seqid = code, gemmi.SeqId(index, ' ')
            atom = gemmi.Atom()
            atom.name, atom.element = 'CA', gemmi.Element('C')
            atom.pos, atom.occ = gemmi.Position(x, 0, 0), 1
            residue.add_atom(atom)
            output.append(residue)
        return output
    target, binder = residues('ALA', [0, 10]), residues('CYS', [1, 11])
    exact = module.contact_recovery(target, binder, target, binder)
    assert exact['reference_residue_contacts'] == 2
    assert exact['contact_recall'] == exact['contact_precision'] == exact['contact_jaccard'] == 1
    assert exact['reference_interface_CA_jointfit_RMSD_A'] < 1e-10
    moved = module.contact_recovery(target, residues('CYS', [1, 4.9]), target, binder)
    assert moved['recovered_residue_contacts'] == 1
    assert moved['contact_recall'] == moved['contact_precision'] == .5
    assert abs(moved['contact_jaccard'] - 1/3) < 1e-10
    boundary = module.contact_recovery(target, residues('CYS', [1, 5]), target, binder)
    assert boundary['predicted_residue_contacts'] == 1  # exactly5Å excluded in both directions
