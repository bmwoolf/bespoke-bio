import importlib.util
import sys
from pathlib import Path

import gemmi
import pytest

sys.path.insert(0, str(Path(__file__).parents[1] / 'scripts'))
SPEC = importlib.util.spec_from_file_location(
    'prepare_candidates', Path(__file__).parents[1] / 'scripts/prepare_candidate_validation.py')
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def assembly_fixture(tmp_path):
    structure, model = gemmi.Structure(), gemmi.Model('1')
    specs = []
    for name, positions, offset in [('A', [1, 2, 3], 0), ('B', [2, 3], 20)]:
        chain = gemmi.Chain(name)
        for position in positions:
            residue = gemmi.Residue()
            residue.name = ['ALA', 'CYS', 'ASP'][position - 1]
            residue.seqid = gemmi.SeqId(position, ' ')
            atom = gemmi.Atom()
            atom.name, atom.element = 'CA', gemmi.Element('C')
            atom.pos = gemmi.Position(position, offset, 0)
            residue.add_atom(atom)
            chain.add_residue(residue)
        model.add_chain(chain)
        specs.append({'chain': name, 'sequence': 'ACD',
                      'mapping': [{'input_position': p} for p in positions]})
    structure.add_model(model)
    source = tmp_path / 'source.cif'
    structure.make_mmcif_document().write_file(str(source))
    return source, specs


def test_normalized_assembly_retains_relative_coordinates_and_missing_masks(tmp_path):
    source, specs = assembly_fixture(tmp_path)
    destination = tmp_path / 'normalized.cif'
    masks = module.normalize_assembly(source, specs, destination)
    assert masks == {'A': [1, 2, 3], 'B': [2, 3]}
    parsed = gemmi.read_structure(str(destination))
    assert len(parsed[0]) == 2
    assert all(e.full_sequence == ['ALA', 'CYS', 'ASP'] for e in parsed.entities)
    assert parsed[0]['A'][1][0].pos.dist(parsed[0]['B'][0][0].pos) == pytest.approx(20)
    assert parsed[0]['B'][0].label_seq == 2
    assert sum(len(c) for c in parsed[0]) == 5  # no missing coordinates fabricated
    with pytest.raises(ValueError, match='target chains only'):
        module.normalize_assembly(source, specs[:1], destination)


def test_template_mapping_mismatch_rejected(tmp_path):
    source, specs = assembly_fixture(tmp_path)
    specs[1]['mapping'] = [{'input_position': 1}, {'input_position': 2}]
    with pytest.raises(ValueError, match='audited mapping'):
        module.normalize_assembly(source, specs, tmp_path / 'bad.cif')
