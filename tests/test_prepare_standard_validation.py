import importlib.util
import json
from pathlib import Path

import gemmi
import pytest

SPEC = importlib.util.spec_from_file_location(
    'prepare_standard', Path(__file__).parents[1] / 'scripts/prepare_standard_validation.py')
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def template_fixture(tmp_path):
    structure = gemmi.Structure()
    model, chain = gemmi.Model('1'), gemmi.Chain('A')
    for position, code in [(2, 'CYS'), (3, 'ASP')]:
        residue = gemmi.Residue()
        residue.name = code
        residue.seqid = gemmi.SeqId(position, ' ')
        atom = gemmi.Atom()
        atom.name, atom.element = 'CA', gemmi.Element('C')
        atom.pos = gemmi.Position(position, 0, 0)
        residue.add_atom(atom)
        chain.add_residue(residue)
    model.add_chain(chain)
    structure.add_model(model)
    source = tmp_path / 'source.cif'
    structure.make_mmcif_document().write_file(str(source))
    ref = {'role': 'template', 'matched_structure': str(source), 'matched_chain': 'A',
           'mapping': [{'input_position': 2}, {'input_position': 3}]}
    controls = tmp_path / 'controls.json'
    controls.write_text(json.dumps({'targets': [
        {'gene_symbol': 'X', 'sequence': 'ACD', 'references': [ref]}]}))
    return controls, ref


def test_partial_template_preserves_unresolved_position_and_full_sequence(tmp_path):
    _, reference = template_fixture(tmp_path)
    destination = tmp_path / 'normalized.cif'
    audit = module.normalized_template(reference, 'ACD', destination, tmp_path)
    parsed = gemmi.read_structure(str(destination))
    assert parsed.entities[0].full_sequence == ['ALA', 'CYS', 'ASP']
    assert [(r.seqid.num, r.label_seq) for r in parsed[0][0]] == [(2, 2), (3, 3)]
    assert audit['unobserved_positions'] == [1]
    assert sum(len(r) for r in parsed[0][0]) == 2
    with pytest.raises(ValueError, match='disagrees'):
        module.normalized_template(reference, 'AAD', destination, tmp_path)


def test_factorial_matched_seeds_and_target_only_complex_template(tmp_path):
    controls, _ = template_fixture(tmp_path)
    output = tmp_path / 'prepared'
    manifest = module.prepare(controls, output, seeds=[11, 22],
                              binders=[{'id': 'candidate', 'target_id': 'X', 'sequence': 'EFG'}])
    assert len(manifest['inputs']) == 8
    assert manifest['prediction_count'] == 16
    assert len(manifest['commands']) == 8  # conditions batch target + binder case
    for row in manifest['inputs']:
        schema = json.loads((output / row['input']).read_text())
        assert schema['sequences'][0]['protein']['sequence'] == 'ACD'
        assert ('msa' in schema['sequences'][0]['protein']) == (row['msa'] == 'empty')
        if row['binder_sequence']:
            assert schema['sequences'][1]['protein'] == {'id': 'B', 'sequence': 'EFG', 'msa': 'empty'}
        if row['target_template']:
            assert schema['templates'] == [{'cif': 'templates/X.cif', 'chain_id': 'A',
                                            'template_id': 'A', 'force': False}]
        else:
            assert 'templates' not in schema
    for command in manifest['commands']:
        assert '--write_full_pae' in command['argv']
        assert '--use_potentials' not in command['argv']
    with pytest.raises(ValueError, match='new or empty'):
        module.prepare(controls, output)


def test_frozen_msa_reused_across_arms_and_seeds_without_gpu(tmp_path, monkeypatch):
    import sys
    import types

    controls, _ = template_fixture(tmp_path)
    output = tmp_path / 'prepared'
    module.prepare(controls, output, seeds=[11, 22])
    calls = []

    def fake_msa(sequences, prefix, **kwargs):
        calls.append(sequences)
        return ['>query\nACD\n>homolog\nA-D\n']

    monkeypatch.setitem(sys.modules, 'boltz.data.msa.mmseqs2',
                        types.SimpleNamespace(run_mmseqs2=fake_msa))
    receipt = module.freeze_msas(output)
    assert calls == [['ACD']]
    assert receipt['X']['sequence_count'] == 2
    manifest = json.loads((output / 'manifest.json').read_text())
    assert all('--use_msa_server' not in c['argv'] for c in manifest['commands'])
    for row in manifest['inputs']:
        protein = json.loads((output / row['input']).read_text())['sequences'][0]['protein']
        assert protein['msa'] == ('empty' if row['msa'] == 'empty' else 'msas/X.a3m')
    module.freeze_msas(output)
    assert calls == [['ACD']]  # repeated freezing does not call network again
    (output / 'results').mkdir()
    with pytest.raises(ValueError, match='after inference'):
        module.freeze_msas(output)
