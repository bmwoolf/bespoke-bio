import importlib.util
import json
import sys
import types
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parents[1] / 'scripts'))
SPEC = importlib.util.spec_from_file_location(
    'positive_msa', Path(__file__).parents[1] / 'scripts/prepare_positive_msa_control.py')
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_native_pairing_replaces_both_msas_without_binder_template(tmp_path, monkeypatch):
    output = tmp_path / 'prepared'
    path = output / 'inputs/old/positive.yaml'
    path.parent.mkdir(parents=True)
    schema = {'version': 1, 'sequences': [
        {'protein': {'id': 'A', 'sequence': 'ACD', 'msa': 'msas/oldtarget.a3m'}},
        {'protein': {'id': 'Z', 'sequence': 'EFG', 'msa': 'msas/oldbinder.a3m'}}],
        'templates': [{'cif': 'templates/target.cif', 'chain_id': ['A'], 'template_id': ['A'], 'force': False}]}
    path.write_text(json.dumps(schema))
    (output / 'manifest.json').write_text(json.dumps({
        'inputs': [{'id': 'positive', 'input': 'inputs/old/positive.yaml', 'condition': 'old', 'binder_chain': 'Z'}],
        'commands': [{'condition': 'old', 'argv': ['boltz', 'predict', 'inputs/old'],
                      'expected_prediction_directory': 'results/old'}],
        'frozen_msas': {'old': {}}, 'limitations': ['old']}))
    calls = []
    def compute_msa(**kwargs):
        calls.append(kwargs)
        for name, seq in kwargs['data'].items():
            (kwargs['msa_dir'] / f'{name}.csv').write_text('key,sequence\n0,' + seq + '\n')
    def parse_csv(path, **kwargs):
        return types.SimpleNamespace(sequences=np.array([(0, 3)], dtype=[('res_start', int), ('res_end', int)]))
    monkeypatch.setitem(sys.modules, 'boltz.main', types.SimpleNamespace(compute_msa=compute_msa))
    monkeypatch.setitem(sys.modules, 'boltz.data.parse.csv', types.SimpleNamespace(parse_csv=parse_csv))
    monkeypatch.setattr(module, 'validate_native', lambda *_: [])
    receipt = module.freeze_paired(output, tmp_path)
    assert len(calls) == 1 and calls[0]['data'] == {'natural_A': 'ACD', 'natural_Z': 'EFG'}
    assert calls[0]['msa_pairing_strategy'] == 'greedy'
    manifest = json.loads((output / 'manifest.json').read_text())
    new_schema = json.loads((output / manifest['inputs'][0]['input']).read_text())
    assert new_schema['templates'] == schema['templates']
    assert all(p['protein']['msa'].endswith('.csv') for p in new_schema['sequences'])
    assert len(receipt) == 2
    assert not path.exists()
