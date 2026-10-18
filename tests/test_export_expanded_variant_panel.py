import importlib.util
import json
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    'export_panel', Path(__file__).parents[1] / 'scripts/export_expanded_variant_panel.py')
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def fixture_campaign(root, status='completed'):
    root.mkdir()
    bundle = root / 'generation/X_reference/live/bundle'
    bundle.mkdir(parents=True)
    (root / 'generation_manifest.json').write_text(json.dumps({'cohorts': [
        {'target_id': 'X', 'state': 'reference', 'sequence': 'ACD', 'status': status,
         'directory': 'generation/X_reference', 'num_designs': 100,
         'private_api_key': 'SECRET_SENTINEL'}]}))
    (bundle / 'bundle.json').write_text(json.dumps({
        'components': [{'id': 'binder_1', 'sequence': 'EFG',
                        'role': 'computational_binder_candidate',
                        'evidence': {'private_path': '/Users/private/SECRET_SENTINEL'}}],
        'evaluations': [{'subject_id': 'binder_1', 'method': 'native filters',
                         'status': 'fail', 'metrics': {'design_to_target_iptm': .2}}]}))
    (root / 'connection.json').write_text('{"token":"SECRET_SENTINEL"}')


def test_export_counts_actual_candidates_and_excludes_operational_data(tmp_path):
    root, destination = tmp_path / 'source', tmp_path / 'public'
    fixture_campaign(root)
    summary = module.export(root, destination)
    assert summary['actual_candidates'] == 1
    assert summary['native_fail'] == 1
    assert summary['native_pass'] == 0
    assert not (destination / 'connection.json').exists()
    combined = ''.join(p.read_text() for p in destination.iterdir())
    assert 'SECRET_SENTINEL' not in combined
    assert '/Users/' not in combined
    assert str(root) not in combined
    assert 'not affinity measurements' in combined
    assert 'native_fail' in (destination / 'binders.fasta').read_text()


def test_unfinished_export_requires_explicit_partial(tmp_path):
    root, destination = tmp_path / 'source', tmp_path / 'public'
    fixture_campaign(root, 'running')
    with pytest.raises(ValueError, match='unfinished'):
        module.export(root, destination)
    assert not destination.exists()
    assert module.export(root, destination, allow_partial=True)['partial']


def test_crossfold_receipt_without_actual_files_is_incomplete(tmp_path):
    root = tmp_path / 'source'
    fixture_campaign(root)
    (root / 'crossfold_tasks.json').write_text(json.dumps({'tasks': [{'id': 't1'}]}))
    (root / 'crossfold').mkdir()
    (root / 'crossfold/receipt.json').write_text(json.dumps({
        'seeds': [11], 'results': [{'id': 't1', 'seed': 11}]}))
    with pytest.raises(ValueError, match='unfinished'):
        module.export(root, tmp_path / 'public')


def test_assembly_labels_match_fasta_csv_and_export_sanitized_audits(tmp_path):
    import csv

    root, destination = tmp_path / 'source', tmp_path / 'public'
    fixture_campaign(root)
    manifest_path = root / 'generation_manifest.json'
    manifest = json.loads(manifest_path.read_text())
    manifest['cohorts'][0]['target_id'] = 'PF4'
    manifest_path.write_text(json.dumps(manifest))
    (root / 'pf4_assembly_audit.json').write_text(json.dumps({
        'candidate_id': 'binder_1', 'placements': []}))
    (root / 'assembly_audits.json').write_text(json.dumps({
        'method': 'Example geometry', 'private_api_key': 'SECRET_SENTINEL',
        'targets': [{'target': 'PF4', 'status': 'assessed', 'pdb': '1F9Q',
                     'candidates': [{'candidate_id': 'binder_1',
                                     'assembly_assessment': 'substantial_overlap_all_placements',
                                     'complex_path': '/Users/private/SECRET_SENTINEL'}]}]}))
    (root / 'report').mkdir()
    (root / 'report/native_passes_test.png').write_bytes(b'fixture-image')
    module.export(root, destination)
    with (destination / 'candidates.csv').open() as handle:
        row = next(csv.DictReader(handle))
    assert row['assembly_assessment'] == 'intact_tetramer_incompatible_rigid_pose'
    assert '|assembly_' + row['assembly_assessment'] in (destination / 'binders.fasta').read_text()
    assert row['native_filter'] == 'fail'  # assembly annotation never overwrites the native result
    assert (destination / 'native_passes_test.png').is_file()
    assert (destination / 'assembly_audits.md').is_file()
    exported = (destination / 'assembly_audits.json').read_text()
    assert 'SECRET_SENTINEL' not in exported and '/Users/' not in exported
    assert 'complex_path' not in exported
