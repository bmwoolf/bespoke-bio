"""Public exporter must fail closed and never copy operational trees."""
import importlib.util
import json
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    'folding_export', Path(__file__).parents[1] / 'scripts/export_folding_validation.py')
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))


def complete(root):
    write(root / 'design_status.json', {'cohorts': [
        {'target': gene, 'status': 'completed'} for gene in MODULE.GENES]})
    write(root / 'focused_validation_status.json', {
        'state': 'fixed_validation_completed_pending_scientific_review'})
    for gene in MODULE.GENES:
        write(root / f'focused_design/{gene}/live/bundle/bundle.json', {'components': []})
    for (path, _, key) in MODULE.REPORTS.values():
        write(root / path, {key: [{'id': 'test', 'status': 'analyzed',
                                  'private': 'DO_NOT_EXPORT'}],
                           'status': 'complete', 'expected_predictions': 1,
                           'completed_predictions': 1, 'analyzed_predictions': 1,
                           'private': 'DO_NOT_EXPORT', 'summary': []})


def test_unfinished_refused_before_destination_creation(tmp_path):
    root, destination = tmp_path / 'campaign', tmp_path / 'public'
    root.mkdir()
    with pytest.raises(ValueError, match='unfinished'):
        MODULE.export(root, destination)
    assert not destination.exists()
    result = MODULE.export(root, destination, allow_partial=True)
    assert result['partial'] and result['incomplete_reasons']
    assert 'PARTIAL' in (destination / 'README.md').read_text()


def test_allowlist_excludes_operational_files_and_extra_fields(tmp_path):
    root, destination = tmp_path / 'campaign', tmp_path / 'public'
    complete(root)
    write(root / 'private/credentials.json', {'token': 'DO_NOT_EXPORT'})
    write(root / 'runtime/cache.json', {'token': 'DO_NOT_EXPORT'})
    (root / 'run.log').write_text('DO_NOT_EXPORT')
    result = MODULE.export(root, destination)
    assert not result['partial']
    assert {p.name for p in destination.iterdir()} == {'README.md', 'science.json'}
    assert 'DO_NOT_EXPORT' not in (destination / 'science.json').read_text()


def test_artifact_symlink_escape_rejected(tmp_path):
    root, destination = tmp_path / 'campaign', tmp_path / 'public'
    complete(root)
    outside = tmp_path / 'secret.cif'
    outside.write_text('private')
    linked = root / 'download/candidate_controls/escape.cif'
    linked.parent.mkdir(parents=True, exist_ok=True)
    linked.symlink_to(outside)
    write(root / 'candidate_report/analysis.json', {
        'expected_predictions': 1, 'analyzed_predictions': 1,
        'cases': [{'id': 'test', 'files': {'structure': 'escape.cif'}}]})
    with pytest.raises(ValueError, match='escapes'):
        MODULE.export(root, destination)
    assert not destination.exists()


def test_only_referenced_scientific_artifacts_are_copied(tmp_path):
    root, destination = tmp_path / 'campaign', tmp_path / 'public'
    complete(root)
    base = root / 'download/candidate_controls/results'
    base.mkdir(parents=True)
    (base / 'model.cif').write_text('data_test\n')
    (base / 'worker.log').write_text('DO_NOT_EXPORT')
    write(root / 'candidate_report/analysis.json', {
        'expected_predictions': 1, 'analyzed_predictions': 1,
        'cases': [{'id': 'test', 'status': 'analyzed',
                   'files': {'structure': 'results/model.cif', 'worker_log': 'results/worker.log'},
                   'binder_target_chainpair_iptm': {'A|Z': 0.7, 'api_token': 'DO_NOT_EXPORT'}}]})
    result = MODULE.export(root, destination)
    row = result['reports']['candidates']['rows'][0]
    assert row['binder_target_chainpair_iptm'] == {'A|Z': 0.7}
    assert (destination / row['structure']['path']).read_text() == 'data_test\n'
    assert len(row['structure']['sha256']) == 64
    assert not list(destination.rglob('*.log'))


def test_workspace_relative_report_paths(tmp_path, monkeypatch):
    monkeypatch.setattr(MODULE, 'WORKSPACE', tmp_path)
    root = tmp_path / 'outputs/folding_validation'
    root.mkdir(parents=True)
    source = root / 'download/model.cif'
    source.parent.mkdir()
    source.write_text('data_test\n')
    result = MODULE.report_source(root, '.', 'outputs/folding_validation/download/model.cif')
    assert result == source
    assert MODULE.report_source(root, 'download', 'model.cif') == source
    assert MODULE.report_source(root, '.', str(source)) == source


def test_posthoc_geometry_preserves_scores_without_local_reference_paths():
    row = MODULE.public_row({'posthoc_target_reference_geometry': {
        'status': 'evaluated', 'references': [{
            'role': 'held_out', 'pdb_id': '1BFG', 'status': 'assessed',
            'reference': '/private/path/structure.cif', 'global_CA_RMSD_A': 1.2,
            'target_input_coverage': 0.98}]},
        'posthoc_target_heldout_CA_RMSD_A': {'mean': 1.2}})
    ref = row['posthoc_target_reference_geometry']['references'][0]
    assert ref['role'] == 'held_out' and ref['pdb_id'] == '1BFG'
    assert ref['global_CA_RMSD_A'] == 1.2
    assert 'reference' not in ref
