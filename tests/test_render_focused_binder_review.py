import importlib.util
from pathlib import Path

SPEC = importlib.util.spec_from_file_location(
    'focused_review', Path(__file__).parents[1] / 'scripts/render_focused_binder_review.py')
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def native(identifier, status='pass'):
    return {'target_id': 'FGF2', 'candidate_id': identifier, 'native_filter_status': status, 'sequence': 'ACD'}


def cases(identifier, values=(.6, .6, .6)):
    result = []
    for seed, score in zip((11, 22, 33), values):
        base = {'target_id': 'FGF2', 'candidate_id': identifier, 'seed': seed, 'status': 'analyzed',
                'binder_plddt_0_100': {'mean': 90},
                'binder_reference_pose': {'binder_own_fold_CA_RMSD_A': .5,
                                          'binder_CA_RMSD_after_target_alignment_A': 4,
                                          'target_CA_RMSD_A': 1}}
        result.append({**base, 'kind': 'fixed_binder_complex',
                       'binder_target_chainpair_iptm': {'A|Z': score, 'Z|A': score + .02},
                       'binder_target_interfaces': {'A': {'severe_heavy_atom_overlap': {
                           'status': 'assessed', 'atom_pairs_below_cutoff': 0}}},
                       'posthoc_target_reference_geometry': {'references': [{
                           'status': 'assessed', 'role': 'held_out', 'pdb_id': '1BFG', 'global_CA_RMSD_A': 1}]}})
        result.append({**base, 'kind': 'binder_alone'})
    return result


def test_conservative_seed_ranking_and_native_failure_exclusion():
    rows = module.rank_candidates([native('unstable'), native('consistent'), native('failed', 'fail')],
                                  cases('unstable', (.9, .2, .9)) + cases('consistent', (.5, .5, .5)) + cases('failed', (.95, .95, .95)))
    ranked = {r['candidate_id']: r for r in rows}
    assert ranked['consistent']['rank_within_target'] == 1
    assert ranked['unstable']['rank_within_target'] == 2
    assert ranked['failed']['rank_within_target'] is None
    assert 'Native filters did not pass' in ranked['failed']['exclusion_reasons']


def test_missing_failed_or_duplicate_results_never_ranked():
    complete = cases('candidate')
    for bad in (complete[:-1], complete + [complete[0]],
                [{**c, 'status': 'unavailable'} if i == 0 else c for i, c in enumerate(complete)]):
        result = module.rank_candidates([native('candidate')], bad)[0]
        assert result['rank_within_target'] is None
        assert result['ranking_status'] == 'excluded'
    incomplete = cases('candidate')
    incomplete[1]['binder_reference_pose'] = {}
    assert module.rank_candidates([native('candidate')], incomplete)[0]['rank_within_target'] is None


def test_html_escapes_candidate_sequence_and_flags(tmp_path):
    identifier = '<script>alert("bad")</script>'
    record = native(identifier)
    record['sequence'] = '<img src=x onerror=alert(1)>'
    ranked = module.rank_candidates([record], cases(identifier))
    ranked[0]['review_flags'].append('<svg onload=alert(2)>')
    report = {'native_candidates': 1, 'native_passed': 1, 'rankable': 1, 'excluded': 0, 'candidates': ranked}
    output = module.render_html(report, tmp_path, tmp_path)
    assert '<script>' not in output and '<img src=x' not in output and '<svg onload' not in output
    assert '&lt;script&gt;' in output and '&lt;img src=x' in output
    assert module.file_within(tmp_path / 'child', '../secret') is None
