#!/usr/bin/env python3
"""Render a cautious local review of focused designs and standard refolding.

Ranks within targets only, using a transparent lexicographic ordering. Native
filter pass plus complete repeated predictions permits comparison, not a claim
of binding. No molecular JavaScript, inference, or modification of inputs.
"""
from __future__ import annotations

import argparse
import html
import json
import math
import os
import statistics
from collections import Counter
from pathlib import Path
from urllib.parse import quote


def finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def file_within(root, relative):
    if not relative:
        return None
    root = Path(root).resolve()
    path = (root / relative).resolve()
    return path if path.is_relative_to(root) and path.is_file() else None


def collect_native(design_root, prepared):
    records = {}
    # Prefer active design bundle; a prepared copied bundle is an equivalent fallback.
    sources = list(Path(design_root).glob('*/live/bundle/bundle.json'))
    sources += list(Path(prepared).glob('provenance/*/native_bundle/bundle.json'))
    for path in sources:
        gene = path.parents[2].name if path.parent.name == 'bundle' else path.parents[1].name
        bundle = json.loads(path.read_text())
        evaluations = {e['subject_id']: e for e in bundle.get('evaluations', [])
                       if e.get('method') == 'BoltzGen native default filters'}
        for component in bundle.get('components', []):
            identifier = component['id']
            evaluation = evaluations.get(identifier, {})
            record = {'target_id': gene, 'candidate_id': identifier, 'sequence': component.get('sequence', ''),
                      'native_filter_status': evaluation.get('status', 'unavailable'),
                      'native_metrics': evaluation.get('metrics', {}),
                      'native_failed_filters': evaluation.get('details', {}).get('failed_filters', []),
                      'native_complex': str(p) if (p := file_within(path.parent, component.get('evidence', {}).get('complex_structure'))) else None,
                      'native_binder': str(p) if (p := file_within(path.parent, component.get('structure'))) else None,
                      'native_bundle': str(path.resolve())}
            records.setdefault((gene, identifier), record)
    return list(records.values())


def target_geometry(case):
    refs = case.get('posthoc_target_reference_geometry', {}).get('references', [])
    assessed = [r for r in refs if r.get('status') == 'assessed' and finite(r.get('global_CA_RMSD_A'))]
    heldout = [r for r in assessed if r.get('role') == 'held_out']
    if heldout or assessed:
        chosen = (heldout or assessed)[0]
        return chosen['global_CA_RMSD_A'], f"Experimental {chosen.get('pdb_id', '')} ({chosen.get('role', '')}); evaluation only"
    pose = case.get('binder_reference_pose', {})
    value = pose.get('target_anchor_CA_RMSD_A', pose.get('target_CA_RMSD_A'))
    return value, 'Original computational design target: agreement only'


def rank_candidates(native, cases, seeds=(11, 22, 33), selected=None):
    """Missing predictions/metrics and native failures cannot become ranked leads."""
    seeds = list(seeds)
    if not seeds or len(set(seeds)) != len(seeds):
        raise ValueError('Ranking requires distinct nonempty expected seeds')
    output = []
    for source in native:
        row = dict(source)
        matching = [c for c in cases if c.get('target_id') == row['target_id']
                    and c.get('candidate_id') == row['candidate_id']]
        complex_rows = [c for c in matching if c.get('kind') == 'fixed_binder_complex']
        alone_rows = [c for c in matching if c.get('kind') == 'binder_alone']
        reasons = []
        if row['native_filter_status'] != 'pass':
            reasons.append('Native filters did not pass')
        key = (row['target_id'], row['candidate_id'])
        if selected is not None and key not in selected:
            reasons.append('Not selected for standard repeated validation')
        for label, values in [('complex', complex_rows), ('binder-alone', alone_rows)]:
            observed = Counter(c.get('seed') for c in values)
            if observed != Counter(seeds) or any(c.get('status') != 'analyzed' for c in values):
                reasons.append(f'{label}: incomplete or duplicate seed results')
        min_interface, alone_confidence, complex_confidence, fold_rmsds, clashes, geometries = [], [], [], [], [], []
        geometry_labels = set()
        heldout_rmsds, pose_rmsds = [], []
        for case in complex_rows:
            if case.get('status') != 'analyzed':
                continue
            scores = case.get('binder_target_chainpair_iptm', {})
            # Focused protocol has one target and one binder, both directed scores required.
            if len(scores) != 2 or not all(finite(x) for x in scores.values()):
                reasons.append('Directional interface confidence missing/nonfinite')
            else:
                min_interface.append(min(scores.values()))
            confidence = case.get('binder_plddt_0_100', {}).get('mean')
            if finite(confidence):
                complex_confidence.append(confidence)
            else:
                reasons.append('Complex binder pLDDT missing/nonfinite')
            interfaces = case.get('binder_target_interfaces', {})
            overlaps = [i.get('severe_heavy_atom_overlap', {}) for i in interfaces.values()]
            if not overlaps or any(i.get('status') != 'assessed' or not finite(i.get('atom_pairs_below_cutoff')) for i in overlaps):
                reasons.append('Severe-clash assessment missing')
            else:
                clashes.append(sum(i['atom_pairs_below_cutoff'] for i in overlaps))
            for reference in case.get('posthoc_target_reference_geometry', {}).get('references', []):
                if reference.get('role') == 'held_out' and reference.get('status') == 'assessed' and finite(reference.get('global_CA_RMSD_A')):
                    heldout_rmsds.append(reference['global_CA_RMSD_A'])
            pose = case.get('binder_reference_pose', {}).get('binder_CA_RMSD_after_target_alignment_A')
            if finite(pose):
                pose_rmsds.append(pose)
            geometry, label = target_geometry(case)
            if finite(geometry):
                geometries.append(geometry)
                geometry_labels.add(label)
            else:
                reasons.append('Target geometry assessment missing')
        for case in alone_rows:
            if case.get('status') != 'analyzed':
                continue
            confidence = case.get('binder_plddt_0_100', {}).get('mean')
            fold = case.get('binder_reference_pose', {}).get('binder_own_fold_CA_RMSD_A')
            if finite(confidence):
                alone_confidence.append(confidence)
            else:
                reasons.append('Binder-alone pLDDT missing/nonfinite')
            if finite(fold):
                fold_rmsds.append(fold)
            else:
                reasons.append('Binder-alone fold agreement missing/nonfinite')
        flags = []
        if heldout_rmsds and max(heldout_rmsds) > 3:
            flags.append('Target-reference RMSD exceeds provisional 3Å review level; inspect target geometry. This is local triage, not a universal binding cutoff.')
        if not heldout_rmsds:
            flags.append('Experimental held-out target geometry unavailable; design-target agreement is not experimental accuracy.')
        if pose_rmsds:
            flags.append(f'Design-pose displacement across seeds: {min(pose_rmsds):.2f}–{max(pose_rmsds):.2f}Å after target alignment; agreement only, with no binding threshold assigned.')
        else:
            flags.append('Design-pose consistency unavailable.')
        if clashes and max(clashes) > 0:
            flags.append(f'Severe-overlap review: up to {max(clashes):g} binder–target atom pairs <2Å.')
        row.update({'ranking_status': 'excluded' if reasons else 'rankable_computational_candidate',
                    'exclusion_reasons': sorted(set(reasons)), 'rank_within_target': None,
                    'minimum_directional_interface_confidence': min(min_interface) if min_interface else None,
                    'minimum_binder_alone_plddt_0_100': min(alone_confidence) if alone_confidence else None,
                    'minimum_complex_binder_plddt_0_100': min(complex_confidence) if complex_confidence else None,
                    'maximum_binder_alone_design_fold_RMSD_A': max(fold_rmsds) if fold_rmsds else None,
                    'maximum_severe_atom_pairs_lt2A': max(clashes) if clashes else None,
                    'maximum_target_geometry_RMSD_A': max(geometries) if geometries else None,
                    'target_geometry_reference': sorted(geometry_labels), 'review_flags': flags,
                    'heldout_target_RMSD_A': {'n': len(heldout_rmsds), 'max': max(heldout_rmsds) if heldout_rmsds else None,
                                             'median': statistics.median(heldout_rmsds) if heldout_rmsds else None},
                    'design_pose_RMSD_A': {'n': len(pose_rmsds), 'max': max(pose_rmsds) if pose_rmsds else None,
                                          'median': statistics.median(pose_rmsds) if pose_rmsds else None},
                    'complex_results': complex_rows, 'binder_alone_results': alone_rows})
        output.append(row)
    for target in sorted({r['target_id'] for r in output}):
        ranked = [r for r in output if r['target_id'] == target and r['ranking_status'] != 'excluded']
        ranked.sort(key=lambda r: (r['maximum_severe_atom_pairs_lt2A'],
                                  -r['minimum_directional_interface_confidence'],
                                  -r['minimum_binder_alone_plddt_0_100'],
                                  r['maximum_binder_alone_design_fold_RMSD_A'],
                                  r['maximum_target_geometry_RMSD_A'], r['candidate_id']))
        for rank, row in enumerate(ranked, 1):
            row['rank_within_target'] = rank
    return sorted(output, key=lambda r: (r['target_id'], r['rank_within_target'] or math.inf, r['candidate_id']))


POLICY = ('Rank within each target only. Native filter pass and complete distinct seeds in BOTH complex and binder-alone contexts are required. '
          'Lexicographic ordering: fewer worst-seed severe heavy-atom pairs <2Å; higher minimum of both directional interface confidence scores across all seeds; '
          'higher minimum binder-alone pLDDT; lower worst binder-alone RMSD to the original design fold; lower worst target-geometry RMSD. '
          'This is a transparent prioritization rule, not a calibrated binding threshold. Geometry remains visible even when it does not break a tie.')


def render_html(report, output, prepared):
    output, prepared = Path(output), Path(prepared)
    esc = lambda value: html.escape(str(value), quote=True)
    def number(value):
        return f'{value:.3f}' if finite(value) else '—'
    def link(path, label):
        if not path or not Path(path).is_file():
            return 'unavailable'
        relative = os.path.relpath(Path(path).resolve(), output.resolve())
        return f'<a href="{esc(quote(relative, safe="/"))}">{esc(label)}</a>'
    body = ['<!doctype html><html><head><meta charset="utf-8"><title>Focused binder review</title>',
            '<style>body{font:16px system-ui;max-width:1450px;margin:32px auto;padding:0 24px;color:#182633;background:#fafbfd}table{border-collapse:collapse;width:100%;font-size:14px}th,td{border-bottom:1px solid #dce2ea;text-align:left;padding:10px;vertical-align:top}code{overflow-wrap:anywhere}a{color:#165c9c}.notice{background:#fff2d3;padding:18px;border-radius:10px}.muted{color:#5d6976}details{margin:16px 0}summary{cursor:pointer}h2{margin-top:32px}</style></head><body>',
            '<h1>Focused binder review</h1>',
            '<p class="notice"><strong>Computational candidates, not confirmed binders.</strong> Fresh standard Boltz predictions use target MSAs and empty binder MSAs, with no structure templates. Native filters and confident predictions do not establish affinity, selectivity or biological activity.</p>',
            f"<p>{report['native_candidates']} actual native candidates; {report['native_passed']} passed native filters. {report['rankable']} have complete repeated checks and can be compared; {report['excluded']} remain excluded from ranking.</p>",
            '<details><summary>How the ranking works</summary><p>' + esc(POLICY) + '</p><p>Three requested seeds are repeat predictions, not an equilibrium ensemble or a significance test. Counts and exclusions come from actual saved outputs. Rank is within each target; interface scores are not calibrated across target sizes.</p></details>',
            '<table><thead><tr><th>Target / candidate</th><th>Rank</th><th>Native filters</th><th>Worst directional interface confidence</th><th>Worst binder-alone pLDDT /100</th><th>Worst binder-alone design RMSD Å</th><th>Worst severe atom pairs &lt;2Å</th><th>Worst target geometry RMSD Å</th></tr></thead><tbody>']
    for row in report['candidates']:
        body.append('<tr><td>' + esc(row['target_id']) + '<br><code>' + esc(row['candidate_id']) + '</code></td><td>' +
                    esc(row['rank_within_target'] or 'Excluded') + '</td><td>' + esc(row['native_filter_status']) + '</td>' +
                    ''.join('<td>' + number(row[key]) + '</td>' for key in (
                        'minimum_directional_interface_confidence', 'minimum_binder_alone_plddt_0_100',
                        'maximum_binder_alone_design_fold_RMSD_A', 'maximum_severe_atom_pairs_lt2A',
                        'maximum_target_geometry_RMSD_A')) + '</tr>')
        body.append('<tr><td colspan="8" class="muted"><strong>Review notes:</strong> ' + esc(' '.join(row['review_flags'])) +
                    ' Held-out target RMSD median/max: ' + number(row['heldout_target_RMSD_A']['median']) + '/' + number(row['heldout_target_RMSD_A']['max']) + ' Å.</td></tr>')
    body.append('</tbody></table><h2>Structures, sequences and what failed</h2><p>All available standard prediction seeds are linked. For an initial visual inspection use the lowest-numbered seed; it is not selected for the best-looking outcome. Original designs are labeled separately.</p>')
    for row in report['candidates']:
        body.append('<details><summary>' + esc(row['target_id'] + ' · ' + row['candidate_id']) + '</summary>')
        if row['exclusion_reasons']:
            body.append('<p><strong>Excluded:</strong> ' + esc('; '.join(row['exclusion_reasons'])) + '</p>')
        if row.get('native_failed_filters'):
            body.append('<p>Native failures: ' + esc(', '.join(row['native_failed_filters'])) + '</p>')
        body.append('<p>Target geometry comparison: ' + esc('; '.join(row['target_geometry_reference']) or 'Unavailable') + '</p>')
        body.append('<p>' + link(row.get('native_complex'), 'Original computational design complex') + ' · ' +
                    link(row.get('native_binder'), 'Original binder structure') + '</p>')
        for label, key in [('Standard complex', 'complex_results'), ('Standard binder alone', 'binder_alone_results')]:
            for case in sorted(row[key], key=lambda r: r.get('seed', -1)):
                path = file_within(prepared, case.get('files', {}).get('structure'))
                body.append('<p>' + esc(f"{label}, seed {case.get('seed')}: {case.get('status')}") + ' · ' +
                            link(path, 'CIF coordinates') + '</p>')
                if case.get('error') or case.get('missing_files'):
                    body.append('<p class="muted">' + esc(case.get('error') or ', '.join(case['missing_files'])) + '</p>')
        body.append('<p>Binder sequence:</p><code>' + esc(row.get('sequence', '')) + '</code></details>')
    body += ['<p><a href="ranking.json">Full ranking, per-seed metrics and exclusions</a></p>',
             '<p class="muted">Native structural filtering, repeated model confidence and experimental validation are separate evidence levels. No candidate is declared experimentally validated here.</p></body></html>']
    return '\n'.join(body)


def build(analysis_path, prepared, design_root, output):
    analysis_path, prepared, output = Path(analysis_path), Path(prepared), Path(output)
    analysis = json.loads(analysis_path.read_text()) if analysis_path.is_file() else {'cases': []}
    manifest_path = prepared / 'manifest.json'
    manifest = json.loads(manifest_path.read_text()) if manifest_path.is_file() else {}
    selected = {(r['target_id'], r['candidate_id']) for r in manifest.get('selection_audit', []) if r.get('selected')}
    native = collect_native(design_root, prepared)
    candidates = rank_candidates(native, analysis.get('cases', []), manifest.get('seeds', [11, 22, 33]),
                                 selected if 'selection_audit' in manifest else None)
    report = {'policy': POLICY, 'analysis_available': analysis_path.is_file(),
              'native_candidates': len(native), 'native_passed': sum(r['native_filter_status'] == 'pass' for r in native),
              'rankable': sum(r['rank_within_target'] is not None for r in candidates),
              'excluded': sum(r['rank_within_target'] is None for r in candidates), 'candidates': candidates}
    output.mkdir(parents=True, exist_ok=True)
    (output / 'ranking.json').write_text(json.dumps(report, indent=2) + '\n')
    (output / 'review.html').write_text(render_html(report, output, prepared))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--analysis', type=Path, default=Path('outputs/folding_validation/focused_validation_report/analysis.json'))
    parser.add_argument('--prepared', type=Path, default=Path('outputs/folding_validation/download/focused_validation'))
    parser.add_argument('--design-root', type=Path, default=Path('outputs/folding_validation/focused_design'))
    parser.add_argument('--output', type=Path, default=Path('outputs/folding_validation/focused_review'))
    args = parser.parse_args()
    report = build(args.analysis, args.prepared, args.design_root, args.output)
    print(json.dumps({k: report[k] for k in ('native_candidates', 'native_passed', 'rankable', 'excluded')}))


if __name__ == '__main__':
    main()
