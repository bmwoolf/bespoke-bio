#!/usr/bin/env python3
"""Export an allowlisted public handoff; never copy operational files or raw logs.

Run AFTER completion:
  python scripts/export_expanded_variant_panel.py outputs/variant_expansion
Defaults to examples/expanded_variant_panel. Refuses a nonempty destination and an
unfinished queue/comparison plan. --allow-partial explicitly labels partial exports.
Counts come from actual native bundle components, not planned design counts or stale reports.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
from datetime import UTC, datetime
from pathlib import Path

METRICS = ('design_to_target_iptm', 'esmfold2_ipsae_min', 'filter_rmsd',
           'filter_rmsd_design', 'designfolding-filter_rmsd', 'design_ptm')
PLOTS = ('fixed_binder_confidence.png', 'target_structure_variation.png')
TERMINAL = {'completed', 'failed', 'cancelled', 'interrupted', 'skipped', 'skipped_time_budget'}


def load(path, default=None):
    return json.loads(path.read_text()) if path.is_file() else default


def safe_id(value):
    if not isinstance(value, str) or not re.fullmatch(r'[A-Za-z0-9_.-]+', value):
        raise ValueError('Public identifiers must contain only letters, digits, dot, hyphen, underscore')
    return value


def protein_sequence(value):
    if not isinstance(value, str) or not re.fullmatch('[ACDEFGHIKLMNPQRSTVWY]+', value):
        raise ValueError('Expected canonical amino-acid sequence')
    return value


def child(root, relative):
    path = (root / relative).resolve()
    if not path.is_relative_to(root):
        raise ValueError('Input path leaves campaign root')
    return path


def public_assembly_audit(value):
    """Allowlist scientific geometry fields; omit file paths and operational metadata."""
    allowed = {'method', 'classification', 'limitations', 'targets', 'target', 'status', 'pdb',
               'assembly_id', 'assembly_chains', 'source', 'native_sha256', 'candidates',
               'candidate_id', 'native_filter_status', 'assembly_assessment', 'placements',
               'anchor_chain', 'other_chain', 'matched_CA', 'alignment_CA_RMSD_A',
               'minimum_heavy_atom_distance_A', 'atom_pairs_lt2A', 'atom_pairs_lt2_5A',
               'atom_pairs_lt4A', 'binder_residues_lt2A', 'backbone_backbone_pairs_lt2A',
               'complex_sha256'}
    if isinstance(value, dict):
        return {k: public_assembly_audit(v) for k, v in value.items() if k in allowed}
    if isinstance(value, list):
        return [public_assembly_audit(v) for v in value]
    if isinstance(value, str):
        return re.sub(r'(?<![:\w])/(?:Users|home|workspace|opt|private|tmp|root|mnt)/[^\s)]+',
                      '[local path omitted]', value)
    return value


def export(root, destination, allow_partial=False):
    root, destination = Path(root).resolve(), Path(destination).resolve()
    generation = load(root / 'generation_manifest.json')
    if not generation or not generation.get('cohorts'):
        raise ValueError('No generation cohorts')
    unfinished = [f"{c['target_id']}/{c['state']}" for c in generation['cohorts']
                  if c.get('status') not in TERMINAL]
    tasks = load(root / 'crossfold_tasks.json', {}).get('tasks', [])
    receipt = load(root / 'crossfold/receipt.json', {})
    seeds = receipt.get('seeds', [11, 22, 33])
    expected = {(t['id'], seed) for t in tasks for seed in seeds}
    observed = {(r['id'], r['seed']) for r in receipt.get('results', [])
                if (root / 'crossfold' / f"seed_{r['seed']}" / 'refold_cif' /
                    f"{safe_id(r['id'])}.cif").is_file()
                and (root / 'crossfold' / f"seed_{r['seed']}" / 'fold_out_npz' /
                     f"{safe_id(r['id'])}.npz").is_file()}
    missing_comparisons = sorted(expected - observed)
    if (unfinished or missing_comparisons) and not allow_partial:
        raise ValueError('Campaign unfinished; wait or explicitly use --allow-partial')
    if destination.exists() and any(destination.iterdir()):
        raise FileExistsError('Destination is nonempty; choose a new directory')
    inputs, candidates, cohort_counts = [], [], []
    for cohort in generation['cohorts']:
        gene, state = safe_id(cohort['target_id']), safe_id(cohort['state'])
        sequence = protein_sequence(cohort['sequence'])
        directory = child(root, cohort['directory'])
        validation = load(child(root, cohort['input_validation']), {}) if cohort.get('input_validation') else {}
        construct = validation.get('requested_construct', {})
        mapping = validation.get('sequence_mapping', [])
        requested_span = [construct.get('start'), construct.get('end')]
        modeled_span = ([mapping[0]['precursor_position'], mapping[-1]['precursor_position']]
                        if mapping else requested_span)
        inputs.append({'target_id': gene, 'state': state, 'sequence': sequence,
                       'sequence_sha256': hashlib.sha256(sequence.encode()).hexdigest(),
                       'input_mode': cohort.get('input_mode', validation.get('input_mode', 'unknown')),
                       'modeled_length': len(sequence), 'requested_precursor_span': requested_span,
                       'modeled_precursor_span': modeled_span,
                       'pdb_id': validation.get('pdb_id'),
                       'selected_model': validation.get('selected_model'),
                       'source_model_count': validation.get('source_model_count'),
                       'hotspots_specified': False})
        bundle_path = directory / 'live/bundle/bundle.json'
        bundle = load(bundle_path, {})
        rows = []
        for component in bundle.get('components', []):
            if component.get('role') != 'computational_binder_candidate':
                continue
            evaluations = [e for e in bundle.get('evaluations', [])
                           if e.get('subject_id') == component['id']]
            filtering = next((e for e in evaluations if 'filter' in e.get('method', '').lower()), {})
            metrics = filtering.get('metrics', {})
            status = filtering.get('status', 'unknown')
            if status not in {'pass', 'fail'}:
                status = 'unknown'
            row = {'target_id': gene, 'state': state, 'binder_id': safe_id(component['id']),
                   'sequence': protein_sequence(component['sequence']), 'native_filter': status,
                   **{m: metrics.get(m) for m in METRICS}}
            rows.append(row)
        candidates.extend(rows)
        cohort_counts.append({'target_id': gene, 'state': state,
                              'generation_status': cohort.get('status', 'unknown'),
                              'bundle_present': bundle_path.is_file(), 'candidates': len(rows),
                              'native_pass': sum(r['native_filter'] == 'pass' for r in rows),
                              'native_fail': sum(r['native_filter'] == 'fail' for r in rows),
                              'native_unknown': sum(r['native_filter'] == 'unknown' for r in rows),
                              'unique_sequences': len({r['sequence'] for r in rows})})
    missing_bundles = [c for c in cohort_counts if c['generation_status'] == 'completed'
                       and not c['bundle_present']]
    if missing_bundles and not allow_partial:
        raise ValueError('Completed cohorts lack downloaded native bundles')
    partial = bool(unfinished or missing_comparisons or missing_bundles)
    summary = {'export_schema_version': 1, 'exported_utc': datetime.now(UTC).isoformat(),
               'partial': partial, 'cohorts': cohort_counts, 'actual_candidates': len(candidates),
               'native_pass': sum(c['native_pass'] for c in cohort_counts),
               'native_fail': sum(c['native_fail'] for c in cohort_counts),
               'native_unknown': sum(c['native_unknown'] for c in cohort_counts),
               'planned_crossfold_predictions': len(expected),
               'reported_crossfold_predictions': len(observed & expected),
               'missing_crossfold_predictions': len(missing_comparisons),
               'binding_validation': 'No experimental binding, affinity or variant selectivity established'}
    audit = load(root / 'pf4_assembly_audit.json', {})
    if audit:
        summary['PF4_assembly'] = {'binder_id': safe_id(audit['candidate_id']),
            'native_filter_status_unchanged': True,
            'assessment': 'Rigid monomer pose incompatible with intact experimental tetramer',
            'atom_pairs_below_2_angstrom_by_placement': [
                p['total_other_subunit_atom_pairs_lt2A'] for p in audit.get('placements', [])],
            'clashing_binder_residues_by_placement': [
                p['unique_binder_residues_clashing_lt2A'] for p in audit.get('placements', [])]}
    assemblies = public_assembly_audit(load(root / 'assembly_audits.json', {}))
    annotations = {(t['target'], c['candidate_id']): c['assembly_assessment']
                   for t in assemblies.get('targets', []) for c in t.get('candidates', [])}
    if audit:
        annotations[('PF4', audit['candidate_id'])] = 'intact_tetramer_incompatible_rigid_pose'
    for candidate in candidates:
        candidate['assembly_assessment'] = annotations.get(
            (candidate['target_id'], candidate['binder_id']), 'not_assessed')
    destination.mkdir(parents=True, exist_ok=True)
    for filename, data in [('inputs.json', inputs), ('summary.json', summary)]:
        (destination / filename).write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')
    (destination / 'targets.fasta').write_text(''.join(
        f">{r['target_id']}|{r['state']}\n{r['sequence']}\n" for r in inputs))
    (destination / 'binders.fasta').write_text(''.join(
        f">{r['target_id']}|{r['state']}|{r['binder_id']}|native_{r['native_filter']}"
        f"|assembly_{r['assembly_assessment']}\n{r['sequence']}\n"
        for r in candidates))
    fields = ['target_id', 'state', 'binder_id', 'sequence', 'native_filter', 'assembly_assessment', *METRICS]
    with (destination / 'candidates.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(candidates)
    task_by_id = {task['id']: task for task in tasks}
    comparison_rows = []
    for result in receipt.get('results', []):
        if (result['id'], result['seed']) not in observed & expected:
            continue
        task = task_by_id[result['id']]
        comparison_rows.append({
            'task_id': safe_id(result['id']), 'target_id': safe_id(task['target_id']),
            'designed_on': task.get('designed_on'), 'scored_on': task.get('scored_on'),
            'binder_id': task.get('binder_id'), 'seed': result['seed'],
            'metrics': {key: value for key, value in result.get('metrics', {}).items()
                        if key in {'iptm', 'ptm', 'protein_iptm', 'design_ptm', 'target_ptm',
                                   'design_to_target_iptm', 'min_design_to_target_pae'}},
        })
    (destination / 'comparison_scores.json').write_text(
        json.dumps(comparison_rows, indent=2, allow_nan=False) + '\n')
    plots = []
    gallery_names = [p.name for p in sorted((root / 'report').glob('native_passes_*.png'))
                     if re.fullmatch(r'native_passes_[A-Za-z0-9_.-]+\.png', p.name)]
    for filename in [*PLOTS, *gallery_names]:
        if (root / 'report' / filename).is_file():
            shutil.copyfile(root / 'report' / filename, destination / filename)
            plots.append(filename)
    if audit and (root / 'pf4_assembly_comparison.png').is_file():
        filename = 'pf4_assembly_comparison.png'
        shutil.copyfile(root / filename, destination / filename)
        plots.append(filename)
    if assemblies:
        (destination / 'assembly_audits.json').write_text(json.dumps(assemblies, indent=2) + '\n')
        notes = ['# Assembly geometry audits', '',
                 'Native filter status is unchanged. No rigid overlap does not establish binding.', '']
        for target in assemblies.get('targets', []):
            notes.append(f"## {target['target']} / {target.get('pdb', '')}: {target['status']}")
            notes.extend(f"- {c['candidate_id']}: {c['assembly_assessment']}"
                         for c in target.get('candidates', []))
            notes.append('')
        notes.extend([assemblies.get('method', ''), assemblies.get('classification', ''), ''])
        notes.extend('- ' + value for value in assemblies.get('limitations', []))
        (destination / 'assembly_audits.md').write_text('\n'.join(notes) + '\n')
    text = ['# Expanded protein binder research panel', '',
            ('**Partial export.**' if partial else '**Completed export of available outputs.**'), '',
            (f"{len(inputs)} target states; **{len(candidates)} actual generated candidates**; "
            f"{summary['native_pass']} native filter passes, {summary['native_fail']} failures, "
            f"{summary['native_unknown']} unknown. Counts use downloaded native bundles."), '',
            '| Target | State | Generation status | Candidates | Pass | Fail | Unknown |',
            '| --- | --- | --- | ---: | ---: | ---: | ---: |']
    text.extend(f"| {c['target_id']} | {c['state']} | {c['generation_status']} | {c['candidates']} | {c['native_pass']} | "
                f"{c['native_fail']} | {c['native_unknown']} |" for c in cohort_counts)
    text.extend(['', '## What these results mean', '',
        '- Native structural filters are computational gates, not proof of binding or therapeutic activity.',
        '- Interface ipTM and ESMFold2 ipSAE are confidence metrics, not affinity measurements.',
        ('- Fixed-binder reference/variant comparisons are exploratory; three model samples do not '
        'establish significance, selectivity, or an equilibrium ensemble.'),
        ('- Generation seeds were uncontrolled. Comparison Lightning/Torch seeds were 11/22/33; '
        'native feature randomness remains nondeterministic.'),
        ('- Experimental inputs use the exact modeled spans in inputs.json, including terminal truncations. '
        'Oligomers, modifications and cofactors may be omitted.'),
        '- SERPINA1 E288V/E366K are controlled substitutions on canonical background, not full phased S/Z haplotypes.',
        '', (f"Comparison progress: {summary['reported_crossfold_predictions']}/"
        f"{summary['planned_crossfold_predictions']} planned predictions reported.")])
    if audit:
        text.extend(['', ('**PF4 assembly caveat:** the monomer lead retains its native filter pass, '
                     'but its rigid pose clashes severely with the intact experimental tetramer. '
                     'It is not an assembly-compatible lead. See summary.json and the assembly image.')])
    text.extend(['', '## Files', '',
                 '- [Input sequences and construct metadata](inputs.json)',
                 '- [Target FASTA](targets.fasta)', '- [Binder FASTA](binders.fasta)',
                 '- [Candidate scores](candidates.csv)', '- [Repeated comparison scores](comparison_scores.json)',
                 '- [Machine-readable summary](summary.json)'])
    if assemblies:
        text.extend(['- [Assembly audit summary](assembly_audits.md)',
                     '- [Assembly geometry data](assembly_audits.json)'])
    for filename in plots:
        text.extend(['', f'![{filename.removesuffix(".png").replace("_", " ")}]({filename})'])
    (destination / 'README.md').write_text('\n'.join(text) + '\n')
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--output', type=Path, default=Path('examples/expanded_variant_panel'))
    parser.add_argument('--allow-partial', action='store_true')
    args = parser.parse_args()
    summary = export(args.root, args.output, args.allow_partial)
    print(json.dumps({'candidates': summary['actual_candidates'], 'partial': summary['partial']}))


if __name__ == '__main__':
    main()
