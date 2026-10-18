#!/usr/bin/env python3
"""Export a confined, allowlisted scientific handoff; never copy the campaign tree.

Default refuses unfinished focused design/validation. --allow-partial labels gaps.
No cloud operations, publication, or execution of scientific pipelines occurs here.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
from pathlib import Path

WORKSPACE = Path(__file__).resolve().parents[1]
GENES = ('FGF2', 'IL1RN')
REPORTS = {
    'controls': ('analysis/analysis.json', 'download/standard_controls', 'rows'),
    'candidates': ('candidate_report/analysis.json', 'download/candidate_controls', 'cases'),
    'positive_control': ('positive_msa_report/analysis.json', 'download/positive_msa_control', 'cases'),
    'variants': ('variant_report/comparison.json', '.', 'models'),
    'focused_validation': ('focused_validation_report/analysis.json',
                           'download/focused_validation', 'cases'),
}
LABELS = ('id', 'target_id', 'candidate_id', 'name', 'kind', 'condition', 'context',
          'source_native_filter', 'binder_chain', 'chain', 'seed', 'status', 'msa',
          'target_template', 'validation_scope', 'sequence_sha256')
METRICS = ('chains', 'reference_comparisons', 'confidence', 'pae_angstrom',
           'chain_plddt_0_100', 'binder_plddt_0_100', 'binder_target_chainpair_iptm',
           'target_target_chainpair_iptm', 'overall_complex_iptm_context_only', 'overall_ptm',
           'binder_target_interfaces', 'target_assembly_integrity', 'binder_reference_pose',
           'plddt_mean', 'plddt_per_position', 'local_windows', 'completed', 'expected',
           'expected_seeds', 'completed_seeds', 'target_plddt', 'reference_rmsd_angstrom',
           'binder_CA_RMSD_after_target_alignment_A', 'binder_own_fold_CA_RMSD_A',
           'target_global_CA_RMSD_A', 'binder_target_atom_pairs_lt2A',
           'posthoc_target_reference_geometry', 'posthoc_target_heldout_CA_RMSD_A',
           'target_anchor_design_CA_RMSD_A')
NATIVE_METRICS = ('design_to_target_iptm', 'esmfold2_ipsae_min', 'filter_rmsd',
                  'filter_rmsd_design', 'designfolding-filter_rmsd', 'design_ptm')
CAVEATS = [
    'Computational predictions only; no experimental binding, selectivity, safety, or benefit demonstrated.',
    'Native design filters and independent fixed-sequence validation are separate evidence.',
    ('When provided, posthoc target reference geometry is evaluated after prediction; those reference coordinates '
    'are not conditioning inputs to the predictor.'),
    'Template-assisted recovery is not independent fold recovery; historical structures may be training data.',
    'Three seeds measure limited prediction variability, not equilibrium ensembles or statistical significance.',
    ('SERPINA1 primary comparison uses shared experimental core positions, excluding variant sites; '
    'all-residue comparisons can be dominated by the unresolved amino-terminal tail.'),
    'SERPINA1 variants are single substitutions on the canonical background, not full natural haplotypes.',
    ('The experimental positive binder and shuffled decoy calibrate this workflow only; '
    'the decoy is not a validated nonbinder and may differ in folding.'),
    ('Oligomer context and competing biological partners require separate interpretation; '
    'a confidence score does not establish binding.'),
]


def load(path):
    return json.loads(path.read_text()) if path.is_file() else {}


def confined(root, path):
    path = Path(path)
    resolved = (path if path.is_absolute() else root / path).resolve()
    if not resolved.is_relative_to(root.resolve()):
        raise ValueError('Artifact escapes campaign root')
    return resolved


def report_source(root, data_root, source):
    """Resolve documented absolute, workspace-relative, or package-relative paths."""
    source = Path(source)
    if source.is_absolute():
        return confined(root, source)
    if root.is_relative_to(WORKSPACE):
        prefix = root.relative_to(WORKSPACE)
        if source.is_relative_to(prefix):
            return confined(root, WORKSPACE / source)
    return confined(root, root / data_root / source)


def identifier(value):
    if not isinstance(value, str) or not re.fullmatch(r'[A-Za-z0-9_.-]+', value):
        raise ValueError('Invalid scientific identifier')
    return value


def numbers(value):
    """Metric trees contain numbers only; no paths or arbitrary string payloads."""
    if value is None or isinstance(value, (bool, int, float)):
        return value
    if isinstance(value, list):
        return [numbers(v) for v in value if not isinstance(v, str)]
    if isinstance(value, dict):
        return {k: (v if isinstance(v, str) else numbers(v)) for k, v in value.items()
                if re.fullmatch(r'[A-Za-z0-9_. |:()><=+/-]+', k)
                and not re.search(r'path|secret|token|password|credential|command|log|host', k, re.IGNORECASE)
                and (not isinstance(v, str) or (k in {'a', 'b', 'name', 'pdb_id', 'role', 'chain', 'status'}
                     and re.fullmatch(r'[A-Za-z0-9_. -]+', v)))}
    return None


def public_row(row):
    result = {k: row[k] for k in LABELS if k in row}
    result.update({k: numbers(row[k]) for k in METRICS if k in row})
    return result


def export(root, destination, allow_partial=False):
    root, destination = Path(root).resolve(), Path(destination).resolve()
    if destination == root or root.is_relative_to(destination) or destination.is_relative_to(root):
        raise ValueError('Destination must be outside campaign root')
    if destination.exists() and any(destination.iterdir()):
        raise FileExistsError('Destination must be empty')
    design = load(root / 'design_status.json')
    validation = load(root / 'focused_validation_status.json')
    gaps = []
    for gene in GENES:
        rows = [r for r in design.get('cohorts', []) if r.get('target') == gene]
        if len(rows) != 1 or rows[0].get('status') != 'completed':
            gaps.append(f'{gene}: focused design unfinished')
        if not (root / f'focused_design/{gene}/live/bundle/bundle.json').is_file():
            gaps.append(f'{gene}: native bundle unavailable')
    if validation.get('state') != 'fixed_validation_completed_pending_scientific_review':
        gaps.append('Focused standard validation unfinished')
    reports = {name: load(root / spec[0]) for name, spec in REPORTS.items()}
    for name, data in reports.items():
        rows_key = REPORTS[name][2]
        if not data or not data.get(rows_key):
            gaps.append(f'{name}: report missing or empty')
        elif name == 'variants':
            if data.get('status') != 'complete':
                gaps.append('variants: incomplete')
        else:
            expected = data.get('expected_predictions', 0)
            done = data.get('analyzed_predictions', data.get('completed_predictions', 0))
            if not expected or done != expected:
                gaps.append(f'{name}: missing predictions')
    if gaps and not allow_partial:
        raise ValueError('Campaign unfinished: ' + '; '.join(gaps))
    # Build in memory first; copy only explicitly referenced structure/confidence artifacts.
    payload = {'schema_version': 1, 'partial': bool(gaps), 'incomplete_reasons': gaps,
               'limitations': CAVEATS, 'reports': {}, 'native_candidates': []}
    copies = []

    def artifact(source, relative, extensions):
        source = confined(root, source)
        if source.suffix not in extensions:
            raise ValueError('Unsupported scientific artifact extension')
        if not source.is_file():
            if not allow_partial:
                raise ValueError('Referenced scientific artifact missing')
            payload['partial'] = True
            payload['incomplete_reasons'].append(f'Missing artifact: {relative}')
            return None
        copies.append((source, relative))
        return {'path': relative, 'sha256': hashlib.sha256(source.read_bytes()).hexdigest()}

    for name, (report_path, data_root, key) in REPORTS.items():
        data = reports[name]
        out = {'rows': [], 'summary': [public_row(r) for r in data.get('summary', [])]}
        for index, row in enumerate(data.get(key, [])):
            item = public_row(row)
            files = row.get('files', row.get('artifacts', {}))
            source = row.get('path') if name == 'variants' else files.get('structure', files.get('cif'))
            if source:
                item['structure'] = artifact(report_source(root, data_root, source),
                    f'structures/{name}/{index:04d}.cif', {'.cif'})
            # NPZ arrays are scientific confidence evidence, not runtime caches.
            for field in ('pae', 'plddt'):
                if field in files:
                    item[field + '_artifact'] = artifact(report_source(root, data_root, files[field]),
                        f'evidence/{name}/{index:04d}_{field}.npz', {'.npz'})
            out['rows'].append(item)
        if 'matched_positive_decoy_deltas' in data:
            out['matched_positive_decoy_deltas'] = numbers(data['matched_positive_decoy_deltas'])
        if name == 'variants':
            for field in ('primary_fit_positions', 'primary_fit_count', 'rmsd_summaries',
                          'secondary_all_residue_fit', 'within', 'between_reference_variant',
                          'matched_seed'):
                if field in data:
                    out[field] = numbers(data[field])
        payload['reports'][name] = out
    for gene in GENES:
        base = root / f'focused_design/{gene}/live/bundle'
        bundle = load(base / 'bundle.json')
        for component in bundle.get('components', []):
            if component.get('role') != 'computational_binder_candidate':
                continue
            cid = identifier(component['id'])
            sequence = component['sequence']
            if not re.fullmatch('[ACDEFGHIKLMNPQRSTVWY]+', sequence):
                raise ValueError('Invalid protein sequence')
            evaluations = [e for e in bundle.get('evaluations', [])
                           if e.get('subject_id') == cid and 'filter' in e.get('method', '').lower()]
            evaluation = evaluations[0] if evaluations else {}
            if component.get('sequence_sha256') not in (None, hashlib.sha256(sequence.encode()).hexdigest()):
                raise ValueError('Native sequence hash mismatch')
            item = {'target': gene, 'id': cid, 'sequence': sequence,
                    'sequence_sha256': hashlib.sha256(sequence.encode()).hexdigest(),
                    'native_filter': evaluation.get('status', 'unknown'),
                    'metrics': {k: numbers(evaluation.get('metrics', {}).get(k)) for k in NATIVE_METRICS}}
            for label, source in [('binder', component.get('structure')),
                                  ('complex', component.get('evidence', {}).get('complex_structure'))]:
                if source:
                    item[label] = artifact(base / source,
                        f'native/{gene}/{cid}_{label}.cif', {'.cif'})
            payload['native_candidates'].append(item)
    plot = root / 'analysis/control_effects.png'
    if plot.exists():
        artifact(plot, 'control_effects.png', {'.png'})
    # Reviewed prose and fixed scientific figures only; no directory recursion.
    for relative in ('FINAL_FINDINGS.md', 'summary_figures/README.md'):
        if (root / relative).is_file():
            artifact(root / relative, relative, {'.md'})
    for stem in ('01_target_controls', '02_native_filtering',
                 '03_focused_candidate_metrics', '04_experimental_positive_control'):
        for extension in ('.png', '.pdf'):
            relative = f'summary_figures/{stem}{extension}'
            if (root / relative).is_file():
                artifact(root / relative, relative, {extension})
    destination.mkdir(parents=True, exist_ok=True)
    for source, relative in copies:
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
    (destination / 'science.json').write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n')
    (destination / 'README.md').write_text(
        '# Folding validation scientific handoff\n\n'
        + ('**PARTIAL — unfinished evidence.**\n\n' if payload['partial'] else
           'Completed computational evidence; scientific interpretation remains required.\n\n')
        + 'See science.json for scores and structure links. Native filters are preserved separately '
        'from standard fixed-sequence validation. All native bundle candidates are retained, '
        'including failures; counts are not restricted to successful examples.\n\n'
        + '\n'.join('- ' + c for c in CAVEATS) + '\n')
    return payload


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path, nargs='?', default=Path('outputs/folding_validation'))
    parser.add_argument('--destination', type=Path, default=Path('examples/folding_validation'))
    parser.add_argument('--allow-partial', action='store_true')
    args = parser.parse_args()
    export(args.root, args.destination, args.allow_partial)
