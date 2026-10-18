#!/usr/bin/env python3
"""Compare completed positive-control poses with withheld experimental 1WQJ.

No input modifications or inference. Scores are contact/pose agreement, not
biological binding validation. Reads both matched empty-binder-MSA and natural
paired-MSA conditions; reports missing seeds explicitly.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from analyze_candidate_validation import read_chains, reference_pose


def audit(root):
    root = Path(root)
    evaluation = json.loads((root / 'candidates/manifest.json').read_text())
    results = []
    for condition in ('candidate_controls', 'positive_msa_control'):
        manifest = json.loads((root / condition / 'manifest.json').read_text())
        row = next(r for r in manifest['inputs'] if r['kind'] == 'experimental_positive')
        expected = {r['chain']: r['sequence'] for r in row['target_chains']}
        expected[row['binder_chain']] = row['binder_sequence']
        for command in manifest['commands']:
            if command['condition'] != row['condition']:
                continue
            path = (root / 'download' / condition / command['expected_prediction_directory'] /
                    row['id'] / f"{row['id']}_model_0.cif")
            result = {'condition': condition, 'seed': command['seed'], 'prediction': str(path), 'status': 'missing'}
            if path.exists():
                try:
                    chains = read_chains(path, expected)
                    result.update(reference_pose(chains, row, evaluation))
                except (ValueError, RuntimeError, OSError) as exc:
                    result.update(status='unavailable', reason=str(exc))
            results.append(result)
    return {'reference': evaluation['positive_control']['experimental_complex_EVALUATION_ONLY'],
            'cases': results, 'limitations': [
                'Reference binder coordinates are evaluation-only and were never prediction inputs.',
                'Matched empty-binder-MSA condition and natural paired-MSA condition differ in binder MSA and target pairing information.',
                'Contact recall, precision and Jaccard use strict heavy-atom distance <5 angstroms at exact sequence positions.',
                'CA-only interface RMSD jointly fits reference-contact residues from both chains. It is not DockQ or a full-backbone iRMSD.',
                'Historical1WQJ may appear in model training; approximate contact recovery is not experimental validation of new binding.',
                'Three requested seeds per condition do not provide a significance test.']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path('outputs/folding_validation'))
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = audit(args.root)
    output = args.output or args.root / 'positive_pose_audit.json'
    output.write_text(json.dumps(result, indent=2) + '\n')
    for row in result['cases']:
        contact = row.get('reference_contact_recovery', {})
        print(row['condition'], row['seed'], row['status'],
              {k: contact.get(k) for k in ('contact_recall', 'contact_precision', 'contact_jaccard',
                                          'reference_interface_CA_jointfit_RMSD_A')})


if __name__ == '__main__':
    main()
