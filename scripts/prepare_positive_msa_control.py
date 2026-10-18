#!/usr/bin/env python3
"""Prepare a separate natural positive-control condition with unpaired binder MSA.

Copies the fixed positive from candidate_controls; changes only its binder MSA
and case/output identity. Never modifies running inputs or supplies binder pose.
Uses the official Boltz MSA client on CPU/network, then native parser validation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

from prepare_candidate_validation import save, validate_native
from prepare_standard_validation import sha256


def prepare(source, output, mol_dir):
    from boltz.data.msa.mmseqs2 import run_mmseqs2
    from boltz.data.parse.a3m import parse_a3m

    source, output = Path(source), Path(output)
    if output.exists() and any(output.iterdir()):
        raise ValueError('Output must be new or empty')
    original = json.loads((source / 'manifest.json').read_text())
    positives = [r for r in original['inputs'] if r['kind'] == 'experimental_positive']
    if len(positives) != 1 or not original.get('native_parser_validated'):
        raise ValueError('One already native-validated positive control required')
    row = positives[0].copy()
    original_id = row['id']
    schema = json.loads((source / row['input']).read_text())
    proteins = {p['protein']['id']: p['protein'] for p in schema['sequences']}
    if proteins[row['binder_chain']]['msa'] != 'empty':
        raise ValueError('Original matched positive must have empty binder MSA')
    output.mkdir(parents=True)
    (output / 'msas').mkdir()
    (output / 'templates').mkdir()
    receipt = {}
    for protein in proteins.values():
        if protein['id'] == row['binder_chain']:
            continue
        old_path = Path(protein['msa'])
        if old_path.is_absolute() or '..' in old_path.parts:
            raise ValueError('Expected portable frozen target MSA')
        shutil.copyfile(source / old_path, output / old_path)
        digest = hashlib.sha256(protein['sequence'].encode()).hexdigest()
        receipt[digest] = {**original['frozen_msas'][digest], 'role': 'target'}
    for template in schema['templates']:
        if row['binder_chain'] in template['chain_id']:
            raise ValueError('Binder template forbidden')
        path = Path(template['cif'])
        shutil.copyfile(source / path, output / path)
    binder = proteins[row['binder_chain']]
    sequence = binder['sequence']
    digest = hashlib.sha256(sequence.encode()).hexdigest()
    path = Path('msas') / (digest + '.a3m')
    alignments = run_mmseqs2([sequence], str(output / 'msas' / (digest + '_raw')),
                            use_env=True, use_pairing=False, host_url='https://api.colabfold.com')
    if len(alignments) != 1:
        raise ValueError('Unexpected MSA response')
    (output / path).write_text(alignments[0])
    lines = alignments[0].splitlines()
    query = ''.join(lines[1:next((i for i in range(1, len(lines)) if lines[i].startswith('>')), len(lines))])
    if query != sequence:
        raise ValueError('Binder MSA query differs from fixed binder')
    parsed = parse_a3m(output / path, taxonomy=None, max_seqs=16384)
    if any((parsed.sequences['res_end'] - parsed.sequences['res_start']) != len(sequence)):
        raise ValueError('Native alignment lengths differ')
    binder['msa'] = str(path)
    receipt[digest] = {'path': str(path), 'sha256': sha256(output / path), 'role': 'natural_binder',
                       'raw_sequence_count': sum(x.startswith('>') for x in lines),
                       'native_deduplicated_count': len(parsed.sequences), 'pairing': 'unpaired'}
    condition = 'natural_protein_msa_target_template'
    row.update({'id': original_id + '_with_binder_msa', 'condition': condition,
                'context': 'positive_natural_protein_binder_msa',
                'comparison_to': original_id,
                'limitation': 'Natural-protein control with target and binder unpaired MSAs. Separate from matched de novo empty-binder-MSA condition. No binder coordinates or disulfide constraints supplied.'})
    row['input'] = f"inputs/{condition}/{row['id']}.yaml"
    (output / row['input']).parent.mkdir(parents=True)
    (output / row['input']).write_text(json.dumps(schema, indent=2) + '\n')
    row['sha256'] = sha256(output / row['input'])
    commands = []
    for previous in original['commands']:
        if previous['condition'] != 'target_template_complexes':
            continue
        command = json.loads(json.dumps(previous))
        command['condition'] = condition
        command['argv'] = [a.replace('target_template_complexes', condition) for a in command['argv']]
        command['expected_prediction_directory'] = command['expected_prediction_directory'].replace('target_template_complexes', condition)
        commands.append(command)
    manifest = {k: original[k] for k in ('package', 'wheel_sha256', 'checkpoint', 'seeds')}
    manifest.update({'source_sha256': sha256(source / 'manifest.json'), 'inputs': [row],
                     'commands': commands, 'prediction_count': len(commands), 'frozen_msas': receipt,
                     'native_parser_validated': False,
                     'limitations': [row['limitation'],
                        'Template-assisted target recovery is not independent folding validation.',
                        'Positive binder is the deposited80-aa IGFBP4 domain; observed62-aa IGF1 lacks unresolved termini.',
                        'Improved confidence does not establish experimental pose recovery; evaluate target-aligned binder RMSD.',
                        'No affinity or biological activity conclusion; historical1WQJ may appear in model training.']})
    save(output, manifest)
    validate_native(output, mol_dir)
    return manifest


def freeze_paired(output, mol_dir):
    """Use the standard CLI's exact paired+unpaired preprocessing function."""
    import csv

    from boltz.data.parse.csv import parse_csv
    from boltz.main import compute_msa

    output = Path(output)
    if (output / 'results').exists():
        raise ValueError('Cannot alter inputs after inference')
    manifest = json.loads((output / 'manifest.json').read_text())
    row = manifest['inputs'][0]
    old_input = output / row['input']
    schema = json.loads(old_input.read_text())
    data = {f"natural_{p['protein']['id']}": p['protein']['sequence'] for p in schema['sequences']}
    directory = output / 'msas/paired_natural'
    directory.mkdir(parents=True, exist_ok=True)
    compute_msa(data=data, target_id='positive_IGF1_IGFBP4', msa_dir=directory,
                msa_server_url='https://api.colabfold.com', msa_pairing_strategy='greedy')
    receipt = {}
    for item in schema['sequences']:
        protein = item['protein']
        path = directory / f"natural_{protein['id']}.csv"
        with path.open() as handle:
            raw = list(csv.DictReader(handle))
        if not raw or raw[0]['sequence'] != protein['sequence']:
            raise ValueError('Native paired CSV first sequence differs from exact query')
        parsed = parse_csv(path, max_seqs=16384)
        if any(parsed.sequences['res_end'] - parsed.sequences['res_start'] != len(protein['sequence'])):
            raise ValueError('Native paired CSV lengths differ from query')
        protein['msa'] = str(path.relative_to(output))
        digest = hashlib.sha256(protein['sequence'].encode()).hexdigest()
        receipt[digest] = {'path': protein['msa'], 'sha256': sha256(path),
                           'native_deduplicated_count': len(parsed.sequences),
                           'raw_sequence_count': len(raw),
                           'paired_rows_including_query': sum(int(r['key']) >= 0 for r in raw),
                           'role': 'natural_binder' if protein['id'] == row['binder_chain'] else 'target',
                           'pairing': 'official compute_msa paired+unpaired, greedy'}
    condition = 'natural_protein_paired_msas_target_template'
    previous_condition = row['condition']
    row['id'] = 'IGF1_IGFBP4__positive_with_paired_msas'
    row['condition'] = condition
    row['context'] = 'positive_natural_protein_paired_msas'
    row['limitation'] = 'Natural-protein control using official Boltz paired+unpaired MSAs for both chains. Both target and binder MSA preprocessing differ from matched de novo condition. No binder template/coordinates or added bond constraints.'
    row['input'] = f"inputs/{condition}/{row['id']}.yaml"
    (output / row['input']).parent.mkdir(parents=True, exist_ok=True)
    (output / row['input']).write_text(json.dumps(schema, indent=2) + '\n')
    row['sha256'] = sha256(output / row['input'])
    (output / 'preliminary_unpaired_input_NOT_USED.json').write_text(old_input.read_text())
    old_input.unlink()
    for command in manifest['commands']:
        command['condition'] = condition
        command['argv'] = [a.replace(previous_condition, condition) for a in command['argv']]
        command['expected_prediction_directory'] = command['expected_prediction_directory'].replace(previous_condition, condition)
    manifest['preliminary_unpaired_msas_not_used'] = manifest['frozen_msas']
    manifest['frozen_msas'] = receipt
    manifest['limitations'][0] = row['limitation']
    manifest['native_parser_validated'] = False
    save(output, manifest)
    validate_native(output, mol_dir)
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path)
    parser.add_argument('--freeze-paired', action='store_true')
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--mol-dir', type=Path, default=Path.home() / '.boltz/mols')
    args = parser.parse_args()
    if args.freeze_paired:
        print(json.dumps(freeze_paired(args.output, args.mol_dir), indent=2))
        return
    if not args.source:
        parser.error('--source is required for initial preparation')
    manifest = prepare(args.source, args.output, args.mol_dir)
    print(json.dumps({'predictions': manifest['prediction_count'], 'frozen_msas': manifest['frozen_msas']}, indent=2))


if __name__ == '__main__':
    main()
