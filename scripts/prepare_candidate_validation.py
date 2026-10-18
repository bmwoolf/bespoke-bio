#!/usr/bin/env python3
"""Prepare fixed candidates and controls for standard Boltz 2.2.1; never infer.

Use --source candidates/manifest.json --output candidate_controls, then
--output candidate_controls --freeze-msas and --validate-native (pinned Boltz
runtime). Target templates alone preserve the experimental assembly geometry;
binders always have empty MSA and no coordinates. All paths in commands/YAMLs
are relative to the prepared directory, which must be inference working dir.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shlex
from pathlib import Path

import gemmi
from prepare_standard_validation import VERSION, WHEEL_SHA256, check_sequence, safe_id, sha256


def normalize_assembly(source, chain_specs, destination):
    structure = gemmi.read_structure(str(source))
    if len(structure) != 1:
        raise ValueError('Template must have exactly one predetermined model')
    expected = {c['chain']: c for c in chain_specs}
    if {c.name for c in structure[0]} != set(expected):
        raise ValueError('Template must contain target chains only')
    structure.entities.clear()
    masks = {}
    for index, chain in enumerate(structure[0], 1):
        spec = expected[chain.name]
        sequence = check_sequence(spec['sequence'])
        positions = []
        for residue in chain:
            position = residue.seqid.num
            code = gemmi.find_tabulated_residue(residue.name).one_letter_code
            if (position < 1 or position > len(sequence) or
                    sequence[position - 1] != code or residue.seqid.icode.strip()):
                raise ValueError('Target template sequence/position mismatch')
            if residue.find_atom('CA', '*') is None:
                raise ValueError('Target template residue has no CA')
            positions.append(position)
            residue.label_seq, residue.subchain = position, chain.name
            residue.entity_id = str(index)
            residue.entity_type = gemmi.EntityType.Polymer
        if positions != sorted(set(positions)):
            raise ValueError('Template residue positions not unique and ordered')
        if positions != [r['input_position'] for r in spec['mapping']]:
            raise ValueError('Template positions differ from audited mapping')
        masks[chain.name] = positions
        entity = gemmi.Entity(str(index))
        entity.entity_type = gemmi.EntityType.Polymer
        entity.polymer_type = gemmi.PolymerType.PeptideL
        entity.subchains = [chain.name]
        entity.full_sequence = gemmi.expand_one_letter_sequence(sequence, gemmi.ResidueKind.AA)
        structure.entities.append(entity)
    structure.connections.clear()
    structure.assemblies.clear()
    structure.make_mmcif_document().write_file(str(destination))
    return masks


def candidate_cases(source):
    cases, alone = [], []
    for candidate in source['fixed_candidates']:
        candidate_id = safe_id(candidate['candidate_id'])
        target = safe_id(candidate['target'])
        binder = check_sequence(candidate['binder_sequence'])
        for key in ('binder', 'target'):
            if hashlib.sha256(candidate[f'{key}_sequence'].encode()).hexdigest() != candidate[f'{key}_sequence_sha256']:
                raise ValueError('Candidate sequence hash mismatch')
        native_filter = 'unavailable'
        bundle = Path(candidate.get('source_bundle', ''))
        if bundle.is_file():
            evaluations = json.loads(bundle.read_text()).get('evaluations', [])
            native_filter = next((e['status'] for e in evaluations if
                                  e.get('subject_id') == candidate_id and
                                  e.get('method') == 'BoltzGen native default filters'), 'unavailable')
        for condition in candidate['conditions']:
            cases.append({'id': f"{target}__{condition['context']}", 'target_id': target,
                          'candidate_id': candidate_id, 'kind': 'fixed_binder_complex',
                          'binder_sequence': binder, 'binder_chain': 'Z',
                          'target_chains': condition['target_chains'],
                          'context': condition['context'],
                          'template_source': condition['target_only_template'],
                          'source_native_filter': native_filter,
                          'limitation': condition['limitation']})
        alone.append({'id': f'{target}__binder_alone', 'target_id': target,
                      'candidate_id': candidate_id, 'kind': 'binder_alone',
                      'binder_sequence': binder, 'binder_chain': 'Z', 'target_chains': [],
                      'context': 'binder_alone', 'source_native_filter': native_filter})
    positive = source['positive_control']
    for name, binder in [('positive', positive['binder_sequence']),
                         ('shuffled_decoy', positive['decoy']['sequence'])]:
        identifier = f'IGF1_IGFBP4__{name}'
        cases.append({'id': identifier, 'target_id': 'IGF1', 'candidate_id': identifier,
                      'kind': 'experimental_positive' if name == 'positive' else 'putative_negative',
                      'binder_sequence': check_sequence(binder), 'binder_chain': 'Z',
                      'target_chains': [{'chain': 'A', 'sequence': positive['target_sequence'],
                                         'mapping': positive['target_mapping']}],
                      'context': name, 'template_source': positive['target_only_template'],
                      'source_native_filter': 'not_applicable',
                      'limitation': positive['limitation'] + ' ' +
                      ('Binder MSA deliberately empty to match de novo protocol; failure does not establish a model bug.'
                       if name == 'positive' else positive['decoy']['limitation'])})
        alone.append({'id': identifier + '__binder_alone', 'target_id': 'IGF1',
                      'candidate_id': identifier, 'kind': 'binder_alone',
                      'binder_sequence': binder, 'binder_chain': 'Z', 'target_chains': [],
                      'context': name + '_binder_alone', 'source_native_filter': 'not_applicable'})
    return cases + alone


def prepare(source_path, output, seeds=(11, 22, 33)):
    source_path, output = Path(source_path), Path(output)
    if output.exists() and any(output.iterdir()):
        raise ValueError('Output directory must be new or empty')
    source = json.loads(source_path.read_text())
    cases = candidate_cases(source)
    output.mkdir(parents=True, exist_ok=True)
    (output / 'templates').mkdir()
    inputs = []
    for case in cases:
        safe_id(case['id'])
        group = 'binder_alone' if case['kind'] == 'binder_alone' else 'target_template_complexes'
        schema = {'version': 1, 'sequences': []}
        for chain in case['target_chains']:
            if chain['chain'] == 'Z':
                raise ValueError('Binder chain collides with target')
            schema['sequences'].append({'protein': {'id': chain['chain'],
                                                    'sequence': check_sequence(chain['sequence'])}})
        schema['sequences'].append({'protein': {'id': 'Z', 'sequence': case['binder_sequence'], 'msa': 'empty'}})
        row = {**case, 'condition': group, 'target_template': bool(case['target_chains'])}
        if case['target_chains']:
            template = Path('templates') / (case['id'] + '.cif')
            row['observed_template_positions'] = normalize_assembly(
                Path(case['template_source']), case['target_chains'], output / template)
            row['template_sha256'] = sha256(output / template)
            chain_ids = [c['chain'] for c in case['target_chains']]
            schema['templates'] = [{'cif': str(template), 'chain_id': chain_ids,
                                    'template_id': chain_ids, 'force': False}]
        path = Path('inputs') / group / (case['id'] + '.yaml')
        (output / path.parent).mkdir(parents=True, exist_ok=True)
        (output / path).write_text(json.dumps(schema, indent=2) + '\n')
        row.update({'input': str(path), 'sha256': sha256(output / path)})
        inputs.append(row)
    commands = []
    for condition in ('target_template_complexes', 'binder_alone'):
        for seed in seeds:
            argv = ['boltz', 'predict', f'inputs/{condition}', '--model', 'boltz2',
                    '--out_dir', f'results/{condition}/seed_{seed}', '--cache', 'cache',
                    '--seed', str(seed), '--recycling_steps', '3', '--sampling_steps', '200',
                    '--diffusion_samples', '1', '--max_parallel_samples', '1', '--devices', '1',
                    '--accelerator', 'gpu', '--output_format', 'mmcif', '--num_workers', '0',
                    '--no_kernels', '--write_full_pae']
            if condition == 'target_template_complexes':
                argv.append('--use_msa_server')
            commands.append({'condition': condition, 'seed': seed, 'argv': argv,
                             'expected_prediction_directory':
                             f'results/{condition}/seed_{seed}/boltz_results_{condition}/predictions'})
    manifest = {'package': f'boltz=={VERSION}', 'wheel_sha256': WHEEL_SHA256,
                'checkpoint': 'boltz2_conf.ckpt', 'source_sha256': sha256(source_path),
                'seeds': list(seeds), 'inputs': inputs, 'commands': commands,
                'prediction_count': len(inputs) * len(seeds), 'native_parser_validated': False,
                'limitations': ['All complexes are target-template conditioned, not independent target folds.',
                                'Template coordinates never include a binder; assembly can still rearrange during inference.',
                                'No binding affinity, selectivity, or biological activity is established by these confidence scores.',
                                'PF4 prior pose clashes with intact tetramer; prior native filter status remains a separate fact.',
                                'Positive binder MSA empty deliberately matches de novo protocol; failure does not establish a model bug.',
                                'Shuffled decoy is a putative negative, not an experimentally validated nonbinder.',
                                'Three requested matched seeds are not equilibrium samples or a significance test.']}
    save(output, manifest)
    return manifest


def save(output, manifest):
    (output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    (output / 'commands.sh').write_text('#!/usr/bin/env bash\nset -euo pipefail\n' +
        '\n'.join(shlex.join(c['argv']) for c in manifest['commands']) + '\n')


def freeze_msas(output):
    from boltz.data.msa.mmseqs2 import run_mmseqs2
    from boltz.data.parse.a3m import parse_a3m

    output = Path(output)
    if (output / 'results').exists():
        raise ValueError('Cannot change inputs after inference')
    manifest = json.loads((output / 'manifest.json').read_text())
    sequences = {}
    for row in manifest['inputs']:
        for chain in row['target_chains']:
            digest = hashlib.sha256(chain['sequence'].encode()).hexdigest()
            sequences[digest] = chain['sequence']
    receipt = {}
    (output / 'msas').mkdir(exist_ok=True)
    for digest, sequence in sequences.items():
        path = Path('msas') / f'{digest}.a3m'
        if not (output / path).exists():
            alignments = run_mmseqs2([sequence], str(output / 'msas' / (digest + '_raw')),
                                    use_env=True, use_pairing=False,
                                    host_url='https://api.colabfold.com')
            if len(alignments) != 1:
                raise ValueError('Unexpected MSA response count')
            (output / path).write_text(alignments[0])
        lines = (output / path).read_text().splitlines()
        query = ''.join(lines[1:next((i for i in range(1, len(lines))
                                     if lines[i].startswith('>')), len(lines))])
        if query != sequence:
            raise ValueError('MSA query mismatch')
        parsed = parse_a3m(output / path, taxonomy=None, max_seqs=16384)
        lengths = parsed.sequences['res_end'] - parsed.sequences['res_start']
        if not all(int(n) == len(sequence) for n in lengths):
            raise ValueError('Native MSA parser alignment length mismatch')
        receipt[digest] = {'path': str(path), 'sha256': sha256(output / path),
                           'raw_sequence_count': sum(line.startswith('>') for line in lines),
                           'native_deduplicated_count': len(parsed.sequences)}
    for row in manifest['inputs']:
        path = output / row['input']
        schema = json.loads(path.read_text())
        for item in schema['sequences']:
            protein = item['protein']
            if protein['id'] != row['binder_chain']:
                digest = hashlib.sha256(protein['sequence'].encode()).hexdigest()
                protein['msa'] = receipt[digest]['path']
        path.write_text(json.dumps(schema, indent=2) + '\n')
        row['sha256'] = sha256(path)
    for command in manifest['commands']:
        command['argv'] = [a for a in command['argv'] if a != '--use_msa_server']
    manifest['frozen_msas'] = receipt
    manifest['native_parser_validated'] = False
    save(output, manifest)
    return receipt


def validate_native(output, mol_dir):
    from boltz.data.feature.featurizerv2 import process_template_features
    from boltz.data.parse.schema import parse_boltz_schema
    from boltz.data.tokenize.boltz2 import Boltz2Tokenizer
    from boltz.data.types import Input

    output, mol_dir = Path(output).resolve(), Path(mol_dir).resolve()
    manifest = json.loads((output / 'manifest.json').read_text())
    receipt = []
    old_cwd = Path.cwd()
    try:
        os.chdir(output)
        for row in manifest['inputs']:
            schema = json.loads(Path(row['input']).read_text())
            parsed = parse_boltz_schema(row['id'], schema, {}, mol_dir, boltz_2=True)
            expected = {c['chain']: c['sequence'] for c in row['target_chains']}
            expected[row['binder_chain']] = row['binder_sequence']
            actual = {str(c.chain_name): parsed.sequences[c.entity_id] for c in parsed.record.chains}
            if actual != expected:
                raise ValueError('Native chain sequences differ')
            records = parsed.record.templates or []
            if {r.query_chain for r in records} != {c['chain'] for c in row['target_chains']}:
                raise ValueError('Template applied to wrong chain set')
            for record in records:
                if (record.query_chain != record.template_chain or record.force or
                        record.query_st != 0 or record.query_en != len(expected[record.query_chain])):
                    raise ValueError('Native template mapping or forcing mismatch')
            for template in parsed.templates.values():
                for chain in template.chains:
                    residues = template.residues[chain['res_idx']:chain['res_idx'] + chain['res_num']]
                    observed = [int(r['res_idx']) + 1 for r in residues if r['is_present']]
                    if observed != row['observed_template_positions'][str(chain['name'])]:
                        raise ValueError('Native observed positions differ from audited chain mapping')
            feature_audit = {}
            if records:
                tokens = Boltz2Tokenizer().tokenize(Input(
                    parsed.structure, {}, record=parsed.record, templates=parsed.templates))
                features = process_template_features(tokens, len(tokens.tokens))
                target_visibility = []
                for chain in parsed.structure.chains:
                    idx = tokens.tokens['asym_id'] == chain['asym_id']
                    name = str(chain['name'])
                    mask_sum = features['template_mask'][0, idx].sum().item()
                    visibility = features['visibility_ids'][0, idx].unique().tolist()
                    if name == row['binder_chain']:
                        if mask_sum != 0:
                            raise ValueError('Binder has template feature conditioning')
                    else:
                        target_visibility.extend(visibility)
                    feature_audit[name] = {'template_token_mask_sum': mask_sum,
                                           'visibility_ids': visibility}
                if len(set(target_visibility)) != 1:
                    raise ValueError('Target assembly chains do not share template geometry visibility')
            receipt.append({'id': row['id'], 'chain_lengths': {k: len(v) for k, v in actual.items()},
                            'template_feature_audit': feature_audit,
                            'template_chains': [r.query_chain for r in records], 'binder_template': False})
    finally:
        os.chdir(old_cwd)
    manifest['native_parser_validated'] = True
    (output / 'native_validation.json').write_text(json.dumps(receipt, indent=2) + '\n')
    save(output, manifest)
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--freeze-msas', action='store_true')
    parser.add_argument('--validate-native', action='store_true')
    parser.add_argument('--mol-dir', type=Path, default=Path.home() / '.boltz/mols')
    args = parser.parse_args()
    if args.freeze_msas:
        print(json.dumps(freeze_msas(args.output), indent=2))
    elif args.validate_native:
        print(json.dumps({'validated': len(validate_native(args.output, args.mol_dir))}))
    elif args.source:
        print(json.dumps({'predictions': prepare(args.source, args.output)['prediction_count']}))
    else:
        parser.error('--source required to prepare')


if __name__ == '__main__':
    main()
