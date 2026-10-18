#!/usr/bin/env python3
"""Prepare, never execute, standard Boltz 2.2.1 folding controls.

Input: researcher controls JSON with targets[{gene_symbol,sequence,references}].
A template reference has role=template, matched_structure, matched_chain, and
mapping[].input_position. Matched CIF residue numbers are exact input positions.
Optional binders JSON: [{id,target_id,sequence}]; adds fixed-sequence complexes
with empty binder MSA and ONLY a target template. No candidate pose is consumed.
Run commands from output directory; portable paths remain relative to that cwd.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shlex
from pathlib import Path

import gemmi

VERSION = '2.2.1'
WHEEL_SHA256 = 'b8c62bbdede1922931d9203118f62c858f11aa699bf91fd4c05a5ed6a6d8b4fc'
AA = set('ACDEFGHIKLMNPQRSTVWY')


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def safe_id(value):
    if not re.fullmatch(r'[A-Za-z0-9_]+', value):
        raise ValueError(f'Unsafe identifier: {value!r}')
    return value


def check_sequence(sequence):
    if not sequence or set(sequence) - AA:
        raise ValueError('Expected a nonempty exact canonical amino-acid sequence')
    return sequence


def normalized_template(reference, sequence, destination, base):
    """Add full polymer sequence metadata; retain only existing target atoms.

    Missing coordinates are not synthesized. Full entity sequence is essential
    for Boltz's native mmCIF parser to retain unresolved terminal positions.
    """
    source = Path(reference['matched_structure'])
    if not source.is_absolute():
        source = base / source
    structure = gemmi.read_structure(str(source))
    if len(structure) != 1 or len(structure[0]) != 1:
        raise ValueError('Template must contain exactly one target model and chain')
    chain = structure[0][0]
    if chain.name != reference.get('matched_chain', 'A'):
        raise ValueError('Template chain differs from recorded chain')
    positions = []
    chain.name = 'A'
    for residue in chain:
        position = residue.seqid.num
        code = gemmi.find_tabulated_residue(residue.name).one_letter_code
        if (position < 1 or position > len(sequence) or
                sequence[position - 1] != code or residue.seqid.icode.strip()):
            raise ValueError(f'Template residue {residue.seqid} disagrees with sequence')
        if residue.find_atom('CA', '*') is None:
            raise ValueError('Template residue has no CA')
        positions.append(position)
        residue.label_seq = position
        residue.subchain = 'A'
        residue.entity_id = '1'
        residue.entity_type = gemmi.EntityType.Polymer
    if positions != sorted(set(positions)):
        raise ValueError('Template positions must be unique and ordered')
    expected = [m['input_position'] for m in reference['mapping']]
    if positions != expected:
        raise ValueError('Template coordinates disagree with audited residue mapping')
    structure.entities.clear()
    entity = gemmi.Entity('1')
    entity.entity_type = gemmi.EntityType.Polymer
    entity.polymer_type = gemmi.PolymerType.PeptideL
    entity.subchains = ['A']
    entity.full_sequence = gemmi.expand_one_letter_sequence(sequence, gemmi.ResidueKind.AA)
    structure.entities.append(entity)
    structure.connections.clear()
    structure.assemblies.clear()
    structure.make_mmcif_document().write_file(str(destination))
    return {'path': destination.name, 'sha256': sha256(destination),
            'source_sha256': sha256(source), 'observed_positions': positions,
            'unobserved_positions': sorted(set(range(1, len(sequence) + 1)) - set(positions)),
            'chain': 'A', 'coordinates_invented': False}


def input_schema(sequence, msa, template=None, binder=None):
    target = {'id': 'A', 'sequence': check_sequence(sequence)}
    if msa == 'empty':
        target['msa'] = 'empty'
    elif msa != 'server':
        raise ValueError('MSA condition must be empty or server')
    result = {'version': 1, 'sequences': [{'protein': target}]}
    if binder:
        result['sequences'].append({'protein': {
            'id': 'B', 'sequence': check_sequence(binder), 'msa': 'empty'}})
    if template:
        result['templates'] = [{'cif': template, 'chain_id': 'A',
                                'template_id': 'A', 'force': False}]
    return result


def prepare(controls_path, output, seeds=(11, 22, 33), binders=(), base=None):
    controls_path, output = Path(controls_path).resolve(), Path(output)
    base = Path.cwd() if base is None else Path(base)
    if output.exists() and any(output.iterdir()):
        raise ValueError('Output directory must be new or empty')
    controls = json.loads(controls_path.read_text())
    # Validate everything before publishing the manifest; partial prep has no manifest.
    output.mkdir(parents=True, exist_ok=True)
    (output / 'templates').mkdir()
    cases, templates = [], {}
    targets = {safe_id(t['gene_symbol']): t for t in controls['targets']}
    for gene, target in targets.items():
        sequence = check_sequence(target['sequence'])
        recorded = target.get('sequence_sha256')
        if recorded and recorded != hashlib.sha256(sequence.encode()).hexdigest():
            raise ValueError(f'{gene} sequence hash mismatch')
        refs = [r for r in target['references'] if r['role'] == 'template']
        if len(refs) != 1:
            raise ValueError(f'{gene} requires exactly one predetermined template')
        templates[gene] = normalized_template(
            refs[0], sequence, output / 'templates' / f'{gene}.cif', base)
        templates[gene]['pdb_id'] = refs[0].get('pdb_id')
        cases.append({'id': gene, 'target_id': gene, 'sequence': sequence,
                      'kind': 'target_alone', 'binder_sequence': None})
    for binder in binders:
        gene = binder['target_id']
        if gene not in targets:
            raise ValueError(f'Unknown binder target {gene}')
        cases.append({'id': safe_id(gene + '__' + binder['id']), 'target_id': gene,
                      'sequence': targets[gene]['sequence'], 'kind': 'fixed_binder_complex',
                      'binder_sequence': check_sequence(binder['sequence'])})
    if len({c['id'] for c in cases}) != len(cases):
        raise ValueError('Case IDs must be unique')
    inputs, commands = [], []
    for msa in ('empty', 'server'):
        for templated in (False, True):
            condition = f'{msa}_msa__' + ('target_template' if templated else 'no_template')
            directory = Path('inputs') / condition
            (output / directory).mkdir(parents=True)
            for case in cases:
                path = directory / (case['id'] + '.yaml')
                schema = input_schema(case['sequence'], msa,
                                      f"templates/{case['target_id']}.cif" if templated else None,
                                      case['binder_sequence'])
                # JSON is YAML 1.2, accepted by the standard Boltz YAML loader.
                (output / path).write_text(json.dumps(schema, indent=2) + '\n')
                inputs.append({**case, 'condition': condition, 'msa': msa,
                               'target_template': templated, 'input': str(path),
                               'sha256': sha256(output / path)})
            for seed in seeds:
                argv = ['boltz', 'predict', str(directory), '--model', 'boltz2',
                        '--out_dir', f'results/{condition}/seed_{seed}', '--cache', 'cache',
                        '--seed', str(seed), '--recycling_steps', '3', '--sampling_steps', '200',
                        '--diffusion_samples', '1', '--max_parallel_samples', '1',
                        '--devices', '1', '--accelerator', 'gpu', '--output_format', 'mmcif',
                        '--num_workers', '0', '--no_kernels', '--write_full_pae']
                if msa == 'server':
                    argv.append('--use_msa_server')
                commands.append({'condition': condition, 'seed': seed, 'argv': argv,
                                 'expected_prediction_directory':
                                 f'results/{condition}/seed_{seed}/boltz_results_{condition}/predictions'})
    manifest = {
        'package': f'boltz=={VERSION}', 'wheel_sha256': WHEEL_SHA256,
        'checkpoint': 'boltz2_conf.ckpt',
        'checkpoint_url': 'https://huggingface.co/boltz-community/boltz-2/resolve/main/boltz2_conf.ckpt',
        'checkpoint_policy': 'Record downloaded checkpoint SHA256 before inference; do not use BoltzGen checkpoints.',
        'source': 'https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md',
        'controls_sha256': sha256(controls_path), 'seeds': list(seeds),
        'templates': templates, 'inputs': inputs, 'commands': commands,
        'prediction_count': len(inputs) * len(seeds),
        'native_parser_validated': False,
        'limitations': [
            'Template-assisted recovery is conditioning, not independent folding validation.',
            'Historical held-out structures are withheld from these inputs, not necessarily model training.',
            'Matched requested seeds do not establish bitwise reproducibility or independent equilibrium samples.',
            'MSA server arms transmit public target sequences to the standard ColabFold service.',
            'Archive generated MSAs and hashes; service versions/results may change between requests.',
            'Confidence metrics are not binding affinity or evidence of biological activity.']}
    (output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    (output / 'commands.sh').write_text('#!/usr/bin/env bash\nset -euo pipefail\n'
                                       '# Run from this prepared directory. No jobs started by preparation.\n' +
                                       '\n'.join(shlex.join(c['argv']) for c in commands) + '\n')
    return manifest


def freeze_msas(output, server='https://api.colabfold.com'):
    """CPU/network-only precomputation using the pinned package's own MSA client.

    Fetch each unique target once; freeze identical A3M inputs across conditions
    and seeds. Refuse an existing inference output to avoid changing its inputs.
    """
    from boltz.data.msa.mmseqs2 import run_mmseqs2

    output = Path(output)
    if (output / 'results').exists():
        raise ValueError('Cannot change inputs after inference outputs exist')
    manifest = json.loads((output / 'manifest.json').read_text())
    msa_dir = output / 'msas'
    msa_dir.mkdir(exist_ok=True)
    targets = {row['target_id']: row['sequence'] for row in manifest['inputs']}
    receipt = {}
    for gene, sequence in targets.items():
        path = msa_dir / f'{gene}.a3m'
        if not path.exists():
            alignments = run_mmseqs2([sequence], str(msa_dir / f'{gene}_raw'),
                                     use_env=True, use_pairing=False, host_url=server)
            if len(alignments) != 1:
                raise ValueError('MSA service returned unexpected number of alignments')
            path.write_text(alignments[0])
        lines = path.read_text().splitlines()
        query = ''.join(lines[1:next((i for i in range(1, len(lines))
                                     if lines[i].startswith('>')), len(lines))])
        if query != sequence:
            raise ValueError(f'{gene} MSA query differs from exact target sequence')
        receipt[gene] = {'path': f'msas/{gene}.a3m', 'sha256': sha256(path),
                         'sequence_count': sum(line.startswith('>') for line in lines),
                         'server': server}
    for row in manifest['inputs']:
        if row['msa'] != 'server':
            continue
        path = output / row['input']
        schema = json.loads(path.read_text())
        schema['sequences'][0]['protein']['msa'] = receipt[row['target_id']]['path']
        path.write_text(json.dumps(schema, indent=2) + '\n')
        row['sha256'] = sha256(path)
    for command in manifest['commands']:
        command['argv'] = [a for a in command['argv'] if a != '--use_msa_server']
    manifest['frozen_msas'] = receipt
    manifest['native_parser_validated'] = False
    (output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    (output / 'commands.sh').write_text('#!/usr/bin/env bash\nset -euo pipefail\n' +
        '\n'.join(shlex.join(c['argv']) for c in manifest['commands']) + '\n')
    return receipt


def validate_native(output, mol_dir):
    """Parse every YAML with standard Boltz; assert sequence and template mapping."""
    import os

    from boltz.data.parse.schema import parse_boltz_schema

    output = Path(output).resolve()
    mol_dir = Path(mol_dir).resolve()
    manifest = json.loads((output / 'manifest.json').read_text())
    original_cwd = Path.cwd()
    rows = []
    try:
        os.chdir(output)
        for row in manifest['inputs']:
            path = Path(row['input'])
            target = parse_boltz_schema(path.stem, json.loads(path.read_text()), {},
                                       mol_dir, boltz_2=True)
            if target.sequences[0] != row['sequence']:
                raise ValueError('Native parser changed target sequence')
            if row['binder_sequence'] and target.sequences[1] != row['binder_sequence']:
                raise ValueError('Native parser changed binder sequence')
            records = target.record.templates or []
            if bool(records) != row['target_template']:
                raise ValueError('Native template presence differs from condition')
            for record in records:
                if record.query_chain != 'A' or record.template_chain != 'A':
                    raise ValueError('Template applied outside target chain')
                if record.query_st != 0 or record.query_en != len(row['sequence']):
                    raise ValueError('Native template alignment changed target numbering')
                if record.force:
                    raise ValueError('Unexpected forced template potential')
            observed = []
            for template in target.templates.values():
                observed = [int(r['res_idx']) + 1 for r in template.residues if r['is_present']]
                if observed != manifest['templates'][row['target_id']]['observed_positions']:
                    raise ValueError('Native template observed-residue mask differs from audit')
            rows.append({'input': row['input'], 'target_length': len(target.sequences[0]),
                         'template_records': len(records), 'observed_template_positions': observed})
    finally:
        os.chdir(original_cwd)
    manifest['native_parser_validated'] = True
    (output / 'native_validation.json').write_text(json.dumps(rows, indent=2) + '\n')
    (output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--controls', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--binders', type=Path)
    parser.add_argument('--seeds', type=int, nargs='+', default=[11, 22, 33])
    parser.add_argument('--freeze-msas', action='store_true')
    parser.add_argument('--validate-native', action='store_true')
    parser.add_argument('--mol-dir', type=Path, default=Path.home() / '.boltz/mols')
    args = parser.parse_args()
    if args.freeze_msas:
        print(json.dumps(freeze_msas(args.output), indent=2))
        return
    if args.validate_native:
        print(json.dumps({'validated': len(validate_native(args.output, args.mol_dir))}))
        return
    if not args.controls:
        parser.error('--controls is required for preparation')
    binders = json.loads(args.binders.read_text()) if args.binders else []
    manifest = prepare(args.controls, args.output, args.seeds, binders)
    print(json.dumps({'inputs': len(manifest['inputs']), 'predictions': manifest['prediction_count']}))


if __name__ == '__main__':
    main()
