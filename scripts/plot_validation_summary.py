#!/usr/bin/env python3
"""Generate scientific summary figures from actual saved validation results only."""
from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze_candidate_validation import aligned_pose_metrics, read_chains
from matplotlib.lines import Line2D
from render_focused_binder_review import collect_native

COLORS = {11: '#0072B2', 22: '#D55E00', 33: '#009E73'}
CONDITIONS = ['empty_msa__no_template', 'server_msa__no_template',
              'empty_msa__target_template', 'server_msa__target_template']
LABELS = ['MSA −\nTemplate −', 'MSA +\nTemplate −',
          'MSA −\nTemplate +', 'MSA +\nTemplate +']


def style():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.titlesize': 11,
                         'axes.labelsize': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                         'figure.facecolor': 'white', 'savefig.facecolor': 'white',
                         'axes.grid': True, 'grid.alpha': .15, 'grid.linewidth': .5,
                         'pdf.fonttype': 42, 'ps.fonttype': 42})


def finish(fig, output, name, footnote):
    fig.text(.04, .02, footnote, fontsize=9, va='bottom', color='#49515c')
    fig.subplots_adjust(bottom=.25 if fig.get_figheight() < 7 else .16,
                        top=.83 if fig.get_figheight() < 7 else .90, hspace=.55, wspace=.32)
    fig.savefig(output / f'{name}.png', dpi=190)
    fig.savefig(output / f'{name}.pdf')
    plt.close(fig)


def legend_seeds(fig):
    fig.legend([Line2D([], [], marker='o', linestyle='', color=c) for c in COLORS.values()],
               [f'Seed {s}' for s in COLORS], loc='upper right', ncol=3, frameon=False)


def target_controls(data, output):
    genes = ['FGF2', 'IL1RN', 'SERPINA1']
    fig, axes = plt.subplots(3, 3, figsize=(13, 10))
    fig.suptitle(f"Target folding controls: {data['completed_predictions']}/{data['expected_predictions']} predictions", x=.04, ha='left', fontsize=16)
    legend_seeds(fig)
    rows = [r for r in data['rows'] if r['status'] == 'completed']
    points = []
    for gi, gene in enumerate(genes):
        for ci, condition in enumerate(CONDITIONS):
            values = [r for r in rows if r['target_id'] == gene and r['condition'] == condition]
            for row in values:
                heldout = next((r['rmsd_angstrom'] for r in row['reference_comparisons']
                                if r['role'] == 'held_out' and r['status'] == 'compared'), None)
                metrics = [row['chains']['A']['plddt']['mean'], row['chains']['A']['ptm'], heldout]
                points.append({'target': gene, 'condition': condition, 'seed': row['seed'],
                               'plddt_0_100': metrics[0], 'ptm': metrics[1], 'heldout_CA_RMSD_A': metrics[2]})
                offset = (list(COLORS).index(row['seed']) - 1) * .09
                for mi, value in enumerate(metrics):
                    if value is not None:
                        axes[gi, mi].scatter(ci + offset, value, color=COLORS[row['seed']], s=32, zorder=3)
        for mi, title in enumerate(['Local pLDDT ↑', 'Global pTM ↑', 'Experimental CA RMSD ↓']):
            ax = axes[gi, mi]
            ax.set_title(f'{gene} · {title}')
            ax.set_xticks(range(4), LABELS, fontsize=8)
            ax.set_xlim(-.4, 3.4)
            ax.set_ylabel(['pLDDT / 100', 'pTM / 1', 'RMSD (Å)'][mi])
            if mi == 0:
                ax.set_ylim(0, 100)
            elif mi == 1:
                ax.set_ylim(0, 1)
            else:
                ax.set_ylim(bottom=0)
    finish(fig, output, '01_target_controls',
           'MSA + = frozen alignment; MSA − = explicit empty alignment. Template + = target coordinates supplied. Dots are predictions.\nExperimental references were withheld from inputs; historical structures may have been in training. Three seeds are not a significance test.')
    return points


def native_overview(native, output):
    fig, axes = plt.subplots(1, 3, figsize=(13, 5.6))
    fig.suptitle('Native design filters and confidence are different measures', x=.04, ha='left', fontsize=16)
    genes = sorted({r['target_id'] for r in native})
    for i, gene in enumerate(genes):
        rows = [r for r in native if r['target_id'] == gene]
        passed = sum(r['native_filter_status'] == 'pass' for r in rows)
        axes[0].bar(i, len(rows) - passed, color='#bbc1c9', label='Native fail' if i == 0 else None)
        axes[0].bar(i, passed, bottom=len(rows) - passed, color='#0072B2', label='Native pass' if i == 0 else None)
        axes[0].text(i, len(rows) + .4, f'{passed}/{len(rows)} pass', ha='center')
        if i < 2:
            ax = axes[i + 1]
            for row in rows:
                metrics = row['native_metrics']
                x, y = metrics.get('design_to_target_iptm'), metrics.get('esmfold2_ipsae_min')
                if x is not None and y is not None:
                    passed = row['native_filter_status'] == 'pass'
                    ax.scatter(x, y, marker='o' if passed else 'x', s=60 if passed else 28,
                               color='#0072B2' if passed else '#9b9fa5', zorder=3 if passed else 2)
                    if passed:
                        ax.annotate(row['candidate_id'][-6:], (x, y), xytext=(5, 18 if row['candidate_id'].endswith('00c4dc') else 5), textcoords='offset points', fontsize=8)
            ax.set(title=gene, xlabel='Native Boltz design→target ipTM', ylabel='Native ESMFold2 ipSAE minimum', xlim=(-.03, 1), ylim=(-.03, 1))
    axes[0].set_xticks(range(len(genes)), genes)
    axes[0].set_ylabel('Actual generated candidates')
    axes[0].set_title('Native filtering outcome')
    axes[0].set_ylim(0, max([sum(r['target_id'] == g for r in native) for g in genes] or [1]) + 5)
    axes[0].legend(frameon=False)
    finish(fig, output, '02_native_filtering',
           'Native pass checks structural recovery and composition; it does NOT require high interface confidence.\nThe two native model scores are distinct metrics. Neither native pass nor either score establishes experimental binding.')


def seed_consistency(cases, prepared):
    results = []
    complex_rows = [r for r in cases if r['kind'] == 'fixed_binder_complex' and r['status'] == 'analyzed']
    for identifier in sorted({r['candidate_id'] for r in complex_rows}):
        selected = [r for r in complex_rows if r['candidate_id'] == identifier]
        for left, right in combinations(selected, 2):
            try:
                # Exact sequences from already-validated prepared YAMLs, not inferred chain roles.
                manifest = json.loads((prepared / 'manifest.json').read_text())
                spec = next(r for r in manifest['inputs'] if r['id'] == left['id'])
                expected = {c['chain']: c['sequence'] for c in spec['target_chains']}
                expected[spec['binder_chain']] = spec['binder_sequence']
                a = read_chains(prepared / left['files']['structure'], expected)
                b = read_chains(prepared / right['files']['structure'], expected)
                def ca(chain):
                    return [[atom.pos.x, atom.pos.y, atom.pos.z] for residue in chain['residues']
                            if (atom := residue.find_atom('CA', '*')) is not None]
                metrics = aligned_pose_metrics(ca(a['A']), ca(b['A']), ca(a['Z']), ca(b['Z']))
                results.append({'candidate_id': identifier, 'seed_pair': [left['seed'], right['seed']],
                                'status': 'assessed', **metrics})
            except (OSError, ValueError, KeyError, StopIteration) as exc:
                results.append({'candidate_id': identifier, 'seed_pair': [left['seed'], right['seed']],
                                'status': 'unavailable', 'reason': str(exc)})
    return results


def focused(data, prepared, output):
    rows = [r for r in data['cases'] if r['status'] == 'analyzed']
    candidates = sorted({(r['target_id'], r['candidate_id']) for r in data['cases']})
    fig, axes = plt.subplots(2, 3, figsize=(13, 9))
    fig.suptitle(f"Focused candidates: standard refolding {data['analyzed_predictions']}/{data['expected_predictions']} complete", x=.04, ha='left', fontsize=16)
    legend_seeds(fig)
    titles = ['A  Binder local confidence ↑', 'B  Directional interface confidence ↑',
              'C  Binder-alone agreement with design ↓', 'D  Original design pose agreement ↓',
              'E  Target vs experimental reference ↓', 'F  Predicted pose consistency across seeds ↓']
    ys = ['Mean pLDDT / 100', 'Chain-pair ipTM / 1', 'Binder CA RMSD (Å)',
          'Binder CA RMSD after target fit (Å)', 'Held-out target CA RMSD (Å)', 'Binder CA RMSD after target fit (Å)']
    for index, (gene, candidate) in enumerate(candidates):
        for row in [r for r in rows if r['candidate_id'] == candidate and r['target_id'] == gene]:
            offset = (list(COLORS).index(row['seed']) - 1) * .08
            alone = row['kind'] == 'binder_alone'
            x = index + offset + (-.18 if alone else .18)
            axes[0, 0].scatter(x, row['binder_plddt_0_100']['mean'], color=COLORS[row['seed']],
                               marker='s' if alone else 'o', s=40)
            pose = row.get('binder_reference_pose', {})
            if alone:
                value = pose.get('binder_own_fold_CA_RMSD_A')
                if value is not None:
                    axes[0, 2].scatter(index + offset, value, color=COLORS[row['seed']], s=40)
            else:
                for direction, marker, shift in [('A|Z', 'o', -.13), ('Z|A', '^', .13)]:
                    value = row['binder_target_chainpair_iptm'].get(direction)
                    if value is not None:
                        axes[0, 1].scatter(index + offset + shift, value, color=COLORS[row['seed']], marker=marker, s=40)
                value = pose.get('binder_CA_RMSD_after_target_alignment_A')
                if value is not None:
                    axes[1, 0].scatter(index + offset, value, color=COLORS[row['seed']], s=40)
                refs = row.get('posthoc_target_reference_geometry', {}).get('references', [])
                for ref in refs:
                    if ref.get('role') == 'held_out' and ref.get('status') == 'assessed':
                        axes[1, 1].scatter(index + offset, ref['global_CA_RMSD_A'], color=COLORS[row['seed']], s=40)
    consistency = seed_consistency(rows, prepared)
    for i, (_, candidate) in enumerate(candidates):
        for j, row in enumerate(r for r in consistency if r['candidate_id'] == candidate and r['status'] == 'assessed'):
            axes[1, 2].scatter(i + (j - 1) * .10, row['binder_CA_RMSD_after_target_alignment_A'], color='#5f4a88', s=40)
    for ax, title, ylabel in zip(axes.flat, titles, ys):
        ax.set_title(title, loc='left')
        ax.set_xticks(range(len(candidates)), [f'{g}\n{c[-6:]}' for g, c in candidates], fontsize=9)
        ax.set_ylabel(ylabel)
        ax.set_xlim(-.5, max(.5, len(candidates) - .5))
        ax.set_ylim(bottom=0)
    axes[0, 0].set_ylim(0, 100)
    axes[0, 1].set_ylim(0, 1)
    axes[0, 0].legend(handles=[Line2D([], [], marker='s', color='#333', linestyle='', label='Alone'),
                                     Line2D([], [], marker='o', color='#333', linestyle='', label='In complex')], frameon=False, fontsize=8, loc='lower left')
    axes[0, 1].legend(handles=[Line2D([], [], marker='o', color='#333', linestyle='', label='A|Z'),
                                     Line2D([], [], marker='^', color='#333', linestyle='', label='Z|A')], frameon=False, fontsize=8, loc='lower left')
    finish(fig, output, '03_focused_candidate_metrics',
           'No coordinates/templates supplied to standard predictions; target frozen MSA, binder empty MSA. A=target, Z=binder.\nC–D compare to computational designs (agreement, NOT accuracy); F compares three non-independent seed pairs. High confidence is not binding evidence.')
    return consistency


def positive(data, output):
    fig, axes = plt.subplots(1, 3, figsize=(12, 5.6))
    fig.suptitle('Known experimental complex: did the local interface recover?', x=.04, ha='left', fontsize=15)
    legend_seeds(fig)
    for row in data['cases']:
        if row['status'] != 'assessed':
            continue
        x = (0 if row['condition'] == 'candidate_controls' else 1) + (list(COLORS).index(row['seed']) - 1) * .09
        contact = row['reference_contact_recovery']
        for ax, value in zip(axes, [contact['contact_recall'], contact['contact_precision'], contact['reference_interface_CA_jointfit_RMSD_A']]):
            ax.scatter(x, value, color=COLORS[row['seed']], s=45)
    for ax, title, ylabel in zip(axes, ['Experimental contact recall ↑', 'Experimental contact precision ↑', 'Reference-interface CA RMSD ↓'],
                                ['Recovered / reference contacts', 'Recovered / predicted contacts', 'CA-only joint-fit RMSD (Å)']):
        ax.set(title=title, ylabel=ylabel, xlim=(-.4, 1.4), xticks=[0, 1], xticklabels=['Target MSA\nBinder empty', 'Both MSAs\nPaired + unpaired'])
        ax.set_ylim(bottom=0)
    axes[0].set_ylim(0, 1)
    axes[1].set_ylim(0, 1)
    finish(fig, output, '04_experimental_positive_control',
           '1WQJ: target-only template in both conditions; binder coordinates withheld. Contacts: occupied heavy-atom distance <5Å.\nCA-only interface RMSD is not DockQ. Natural-chain MSA success does not establish that de novo designed binders work; training independence is unknown.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path('outputs/folding_validation'))
    args = parser.parse_args()
    root, output = args.root, args.root / 'summary_figures'
    output.mkdir(parents=True, exist_ok=True)
    style()
    controls = json.loads((root / 'analysis/analysis.json').read_text())
    focused_data = json.loads((root / 'focused_validation_report/analysis.json').read_text())
    positive_data = json.loads((root / 'positive_pose_audit.json').read_text())
    prepared = root / 'download/focused_validation'
    native = collect_native(root / 'focused_design', prepared)
    control_points = target_controls(controls, output)
    native_overview(native, output)
    consistency = focused(focused_data, prepared, output)
    positive(positive_data, output)
    (output / 'figure_data.json').write_text(json.dumps({'controls': control_points,
        'focused_seed_consistency': consistency,
        'native_counts': {g: {'total': sum(r['target_id'] == g for r in native),
                             'pass': sum(r['target_id'] == g and r['native_filter_status'] == 'pass' for r in native)}
                          for g in sorted({r['target_id'] for r in native})},
        'source_reports': ['../analysis/analysis.json', '../focused_validation_report/analysis.json', '../positive_pose_audit.json']}, indent=2) + '\n')
    (output / 'README.md').write_text('''# Metrics summary figures

PNG figures are accompanied by vector PDF copies. All values come from saved outputs; no inference was run to make these figures.

1. **Target controls:** 2×2 frozen MSA / explicit empty MSA and template / no-template conditions. pLDDT (0–100) measures local model confidence; pTM (0–1) global model confidence. CA RMSD (Å) compares exact mapped residues with a separate experimental structure. Template-conditioned recovery is distinct from no-template folding. Historical PDBs are not guaranteed absent from training.
2. **Native filtering:** actual generated candidates and original native filter pass/fail. Native filters test structural recovery and composition, without an interface-confidence acceptance threshold. Native Boltz ipTM and ESMFold2 ipSAE are distinct metrics, not affinity values.
3. **Focused standard refolding:** target MSA, binder empty MSA, no coordinates/templates. Individual seeds show binder-alone and complex pLDDT, both raw chain-pair ipTM directions, binder-alone fold RMSD to design, binder pose RMSD to design after target alignment, experimental target RMSD, and seed-pair pose RMSD after target alignment. Lower design RMSD means agreement with a computational design, not experimental accuracy. Low seed-pair RMSD means repeatability, not proof that the pose is real. Candidate labels show the final six identifier characters; full IDs are in source analysis JSON.
4. **Experimental positive control:** 1WQJ IGF1–IGFBP4 domain with target-only template. Empty-binder-MSA and both-natural-chain paired+unpaired MSA arms are different protocols. Contact recall/precision use exact sequence positions at strict heavy-atom distance <5 Å. Interface RMSD jointly fits CAs of the reference-contact residues from both chains; this is not DockQ or full-backbone iRMSD. This calibrates a known natural interaction, not new designed binders.

Three requested seeds are not equilibrium samples or a significance test. The three seed pairs share predictions and are not independent. No figure applies a calibrated binding-success cutoff, and no experimentally validated designed binder is claimed. Severe clashes and full per-seed details remain in the source candidate report and focused review.

Sources: `../analysis/analysis.json`, `../focused_validation_report/analysis.json`, focused native bundles, `../positive_pose_audit.json`. Derived point data and seed-pair computations are in `figure_data.json`.
''')
    print(json.dumps({'output': str(output), 'controls': f"{controls['completed_predictions']}/{controls['expected_predictions']}",
                      'focused': f"{focused_data['analyzed_predictions']}/{focused_data['expected_predictions']}",
                      'seed_pairs': len(consistency)}))


if __name__ == '__main__':
    main()
