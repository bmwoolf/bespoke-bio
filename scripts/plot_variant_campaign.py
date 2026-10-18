#!/usr/bin/env python3
"""Static descriptive confidence and structure plots; no significance or affinity claims."""
import argparse
import json
import math
import os
from collections import defaultdict
from pathlib import Path


def plot(report_path, manifest_path, output):
    os.environ.setdefault('MPLCONFIGDIR', '/tmp/biocad-campaign-matplotlib')
    import matplotlib

    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np

    report = json.loads(Path(report_path).read_text())
    manifest = json.loads(Path(manifest_path).read_text()) if Path(manifest_path).is_file() else {}
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    made = []
    groups = defaultdict(lambda: defaultdict(dict))
    labels = {}
    for row in manifest.get('cross_scores', []):
        value = row['metrics'].get('design_to_target_iptm')
        if value is None or not math.isfinite(value):
            continue
        key = (row['target_id'], row['designed_on'], row['binder_sequence'])
        groups[key][row['scored_on']][row['seed']] = value
        labels[key] = row['binder_id']
    selections_path = Path(report_path).parent.parent / 'crossfold_tasks.json'
    selections = (json.loads(selections_path.read_text()).get('selected', [])
                  if selections_path.is_file() else [])
    original_status = {(r['target_id'], r['designed_on'], r['sequence']):
                       (r.get('filter') or {}).get('status', 'unknown') for r in selections}
    panels = []
    for key, states in groups.items():
        for variant in states:
            if variant != 'reference':
                panels.append((key, variant, states))
    if panels:
        columns = min(3, len(panels))
        rows = math.ceil(len(panels) / columns)
        fig, axes = plt.subplots(rows, columns, figsize=(4.6 * columns, 4.0 * rows),
                                 squeeze=False, sharey=True)
        seed_values = sorted({seed for states in groups.values() for scores in states.values()
                              for seed in scores})
        colors = {seed: plt.get_cmap('tab10')(i % 10) for i, seed in enumerate(seed_values)}
        for axis, (key, variant, states) in zip(axes.flat, panels):
            reference, alternate = states.get('reference', {}), states[variant]
            for seed in sorted(reference.keys() | alternate.keys()):
                xs, ys = [], []
                for x, scores in [(0, reference), (1, alternate)]:
                    if seed in scores:
                        xs.append(x)
                        ys.append(scores[seed])
                axis.plot(xs, ys, marker='o', color=colors[seed], alpha=.85, linewidth=1.3,
                          label=f'Seed {seed}')
            axis.set_xticks([0, 1], ['Reference', variant])
            axis.set_xlim(-.2, 1.2)
            axis.set_ylim(0, 1)
            axis.grid(axis='y', alpha=.2)
            axis.set_ylabel('Interface ipTM (model confidence)')
            status = original_status.get(key, 'unknown').upper()
            suffix = ' · exploratory' if status == 'FAIL' else ''
            axis.set_title(f'{key[0]} · {labels[key][-8:]}\nDesigned on {key[1]}\n'
                           f'Original native filter: {status}{suffix}', fontsize=10,
                           color='#a33a28' if status == 'FAIL' else '#172536')
        for axis in list(axes.flat)[len(panels):]:
            axis.set_visible(False)
        handles = [plt.Line2D([0], [0], color=colors[s], marker='o', label=f'Seed {s}')
                   for s in seed_values]
        fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.5, .96),
                   ncol=min(6, len(handles)), frameon=False)
        fig.suptitle('Same binder, reference versus variant target', fontsize=15, y=.995)
        fig.text(.5, .01, 'Each line joins the same seed. Confidence is not affinity. '
                 '\nFilter pass is not binding validation; few replicates do not establish significance.', ha='center', fontsize=9)
        fig.tight_layout(rect=[0, .035, 1, .91])
        path = output / 'fixed_binder_confidence.png'
        fig.savefig(path, dpi=170, facecolor='white')
        plt.close(fig)
        made.append(str(path))
    geometries = defaultdict(list)
    for row in report.get('target_structural_comparisons', []):
        kind = row.get('comparison_type', 'within_state' if row['state'] == 'reference'
                       else 'reference_vs_variant')
        geometries[(row['target_id'], row['state'], kind)].append(row['ca_rmsd_angstrom'])
    if geometries:
        fig, axis = plt.subplots(figsize=(11, max(4, .5 * len(geometries) + 1.8)))
        ticks = []
        for i, (key, values) in enumerate(sorted(geometries.items())):
            offsets = np.linspace(-.12, .12, len(values)) if len(values) > 1 else [0]
            color = '#267b98' if key[2] == 'within_state' else '#b26927'
            axis.scatter(values, i + np.asarray(offsets), color=color, alpha=.75, s=28)
            axis.plot([min(values), max(values)], [i, i], color=color, alpha=.35)
            ticks.append(f'{key[0]} {key[1]} · ' +
                         ('within state' if key[2] == 'within_state' else 'vs reference') +
                         f' (pairs={len(values)})')
        axis.set_yticks(range(len(ticks)), ticks)
        axis.invert_yaxis()
        axis.set_xlim(left=0)
        axis.set_xlabel('Target CA RMSD after proper rigid alignment (Å)')
        axis.grid(axis='x', alpha=.2)
        axis.set_title('Target-alone replicate structural variation', pad=15)
        fig.text(.5, .02, 'Dots are pairwise comparisons and are not independent observations. '
                 'Model samples are not an equilibrium ensemble.', ha='center', fontsize=9)
        fig.tight_layout(rect=[0, .055, 1, 1])
        path = output / 'target_structure_variation.png'
        fig.savefig(path, dpi=170, facecolor='white')
        plt.close(fig)
        made.append(str(path))
    return made


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report', type=Path)
    parser.add_argument('manifest', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    print(json.dumps(plot(args.report, args.manifest, args.output)))
