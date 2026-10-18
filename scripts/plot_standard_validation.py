#!/usr/bin/env python3
"""Plot complete standard folding controls; ranges are repeats, not confidence intervals."""
import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--analysis', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = json.loads(args.analysis.read_text())
    if report['completed_predictions'] != report['expected_predictions']:
        raise ValueError('Complete the expected control grid before plotting')
    conditions = ['empty_msa__no_template', 'server_msa__no_template',
                  'empty_msa__target_template', 'server_msa__target_template']
    labels = ['Sequence only', 'MSA only', 'Template only', 'MSA + template']
    colors = ['#9ba8b5', '#1b9380', '#9080b3', '#5a559c']
    fig, axes = plt.subplots(3, 2, figsize=(12, 9), constrained_layout=True)
    for row, gene in enumerate(['FGF2', 'IL1RN', 'SERPINA1']):
        entries = {s['condition']: s for s in report['summary'] if s['target_id'] == gene}
        confidence = [entries[c]['target_plddt'] for c in conditions]
        rmsd = [next(x for x in entries[c]['reference_rmsd_angstrom']
                     if x['role'] == 'held_out') for c in conditions]
        for column, (stats, title) in enumerate([
                (confidence, 'Fold confidence · higher is better'),
                (rmsd, 'Experimental structure error · lower is better')]):
            ax = axes[row, column]
            values = np.array([s['median'] for s in stats])
            error = np.array([[s['median'] - s['min'] for s in stats],
                              [s['max'] - s['median'] for s in stats]])
            bars = ax.bar(np.arange(4), values, color=colors, yerr=error, capsize=4)
            ax.bar_label(bars, labels=[f'{v:.1f}' for v in values], padding=5, fontsize=10)
            ax.set_xticks(range(4), labels, fontsize=9)
            ax.spines[['top', 'right']].set_visible(False)
            ax.set_title(f'{gene} — {title}', fontsize=11, loc='left')
            ax.set_ylim(0, 108 if column == 0 else max(values.max() * 1.25, 2.5))
            ax.set_ylabel('Target pLDDT / 100' if column == 0 else 'Matched Cα RMSD (Å)')
    fig.suptitle('MSA and template effects on target folding', fontsize=18)
    fig.supxlabel('Three repeats per condition; bars show medians, whiskers show ranges.\n'
                  'Templates condition the prediction. Historical test structures may be in training.', fontsize=10)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180)
    plt.close(fig)


if __name__ == '__main__':
    main()
