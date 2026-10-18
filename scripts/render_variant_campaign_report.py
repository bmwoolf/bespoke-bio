#!/usr/bin/env python3
"""Render campaign analysis to report/summary.md and report/review.html.

Usage: python scripts/render_variant_campaign_report.py outputs/variant_campaign
Reads report/analysis_report.json, generation_manifest.json and optional
crossfold_tasks.json. Writes only ROOT/report. No scientific inference or cloud use.
"""
import argparse
import html
import importlib.util
import json
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from urllib.parse import quote


def stats(metric):
    if not metric or not metric.get("n"):
        return "—"
    return (f"{metric['mean']:.3f} ± {metric['std']:.3f} "
            f"[{metric['min']:.3f}, {metric['max']:.3f}]; n={metric['n']}")


def render(root):
    root = Path(root).resolve()
    output = root / "report"
    report = json.loads((output / "analysis_report.json").read_text())
    generation = json.loads((root / "generation_manifest.json").read_text())
    task_path = root / "crossfold_tasks.json"
    tasks = json.loads(task_path.read_text()) if task_path.is_file() else {}
    identities = {(t["target_id"], t["designed_on"], t["binder_sequence"]): t
                  for t in tasks.get("tasks", []) if t.get("binder_sequence")}
    selections = {(t['target_id'], t['designed_on'], t['sequence']):
                  (t.get('filter') or {}).get('status', 'unknown')
                  for t in tasks.get('selected', [])}
    for key, identity in identities.items():
        status = selections.get(key, 'unknown')
        identity['original_filter_status'] = (status.upper() +
            (' — exploratory failed candidate' if status == 'fail' else ''))
    md = ["# Variant binder campaign", "",
          ("Computational predictions only. Scores describe model confidence, not binding affinity "
          "or variant selectivity. Replicate folds sample model outputs, not equilibrium ensembles."), "",
          ("Numbers show mean ± population standard deviation [minimum, maximum]. "
          "Small sample counts do not support significance claims."), ""]
    web = ['<!doctype html><meta charset="utf-8"><title>Variant binder campaign</title>',
           ('<style>body{font:15px system-ui;margin:32px;color:#172536}table{border-collapse:collapse;'
           'width:100%;margin:20px 0}th,td{padding:9px;text-align:left;border-bottom:1px solid #ddd;'
           'vertical-align:top}th{background:#eef3f7}a{color:#1263a3}h2{margin-top:35px}'
           '.wrap{overflow:auto}p{max-width:1050px}</style>'),
           '<h1>Variant binder campaign</h1>']
    web.extend(f"<p>{html.escape(p)}</p>" for p in (md[2], md[4]))

    seed_note = ("Generation seeds were not controlled. Only comparison runs set Lightning/Torch "
                 "seeds 11, 22 and 33; native feature RNG remains nondeterministic. Matched-seed "
                 "comparisons are bookkeeping matches, not identical stochastic draws.")
    md.extend([seed_note, ""])
    web.append(f"<p>{html.escape(seed_note)}</p>")

    if report.get('fixed_binder_cross_scoring') or report.get('target_structural_comparisons'):
        runtime = (Path(sys.executable) if importlib.util.find_spec('matplotlib') else
                   Path(__file__).resolve().parents[1] / 'vendor/boltzgen-env/bin/python')
        if runtime.is_file():
            result = subprocess.run([
                str(runtime), str(Path(__file__).with_name('plot_variant_campaign.py')),
                str(output / 'analysis_report.json'), str(output / 'analysis_manifest.json'),
                str(output)], check=True, capture_output=True, text=True)
            for image_path in json.loads(result.stdout):
                relative = quote(os.path.relpath(image_path, output), safe='/')
                caption = ('Fixed-binder confidence, paired by seed; not affinity or significance'
                           if 'confidence' in image_path else
                           'Target-only structural variation; pairs are not independent samples')
                md.extend([f'![{caption}]({relative})', ''])
                web.append(f'<figure><img style="max-width:100%" src="{relative}" '
                           f'alt="{html.escape(caption)}"><figcaption>{html.escape(caption)}'
                           '</figcaption></figure>')

    def table(title, headers, rows):
        md.extend([f"## {title}", "", "| " + " | ".join(headers) + " |",
                   "| " + " | ".join("---" for _ in headers) + " |"])
        web.append(f'<h2>{html.escape(title)}</h2><div class="wrap"><table><tr>' +
                   ''.join(f'<th>{html.escape(h)}</th>' for h in headers) + '</tr>')
        for row in rows:
            strings = [str(v) for v in row]
            md.append("| " + " | ".join(v.replace('|', '\\|') for v in strings) + " |")
            web.append('<tr>' + ''.join(f'<td>{html.escape(v)}</td>' for v in strings) + '</tr>')
        if not rows:
            md.append("| " + " | ".join(["No completed data"] + [""] * (len(headers)-1)) + " |")
        web.append('</table></div>')
        md.append("")

    panel_path = root / 'research/panel.json'
    if panel_path.is_file():
        panel = json.loads(panel_path.read_text())
        short_caveats = {
            'GDF15': 'Native disulfide-linked dimer omitted. Literature construct starts 197; '
                    'UniProt mature annotation starts 195 (H8D). Assay recognition is not bioactivity.',
            'F7': 'Zymogen; activation cleavage, calcium, Gla modifications, glycans and tissue '
                  'factor context omitted. MANE transcript uses different numbering.',
            'APOE': 'Apo monomer omits lipid particles and oligomers; domain orientation is flexible. '
                    'Injury-related biology does not imply an intervention benefit.',
            'VEGFA': 'VEGF165 isoform; native dimer, glycans and receptor/heparin context omitted. '
                     'Reference only; common regulatory variants are not coding substitutions.',
        }
        rows = []
        for target in panel['targets']:
            construct = target['construct']
            variants = '; '.join(v['precursor_change'] + ' → ' + v['construct_change'] +
                                 ' (' + v.get('variant_id', 'no rsID') + ')'
                                 for v in target.get('variants', [])) or 'Reference only'
            rows.append([target['gene_symbol'] + ' / ' + target['accession'], construct['length'],
                         f"{construct['start']}–{construct['end']}", variants,
                         short_caveats.get(target['gene_symbol'], '; '.join(target.get('caveats', [])))])
        table('Target panel', ['Protein / accession', 'Construct aa', 'Precursor span',
                              'Substitution: precursor → modeled chain (rsID)',
                              'Biological context and omissions'], rows)
        notes = root / 'research/SOURCE_NOTES.md'
        if notes.is_file():
            url = quote(os.path.relpath(notes, output), safe='/')
            md.extend([f'[Source notes and numbering evidence]({url})', ''])
            web.append(f'<p><a href="{url}">Source notes and numbering evidence</a></p>')

    table("Generated cohorts", ["Target", "State", "Candidates / unique", "Pass / fail / unknown",
                                "Interface ipTM", "ESMFold2 ipSAE min", "Recovery RMSD (Å)"], [
        [c['target_id'], c['state'], f"{c['n_candidates']} / {c['n_unique_sequences']}",
         f"{c['n_pass']} / {c['n_fail']} / {c['n_unknown_filter']}",
         stats(c['metrics'].get('design_to_target_iptm')),
         stats(c['metrics'].get('esmfold2_ipsae_min')), stats(c['metrics'].get('filter_rmsd'))]
        for c in report['cohorts']])
    from analyze_variant_campaign import distribution

    target_confidence = defaultdict(list)
    for row in report.get("target_alone_scores", []):
        value = row["metrics"].get("ptm")
        if isinstance(value, (int, float)):
            target_confidence[(row["target_id"], row["state"])].append(value)
    table("Target-alone confidence (pTM)", ["Target", "State", "pTM"],
          [[*key, stats(distribution(values))]
           for key, values in sorted(target_confidence.items())])

    geometries = defaultdict(list)
    for r in report.get('target_structural_comparisons', []):
        kind = r.get('comparison_type', 'within_state' if r['state'] == 'reference'
                     else 'reference_vs_variant')
        geometries[(r['target_id'], r['state'], kind)].append(r['ca_rmsd_angstrom'])
    table("Target-alone structure variation", ["Target", "State", "Comparison", "CA RMSD (Å)"],
          [[*key, stats(distribution(values))] for key, values in sorted(geometries.items())])
    fixed = []
    for row in report.get('fixed_binder_cross_scoring', []):
        identity = identities.get((row['target_id'], row['designed_on'], row['binder_sequence']), {})
        for state, metrics in row['scores_by_target_state'].items():
            fixed.append([row['target_id'], identity.get('binder_id', 'unmapped'), row['designed_on'],
                          identity.get('selection_reason', 'unspecified'),
                          identity.get('original_filter_status', 'UNKNOWN'), state,
                          stats(metrics.get('design_to_target_iptm'))])
    table("Fixed-binder cross-fold scores", ["Target", "Binder", "Designed on", "Selection",
                                             "Original native filter", "Scored on", "Interface ipTM"], fixed)
    deltas = []
    for row in report.get('fixed_binder_matched_seed_deltas', []):
        identity = identities.get((row['target_id'], row['designed_on'], row['binder_sequence']), {})
        deltas.append([row['target_id'], identity.get('binder_id', 'unmapped'), row['designed_on'],
                       identity.get('original_filter_status', 'UNKNOWN'),
                       row['scored_on'], row['n_matched_seeds'],
                       stats(row['delta_distributions'].get('design_to_target_iptm'))])
    table("Matched-seed differences: variant minus reference", ["Target", "Binder", "Designed on",
          "Original native filter", "Variant", "Paired seeds", "Interface ipTM difference"], deltas)
    contacts = []
    for row in report.get('variant_residue_contacts', []):
        residues = row.get('residues', [])
        detail = '; '.join(f"{r['residue']}{r['sequence_position']}: " +
                          (f"{r['min_heavy_atom_distance_angstrom']:.2f} Å; "
                           f"within 5 Å: {r['within_cutoff']}" if
                           r['min_heavy_atom_distance_angstrom'] is not None else 'unavailable')
                          for r in residues) or row.get('reason', 'no substitutions')
        contacts.append([row['target_id'], row['source'], row['binder_id'], row['scored_on'],
                         row['compared_variant'], row.get('seed', '—'), detail])
    table("Mutation-to-binder contacts (arbitrary geometric cutoff)", ["Target", "Source", "Binder",
          "Scored on", "Compared variant", "Seed", "Minimum heavy-atom distance"], contacts)
    table("Completeness", ["Target", "State", "Generation status", "Bundle downloaded"],
          [[r['target_id'], r['state'], r['generation_status'], r['bundle_available']]
           for r in report.get('availability', [])])
    missing = report.get('missing_refold_outputs', [])
    table("Missing refold outputs", ["Output"], [[p] for p in missing])
    files = {root / name for name in ("binders.fasta", "candidates.csv",
                                     "download_verification.json", "compute_cleanup.json")
             if (root / name).is_file()}
    for c in generation['cohorts']:
        directory = root / c['directory'] / 'live'
        for pattern in ('bundle/**/*.cif', 'bundle/bundle.json', '**/*metrics*.csv'):
            files.update(p for p in directory.glob(pattern) if p.is_file())
    for row in report.get('target_alone_scores', []) + report.get('variant_residue_contacts', []):
        if row.get('structure') and Path(row['structure']).is_file():
            files.add(Path(row['structure']))
    md.extend(['## Downloaded structures and native files', ''])
    web.append('<h2>Downloaded structures and native files</h2><ul>')
    for path in sorted(files):
        relative = os.path.relpath(path, output)
        url = quote(relative, safe='/')
        md.append(f'- [{relative}]({url})')
        web.append(f'<li><a href="{html.escape(url, quote=True)}">{html.escape(relative)}</a></li>')
    web.append('</ul><h2>Limitations</h2><ul>')
    md.extend(['', '## Limitations', ''])
    for limitation in report.get('limitations', []):
        md.append('- ' + limitation)
        web.append('<li>' + html.escape(limitation) + '</li>')
    web.append('</ul>')
    (output / 'summary.md').write_text('\n'.join(md) + '\n')
    (output / 'review.html').write_text('\n'.join(web) + '\n')
    return {'markdown': str(output / 'summary.md'), 'html': str(output / 'review.html'),
            'linked_files': len(files)}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    print(json.dumps(render(parser.parse_args().root)))
