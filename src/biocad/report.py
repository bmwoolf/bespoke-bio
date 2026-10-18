"""Offline orthographic geometry review, no external viewer dependencies."""

import html
from pathlib import Path

import gemmi
import numpy as np

from .structures import atoms


def render(out, response):
    out = Path(out)
    parts = [
        '<!doctype html><meta charset="utf-8"><title>BioCAD geometry review</title><style>body{font:15px system-ui;margin:32px;max-width:1200px}svg{background:#f2f5f8;border:1px solid #ddd}table{border-collapse:collapse}td,th{padding:8px;border:1px solid #ddd}code{font-size:12px}</style>',
        "<h1>PD-L1 molecular geometry review</h1><p>Experimental coordinate replay and transferred hypotheses. No generated binder or affinity prediction. Projections show observed Cα atoms; use mmCIF for all atoms.</p>",
    ]
    colors = ["#1685aa", "#e07528", "#7f54c6", "#479950"]
    for file in [
        "component.cif",
        "transferred_complex.cif",
        "proposed.cif",
        "deliberate_overlap.cif",
    ]:
        st = gemmi.read_structure(str(out / file))
        chains = [(c.name, atoms(c, ca=True)[0]) for c in st[0]]
        allxyz = np.concatenate([p for _, p in chains])
        center = allxyz.mean(0)
        scale = 300 / max(np.ptp(allxyz, axis=0).max(), 1)
        parts.append(
            f'<h2>{file}</h2><a href="{file}">mmCIF coordinates</a><p>'
            + ", ".join(
                f'<span style="color:{colors[i]}">chain {html.escape(n)}</span>'
                for i, (n, _) in enumerate(chains)
            )
            + "</p>"
        )
        for axes, label in [([0, 1], "XY"), ([0, 2], "XZ"), ([1, 2], "YZ")]:
            parts.append(
                f'<svg width="360" height="360" viewBox="0 0 360 360"><text x="12" y="20">{label}, Å</text>'
            )
            for i, (name, xyz) in enumerate(chains):
                projected = (xyz - center)[:, axes] * scale + 180
                for x, y in projected:
                    parts.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="2" fill="{colors[i]}"/>')
            parts.append("</svg>")
    parts.append(
        "<h2>Evaluations</h2><table><tr><th>Subject</th><th>Method</th><th>Status</th><th>Metrics</th></tr>"
    )
    for e in response.evaluations:
        parts.append(
            "<tr>"
            + "".join(
                f"<td>{html.escape(str(v))}</td>"
                for v in [e.subject_id, e.method, e.status, e.metrics]
            )
            + "</tr>"
        )
    parts.append(
        "</table><p>Proposed fusion sequences omit unresolved native terminal residues. A linker reach bound is necessary only; the connected construct is unmodeled. “Assessed” does not mean experimentally functional.</p>"
    )
    (out / "review.html").write_text("\n".join(parts))
