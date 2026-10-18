"""Check downloaded RFD3 coordinates locally; this is not a folding/binding score."""

import argparse
import json
from itertools import pairwise
from pathlib import Path

import gemmi

from biocad.motifs import check_motif_geometry
from biocad.structures import protein_residues

parser = argparse.ArgumentParser()
parser.add_argument("candidate")
parser.add_argument("topology")
parser.add_argument("--plan", default="examples/gpu_smoke")
parser.add_argument("--output", required=True)
args = parser.parse_args()
plan = json.loads((Path(args.plan) / "motif_plan.json").read_text())
topology = next(t for t in plan["topologies"] if t["id"] == args.topology)
report = check_motif_geometry(args.plan, args.candidate, args.topology)
structure = gemmi.read_structure(args.candidate)
residues = protein_residues(structure[0]["A"])
report["scaffold_residues"] = len(residues)
report["expected_scaffold_residues"] = topology["scaffold_length"]
report["backbone_gaps"] = []
for first, second in pairwise(residues):
    carbon, nitrogen = first.find_atom("C", "*"), second.find_atom("N", "*")
    distance = None if carbon is None or nitrogen is None else carbon.pos.dist(nitrogen.pos)
    if distance is None or not 1.1 <= distance <= 1.8:
        report["backbone_gaps"].append({"after": str(first.seqid), "distance_angstrom": distance})
report["backbone_screen"] = (
    "passed"
    if len(residues) == topology["scaffold_length"] and not report["backbone_gaps"]
    else "failed"
)
report["limits"] = (
    "Geometry triage only; sequence design, independent refolding, packing, disulfide and functional validation remain required."
)
Path(args.output).write_text(json.dumps(report, indent=2) + "\n")
