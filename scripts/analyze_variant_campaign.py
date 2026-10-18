#!/usr/bin/env python3
"""Summarize computational cohorts without interpreting confidence as affinity.

Usage: python scripts/analyze_variant_campaign.py manifest.json --output report.json
Paths resolve relative to the manifest. Minimal manifest:
{"cohorts": [{"target_id":"X", "state":"reference", "sequence":"ACD...",
              "bundle":"run/bundle.json", "structures":[
                 {"path":"target.cif", "chain":"B"}]}],
 "cross_scores": [{"target_id":"X", "binder_id":"b1", "binder_sequence":"...",
                   "designed_on":"reference", "scored_on":"variant",
                   "metrics":{"iptm":0.5}}]}

Use multiple cohort entries for replicate runs. A structure's residue numbering must
be 1-based positions in its cohort sequence, or provide residue_positions, a list of
sequence positions in polymer-residue order (including residues lacking CA). These
positions are verified against residue identity. Only equal-length substitution
variants are supported for structural comparison. Precursor/mature/domain numbering
must be explicitly mapped by the caller. Structures are optional; ensemble RMSDs are
all pairwise comparisons to reference structures, including within-reference pairs.
Cross scores must refer to the SAME fixed binder sequence to compute reference-versus-
variant deltas; generated cohorts alone cannot establish variant selectivity.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np


def proper_rmsd(mobile, reference):
    """Least-squares CA RMSD using a proper rotation (never a reflection)."""
    x, y = np.asarray(mobile, dtype=float), np.asarray(reference, dtype=float)
    if x.shape != y.shape or x.ndim != 2 or x.shape[1] != 3 or len(x) < 3:
        raise ValueError("RMSD needs matching N x 3 arrays with at least three positions")
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError("Nonfinite coordinates")
    x, y = x - x.mean(axis=0), y - y.mean(axis=0)
    u, _, vt = np.linalg.svd(x.T @ y)
    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(u @ vt)
    rotation = u @ correction @ vt
    return float(np.sqrt(np.mean(np.sum((x @ rotation - y) ** 2, axis=1))))


def distribution(values):
    values = [float(v) for v in values if isinstance(v, (int, float))
              and not isinstance(v, bool) and math.isfinite(v)]
    if not values:
        return {"n": 0}
    return {"n": len(values), "min": min(values), "max": max(values),
            "mean": float(np.mean(values)), "median": float(np.median(values)),
            "std": float(np.std(values)), "values": values}


def candidate_rows(bundle):
    rows = []
    for component in bundle.get("components", []):
        if component.get("role") != "computational_binder_candidate":
            continue
        evaluations = [e for e in bundle.get("evaluations", [])
                       if e.get("subject_id") == component["id"]]
        metrics, filter_status = {}, None
        for evaluation in evaluations:
            metrics.update(evaluation.get("metrics", {}))
            if "filter" in evaluation.get("method", "").lower():
                filter_status = evaluation.get("status")
        rows.append({"id": component["id"], "sequence": component.get("sequence", ""),
                     "filter_status": filter_status, "metrics": metrics})
    return rows


def read_ca(spec, sequence, root):
    import gemmi

    structure = gemmi.read_structure(str(root / spec["path"]))
    chain = structure[0][spec["chain"]]
    residues = [r for r in chain if gemmi.find_tabulated_residue(r.name).is_amino_acid()]
    positions = spec.get("residue_positions")
    if positions is None:
        if any(r.seqid.icode.strip() for r in residues):
            raise ValueError("Insertion codes require explicit residue_positions")
        positions = [r.seqid.num for r in residues]
    if len(positions) != len(residues) or len(set(positions)) != len(positions):
        raise ValueError("residue_positions must uniquely map every polymer residue")
    result = {}
    for position, residue in zip(positions, residues):
        if not isinstance(position, int) or not 1 <= position <= len(sequence):
            raise ValueError("Structure residue position outside supplied sequence")
        aa = gemmi.find_tabulated_residue(residue.name).one_letter_code.upper()
        if aa != sequence[position - 1]:
            raise ValueError(f"Residue identity mismatch at sequence position {position}")
        atoms = [a for a in residue if a.name == "CA"]
        if atoms:
            atom = max(atoms, key=lambda a: a.occ)
            result[position] = [atom.pos.x, atom.pos.y, atom.pos.z]
    return result


def variant_contacts(path, target_sequence, binder_sequence, positions, cutoff=5.0):
    """Descriptive mutation-to-binder geometry; identify chains by exact sequences."""
    import gemmi

    if not Path(path).is_file():
        return {"status": "unavailable", "reason": "complex_structure_missing"}
    structure = gemmi.read_structure(str(path))
    chains = []
    for chain in structure[0]:
        residues = [r for r in chain if gemmi.find_tabulated_residue(r.name).is_amino_acid()]
        sequence = "".join(gemmi.find_tabulated_residue(r.name).one_letter_code.upper()
                           for r in residues)
        chains.append((chain.name, sequence, residues))
    targets = [c for c in chains if c[1] == target_sequence]
    binders = [c for c in chains if c[1] == binder_sequence]
    if len(targets) != 1 or len(binders) != 1 or targets[0] is binders[0]:
        return {"status": "unavailable", "reason": "exact_sequence_chain_mapping_ambiguous_or_missing",
                "target_matches": len(targets), "binder_matches": len(binders)}

    def heavy(residues):
        coords = []
        for residue in residues:
            atoms = {}
            for atom in residue:
                if atom.element.atomic_number <= 1 or atom.occ <= 0:
                    continue
                if atom.name not in atoms or atom.occ > atoms[atom.name].occ:
                    atoms[atom.name] = atom
            coords.extend([[a.pos.x, a.pos.y, a.pos.z] for a in atoms.values()])
        return np.asarray(coords, dtype=float).reshape(-1, 3)

    binder_atoms = heavy(binders[0][2])
    rows = []
    for position in positions:
        if not isinstance(position, int) or not 1 <= position <= len(target_sequence):
            raise ValueError("Mutation position outside target sequence")
        residue = targets[0][2][position - 1]
        atoms = heavy([residue])
        distance = None
        if len(atoms) and len(binder_atoms):
            distance = float(np.linalg.norm(atoms[:, None] - binder_atoms[None, :], axis=-1).min())
        rows.append({"sequence_position": position, "residue": target_sequence[position - 1],
                     "structure_residue_id": str(residue.seqid),
                     "min_heavy_atom_distance_angstrom": distance,
                     "within_cutoff": None if distance is None else distance <= cutoff})
    return {"status": "assessed", "target_chain": targets[0][0], "binder_chain": binders[0][0],
            "cutoff_angstrom": cutoff, "cutoff_scope": "arbitrary descriptive geometric cutoff",
            "alternate_atoms": "highest occupancy per atom name; first on ties",
            "residues": rows}


def exact_sequence_ca(path, sequence):
    """Map CA coordinates by ordered sequence positions after exact chain matching."""
    import gemmi

    structure = gemmi.read_structure(str(path))
    matches = []
    for chain in structure[0]:
        residues = [r for r in chain if gemmi.find_tabulated_residue(r.name).is_amino_acid()]
        observed = "".join(gemmi.find_tabulated_residue(r.name).one_letter_code.upper()
                           for r in residues)
        if observed == sequence:
            matches.append(residues)
    if len(matches) != 1:
        raise ValueError("Exact sequence chain mapping missing or ambiguous")
    result = {}
    for index, residue in enumerate(matches[0], 1):
        atom = residue.find_atom("CA", "*")
        if atom is not None:
            result[index] = [atom.pos.x, atom.pos.y, atom.pos.z]
    return result


def native_fold_confidence(path, target_sequence, binder_sequence=None, positions=()):
    """Read native FoldingWriter CA B-factors as normalized pLDDT, never crystallographic B."""
    import gemmi

    if not Path(path).is_file():
        return {"status": "unavailable", "reason": "structure_missing"}
    st = gemmi.read_structure(str(path))
    chains = []
    for chain in st[0]:
        residues = [r for r in chain if gemmi.find_tabulated_residue(r.name).is_amino_acid()]
        seq = "".join(gemmi.find_tabulated_residue(r.name).one_letter_code.upper()
                      for r in residues)
        chains.append((chain.name, seq, residues))
    summaries = {}
    mutated = []
    for role, sequence in [("target", target_sequence), ("binder", binder_sequence)]:
        if sequence is None:
            continue
        matched = [c for c in chains if c[1] == sequence]
        if len(matched) != 1:
            return {"status": "unavailable", "reason": "exact_sequence_chain_mapping_ambiguous_or_missing"}
        name, _, residues = matched[0]
        values = []
        for index, residue in enumerate(residues, 1):
            atom = residue.find_atom("CA", "*")
            value = None if atom is None else float(atom.b_iso)
            if value is not None and (not math.isfinite(value) or not 0 <= value <= 1):
                return {"status": "unavailable", "reason": "B_factor_not_normalized_native_confidence"}
            if value is not None:
                values.append(value)
            if role == "target" and index in positions:
                mutated.append({"sequence_position": index, "residue": sequence[index - 1],
                                "normalized_plddt": value})
        stats = distribution(values)
        stats.pop("values", None)
        summaries[role] = {"chain": name, **stats}
    return {"status": "assessed", "source": "native FoldingWriter CA B-factor; normalized 0..1",
            "chains": summaries, "variant_positions": mutated,
            "interpretation": "Local prediction confidence, not interface confidence or affinity"}


def analyze(manifest, root=Path(".")):
    root = Path(root)
    grouped, sequences, structures = defaultdict(list), {}, defaultdict(list)
    for cohort in manifest.get("cohorts", []):
        key = (cohort["target_id"], cohort["state"])
        sequence = cohort["sequence"].upper()
        if key in sequences and sequences[key] != sequence:
            raise ValueError(f"Conflicting sequences for {key}")
        sequences[key] = sequence
        grouped[key]  # retain empty/structure-only cohorts
        if cohort.get("bundle"):
            bundle = json.loads((root / cohort["bundle"]).read_text())
            grouped[key].extend(candidate_rows(bundle))
        for spec in cohort.get("structures", []):
            structures[key].append((spec["path"], read_ca(spec, sequence, root)))
    summary, geometry = [], []
    for key, rows in sorted(grouped.items()):
        target, state = key
        ref = sequences.get((target, "reference"))
        seq = sequences[key]
        mutations = None
        if ref is not None and len(ref) == len(seq):
            mutations = [{"position": i, "reference": a, "variant": b}
                         for i, (a, b) in enumerate(zip(ref, seq), 1) if a != b]
        values = defaultdict(list)
        for row in rows:
            for metric, value in row["metrics"].items():
                values[metric].append(value)
        unique = {r["sequence"] for r in rows if r["sequence"]}
        summary.append({"target_id": target, "state": state, "n_candidates": len(rows),
                        "n_unique_sequences": len(unique), "sequences": sorted(unique),
                        "n_pass": sum(r["filter_status"] == "pass" for r in rows),
                        "n_fail": sum(r["filter_status"] == "fail" for r in rows),
                        "n_unknown_filter": sum(r["filter_status"] not in ("pass", "fail")
                                                for r in rows),
                        "mutations": mutations, "metrics": {k: distribution(v)
                                                             for k, v in values.items()},
                        "candidates": rows})
        if ref is None or len(ref) != len(seq):
            continue
        refs = structures.get((target, "reference"), [])
        pairs = [("within_state", left, right)
                 for left, right in itertools.combinations(structures[key], 2)]
        if state != "reference":
            pairs.extend(("reference_vs_variant", left, right)
                         for left, right in itertools.product(refs, structures[key]))
        for comparison_type, (ref_path, ref_ca), (path, ca) in pairs:
            positions = sorted(ref_ca.keys() & ca.keys())
            if len(positions) < 3:
                raise ValueError("Insufficient matched CA positions for structural comparison")
            geometry.append({"target_id": target, "state": state,
                             "comparison_type": comparison_type, "reference": ref_path,
                             "structure": path, "matched_positions": positions,
                             "n_matched": len(positions),
                             "ca_rmsd_angstrom": proper_rmsd([ca[p] for p in positions],
                                                             [ref_ca[p] for p in positions])})
    fixed = defaultdict(lambda: defaultdict(list))
    for row in manifest.get("cross_scores", []):
        if not row.get("binder_sequence"):
            raise ValueError("Fixed-binder cross scoring requires binder_sequence")
        key = (row["target_id"], row["binder_sequence"], row["designed_on"])
        fixed[key][row["scored_on"]].append(row.get("metrics", {}))
    cross = []
    for (target, binder_sequence, designed_on), states in fixed.items():
        summaries = {}
        for state, metrics in states.items():
            names = set().union(*(m.keys() for m in metrics))
            summaries[state] = {k: distribution([m.get(k) for m in metrics]) for k in names}
        deltas = {}
        for state, metrics in summaries.items():
            if state == "reference":
                continue
            deltas[state] = {k: v["mean"] - summaries["reference"][k]["mean"]
                             for k, v in metrics.items() if v["n"]
                             and summaries.get("reference", {}).get(k, {}).get("n")}
        cross.append({"target_id": target, "binder_sequence": binder_sequence,
                      "designed_on": designed_on, "scores_by_target_state": summaries,
                      "mean_metric_delta_vs_reference": deltas})
    return {"cohorts": summary, "target_structural_comparisons": geometry,
            "fixed_binder_cross_scoring": cross,
            "interpretation": "Confidence and recovery metrics are not affinity. Generated cohort "
            "differences do not establish selectivity. Fixed-binder cross scores are computational "
            "hypotheses requiring experimental binding and selectivity validation."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = analyze(json.loads(args.manifest.read_text()), args.manifest.parent)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
