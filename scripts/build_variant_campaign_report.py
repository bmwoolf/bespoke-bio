#!/usr/bin/env python3
"""Build an incomplete-safe campaign analysis from downloaded outputs.

python scripts/build_variant_campaign_report.py outputs/variant_campaign \
  --refold outputs/variant_campaign/crossfold --output outputs/variant_campaign/report

Generation inputs: ROOT/generation_manifest.json; each cohort directory contains
live/bundle/bundle.json. Optional --refold directories contain receipt.json and
seed_N/refold_cif/TASK.cif from refold_fixed_sequences.py (remote absolute archive
paths in receipts are ignored). Task metadata should include target_id, scored_on,
and for complexes binder_id + designed_on. Unique target sequences may resolve
missing target_id/scored_on, but binder origin is never guessed.
Generated-complex targets are NOT mixed with target-alone structural ensembles.
Matched seed differences are descriptive predictions, not affinity/selectivity.
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

from analyze_variant_campaign import analyze, distribution, variant_contacts


def matched_seed_deltas(rows):
    groups = defaultdict(lambda: defaultdict(dict))
    for row in rows:
        key = (row["target_id"], row["binder_sequence"], row["designed_on"])
        slot = groups[key][row["scored_on"]]
        seed = row["seed"]
        if seed in slot:
            raise ValueError(f"Duplicate fixed-binder score for {key}, seed {seed}")
        slot[seed] = row["metrics"]
    output = []
    for (target, sequence, origin), states in groups.items():
        reference = states.get("reference", {})
        for state, scores in states.items():
            if state == "reference":
                continue
            pairs, delta_values = [], defaultdict(list)
            for seed in sorted(reference.keys() & scores.keys()):
                metrics = {}
                for name in reference[seed].keys() & scores[seed].keys():
                    a, b = reference[seed][name], scores[seed][name]
                    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
                        metrics[name] = b - a
                        delta_values[name].append(b - a)
                pairs.append({"seed": seed, "variant_minus_reference": metrics})
            output.append({"target_id": target, "binder_sequence": sequence,
                           "designed_on": origin, "scored_on": state,
                           "n_matched_seeds": len(pairs), "paired_deltas": pairs,
                           "delta_distributions": {k: distribution(v)
                                                   for k, v in delta_values.items()},
                           "unmatched_reference_seeds": sorted(reference.keys() - scores.keys()),
                           "unmatched_variant_seeds": sorted(scores.keys() - reference.keys())})
    return output


def collect_contacts(manifest):
    cohorts = {(c["target_id"], c["state"]): c for c in manifest["cohorts"]}
    results = []

    def append_contacts(target, state, binder_sequence, path, metadata):
        ref = cohorts.get((target, "reference"))
        if ref is None:
            return
        variants = [c for (gene, variant), c in cohorts.items()
                    if gene == target and variant != "reference"
                    and (state == "reference" or state == variant)]
        sequence = cohorts[(target, state)]["sequence"]
        for variant in variants:
            info = {**metadata, "target_id": target, "scored_on": state,
                    "compared_variant": variant["state"], "structure": str(path)}
            if len(ref["sequence"]) != len(variant["sequence"]):
                info.update(status="unavailable", reason="non_substitution_sequence_mapping")
            else:
                positions = [i for i, (a, b) in enumerate(
                    zip(ref["sequence"], variant["sequence"]), 1) if a != b]
                info.update(variant_contacts(path, sequence, binder_sequence, positions))
            results.append(info)

    for cohort in cohorts.values():
        if not cohort.get("bundle"):
            continue
        bundle_path = Path(cohort["bundle"])
        bundle = json.loads(bundle_path.read_text())
        for component in bundle.get("components", []):
            if component.get("role") != "computational_binder_candidate":
                continue
            path = component.get("evidence", {}).get("complex_structure")
            if not path:
                continue
            append_contacts(cohort["target_id"], cohort["state"], component["sequence"],
                            bundle_path.parent / path,
                            {"source": "generated_complex", "binder_id": component["id"],
                             "designed_on": cohort["state"]})
    for row in manifest["cross_scores"]:
        append_contacts(row["target_id"], row["scored_on"], row["binder_sequence"],
                        row["structure"], {"source": "fixed_binder_crossfold",
                        "binder_id": row["binder_id"], "seed": row["seed"],
                        "designed_on": row["designed_on"]})
    return results


def build(root, refold_dirs=()):
    root = Path(root).resolve()
    original = json.loads((root / "generation_manifest.json").read_text())
    manifest = {"cohorts": [], "cross_scores": []}
    availability, missing, target_scores = [], [], []
    by_key, by_sequence = {}, defaultdict(list)
    for row in original["cohorts"]:
        key = (row["target_id"], row["state"])
        if key in by_key:
            raise ValueError(f"Repeated generation cohort {key}")
        cohort = {"target_id": key[0], "state": key[1], "sequence": row["sequence"],
                  "structures": []}
        bundle_path = root / row["directory"] / "live/bundle/bundle.json"
        if bundle_path.is_file():
            cohort["bundle"] = str(bundle_path)
        availability.append({"target_id": key[0], "state": key[1],
                             "generation_status": row.get("status", "unknown"),
                             "bundle_available": bundle_path.is_file()})
        manifest["cohorts"].append(cohort)
        by_key[key] = cohort
        by_sequence[row["sequence"]].append(key)
    for supplied in refold_dirs:
        refold = Path(supplied).resolve()
        receipt_path = refold / "receipt.json"
        if not receipt_path.is_file():
            missing.append(str(receipt_path))
            continue
        receipt = json.loads(receipt_path.read_text())
        tasks = {t["id"]: t for t in receipt["tasks"]}
        result_keys = set()
        for result in receipt.get("results", []):
            task = tasks[result["id"]]
            possible = by_sequence.get(task["target_sequence"], [])
            if task.get("target_id") and task.get("scored_on"):
                key = (task["target_id"], task["scored_on"])
            elif len(possible) == 1:
                key = possible[0]
            else:
                raise ValueError(f"Cannot identify target state for {task['id']}")
            if key not in by_key or by_key[key]["sequence"] != task["target_sequence"]:
                raise ValueError(f"Target metadata/sequence mismatch for {task['id']}")
            seed = result["seed"]
            result_keys.add((task["id"], seed))
            cif = refold / f"seed_{seed}" / "refold_cif" / f"{task['id']}.cif"
            if not cif.is_file():
                missing.append(str(cif))
            if task.get("binder_sequence"):
                if not task.get("designed_on"):
                    raise ValueError(f"Complex {task['id']} needs explicit designed_on metadata")
                manifest["cross_scores"].append({"target_id": key[0], "scored_on": key[1],
                    "binder_id": task.get("binder_id", task["id"]),
                    "binder_sequence": task["binder_sequence"], "designed_on": task["designed_on"],
                    "seed": seed, "metrics": result["metrics"], "structure": str(cif)})
            else:
                target_scores.append({"target_id": key[0], "state": key[1], "seed": seed,
                                      "metrics": result["metrics"], "structure": str(cif)})
                if cif.is_file():
                    by_key[key]["structures"].append({"path": str(cif), "chain": "A"})
        for task_id in tasks:
            for seed in receipt.get("seeds", []):
                if (task_id, seed) not in result_keys:
                    missing.append(f"prediction:{refold.name}/{task_id}/seed_{seed}")
    report = analyze(manifest)
    report["variant_residue_contacts"] = collect_contacts(manifest)
    report["availability"] = availability
    report["missing_refold_outputs"] = missing
    report["target_alone_scores"] = target_scores
    report["fixed_binder_matched_seed_deltas"] = matched_seed_deltas(manifest["cross_scores"])
    report["complete_generation"] = all(a["bundle_available"] for a in availability)
    report["refold_evidence_supplied"] = bool(refold_dirs)
    report["limitations"] = [
        "Monomeric constructs may omit native oligomers, modifications and cofactors.",
        ("Native design evaluation uses target templates; these fixed-sequence comparisons do not. "
         "Absolute scores across those protocols are not directly interchangeable."),
        "Low-confidence target-only folds cannot establish biological conformational effects of a variant.",
        ("Generated cohorts contain independently sampled binders; cohort metric shifts are not "
         "fixed-binder variant effects."),
        "Matched seed cross-fold deltas are computational confidence differences, not affinity.",
        "Three seeds provide exploratory variability only; no statistical significance is claimed.",
    ]
    return manifest, report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--refold", type=Path, action="append", default=[])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    manifest, report = build(args.root, args.refold)
    args.output.mkdir(parents=True, exist_ok=True)
    for name, value in [("analysis_manifest", manifest), ("analysis_report", report)]:
        (args.output / f"{name}.json").write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"cohorts": len(report["cohorts"]),
                      "candidates": sum(c["n_candidates"] for c in report["cohorts"]),
                      "structural_comparisons": len(report["target_structural_comparisons"]),
                      "fixed_binder_comparisons": len(report["fixed_binder_matched_seed_deltas"]),
                      "missing_refold_outputs": len(report["missing_refold_outputs"])}))


if __name__ == "__main__":
    main()
