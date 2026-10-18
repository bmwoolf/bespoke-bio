#!/usr/bin/env python3
"""Analyze standard Boltz controls without inference or changing native artifacts."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import statistics
from collections import defaultdict
from pathlib import Path

import gemmi
import numpy as np

from biocad.geometry import fit


def chain_data(path):
    structure = gemmi.read_structure(str(path))
    if len(structure) != 1:
        raise ValueError("Expected exactly one prediction model")
    chains = {}
    for chain in structure[0]:
        residues = []
        for residue in chain:
            info = gemmi.find_tabulated_residue(residue.name)
            if not info.is_amino_acid():
                continue
            ca = residue.find_atom("CA", "*")
            if ca is None:
                raise ValueError(f"Missing CA in chain {chain.name} residue {residue.seqid}")
            if residue.seqid.icode.strip():
                raise ValueError("Expected input-position numbering without insertion codes")
            residues.append(
                (
                    residue.seqid.num,
                    info.one_letter_code,
                    np.array([ca.pos.x, ca.pos.y, ca.pos.z], dtype=float),
                    ca.b_iso,
                )
            )
        if not residues:
            continue
        positions = [r[0] for r in residues]
        if len(set(positions)) != len(positions):
            raise ValueError("Duplicate CA input positions")
        if chain.name in chains:
            raise ValueError("Duplicate chain identifier")
        chains[chain.name] = {
            "sequence": "".join(r[1] for r in residues),
            "positions": positions,
            "ca": {r[0]: r[2] for r in residues},
            "letters": {r[0]: r[1] for r in residues},
            "bfactor": np.array([r[3] for r in residues]),
        }
    return chains


def distribution(values):
    values = [float(v) for v in values if v is not None and np.isfinite(v)]
    return {
        "n": len(values),
        "mean": statistics.mean(values) if values else None,
        "median": statistics.median(values) if values else None,
        "min": min(values) if values else None,
        "max": max(values) if values else None,
    }


def compare_reference(predicted, reference, sequence, base):
    path = Path(reference["matched_structure"])
    if not path.is_absolute():
        path = Path(base) / path
    result = {
        "role": reference["role"],
        "pdb_id": reference.get("pdb_id"),
        "matched_structure": str(path),
        "rmsd_angstrom": None,
    }
    if not path.is_file():
        return {**result, "status": "missing_reference"}
    expected = [m["input_position"] for m in reference["mapping"]]
    if len(set(expected)) != len(expected):
        raise ValueError("Duplicate reference mapping positions")
    ref = chain_data(path)[reference.get("matched_chain", "A")]
    if ref["positions"] != expected:
        raise ValueError("Reference coordinates differ from audited mapping")
    for position in expected:
        if not 1 <= position <= len(sequence) or ref["letters"][position] != sequence[position - 1]:
            raise ValueError("Reference sequence differs at an audited position")
    common = [p for p in expected if p in predicted["ca"]]
    result.update(
        mapped_positions=len(expected),
        compared_positions=len(common),
        excluded_prediction_positions=sorted(set(expected) - set(common)),
        input_coverage=len(common) / len(sequence),
        compared_input_positions=common,
    )
    if len(common) < 3:
        return {**result, "status": "insufficient_coordinates"}
    # One proper rigid fit; no residue correspondence inferred from geometric proximity.
    _, rmsd = fit([predicted["ca"][p] for p in common], [ref["ca"][p] for p in common])
    return {
        **result,
        "status": "compared",
        "rmsd_angstrom": rmsd,
        "scope": "template-conditioned recovery"
        if reference["role"] == "template"
        else "held out from supplied inputs; training independence unknown",
    }


def parse_prediction(directory, case, target, base):
    directory = Path(directory)
    name = case["id"] + "_model_0"
    files = {
        "cif": directory / (name + ".cif"),
        "confidence": directory / ("confidence_" + name + ".json"),
        "plddt": directory / ("plddt_" + name + ".npz"),
    }
    result = {
        "status": "completed",
        "artifacts": {k: str(v) for k, v in files.items()},
        "missing_artifacts": [k for k, p in files.items() if not p.is_file()],
        "chains": {},
        "reference_comparisons": [],
        "warnings": [],
    }
    if result["missing_artifacts"]:
        result["status"] = "incomplete_results"
    if not files["cif"].is_file():
        result["status"] = "missing_results"
        return result
    try:
        chains = chain_data(files["cif"])
        expected = {"A": case["sequence"]}
        if case.get("binder_sequence"):
            expected["B"] = case["binder_sequence"]
        if list(chains) != list(expected):
            raise ValueError("Prediction chain set/order does not match the input")
        for name, sequence in expected.items():
            if chains[name]["sequence"] != sequence:
                raise ValueError(f"Exact observed chain sequence mismatch: {name}")
            if chains[name]["positions"] != list(range(1, len(sequence) + 1)):
                raise ValueError(
                    f"Prediction chain is not numbered by exact input position: {name}"
                )

        def reject_nonfinite(value):
            raise ValueError(f"Nonfinite native confidence value: {value}")

        confidence = (
            json.loads(files["confidence"].read_text(), parse_constant=reject_nonfinite)
            if files["confidence"].is_file()
            else {}
        )
        if not isinstance(confidence, dict):
            raise TypeError("Native confidence JSON must be an object")
        result["confidence"] = confidence
        arrays = None
        if files["plddt"].is_file():
            with np.load(files["plddt"], allow_pickle=False) as data:
                arrays = np.asarray(data["plddt"], dtype=float).reshape(-1)
            if len(arrays) != sum(len(c["sequence"]) for c in chains.values()):
                raise ValueError("Native pLDDT token count differs from protein residue count")
            if not np.all(np.isfinite(arrays)) or np.any(arrays < 0) or np.any(arrays > 1.00001):
                raise ValueError("Native Boltz pLDDT NPZ must use the finite 0..1 scale")
        offset = 0
        for index, (name, chain) in enumerate(chains.items()):
            length = len(chain["sequence"])
            values = arrays[offset : offset + length] * 100 if arrays is not None else None
            offset += length
            # B factors are 100 if writer had no confidence; never silently label them pLDDT.
            source = (
                "native token pLDDT NPZ (0..1 scaled to 0..100)" if values is not None else None
            )
            if values is None and confidence:
                values, source = chain["bfactor"], "native CIF CA B factors; NPZ absent"
            if (
                values is not None
                and arrays is not None
                and np.max(np.abs(values - chain["bfactor"])) > 0.05
            ):
                result["warnings"].append(f"Chain {name}: CIF and NPZ pLDDT disagree")
            result["chains"][name] = {
                "sequence_exact": True,
                "residues": length,
                "plddt_scale": "0..100",
                "plddt_source": source,
                "plddt": distribution([] if values is None else values),
                "fraction_plddt_ge70": float(np.mean(values >= 70)) if values is not None else None,
                "fraction_plddt_ge90": float(np.mean(values >= 90)) if values is not None else None,
                "ptm": confidence.get("chains_ptm", {}).get(str(index)),
            }
        result["reference_comparisons"] = [
            compare_reference(chains["A"], ref, target["sequence"], base)
            for ref in target["references"]
        ]
        pae_path = directory / ("pae_" + case["id"] + "_model_0.npz")
        if pae_path.is_file():
            with np.load(pae_path, allow_pickle=False) as data:
                pae = np.asarray(data["pae"], dtype=float)
            total = sum(len(c["sequence"]) for c in chains.values())
            if pae.shape != (total, total) or not np.isfinite(pae).all():
                raise ValueError("PAE must be a finite residue-by-residue matrix")
            length = len(case["sequence"])
            result["pae_angstrom"] = {"target_internal_mean": float(pae[:length, :length].mean())}
            if case.get("binder_sequence"):
                result["pae_angstrom"].update(
                    target_to_binder_mean=float(pae[:length, length:].mean()),
                    binder_to_target_mean=float(pae[length:, :length].mean()),
                )
    except (ValueError, TypeError, KeyError, OSError, RuntimeError) as exc:
        result.update(status="invalid_artifact", error=str(exc))
    return result


def generation_context(prior_root, targets, base):
    """Contextual prior complex refolds, not matched or blinded validation controls."""
    rows = []
    if prior_root is None:
        return rows
    for gene, target in targets.items():
        directory = (
            Path(prior_root)
            / (gene + "_reference")
            / "live/downloaded/results/intermediate_designs_inverse_folded/refold_cif"
        )
        for path in sorted(directory.glob("*.cif")):
            row = {
                "target_id": gene,
                "prediction": str(path),
                "context": "Native BoltzGen generation-complex refold with target_templates=true; not a custom fixed-sequence comparison baseline",
                "protocol": "BoltzGen FromGeneratedDataModule; target template conditioned",
            }
            config_path = directory.parents[1] / "config/folding.yaml"
            config_text = config_path.read_text() if config_path.is_file() else ""
            row["config_path"] = str(config_path)
            row["protocol_verified"] = (
                "FromGeneratedDataModule" in config_text and "target_templates: true" in config_text
            )
            if not row["protocol_verified"]:
                row.update(
                    status="unverified_protocol",
                    error="Native generation conditioning config missing or unexpected",
                )
                rows.append(row)
                continue
            row["config_sha256"] = hashlib.sha256(config_path.read_bytes()).hexdigest()
            try:
                candidates = [
                    c for c in chain_data(path).values() if c["sequence"] == target["sequence"]
                ]
                if len(candidates) != 1:
                    raise ValueError("Prior target sequence is not uniquely exact")
                chain = candidates[0]
                if (
                    not np.isfinite(chain["bfactor"]).all()
                    or np.any(chain["bfactor"] < 0)
                    or np.any(chain["bfactor"] > 1.00001)
                ):
                    raise ValueError("Pinned custom BoltzGen refold B factors must use 0..1 pLDDT")
                row.update(
                    status="compared",
                    target_CA_bfactor=distribution(chain["bfactor"]),
                    target_CA_bfactor_native_scale="0..1 (pinned BoltzGen predict writer)",
                    target_plddt_0_100=distribution(chain["bfactor"] * 100),
                    reference_comparisons=[
                        compare_reference(chain, ref, target["sequence"], base)
                        for ref in target["references"]
                    ],
                )
            except (ValueError, TypeError, KeyError, OSError, RuntimeError) as exc:
                row.update(status="invalid_artifact", error=str(exc))
            rows.append(row)
    return rows


def crossfold_context(root, targets, base):
    """The actual previous custom empty-MSA, template-free repeated comparisons."""
    if root is None or not (Path(root) / "receipt.json").is_file():
        return []
    root = Path(root)
    receipt_path = root / "receipt.json"
    receipt = json.loads(receipt_path.read_text())
    if (
        receipt.get("mode") != "fixed_sequence_refolding"
        or receipt.get("target_templates") is not False
        or receipt.get("sequence_redesign") is not False
    ):
        raise ValueError(
            "Prior comparison receipt does not verify the expected fixed-sequence protocol"
        )
    tasks = {t["id"]: t for t in receipt["tasks"]}
    rows = []
    for result in receipt["results"]:
        task = tasks[result["id"]]
        gene = task["target_id"]
        if gene not in targets or task["target_sequence"] != targets[gene]["sequence"]:
            continue
        seed = result["seed"]
        directory = root / f"seed_{seed}"
        cif = directory / "refold_cif" / (task["id"] + ".cif")
        archive = directory / "fold_out_npz" / (task["id"] + ".npz")
        row = {
            "target_id": gene,
            "id": task["id"],
            "seed": seed,
            "kind": "custom_complex_comparison"
            if task.get("binder_sequence")
            else "custom_target_alone",
            "designed_on": task.get("designed_on"),
            "scored_on": task.get("scored_on"),
            "prediction": str(cif),
            "archive": str(archive),
            "protocol": "Custom BoltzGen FromYamlDataModule + fixed-sequence wrapper; not standard Boltz CLI",
            "target_templates": False,
            "sequence_redesign": False,
            "checkpoint": receipt["checkpoint"],
            "receipt_sha256": hashlib.sha256(receipt_path.read_bytes()).hexdigest(),
            "context": "Actual prior custom comparison; empty MSA, no target templates, fixed sequences",
        }
        config = directory / "config.yaml"
        native_input = root / "inputs" / (task["id"] + ".yaml")
        row["config_sha256"] = (
            hashlib.sha256(config.read_bytes()).hexdigest() if config.is_file() else None
        )
        row["input_sha256"] = (
            hashlib.sha256(native_input.read_bytes()).hexdigest()
            if native_input.is_file()
            else None
        )
        if not cif.is_file() or not archive.is_file():
            row.update(status="missing_results")
            rows.append(row)
            continue
        try:
            chains = chain_data(cif)
            expected = {"A": task["target_sequence"]}
            if task.get("binder_sequence"):
                expected["B"] = task["binder_sequence"]
            if list(chains) != list(expected) or any(
                chains[k]["sequence"] != v for k, v in expected.items()
            ):
                raise ValueError(
                    "Prior custom comparison chain sequences are not exactly the archived inputs"
                )
            chain = chains["A"]
            if (
                not np.isfinite(chain["bfactor"]).all()
                or np.any(chain["bfactor"] < 0)
                or np.any(chain["bfactor"] > 1.00001)
            ):
                raise ValueError("Custom comparison pLDDT must use native 0..1 B-factor units")
            metrics = {}
            with np.load(archive, allow_pickle=False) as data:
                for key in [
                    "ptm",
                    "iptm",
                    "design_ptm",
                    "target_ptm",
                    "complex_plddt",
                    "design_to_target_iptm",
                ]:
                    if key in data and data[key].size == 1:
                        value = float(data[key].reshape(-1)[0])
                        metrics[key] = value if np.isfinite(value) else None
            row.update(
                status="compared",
                target_plddt_0_100=distribution(chain["bfactor"] * 100),
                metrics=metrics,
                reference_comparisons=[
                    compare_reference(chain, ref, targets[gene]["sequence"], base)
                    for ref in targets[gene]["references"]
                ],
            )
        except (ValueError, TypeError, KeyError, OSError, RuntimeError) as exc:
            row.update(status="invalid_artifact", error=str(exc))
        rows.append(row)
    return rows


def analyze(
    controls_path,
    prepared,
    output,
    *,
    results_root=None,
    base=None,
    prior_root=None,
    prior_crossfold=None,
):
    prepared, output = Path(prepared), Path(output)
    controls_path = Path(controls_path)
    controls = json.loads(controls_path.read_text())
    manifest = json.loads((prepared / "manifest.json").read_text())
    if manifest.get("controls_sha256") != hashlib.sha256(controls_path.read_bytes()).hexdigest():
        raise ValueError("Controls changed after input preparation")
    base = Path.cwd() if base is None else Path(base)
    results_root = prepared if results_root is None else Path(results_root)
    targets = {t["gene_symbol"]: t for t in controls["targets"]}
    rows = []
    for command in manifest["commands"]:
        for case in manifest["inputs"]:
            if case["condition"] != command["condition"]:
                continue
            row = {
                k: case[k]
                for k in ["id", "target_id", "kind", "condition", "msa", "target_template"]
            }
            row["seed"] = command["seed"]
            row["validation_scope"] = (
                "template-assisted recovery; not independent folding validation"
                if case["target_template"]
                else "no structure template supplied"
            )
            directory = results_root / command["expected_prediction_directory"] / case["id"]
            row.update(parse_prediction(directory, case, targets[case["target_id"]], base))
            rows.append(row)
    groups = defaultdict(list)
    for row in rows:
        groups[(row["id"], row["condition"])].append(row)
    summary = []
    for (case_id, condition), group in sorted(groups.items()):
        complete = [r for r in group if r["status"] == "completed"]
        references = defaultdict(list)
        for row in complete:
            for ref in row["reference_comparisons"]:
                if ref["rmsd_angstrom"] is not None:
                    references[(ref["role"], ref["pdb_id"])].append(ref["rmsd_angstrom"])
        summary.append(
            {
                "id": case_id,
                "target_id": group[0]["target_id"],
                "condition": condition,
                "validation_scope": group[0]["validation_scope"],
                "expected_seeds": [r["seed"] for r in group],
                "completed_seeds": len(complete),
                "statuses": {r["seed"]: r["status"] for r in group},
                "target_plddt": distribution(r["chains"]["A"]["plddt"]["mean"] for r in complete),
                "reference_rmsd_angstrom": [
                    {"role": role, "pdb_id": pdb, **distribution(values)}
                    for (role, pdb), values in references.items()
                ],
            }
        )
    report = {
        "package": manifest["package"],
        "expected_predictions": len(rows),
        "completed_predictions": sum(r["status"] == "completed" for r in rows),
        "rows": rows,
        "summary": summary,
        "prior_generation_context": generation_context(prior_root, targets, base),
        "prior_custom_context": crossfold_context(
            prior_crossfold or (Path(prior_root).parent / "crossfold" if prior_root else None),
            targets,
            base,
        ),
        "limitations": manifest.get("limitations", [])
        + [
            "Three requested seeds are replicate predictions, not independent equilibrium samples.",
            "Low template RMSD in template-conditioned arms is not independent folding validation.",
            "No success threshold or binding-affinity claim is inferred from confidence or RMSD.",
        ],
        "sources": [
            "https://github.com/jwohlwend/boltz/blob/main/src/boltz/data/write/writer.py",
            "https://github.com/jwohlwend/boltz/blob/main/src/boltz/data/write/mmcif.py",
        ],
    }
    prior_groups = defaultdict(list)
    for row in report["prior_custom_context"]:
        if row["status"] == "compared":
            prior_groups[(row["target_id"], row["kind"])].append(row)
    report["prior_custom_summary"] = [
        {
            "target_id": gene,
            "kind": kind,
            "predictions": len(group),
            "target_plddt_0_100": distribution(r["target_plddt_0_100"]["mean"] for r in group),
            "held_out_rmsd_angstrom": distribution(
                ref["rmsd_angstrom"]
                for row in group
                for ref in row["reference_comparisons"]
                if ref["role"] == "held_out"
            ),
        }
        for (gene, kind), group in sorted(prior_groups.items())
    ]
    output.mkdir(parents=True, exist_ok=True)
    (output / "analysis.json").write_text(json.dumps(report, indent=2, allow_nan=False))
    columns = [
        "id",
        "target_id",
        "condition",
        "seed",
        "status",
        "validation_scope",
        "target_plddt_mean",
        "binder_plddt_mean",
        "ptm",
        "iptm",
        "template_rmsd_angstrom",
        "held_out_rmsd_angstrom",
        "missing_artifacts",
    ]
    with (output / "predictions.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            flat = {k: row[k] for k in columns if k in row}
            flat.update(
                target_plddt_mean=row["chains"].get("A", {}).get("plddt", {}).get("mean"),
                binder_plddt_mean=row["chains"].get("B", {}).get("plddt", {}).get("mean"),
                ptm=row.get("confidence", {}).get("ptm"),
                iptm=row.get("confidence", {}).get("iptm"),
                missing_artifacts=",".join(row["missing_artifacts"]),
            )
            for role in ["template", "held_out"]:
                values = [
                    r["rmsd_angstrom"]
                    for r in row["reference_comparisons"]
                    if r["role"] == role and r["rmsd_angstrom"] is not None
                ]
                flat[role + "_rmsd_angstrom"] = min(values) if values else None
            writer.writerow(flat)
    lines = [
        "# Standard Boltz folding controls",
        "",
        f"Completed {report['completed_predictions']} of {len(rows)} expected predictions.",
        "",
        "| Case | Condition | Seeds complete | Target pLDDT median | Held-out CA RMSD median (Å) |",
        "|---|---|---:|---:|---:|",
    ]
    for item in summary:
        held = [r for r in item["reference_rmsd_angstrom"] if r["role"] == "held_out"]
        rmsd = ", ".join(f"{r['pdb_id']}: {r['median']:.3f}" for r in held) or "—"
        confidence = item["target_plddt"]["median"]
        lines.append(
            f"| {item['id']} | {item['condition']} | {item['completed_seeds']}/{len(item['expected_seeds'])} | "
            + (f"{confidence:.1f}" if confidence is not None else "—")
            + f" | {rmsd} |"
        )
    if report["prior_custom_summary"]:
        lines += [
            "",
            "Actual prior custom comparisons (empty MSA, no target templates; separate from native generation refolds):",
            "",
            "| Target / prior mode | Predictions | Target pLDDT median (0–100) | Held-out CA RMSD median (Å) |",
            "|---|---:|---:|---:|",
        ]
        for prior in report["prior_custom_summary"]:
            confidence = prior["target_plddt_0_100"]["median"]
            rmsd = prior["held_out_rmsd_angstrom"]["median"]
            lines.append(
                f"| {prior['target_id']} / {prior['kind']} | {prior['predictions']} | {confidence:.1f} | "
                + (f"{rmsd:.3f}" if rmsd is not None else "—")
                + " |"
            )
    lines += [
        "",
        "Template-assisted arms measure conditioned recovery, not independent folding validation. Held-out structures are withheld from supplied inputs; model-training independence is unknown.",
        "pLDDT is reported on a 0–100 scale. RMSD uses one proper CA rigid alignment on audited common input positions. Missing or invalid predictions remain explicit.",
        "Confidence and geometry do not establish binding or biological activity.",
        "",
    ]
    (output / "README.md").write_text("\n".join(lines))
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--controls", type=Path, required=True)
    p.add_argument("--prepared", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--results-root", type=Path)
    p.add_argument("--base", type=Path)
    p.add_argument(
        "--prior-root",
        type=Path,
        help="Native generation directory; sibling crossfold is auto-included separately",
    )
    p.add_argument(
        "--prior-crossfold", type=Path, help="Explicit actual prior custom comparison directory"
    )
    args = p.parse_args()
    result = analyze(
        args.controls,
        args.prepared,
        args.output,
        results_root=args.results_root,
        base=args.base,
        prior_root=args.prior_root,
        prior_crossfold=args.prior_crossfold,
    )
    print(
        json.dumps(
            {
                "expected": result["expected_predictions"],
                "completed": result["completed_predictions"],
            }
        )
    )


if __name__ == "__main__":
    main()
