#!/usr/bin/env python3
"""Prepare standard Boltz validation of passing focused designs; never run inference.

Targets use verified existing frozen MSAs; de novo binders use empty MSAs.
No coordinates/templates are supplied for any chain. All native source candidates
remain in provenance, including failures and passing candidates not selected.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import shlex
import shutil
from importlib.metadata import version
from pathlib import Path

VERSION = "2.2.1"
WHEEL_SHA256 = "b8c62bbdede1922931d9203118f62c858f11aa699bf91fd4c05a5ed6a6d8b4fc"
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def experimental_references(target):
    """Audit experimental coordinates for posthoc analysis, never predictor input."""
    import gemmi

    references = []
    for ref in target.get("references", []):
        source = Path(ref["source"])
        matched = Path(ref["matched_structure"])
        if sha(source) != ref["source_sha256"]:
            raise ValueError("Experimental reference source SHA256 mismatch")
        structure = gemmi.read_structure(str(matched))
        if len(structure) != 1:
            raise ValueError("Experimental matched reference requires one model")
        residues = list(structure[0][ref["matched_chain"]])
        positions = [int(m["input_position"]) for m in ref["mapping"]]
        if len(set(positions)) != len(positions) or [r.seqid.num for r in residues] != positions:
            raise ValueError("Experimental reference mapping/numbering mismatch")
        for residue, position, mapping in zip(residues, positions, ref["mapping"]):
            observed = gemmi.find_tabulated_residue(residue.name).one_letter_code
            if (not 1 <= position <= len(target["sequence"])
                    or observed != target["sequence"][position - 1]
                    or mapping["observed_residue"] != observed
                    or mapping["reference_residue"] != observed
                    or residue.find_atom("CA", "*") is None):
                raise ValueError("Experimental reference exact sequence/CA mapping mismatch")
        if not re.fullmatch(r"[A-Za-z0-9]+", ref["pdb_id"]):
            raise ValueError("Unsafe experimental PDB identifier")
        references.append((matched, {"role": ref["role"], "pdb_id": ref["pdb_id"],
            "reference_chain": ref["matched_chain"], "prediction_chain": "A",
            "input_positions": positions, "sha256": sha(matched),
            "source_sha256": ref["source_sha256"],
            "historical_role_only": "Coordinates reserved for posthoc evaluation; no predictive templates"}))
    return references


def sequence_hash(sequence):
    if not sequence or set(sequence) - AMINO_ACIDS:
        raise ValueError("Exact canonical protein sequence required")
    return hashlib.sha256(sequence.encode()).hexdigest()


def contained(root, relative):
    root = Path(root).resolve()
    path = (root / relative).resolve()
    if not path.is_relative_to(root) or not path.is_file():
        raise ValueError("Source artifact must be a contained existing file")
    return path


def msa_query(path):
    query, started = [], False
    for line in Path(path).read_text().splitlines():
        if line.startswith("#"):
            continue
        if line.startswith(">"):
            if started:
                break
            started = True
        elif started:
            query.append(line.strip())
    return "".join(query)


def select_candidates(bundle, target_sequence, limit):
    states = bundle.get("target_states", [])
    if len(states) != 1 or list(states[0]["chains"].values()) != [target_sequence]:
        raise ValueError("Native bundle target sequence differs from experimental controls")
    target_id = states[0]["id"]
    evaluations = {
        e["subject_id"]: e
        for e in bundle["evaluations"]
        if e["method"] == "BoltzGen native default filters"
    }
    audit = []
    for component in bundle["components"]:
        digest = sequence_hash(component["sequence"])
        if digest != component["sequence_sha256"]:
            raise ValueError("Native binder sequence hash mismatch")
        if not re.fullmatch(r"[A-Za-z0-9_]+", component["id"]):
            raise ValueError("Unsafe native candidate identifier")
        evaluation = evaluations.get(component["id"])
        if not evaluation or evaluation["target_state_id"] != target_id:
            raise ValueError("Native filter evaluation missing or linked to another target")
        score = evaluation.get("metrics", {}).get("quality_score")
        score = float(score) if score is not None else None
        if score is not None and not math.isfinite(score):
            raise ValueError("Nonfinite native quality score")
        audit.append(
            {
                "candidate_id": component["id"],
                "sequence_sha256": digest,
                "native_filter_status": evaluation["status"],
                "native_quality_score": score,
                "native_evaluation": evaluation,
                "component": component,
                "selected": False,
            }
        )
    passed = [r for r in audit if r["native_filter_status"] == "pass"]
    if not passed:
        raise ValueError("No native-filter-passing binders; refusing to promote failed candidates")
    passed.sort(
        key=lambda r: (
            -(r["native_quality_score"] if r["native_quality_score"] is not None else -math.inf),
            r["candidate_id"],
        )
    )
    selected, seen = [], set()
    for row in passed:
        if row["sequence_sha256"] in seen or len(selected) >= limit:
            continue
        row["selected"] = True
        selected.append(row)
        seen.add(row["sequence_sha256"])
    return selected, audit


def positive_context(path):
    if path is None or not Path(path).is_file():
        return {"status": "not_provided", "interpretation": "No positive-control result invented"}
    source = json.loads(Path(path).read_text())
    rows = [
        r
        for r in source.get("cases", [])
        if r.get("condition") == "positive_msa_control" and r.get("status") == "assessed"
    ]
    return {
        "source": str(path),
        "sha256": sha(path),
        "pdb_id": "1WQJ",
        "scope": "Previous experimental positive control with paired both-chain MSAs and target-only template; not proof of de novo binder recovery",
        "results": [
            {
                "seed": r["seed"],
                "reference_contact_recall": r["reference_contact_recovery"]["contact_recall"],
                "interface_CA_jointfit_RMSD_A": r["reference_contact_recovery"][
                    "reference_interface_CA_jointfit_RMSD_A"
                ],
            }
            for r in rows
        ],
    }


def prepare(
    design_root,
    controls_path,
    standard_controls,
    output,
    *,
    genes=("FGF2", "IL1RN"),
    per_gene=4,
    seeds=(11, 22, 33),
    positive_audit=None,
):
    design_root, controls_path, standard_controls, output = map(
        Path, (design_root, controls_path, standard_controls, output)
    )
    if not 1 <= per_gene <= 4 or len(set(genes)) != len(genes) or not genes:
        raise ValueError("Select one to four binders per distinct requested target")
    if (
        not seeds
        or any(isinstance(s, bool) or not isinstance(s, int) or s < 0 for s in seeds)
        or len(set(seeds)) != len(seeds)
    ):
        raise ValueError("Seeds must be distinct nonnegative integers")
    if output.exists() and any(output.iterdir()):
        raise ValueError("Output must be new or empty")
    controls = json.loads(controls_path.read_text())
    standard = json.loads((standard_controls / "manifest.json").read_text())
    if standard.get("controls_sha256") != sha(controls_path):
        raise ValueError("Experimental control sequence/source hash changed")
    targets = {t["gene_symbol"]: t for t in controls["targets"]}
    sources = []
    # Validate every source and selection before creating outputs, including no-pass failures.
    for gene in genes:
        if not re.fullmatch(r"[A-Za-z0-9_]+", gene) or gene not in targets:
            raise ValueError("Unknown or unsafe target identifier")
        sequence = targets[gene]["sequence"]
        digest = sequence_hash(sequence)
        if targets[gene].get("sequence_sha256", digest) != digest:
            raise ValueError("Experimental target sequence hash mismatch")
        bundle_root = design_root / gene / "live/bundle"
        bundle_path = bundle_root / "bundle.json"
        bundle = json.loads(bundle_path.read_text())
        try:
            selected, audit = select_candidates(bundle, sequence, per_gene)
        except ValueError as exc:
            raise ValueError(f"{gene}: {exc}") from exc
        plan_path = design_root / gene / "prepared/boltzgen_plan.json"
        plan = json.loads(plan_path.read_text())
        if plan.get("target_sequence_sha256") != digest:
            raise ValueError("Focused design target hash differs from controls")
        frozen = standard["frozen_msas"][gene]
        msa = contained(standard_controls, frozen["path"])
        if sha(msa) != frozen["sha256"] or msa_query(msa) != sequence:
            raise ValueError("Frozen MSA hash or exact query sequence mismatch")
        for row in audit:
            component = row["component"]
            contained(bundle_root, component["structure"])
            complex_path = contained(bundle_root, component["evidence"]["complex_structure"])
            if component["evidence"].get("source_sha256") != sha(complex_path):
                raise ValueError("Native design complex hash mismatch")
        sources.append((gene, sequence, bundle_root, plan_path, frozen, msa, selected, audit, experimental_references(targets[gene])))
    output.mkdir(parents=True, exist_ok=True)
    inputs, selections, provenance, frozen_receipt = [], [], {}, {}
    for gene, sequence, bundle_root, plan_path, frozen, msa, selected, audit, references in sources:
        provenance_dir = output / "provenance" / gene
        shutil.copytree(bundle_root, provenance_dir / "native_bundle")
        shutil.copy2(plan_path, provenance_dir / "design_plan.json")
        metrics = (
            design_root
            / gene
            / "live/downloaded/results/final_ranked_designs/all_designs_metrics.csv"
        )
        if metrics.is_file():
            shutil.copy2(metrics, provenance_dir / metrics.name)
        analysis_references = []
        for reference_path, metadata in references:
            relative = Path("provenance") / gene / "experimental_refs" / (metadata["pdb_id"] + "_matched.cif")
            (output / relative.parent).mkdir(parents=True, exist_ok=True)
            shutil.copy2(reference_path, output / relative)
            analysis_references.append({**metadata, "path": str(relative)})
        files = {
            str(p.relative_to(provenance_dir)): sha(p)
            for p in sorted(provenance_dir.rglob("*"))
            if p.is_file()
        }
        provenance[gene] = {
            "source_bundle": str(bundle_root / "bundle.json"),
            "source_bundle_sha256": sha(bundle_root / "bundle.json"),
            "all_candidates": len(audit),
            "native_passed": sum(r["native_filter_status"] == "pass" for r in audit),
            "selected": len(selected),
            "raw_metrics_copied": metrics.is_file(),
            "copied_files_sha256": files,
        }
        msa_relative = Path("msas") / (gene + ".a3m")
        (output / msa_relative.parent).mkdir(exist_ok=True)
        shutil.copy2(msa, output / msa_relative)
        frozen_receipt[gene] = {
            **frozen,
            "path": str(msa_relative),
            "sequence_sha256": sequence_hash(sequence),
        }
        for row in audit:
            selections.append(
                {"target_id": gene, **{k: v for k, v in row.items() if k != "component"}}
            )
        for row in selected:
            candidate = row["component"]
            identifier = gene + "__" + candidate["id"]
            source_pose = str(
                Path("provenance")
                / gene
                / "native_bundle"
                / candidate["evidence"]["complex_structure"]
            )
            common = {
                "target_id": gene,
                "candidate_id": candidate["id"],
                "binder_sequence": candidate["sequence"],
                "binder_sequence_sha256": candidate["sequence_sha256"],
                "binder_chain": "Z",
                "source_native_filter": "pass",
                "target_template": False,
                "source_native_evaluation": row["native_evaluation"],
                "source_design_complex_pose_EVALUATION_ONLY": source_pose,
                "source_design_complex_pose_sha256": candidate["evidence"]["source_sha256"],
                "pose_interpretation": "Agreement with a computational design pose only; not experimental accuracy",
                "target_integrity": "not_applicable: no target template supplied",
                "target_analysis_references_EVALUATION_ONLY": analysis_references,
            }
            for alone in (False, True):
                condition = "binder_alone" if alone else "target_msa_only_complexes"
                target_chains = (
                    []
                    if alone
                    else [
                        {
                            "chain": "A",
                            "sequence": sequence,
                            "sequence_sha256": sequence_hash(sequence),
                            "mapping": [{"input_position": p} for p in range(1, len(sequence) + 1)],
                        }
                    ]
                )
                schema = {"version": 1, "sequences": []}
                if not alone:
                    schema["sequences"].append(
                        {"protein": {"id": "A", "sequence": sequence, "msa": str(msa_relative)}}
                    )
                schema["sequences"].append(
                    {"protein": {"id": "Z", "sequence": candidate["sequence"], "msa": "empty"}}
                )
                name = identifier + ("__binder_alone" if alone else "__msa_only")
                path = Path("inputs") / condition / (name + ".yaml")
                (output / path.parent).mkdir(parents=True, exist_ok=True)
                (output / path).write_text(json.dumps(schema, indent=2) + "\n")
                inputs.append(
                    {
                        **common,
                        "id": name,
                        "kind": "binder_alone" if alone else "fixed_binder_complex",
                        "context": condition,
                        "condition": condition,
                        "target_chains": target_chains,
                        "input": str(path),
                        "sha256": sha(output / path),
                    }
                )
    commands = []
    for condition in ("target_msa_only_complexes", "binder_alone"):
        for seed in seeds:
            argv = [
                "boltz",
                "predict",
                f"inputs/{condition}",
                "--model",
                "boltz2",
                "--out_dir",
                f"results/{condition}/seed_{seed}",
                "--cache",
                "cache",
                "--seed",
                str(seed),
                "--recycling_steps",
                "3",
                "--sampling_steps",
                "200",
                "--diffusion_samples",
                "1",
                "--max_parallel_samples",
                "1",
                "--devices",
                "1",
                "--accelerator",
                "gpu",
                "--output_format",
                "mmcif",
                "--num_workers",
                "0",
                "--no_kernels",
                "--write_full_pae",
            ]
            commands.append(
                {
                    "condition": condition,
                    "seed": seed,
                    "argv": argv,
                    "expected_prediction_directory": f"results/{condition}/seed_{seed}/boltz_results_{condition}/predictions",
                }
            )
    manifest = {
        "package": f"boltz=={VERSION}",
        "wheel_sha256": WHEEL_SHA256,
        "checkpoint": "boltz2_conf.ckpt",
        "controls_sha256": sha(controls_path),
        "standard_control_manifest_sha256": sha(standard_controls / "manifest.json"),
        "seeds": list(seeds),
        "inputs": inputs,
        "commands": commands,
        "prediction_count": len(inputs) * len(seeds),
        "native_parser_validated": False,
        "frozen_msas": frozen_receipt,
        "source_provenance": provenance,
        "selection_policy": "Native filter pass required; quality_score descending then candidate ID; unique binder sequences; up to four per gene. No promotion of failed candidates.",
        "selection_audit": selections,
        "positive_control_context": positive_context(positive_audit),
        "limitations": [
            "No template coordinates supplied for ANY chain; no design pose enters prediction.",
            "Target MSA reused exactly from successful target-folding controls; binder MSA empty.",
            "Native filter pass and positive-control success do not establish designed-binder binding.",
            "No-template target integrity against an input template is not applicable.",
            "Original design pose may be used only for agreement analysis, never experimental accuracy.",
            "Three requested seeds are replicate predictions, not equilibrium or significance samples.",
        ],
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (output / "commands.sh").write_text(
        "#!/usr/bin/env bash\nset -euo pipefail\n"
        + "\n".join(shlex.join(c["argv"]) for c in commands)
        + "\n"
    )
    return manifest


def validate_native(output, mol_dir):
    if version("boltz") != VERSION:
        raise ValueError("Native validation requires standard boltz==2.2.1")
    from boltz.data.parse.a3m import parse_a3m
    from boltz.data.parse.schema import parse_boltz_schema

    output, mol_dir = Path(output).resolve(), Path(mol_dir).resolve()
    manifest = json.loads((output / "manifest.json").read_text())
    receipt = []
    original = Path.cwd()
    try:
        os.chdir(output)
        for row in manifest["inputs"]:
            path = Path(row["input"])
            if sha(path) != row["sha256"]:
                raise ValueError("Prepared input changed before native validation")
            schema = json.loads(path.read_text())
            binder = next(p["protein"] for p in schema["sequences"] if p["protein"]["id"] == "Z")
            if binder.get("msa") != "empty":
                raise ValueError("De novo binder MSA must remain explicitly empty")
            if schema.get("templates"):
                raise ValueError("No template conditioning is permitted")
            parsed = parse_boltz_schema(row["id"], schema, {}, mol_dir, boltz_2=True)
            expected = {c["chain"]: c["sequence"] for c in row["target_chains"]}
            expected["Z"] = row["binder_sequence"]
            actual = {
                str(c.chain_name): parsed.sequences[c.entity_id] for c in parsed.record.chains
            }
            if actual != expected or parsed.record.templates or parsed.templates:
                raise ValueError("Native chain sequences or no-template policy differ")
            receipt.append(
                {
                    "id": row["id"],
                    "chain_lengths": {k: len(v) for k, v in actual.items()},
                    "template_records": 0,
                    "binder_msa": "empty",
                }
            )
        for gene, msa in manifest["frozen_msas"].items():
            if (
                sha(msa["path"]) != msa["sha256"]
                or sequence_hash(msa_query(msa["path"])) != msa["sequence_sha256"]
            ):
                raise ValueError("Frozen target MSA changed before native validation")
            parsed = parse_a3m(Path(msa["path"]), taxonomy=None, max_seqs=16384)
            receipt.append(
                {"target_id": gene, "native_deduplicated_msa_sequences": len(parsed.sequences)}
            )
    finally:
        os.chdir(original)
    manifest["native_parser_validated"] = True
    (output / "native_validation.json").write_text(json.dumps(receipt, indent=2) + "\n")
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--design-root", type=Path, default=Path("outputs/folding_validation/focused_design")
    )
    parser.add_argument(
        "--controls", type=Path, default=Path("outputs/folding_validation/research/controls.json")
    )
    parser.add_argument(
        "--standard-controls",
        type=Path,
        default=Path("outputs/folding_validation/standard_controls"),
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--genes", nargs="+", default=["FGF2", "IL1RN"])
    parser.add_argument("--per-gene", type=int, default=4)
    parser.add_argument("--seeds", nargs="+", type=int, default=[11, 22, 33])
    parser.add_argument(
        "--positive-audit",
        type=Path,
        default=Path("outputs/folding_validation/positive_pose_audit.json"),
    )
    parser.add_argument("--validate-native", action="store_true")
    parser.add_argument("--mol-dir", type=Path, default=Path.home() / ".boltz/mols")
    args = parser.parse_args()
    if args.validate_native:
        print(json.dumps({"validated": len(validate_native(args.output, args.mol_dir))}))
    else:
        result = prepare(
            args.design_root,
            args.controls,
            args.standard_controls,
            args.output,
            genes=args.genes,
            per_gene=args.per_gene,
            seeds=args.seeds,
            positive_audit=args.positive_audit,
        )
        print(
            json.dumps(
                {
                    "selected": sum(v["selected"] for v in result["source_provenance"].values()),
                    "predictions": result["prediction_count"],
                }
            )
        )


if __name__ == "__main__":
    main()
