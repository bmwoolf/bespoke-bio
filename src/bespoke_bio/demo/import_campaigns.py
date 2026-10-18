"""Validate and publish Jacob's downloaded campaign bundles without running inference."""

import argparse
import csv
import json
import math
import re
import shutil
import tempfile
from hashlib import sha256
from io import BytesIO
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

import numpy as np

from .import_science import chain_file, read_structure
from .published import PublishedSnapshot

RELEASE = "https://github.com/bmwoolf/bespoke-bio/releases/tag/protein-results-2026-10-18"
METRICS = (
    "iptm",
    "ptm",
    "protein_iptm",
    "design_ptm",
    "target_ptm",
    "design_to_target_iptm",
    "min_design_to_target_pae",
)


def read_json(path):
    return json.loads(Path(path).read_text())


def within(root, relative):
    root = Path(root).resolve()
    path = root / relative
    if path.is_symlink() or not path.resolve().is_relative_to(root):
        raise ValueError("Artifact path escapes the source directory")
    return path


def verified_bundle(root):
    root = Path(root)
    hashes = read_json(root / "manifest.json")
    for relative, digest in hashes.items():
        if sha256(within(root, relative).read_bytes()).hexdigest() != digest:
            raise ValueError(f"Bundle checksum mismatch: {relative}")
    bundle = read_json(root / "bundle.json")
    if len(bundle["target_states"]) != 1:
        raise ValueError("Expected one explicit target state")
    target = bundle["target_states"][0]
    if len(target["chains"]) != 1:
        raise ValueError("Expected one target protein chain")
    ids = [c["id"] for c in bundle["components"]]
    if len(ids) != len(set(ids)) or any(not re.fullmatch(r"binder_[a-f0-9]+", c) for c in ids):
        raise ValueError("Duplicate or invalid candidate identity")
    for component in bundle["components"]:
        if sha256(component["sequence"].encode()).hexdigest() != component["sequence_sha256"]:
            raise ValueError("Candidate sequence checksum mismatch")
        for path in (component["structure"], component["evidence"]["complex_structure"]):
            if path not in hashes:
                raise ValueError("Candidate structure is absent from the checksum manifest")
    return bundle, next(iter(target["chains"].values()))


def complex_roles(data, target_sequence, binder_sequence):
    model, chains = read_structure(data)
    keys = [chain for chain, sequence in chains.items() if sequence == binder_sequence]
    targets = [chain for chain, sequence in chains.items() if sequence == target_sequence]
    if len(chains) != 2 or len(keys) != 1 or len(targets) != 1 or keys == targets:
        raise ValueError("Complex target/key sequences or chain roles do not match")
    return model, {"target": targets, "binder": keys}


def metric_units(values):
    return {
        key: "Å" if any(s in key.lower() for s in ("rmsd", "angstrom", "pae")) else "0–1"
        for key in values
    }


def import_campaigns(source, archives, snapshot, destination, revision, shared_bundle=None):
    source, archives, snapshot, destination = map(Path, (source, archives, snapshot, destination))
    if destination.exists():
        raise FileExistsError(destination)
    manifest = read_json(snapshot / "manifest.json")
    PublishedSnapshot(snapshot, manifest["catalog_sha256"])
    indexed = list(csv.DictReader((source / "docs/TARGET_DATA_INDEX.csv").open()))
    identities = [(r["campaign"], r["target"], r["state"]) for r in indexed]
    if len(set(identities)) != len(identities):
        raise ValueError("Ambiguous target index")
    expanded = {
        (r["target_id"], r["state"]): r
        for r in read_json(source / "examples/expanded_variant_panel/inputs.json")
    }
    panel = {
        r["gene_symbol"]: r
        for r in read_json(source / "examples/variant_panel/panel.json")["targets"]
    }
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=".campaigns-", dir=destination.parent))
    counts = {
        "imported_jobs": 0,
        "native_candidates": 0,
        "comparison_complexes": 0,
        "target_controls": 0,
    }
    try:
        shutil.copytree(snapshot, temporary, dirs_exist_ok=True)

        # Preserve the richer independent refolds already published for the focused cohorts.
        for row in (r for r in indexed if r["campaign"] == "folding_validation"):
            focused, sequence = verified_bundle(
                archives
                / "scientific_campaign/folding_validation/focused_design"
                / row["target"]
                / "live/bundle"
            )
            expected = {c["id"]: c["sequence"] for c in focused["components"]}
            matches = [
                entry
                for entry in manifest["jobs"].values()
                if entry.get("design", {}).get("name") == row["target"]
                and entry.get("design", {}).get("amino_acid_sequence") == sequence
                and {c["id"]: c["sequence"] for c in entry["results"]["candidates"]} == expected
            ]
            if len(matches) != 1 or len(expected) != int(row["candidates"]):
                raise ValueError("Existing focused study differs from the complete source archive")

        def save(job_id, filename, data):
            digest = sha256(data).hexdigest()
            (temporary / "files" / digest).write_bytes(data)
            url = f"/api/jobs/{job_id}/artifacts/{filename}"
            if url in manifest["files"] and manifest["files"][url]["sha256"] != digest:
                raise ValueError("Artifact identity collision")
            manifest["files"][url] = {"sha256": digest, "size": len(data)}
            return url

        # Older API snapshots kept the target only inside the complex. Resolve roles
        # against their original bundle, never by assuming a chain letter.
        for job_id, entry in manifest["jobs"].items():
            candidates = entry["results"]["candidates"]
            if all(c.get("complex_chains") and c.get("target_structure_url") for c in candidates):
                continue
            archive = (
                temporary / "files" / manifest["files"][entry["results"]["download_url"]]["sha256"]
            )
            with ZipFile(archive) as saved:
                paths = [p for p in saved.namelist() if p.endswith("bundle/bundle.json")]
                if len(paths) != 1:
                    raise ValueError("Ambiguous original saved bundle")
                bundle = json.loads(saved.read(paths[0]))
            states = {s["id"]: s for s in bundle["target_states"]}
            for candidate in candidates:
                state = states[candidate["target_state_id"]]
                if len(state["chains"]) != 1:
                    raise ValueError("Ambiguous saved target chains")
                sequence = next(iter(state["chains"].values()))
                component = next(c for c in bundle["components"] if c["id"] == candidate["id"])
                if component["sequence"] != candidate["sequence"]:
                    raise ValueError("Saved key differs from its source bundle")
                digest = manifest["files"][candidate["complex_structure_url"]]["sha256"]
                model, roles = complex_roles(
                    (temporary / "files" / digest).read_bytes(), sequence, candidate["sequence"]
                )
                candidate["complex_chains"] = roles
                candidate["target_structure_url"] = save(
                    job_id, candidate["id"] + "/target.cif", chain_file(model, roles["target"][0])
                )
                for metric, value in candidate.get("native_filter_metrics", {}).items():
                    if (
                        metric in candidate["confidence_metrics"]
                        and candidate["confidence_metrics"][metric] != value
                    ):
                        raise ValueError("Saved native and interface metrics disagree")
                    candidate["confidence_metrics"][metric] = value
                candidate["confidence_units"] = metric_units(candidate["confidence_metrics"])

        def add_bundle(bundle_root, identity, name, state, description, provenance, comparisons=()):
            bundle, target_sequence = verified_bundle(bundle_root)
            expected = provenance.get("sequence", target_sequence)
            if expected != target_sequence:
                raise ValueError("Bundle target differs from declared generation input")
            job_id = sha256(f"jacob-campaign:{identity}".encode()).hexdigest()[:32]
            if job_id in manifest["jobs"]:
                raise ValueError("This source job was already imported")
            target_id = f"jacob:{identity}"
            normalized = []
            extra_files = {}
            for component in bundle["components"]:
                cid, sequence = component["id"], component["sequence"]
                binder = within(bundle_root, component["structure"]).read_bytes()
                complex_data = within(
                    bundle_root, component["evidence"]["complex_structure"]
                ).read_bytes()
                _, chains = read_structure(binder)
                if list(chains.values()) != [sequence]:
                    raise ValueError("Independent binder sequence mismatch")
                model, roles = complex_roles(complex_data, target_sequence, sequence)
                evaluations = [e for e in bundle["evaluations"] if e["subject_id"] == cid]
                filters = [
                    e for e in evaluations if e["method"] == "BoltzGen native default filters"
                ]
                if len(filters) != 1 or filters[0]["status"] not in {"pass", "fail", "unknown"}:
                    raise ValueError("Missing or ambiguous native filter result")
                native = filters[0]
                if native["target_state_id"] != bundle["target_states"][0]["id"]:
                    raise ValueError("Filter belongs to another target")
                metrics = dict(native["metrics"])
                for e in evaluations:
                    for key, value in e["metrics"].items():
                        if key in metrics and metrics[key] != value:
                            raise ValueError("Conflicting native metrics")
                        metrics[key] = value
                if any(v is not None and not math.isfinite(v) for v in metrics.values()):
                    raise ValueError("Nonfinite native metrics")
                candidate = {
                    "id": cid,
                    "sequence": sequence,
                    "native_filter_passed": {"pass": True, "fail": False, "unknown": None}[
                        native["status"]
                    ],
                    "failed_filters": native["details"].get("failed_filters", []),
                    "confidence_metrics": metrics,
                    "confidence_units": metric_units(metrics),
                    "binder_structure_url": save(job_id, f"{cid}/binder.cif", binder),
                    "complex_structure_url": save(job_id, f"{cid}/complex.cif", complex_data),
                    "target_structure_url": save(
                        job_id, f"{cid}/target.cif", chain_file(model, roles["target"][0])
                    ),
                    "complex_chains": roles,
                    "comparison_scores": [],
                }
                for row in comparisons:
                    task = row["task"]
                    if task.get("binder_id") != cid:
                        continue
                    if task["binder_sequence"] != sequence:
                        raise ValueError("Comparison binder sequence mismatch")
                    data = row["structure"].read_bytes()
                    comparison_model, comparison_roles = complex_roles(
                        data, task["target_sequence"], sequence
                    )
                    prefix = f"comparisons/{task['id']}/seed-{row['seed']}"
                    # Both isolated views are explicitly extracted from this comparison pose.
                    record = {
                        "id": f"{task['id']}-{row['seed']}",
                        "seed": row["seed"],
                        "scored_on": task["scored_on"],
                        "metrics": row["metrics"],
                        "complex_chains": comparison_roles,
                        "complex_structure_url": save(job_id, prefix + "/complex.cif", data),
                        "binder_structure_url": save(
                            job_id,
                            prefix + "/binder.cif",
                            chain_file(comparison_model, comparison_roles["binder"][0]),
                        ),
                        "target_structure_url": save(
                            job_id,
                            prefix + "/target.cif",
                            chain_file(comparison_model, comparison_roles["target"][0]),
                        ),
                    }
                    candidate["comparison_scores"].append(record)
                    extra_files[prefix + "/complex.cif"] = data
                    extra_files[prefix + "/evidence.json"] = json.dumps(
                        {"task": task, "seed": row["seed"], "metrics": row["metrics"]}
                    ).encode()
                    counts["comparison_complexes"] += 1
                normalized.append(candidate)
            archive = BytesIO()
            with ZipFile(archive, "w", ZIP_DEFLATED) as output:
                for path in sorted(Path(bundle_root).rglob("*")):
                    if path.is_file() and not path.is_symlink():
                        output.write(path, "bundle/" + str(path.relative_to(bundle_root)))
                for path, data in extra_files.items():
                    output.writestr(path, data)
                output.writestr("provenance.json", json.dumps(provenance, allow_nan=False))
            zip_url = save(job_id, "results.zip", archive.getvalue())
            download = f"/api/jobs/{job_id}/download"
            manifest["files"][download] = manifest["files"][zip_url]
            entry = {
                "job": {
                    "job_id": job_id,
                    "target_id": target_id,
                    "status": "completed",
                    "mode": "published_native_batch",
                    "execution_performed": True,
                    "dispatch_attempted": None,
                    "progress": None,
                    "error": None,
                },
                "results": {
                    "run_id": job_id,
                    "target_id": target_id,
                    "candidates": normalized,
                    "download_url": download,
                },
                "design": {
                    "target_id": target_id,
                    "name": name,
                    "description": description,
                    "source": "Jacob’s campaigns",
                    "state": state,
                    "source_revision": revision,
                    "amino_acid_sequence": target_sequence,
                    "sequence_sha256": sha256(target_sequence.encode()).hexdigest(),
                    "input": provenance,
                    "complete_study_url": RELEASE,
                    "target_view": "Target and Key are extracted from the selected predicted complex. Original designs and fixed-sequence comparisons use different protocols.",
                    "limitations": [
                        "Computational predictions; binding, affinity and selectivity are unmeasured.",
                        "Constructs can omit residues, modifications, cofactors and native assembly partners.",
                    ],
                },
            }
            manifest["jobs"][job_id] = entry
            counts["imported_jobs"] += 1
            counts["native_candidates"] += len(normalized)
            return entry

        for campaign, dirname in (
            ("variant_campaign", "variant_campaign_results"),
            ("variant_expansion", "variant_expansion_results"),
        ):
            root = archives / dirname / campaign
            cohorts = read_json(root / "generation_manifest.json")["cohorts"]
            rows = [r for r in indexed if r["campaign"] == campaign]
            if {(r["target"], r["state"]) for r in rows} != {
                (c["target_id"], c["state"]) for c in cohorts
            }:
                raise ValueError("Archive and target index disagree")
            tasks = read_json(root / "crossfold_tasks.json")["tasks"]
            sequences = {(c["target_id"], c["state"]): c["sequence"] for c in cohorts}
            if len(sequences) != len(cohorts):
                raise ValueError("Duplicate generation state")
            comparisons = []
            for task in tasks:
                if not re.fullmatch(r"[A-Za-z0-9_]+", task["id"]):
                    raise ValueError("Invalid comparison identity")
                if task["target_sequence"] != sequences[(task["target_id"], task["scored_on"])]:
                    raise ValueError("Comparison state differs from the exact source construct")
                for seed in (11, 22, 33):
                    path = root / f"crossfold/seed_{seed}/refold_cif/{task['id']}.cif"
                    if "binder_id" not in task:
                        _, chains = read_structure(path.read_bytes())
                        if list(chains.values()) != [task["target_sequence"]]:
                            raise ValueError("Target-only control sequence mismatch")
                        counts["target_controls"] += 1
                        continue
                    npz = root / f"crossfold/seed_{seed}/fold_out_npz/{task['id']}.npz"
                    with np.load(npz, allow_pickle=False) as values:
                        metrics = {k: float(values[k].item()) for k in METRICS}
                    if any(not math.isfinite(v) for v in metrics.values()):
                        raise ValueError("Nonfinite comparison scores")
                    comparisons.append(
                        {"task": task, "seed": seed, "structure": path, "metrics": metrics}
                    )
            consumed = 0
            for row in rows:
                name, state = row["target"], row["state"]
                cohort = next(c for c in cohorts if (c["target_id"], c["state"]) == (name, state))
                if cohort["status"] != "completed":
                    raise ValueError("Unfinished campaign cohort")
                selected = [
                    r
                    for r in comparisons
                    if (r["task"]["target_id"], r["task"]["designed_on"]) == (name, state)
                ]
                provenance = {"campaign": campaign, "sequence": cohort["sequence"], "state": state}
                if campaign == "variant_expansion":
                    provenance.update(expanded[(name, state)])
                else:
                    target = panel[name]
                    provenance.update(
                        {
                            k: target[k]
                            for k in ("accession", "construct", "caveats", "sequence_source")
                        }
                    )
                if provenance["sequence"] != cohort["sequence"]:
                    raise ValueError("Construct export and generation sequence disagree")
                entry = add_bundle(
                    within(root, cohort["directory"]) / "live/bundle",
                    f"{campaign}:{name}:{state}",
                    name,
                    state,
                    f"{'Reference' if state == 'reference' else state} · {'Original panel' if campaign == 'variant_campaign' else 'Expanded panel'}",
                    provenance,
                    selected,
                )
                candidates = entry["results"]["candidates"]
                if len(candidates) != int(row["candidates"]):
                    raise ValueError("Candidate count differs from target index")
                consumed += sum(len(c["comparison_scores"]) for c in candidates)
                folder = (
                    "expanded_variant_panel" if campaign == "variant_expansion" else "variant_panel"
                )
                summaries = list(
                    csv.DictReader((source / f"examples/{folder}/candidates.csv").open())
                )
                for c in candidates:
                    match = [
                        r
                        for r in summaries
                        if (
                            r["target_id"],
                            r["state"],
                            r.get("binder_id", r.get("id", "").removeprefix(f"{name}_{state}_")),
                        )
                        == (name, state, c["id"])
                    ]
                    if len(match) != 1 or match[0]["sequence"] != c["sequence"]:
                        raise ValueError("Git candidate index differs from archive")
                    exported = match[0]
                    status = exported.get("native_filter", exported.get("native_filter_status"))
                    if (
                        c["native_filter_passed"]
                        != {"pass": True, "fail": False, "unknown": None}[status]
                    ):
                        raise ValueError("Git native filter status differs from archive")
                    for key, value in c["confidence_metrics"].items():
                        if (
                            exported.get(key)
                            and value is not None
                            and not math.isclose(float(exported[key]), value, abs_tol=1e-9)
                        ):
                            raise ValueError("Git native metrics differ from archive")
                    assessment = exported.get("assembly_assessment")
                    if assessment and assessment not in {"not_assessed", "not assessed"}:
                        c["assembly_assessment"] = assessment
            if consumed != len(comparisons):
                raise ValueError("Some comparisons could not be joined to a candidate")

        add_bundle(
            source / "examples/boltzgen_api_results/bundle",
            "pdl1-api-smoke",
            "PD-L1",
            "reference",
            "Sequence input · API test",
            {"input_mode": "sequence_only"},
        )
        if shared_bundle:
            add_bundle(
                Path(shared_bundle),
                "666dd196f6324942b258a3e28ca38fc7",
                "PD-L1",
                "reference",
                "Sequence input · Shared service test",
                {"input_mode": "sequence_only"},
            )
        (temporary / "manifest.json").write_text(json.dumps(manifest, allow_nan=False))
        PublishedSnapshot(temporary, manifest["catalog_sha256"])
        counts.update(
            total_jobs=len(manifest["jobs"]),
            total_candidates=sum(
                len(j["results"]["candidates"]) for j in manifest["jobs"].values()
            ),
        )
        (temporary / "import-report.json").write_text(json.dumps(counts, indent=2))
        temporary.rename(destination)
        return counts
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("source", "archives", "snapshot", "output", "revision"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--shared-bundle")
    args = parser.parse_args()
    print(
        json.dumps(
            import_campaigns(
                args.source,
                args.archives,
                args.snapshot,
                args.output,
                args.revision,
                args.shared_bundle,
            )
        )
    )


if __name__ == "__main__":
    main()
