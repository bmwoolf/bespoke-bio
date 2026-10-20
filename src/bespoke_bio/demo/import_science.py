"""Import Jacob's explicit public scientific export into an immutable demo snapshot."""

import argparse
import csv
import json
import math
import shutil
import tempfile
from collections import defaultdict
from hashlib import sha256
from io import BytesIO, StringIO
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

from Bio.Data.PDBData import protein_letters_3to1
from Bio.PDB import MMCIFIO, MMCIFParser, Select

from .published import PublishedSnapshot


def read_structure(data):
    structure = MMCIFParser(QUIET=True).get_structure("prediction", StringIO(data.decode()))
    if len(structure) != 1:
        raise ValueError("Expected one coordinate model")
    chains = {}
    for chain in structure[0]:
        residues = list(chain)
        sequence = "".join(protein_letters_3to1[r.resname] for r in residues)
        if not sequence or any(not {"N", "CA", "C", "O"} <= set(r.child_dict) for r in residues):
            raise ValueError("Incomplete protein backbone")
        if not all(math.isfinite(float(v)) for atom in chain.get_atoms() for v in atom.coord):
            raise ValueError("Nonfinite coordinates")
        chains[chain.id] = sequence
    return structure, chains


def chain_file(structure, chain_id):
    class Chain(Select):
        def accept_chain(self, chain):
            return chain.id == chain_id

    stream = StringIO()
    writer = MMCIFIO()
    writer.set_structure(structure)
    writer.save(stream, select=Chain())
    return stream.getvalue().encode()


def import_science(source, snapshot, destination, revision):
    source, snapshot, destination = map(Path, (source, snapshot, destination))
    if destination.exists():
        raise FileExistsError(destination)
    manifest = json.loads((snapshot / "manifest.json").read_text())
    PublishedSnapshot(snapshot, manifest["catalog_sha256"])
    science = json.loads((source / "science.json").read_text())
    if science["schema_version"] != 1 or science["partial"]:
        raise ValueError("Use a complete public scientific export")

    def verified(entry):
        path = (source / entry["path"]).resolve()
        if not path.is_relative_to(source.resolve()) or path.is_symlink():
            raise ValueError("Invalid source artifact path")
        data = path.read_bytes()
        if sha256(data).hexdigest() != entry["sha256"]:
            raise ValueError(f"Source checksum mismatch: {entry['path']}")
        return data

    # Validate every exported model, including controls and failed designs.
    indexed = list(csv.DictReader((source / "STRUCTURE_INDEX.csv").open()))
    for row in indexed:
        read_structure(verified({"path": row["file"], "sha256": row["sha256"]}))

    def verify_evidence(value):
        if isinstance(value, dict):
            if "path" in value and "sha256" in value:
                verified(value)
            for child in value.values():
                verify_evidence(child)
        elif isinstance(value, list):
            for child in value:
                verify_evidence(child)

    verify_evidence(science)
    evaluations = science["reports"]["focused_validation"]["rows"]
    controls = science["reports"]["controls"]["rows"]
    consumed_evaluations = set()

    def exactly_one(rows):
        if len(rows) != 1:
            raise ValueError("Missing or ambiguous target/candidate/context/seed mapping")
        return rows[0]

    groups = defaultdict(list)
    for candidate in science["native_candidates"]:
        groups[candidate["target"]].append(candidate)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=".science-", dir=destination.parent))
    try:
        shutil.copytree(snapshot, temporary, dirs_exist_ok=True)
        imported = {}

        def save(job_id, filename, data):
            digest = sha256(data).hexdigest()
            (temporary / "files" / digest).write_bytes(data)
            url = f"/api/jobs/{job_id}/artifacts/{filename}"
            manifest["files"][url] = {"sha256": digest, "size": len(data)}
            return url

        archive = BytesIO()
        with ZipFile(archive, "w", ZIP_DEFLATED) as output:
            for path in sorted(source.rglob("*")):
                if path.is_file() and not path.is_symlink():
                    output.write(path, path.relative_to(source))
        for name, candidates in groups.items():
            job_id = sha256(f"jacob-public:{revision}:{name}".encode()).hexdigest()[:32]
            target_id = f"jacob:{revision[:12]}:{name}"
            target_sequence = None
            normalized = []
            seen = set()
            for candidate in candidates:
                cid, sequence = candidate["id"], candidate["sequence"]
                if cid in seen or "/" in cid or ".." in cid:
                    raise ValueError("Duplicate or invalid candidate identity")
                seen.add(cid)
                if sha256(sequence.encode()).hexdigest() != candidate["sequence_sha256"]:
                    raise ValueError("Candidate sequence checksum mismatch")
                binder_data, complex_data = (
                    verified(candidate["binder"]),
                    verified(candidate["complex"]),
                )
                _, binder_chains = read_structure(binder_data)
                complex_model, chains = read_structure(complex_data)
                if list(binder_chains.values()) != [sequence] or len(chains) != 2:
                    raise ValueError("Unexpected native binder or complex composition")
                binder_ids = [chain for chain, seq in chains.items() if seq == sequence]
                if len(binder_ids) != 1:
                    raise ValueError("Ambiguous binder identity in complex")
                target_chain = next(chain for chain in chains if chain != binder_ids[0])
                if target_sequence is not None and chains[target_chain] != target_sequence:
                    raise ValueError("Different target constructs in one native batch")
                target_sequence = chains[target_chain]
                native = candidate["native_filter"]
                if native not in {"pass", "fail", "unknown"}:
                    raise ValueError("Unknown native filter status")
                metrics = candidate["metrics"]
                if any(v is not None and not math.isfinite(v) for v in metrics.values()):
                    raise ValueError("Nonfinite candidate metrics")
                normalized.append(
                    {
                        "id": cid,
                        "sequence": sequence,
                        "native_filter_passed": {"pass": True, "fail": False, "unknown": None}[
                            native
                        ],
                        # The public export does not include named failure reasons. Do not invent them.
                        "failed_filters": [],
                        "confidence_metrics": metrics,
                        "confidence_units": {k: "Å" if "rmsd" in k else "0–1" for k in metrics},
                        "binder_structure_url": save(job_id, f"{cid}/binder.cif", binder_data),
                        "complex_structure_url": save(job_id, f"{cid}/complex.cif", complex_data),
                        "target_structure_url": save(
                            job_id, f"{cid}/target.cif", chain_file(complex_model, target_chain)
                        ),
                    }
                )
                item = normalized[-1]
                item["complex_chains"] = {"target": [target_chain], "binder": binder_ids}
                item["predictions"] = []
                for row in evaluations:
                    if (row.get("target_id"), row.get("candidate_id"), row.get("kind")) != (
                        name,
                        cid,
                        "fixed_binder_complex",
                    ):
                        continue
                    if row["context"] != "target_msa_only_complexes":
                        raise ValueError("Unsupported evaluation context")
                    seed = row["seed"]
                    alone = exactly_one(
                        [
                            r
                            for r in evaluations
                            if (
                                r.get("target_id"),
                                r.get("candidate_id"),
                                r.get("kind"),
                                r.get("context"),
                                r.get("seed"),
                            )
                            == (name, cid, "binder_alone", "binder_alone", seed)
                        ]
                    )
                    control = exactly_one(
                        [
                            r
                            for r in controls
                            if (r["target_id"], r["condition"], r["seed"])
                            == (name, "server_msa__no_template", seed)
                        ]
                    )
                    complex_bytes = verified(row["structure"])
                    key_bytes = verified(alone["structure"])
                    target_bytes = verified(control["structure"])
                    _, evaluated_chains = read_structure(complex_bytes)
                    _, key_chains = read_structure(key_bytes)
                    _, target_chains = read_structure(target_bytes)
                    binder_chain = row["binder_chain"]
                    if (
                        evaluated_chains != {"A": target_sequence, binder_chain: sequence}
                        or list(key_chains.values()) != [sequence]
                        or list(target_chains.values()) != [target_sequence]
                    ):
                        raise ValueError(
                            "Evaluated sequence or chain role differs from native input"
                        )
                    if row["source_native_filter"] != native:
                        raise ValueError("Native filter differs in evaluation")
                    pair = row["binder_target_chainpair_iptm"]
                    confidence = {
                        "target_plddt": row["chain_plddt_0_100"]["A"]["mean"],
                        "binder_plddt": row["binder_plddt_0_100"]["mean"],
                        "target_to_binder_iptm": pair[f"A|{binder_chain}"],
                        "binder_to_target_iptm": pair[f"{binder_chain}|A"],
                        "close_atomic_overlaps": row["binder_target_interfaces"]["A"][
                            "severe_heavy_atom_overlap"
                        ]["atom_pairs_below_cutoff"],
                    }
                    evidence = {"complex": row, "binder_alone": alone, "target_alone": control}
                    consumed_evaluations.update(
                        (row["structure"]["path"], alone["structure"]["path"])
                    )
                    item["predictions"].append(
                        {
                            "id": f"standard-{seed}",
                            "seed": seed,
                            "binder_alone_plddt": alone["binder_plddt_0_100"]["mean"],
                            "target_alone_plddt": control["chains"]["A"]["plddt"]["mean"],
                            "context": "Target MSA · no templates",
                            "confidence_metrics": confidence,
                            "complex_chains": {"target": ["A"], "binder": [binder_chain]},
                            "confidence_units": {
                                "target_plddt": "0–100",
                                "binder_plddt": "0–100",
                                "target_to_binder_iptm": "0–1",
                                "binder_to_target_iptm": "0–1",
                                "close_atomic_overlaps": "atom pairs < 2 Å",
                            },
                            "complex_structure_url": save(
                                job_id, row["structure"]["path"], complex_bytes
                            ),
                            "binder_structure_url": save(
                                job_id, alone["structure"]["path"], key_bytes
                            ),
                            "target_structure_url": save(
                                job_id, control["structure"]["path"], target_bytes
                            ),
                            "evidence_url": save(
                                job_id,
                                f"{cid}/seed-{seed}.json",
                                json.dumps(evidence, allow_nan=False).encode(),
                            ),
                        }
                    )
                seeds = [p["seed"] for p in item["predictions"]]
                if seeds and sorted(seeds) != [11, 22, 33]:
                    raise ValueError("Expected exactly the three declared evaluation seeds")
            if target_sequence is None:
                raise ValueError("Empty native candidate batch")
            zip_url = save(job_id, "scientific-handoff.zip", archive.getvalue())
            download = f"/api/jobs/{job_id}/download"
            manifest["files"][download] = manifest["files"][zip_url]
            imported[job_id] = {
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
                    "description": "Focused design · experimental construct",
                    "source": "Evaluated examples",
                    "source_revision": revision,
                    "amino_acid_sequence": target_sequence,
                    "sequence_sha256": sha256(target_sequence.encode()).hexdigest(),
                    "target_view": "Evaluated: independent MSA fold without templates. Native: target chain extracted from the selected design complex.",
                    "filter_reasons": "Named failure reasons are absent from this public export",
                    "validation": "Standard refolds differ from native design poses. See the complete study in Download.",
                    "limitations": science["limitations"],
                },
            }
        if consumed_evaluations != {r["structure"]["path"] for r in evaluations}:
            raise ValueError(
                "Some evaluated predictions could not be joined to their native candidates"
            )
        manifest["jobs"] = {**imported, **manifest["jobs"]}
        (temporary / "manifest.json").write_text(json.dumps(manifest, allow_nan=False))
        PublishedSnapshot(temporary, manifest["catalog_sha256"])
        temporary.rename(destination)
        return {
            "imported_jobs": len(imported),
            "native_candidates": sum(map(len, groups.values())),
            "validated_source_structures": len(indexed),
            "output": str(destination),
        }
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("source", "snapshot", "output", "revision"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    print(json.dumps(import_science(args.source, args.snapshot, args.output, args.revision)))


if __name__ == "__main__":
    main()
