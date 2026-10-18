"""Offline import of the upstream targets.jsonl/provenance.json package."""

import json
import shutil
from pathlib import Path

import gemmi
from pydantic import BaseModel, Field

from bespoke_bio.api import Protein

from .contracts import TargetState
from .pipeline import save
from .structures import prepare, protein_residues, sequence, sha


class ImportedTargets(BaseModel):
    schema_version: str = "0.1.0"
    target_states: list[TargetState] = Field(default_factory=list)
    records: list[dict] = Field(default_factory=list)
    provenance: dict = Field(default_factory=dict)


def contained_path(root, relative):
    path = (root / relative).resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError("Artifact path escapes source package")
    return path


def declared_sequence(structure, chain):
    subchains = {r.subchain for r in protein_residues(chain)}
    entities = [
        e for e in structure.entities if subchains.intersection(e.subchains) and e.full_sequence
    ]
    if len(entities) == 1:
        return gemmi.one_letter_code(entities[0].full_sequence), "deposited_polymer"
    return sequence(chain), "observed_only"


def import_targets(bundle, output, selections=None):
    root, out = Path(bundle).resolve(), Path(output)
    if out.exists() and any(out.iterdir()):
        raise ValueError("Use an empty output directory")
    out.mkdir(parents=True, exist_ok=True)
    (out / "sources").mkdir()
    (out / "prepared").mkdir()
    provenance = json.loads((root / "provenance.json").read_text())
    records = [
        Protein.model_validate_json(line)
        for line in (root / "targets.jsonl").read_text().splitlines()
        if line.strip()
    ]
    result = ImportedTargets(provenance=provenance)
    selections = selections or {}
    by_id = {}
    seen_ids = set()
    for protein in records:
        if protein.target_id in seen_ids:
            raise ValueError("Duplicate upstream target ID")
        seen_ids.add(protein.target_id)
        item = {
            "upstream_target_id": protein.target_id,
            "sequence_hash": protein.sequence_hash,
            "phase_scenario": protein.phase_scenario,
            "source_variant_ids": protein.source_variant_ids,
            "status": "unresolved",
            "target_state_id": None,
        }
        if protein.sequence_hash != "sha256:" + sha(protein.amino_acid_sequence.encode()):
            item.update(status="failed", reason="Upstream sequence hash does not match sequence")
            result.records.append(item)
            continue
        if (
            protein.structure is None
            or protein.structure.status != "complete"
            or not protein.structure.file
        ):
            item["reason"] = "No resolved structure in upstream package"
            result.records.append(item)
            continue
        try:
            path = contained_path(root, protein.structure.file)
            st = gemmi.read_structure(str(path))
            st.setup_entities()
            if len(st) != 1:
                raise ValueError("Explicit model selection required for multi-model structures")
            choices = {c.name: declared_sequence(st, c) for c in st[0] if protein_residues(c)}
            if protein.target_id in selections:
                selected = selections[protein.target_id]
                if len(selected) != 1:
                    raise ValueError(
                        "One polypeptide chain per upstream protein is required; partner context is separate"
                    )
                chain = selected[0]
                if chain not in choices:
                    raise ValueError("Selected chain absent")
            else:
                matches = [c for c, (s, _) in choices.items() if s == protein.amino_acid_sequence]
                if len(matches) != 1:
                    raise ValueError(
                        f"Expected one exact sequence-matching chain; found {len(matches)}. Supply explicit author-chain selection."
                    )
                chain = matches[0]
            if choices[chain][0] != protein.amino_acid_sequence:
                raise ValueError(
                    "Selected structure chain sequence does not match upstream protein"
                )
            sid = (
                "state_"
                + sha((sha(path.read_bytes()) + chain + protein.sequence_hash).encode())[:16]
            )
            if sid not in by_id:
                suffix = path.suffix or ".cif"
                source = out / "sources" / (sha(path.read_bytes()) + suffix)
                if not source.exists():
                    shutil.copyfile(path, source)
                evidence = "experimental" if protein.structure.source == "RCSB PDB" else "predicted"
                state = prepare(
                    source,
                    out / "prepared" / (sid + ".cif"),
                    [chain],
                    evidence=evidence,
                    context={
                        "upstream_provenance": provenance,
                        "upstream_structure": protein.structure.model_dump(mode="json"),
                        "exact_protein_sequence": protein.amino_acid_sequence,
                        "exact_protein_sequence_sha256": protein.sequence_hash,
                        "sequence_match_basis": choices[chain][1],
                        "reference_sequences": {},
                        "upstream_records": [],
                    },
                )
                state.id = sid
                state.source = str(source.relative_to(out))
                state.prepared = str(Path(state.prepared).relative_to(out))
                # label_seq_id maps deposited polymer positions, not author residue numbering.
                for mapping in state.residue_mapping:
                    position = mapping["label_seq_id"]
                    if (
                        choices[chain][1] == "observed_only"
                        and state.chains[chain] == protein.amino_acid_sequence
                    ):
                        position = mapping["observed_sequence_position"]
                    mapping["protein_sequence_position"] = position
                    if position is not None:
                        letter = gemmi.find_tabulated_residue(
                            mapping["residue_name"]
                        ).one_letter_code.upper()
                        if (
                            not 1 <= position <= len(protein.amino_acid_sequence)
                            or protein.amino_acid_sequence[position - 1] != letter
                        ):
                            raise ValueError(
                                "Deposited residue mapping disagrees with upstream sequence"
                            )
                by_id[sid] = state
                result.target_states.append(state)
            state = by_id[sid]
            state.variants = sorted(set(state.variants + protein.source_variant_ids))
            state.context["upstream_records"].append(protein.model_dump(mode="json"))
            state.context["reference_sequences"][protein.target_id] = (
                protein.reference_amino_acid_sequence
            )
            state.measurements.append(
                {
                    "kind": "measured_RNA",
                    "units": "TPM",
                    "value": protein.expression.tpm,
                    "raw_count": protein.expression.raw_count,
                    "transcript_id": protein.transcript_id,
                    "upstream_target_id": protein.target_id,
                    "protein_abundance": None,
                }
            )
            item.update(status="imported", target_state_id=sid, author_chain=chain)
        except (ValueError, RuntimeError, OSError, KeyError) as exc:
            item.update(status="failed", reason=str(exc))
        result.records.append(item)
    save(out / "targets.json", result)
    (out / "chains.fasta").write_text(
        "".join(
            f">{s.id}:{c}\n{seq}\n" for s in result.target_states for c, seq in s.chains.items()
        )
    )
    save(
        out / "manifest.json",
        {
            "input_hashes": {
                n: sha((root / n).read_bytes()) for n in ["targets.jsonl", "provenance.json"]
            },
            "files": {
                str(p.relative_to(out)): sha(p.read_bytes()) for p in out.rglob("*") if p.is_file()
            },
            "imported_records": sum(r["status"] == "imported" for r in result.records),
            "failed_records": sum(r["status"] == "failed" for r in result.records),
            "unresolved_records": sum(r["status"] == "unresolved" for r in result.records),
        },
    )
    return result
