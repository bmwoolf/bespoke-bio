"""Thin input adapter for the upstream BoltzGen end-to-end protein pipeline."""

import json
import shutil
from pathlib import Path

import gemmi

from .contracts import TargetState
from .pipeline import save
from .structures import protein_residues, sequence, sha, write_chains

COMMIT = "a07fea1153111ef453aaea7a15697ee05de9291f"


def prepare(fixtures, output):
    """Backward-compatible public PD-L1 fixture."""
    return prepare_target(
        Path(fixtures) / "5JDS.cif", output, chain="A", name="pdl1", evidence="experimental"
    )


def prepare_target(
    source,
    output,
    *,
    chain=None,
    name="design",
    binder_length=75,
    num_designs=2,
    evidence="predicted",
):
    """Prepare one observed structural chain or one sequence-only FASTA; no epitope."""
    source, output = Path(source), Path(output)
    if not name or any(
        c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-" for c in name
    ):
        raise ValueError("Name must be a safe identifier")
    if not 25 <= binder_length <= 300 or not 1 <= num_designs <= 10000:
        raise ValueError("Invalid binder length or batch size")
    output.mkdir(parents=True, exist_ok=True)
    saved = output / ("source" + "".join(source.suffixes))
    shutil.copy2(source, saved)
    is_fasta = source.suffix.lower() in {".fa", ".fasta", ".faa"}
    mapping = []
    if is_fasta:
        lines = source.read_text().splitlines()
        if sum(line.startswith(">") for line in lines) != 1:
            raise ValueError("Exactly one FASTA record is required")
        target_seq = "".join(line.strip() for line in lines if not line.startswith(">")).upper()
        if not target_seq or set(target_seq) - set("ARNDCQEGHILKMFPSTWYV"):
            raise ValueError("Target must contain canonical amino acids")
        target_entity = {"protein": {"id": "A", "sequence": target_seq}}
        target_path = ""
        source_chain = "A"
    else:
        structure = gemmi.read_structure(str(source))
        if len(structure) != 1:
            raise ValueError("Explicit single-model input is required")
        available = [c.name for c in structure[0] if protein_residues(c)]
        if chain is None and len(available) != 1:
            raise ValueError("Select an author chain explicitly for multichain structures")
        source_chain = chain or available[0]
        if source_chain not in available:
            raise ValueError("Requested target chain is absent")
        structure.remove_alternative_conformations()
        selected = structure[0][source_chain].clone()
        normalized = gemmi.Chain("A")
        for index, residue in enumerate(protein_residues(selected), 1):
            copy = residue.clone()
            copy.seqid = gemmi.SeqId(index, " ")
            copy.label_seq = index
            copy.subchain = "A"
            normalized.add_residue(copy)
            mapping.append(
                {
                    "input_number": index,
                    "source_chain": source_chain,
                    "source_author_number": residue.seqid.num,
                    "source_insertion": residue.seqid.icode.strip(),
                    "source_label_seq_id": residue.label_seq,
                    "residue_name": residue.name,
                }
            )
        target_seq = sequence(normalized)
        if not target_seq or set(target_seq) - set("ARNDCQEGHILKMFPSTWYV"):
            raise ValueError("Observed target contains unsupported amino acids")
        target = write_chains([normalized], output / "target.cif")
        for entity in target.entities:
            entity.full_sequence = [r.name for r in normalized]
        target.make_mmcif_document().write_file(str(output / "target.cif"))
        target_path = "target.cif"
        target_entity = {"file": {"path": target_path, "include": [{"chain": {"id": "A"}}]}}
    state = TargetState(
        id="target_" + sha((sha(source.read_bytes()) + source_chain).encode())[:16],
        source=saved.name,
        prepared=target_path,
        source_sha256=sha(source.read_bytes()),
        chains={"A": target_seq},
        residue_mapping=mapping,
        unresolved=[],
        evidence=evidence,
        context={
            "input_mode": "sequence_only" if is_fasta else "observed_structure",
            "source_author_chain": source_chain,
            "coordinate_normalization": "observed residues only; source retained",
            "sequence_only_prediction_performed": False,
        },
    )
    save(output / "target_state.json", state.model_dump())
    save(output / "target_residue_mapping.json", mapping)
    spec = {"entities": [{"protein": {"id": "B", "sequence": binder_length}}, target_entity]}
    # JSON is a strict YAML subset and avoids another project dependency.
    (output / f"{name}.yaml").write_text(json.dumps(spec, indent=2) + "\n")
    plan = {
        "engine": "boltzgen",
        "commit": COMMIT,
        "protocol": "protein-anything",
        "status": "prepared_not_executed",
        "target": str(source),
        "target_state_id": state.id,
        "target_sequence_sha256": sha(target_seq.encode()),
        "target_input_mode": state.context["input_mode"],
        "target_residues": len(target_seq),
        "site_selection": {
            "method": "unconstrained",
            "specified_epitope": None,
            "automatic_hotspot_detection": False,
        },
        "binder_residues": binder_length,
        "num_designs": num_designs,
        "inverse_fold_num_sequences": 1,
        "budget": num_designs,
        "note": "budget counts final candidates, not dollars; cloud spend cap is separate",
        "cloud_spend_cap_usd": 100,
        "interleaving": False,
        "stages": [
            "design",
            "inverse_folding",
            "folding",
            "design_folding",
            "esmfold2_scoring",
            "analysis",
            "filtering",
        ],
        "target_sha256": sha((output / target_path).read_bytes())
        if target_path
        else state.source_sha256,
        "command": [
            "boltzgen",
            "run",
            f"{name}.yaml",
            "--output",
            "results",
            "--protocol",
            "protein-anything",
            "--num_designs",
            str(num_designs),
            "--budget",
            str(num_designs),
            "--diffusion_batch_size",
            "1",
            "--inverse_fold_num_sequences",
            "1",
            "--num_workers",
            "0",
            "--devices",
            "1",
            "--use_kernels",
            "false",
        ],
        "acceptance": "Pipeline completion is distinct from passing candidate filters or experimental binding.",
    }
    save(output / "boltzgen_plan.json", plan)
    return plan


def summarize(results, output):
    """Retain upstream pass/fail flags; never substitute ranking for acceptance."""
    import csv

    results = Path(results)
    metrics = results / "final_ranked_designs/all_designs_metrics.csv"
    if not metrics.exists():
        raise FileNotFoundError(
            "Full filtering output is missing; pipeline cannot be reported complete"
        )
    with metrics.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    candidates = []
    selected_metrics = [
        "design_ptm",
        "design_to_target_iptm",
        "filter_rmsd",
        "filter_rmsd_design",
        "designfolding-filter_rmsd",
        "esmfold2_ipsae_min",
        "quality_score",
    ]
    for row in rows:
        pass_columns = {
            k: v for k, v in row.items() if k.startswith("pass_") and k.endswith("_filter")
        }
        overall = row.get("pass_filters", "").lower()
        if overall not in {"true", "false"}:
            raise ValueError("Missing explicit overall filter result")
        candidates.append(
            {
                "id": row["id"],
                "file_name": row["file_name"],
                "passed": overall == "true",
                "failed_filters": [k for k, v in pass_columns.items() if v.lower() != "true"],
                "upstream_filter_flags": pass_columns,
                "sequence": row.get("designed_chain_sequence") or row.get("designed_sequence"),
                "metrics": {k: row[k] for k in selected_metrics if k in row},
            }
        )
    report = {
        "status": "pipeline_completed",
        "total": len(candidates),
        "passed": sum(c["passed"] for c in candidates),
        "candidates": candidates,
        "evidence": "computational pipeline filters only; no measured binding",
    }
    save(output, report)
    return report
