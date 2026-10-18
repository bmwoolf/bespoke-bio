"""Prepare the ordinary single-target folding/design baseline; no inference implied."""

from pathlib import Path

import gemmi

from .pipeline import save
from .structures import protein_residues, sequence, sha, write_chains


def prepare(fixtures, output):
    fixtures, output = Path(fixtures), Path(output)
    output.mkdir(parents=True, exist_ok=True)
    source = fixtures / "5JDS.cif"
    structure = gemmi.read_structure(str(source))
    target, binder = structure[0]["A"], structure[0]["B"]
    target_seq, binder_seq = sequence(target), sequence(binder)
    normalized = gemmi.Chain("A")
    mapping = []
    for index, residue in enumerate(protein_residues(target), 1):
        copy = residue.clone()
        copy.seqid = gemmi.SeqId(index, " ")
        copy.label_seq = index
        copy.subchain = "A"
        normalized.add_residue(copy)
        mapping.append(
            {
                "input_number": index,
                "source_author_number": residue.seqid.num,
                "source_insertion": residue.seqid.icode.strip(),
            }
        )
    write_chains([normalized], output / "target.cif")
    write_chains([binder], output / "fold_reference.cif")
    save(output / "target_residue_mapping.json", mapping)
    fold = [
        {"name": "kn035_observed_construct", "components": [{"seq": binder_seq, "chain_id": "A"}]}
    ]
    save(output / "fold_inputs.json", fold)
    # Canonical, motif-free binder generation. Only the target is fixed.
    design = {
        "pdl1_baseline": {
            "input": "target.cif",
            "contig": f"75,/0,A1-{len(target_seq)}",
            "dialect": 2,
            "select_fixed_atoms": True,
            "select_hotspots": "A37,A41,A44",
            "infer_ori_strategy": "hotspots",
        }
    }
    save(output / "binder_inputs.json", design)
    plan = {
        "status": "prepared_not_executed",
        "interleaving": False,
        "assembly_redesign": False,
        "benchmark": "5JDS observed KN035/PD-L1 constructs; observed-only truncations",
        "source_sha256": sha(source.read_bytes()),
        "fold_control": {
            "sequence": binder_seq,
            "length": len(binder_seq),
            "sequence_sha256": sha(binder_seq.encode()),
            "template_free": True,
            "msa": "single_sequence",
            "reference_chain": "B",
        },
        "target": {
            "sequence": target_seq,
            "length": len(target_seq),
            "author_chain": "A",
            "hotspots": [37, 41, 44],
        },
        "stages": [
            "RF3 sequence-only folding control",
            "RFD3 ordinary single-target binder generation",
            "ProteinMPNN sequence design with target fixed",
            "RF3 sequence-only binder refolding",
            "RF3 target+binder cofolding",
            "geometry/confidence review and bounded next batch",
        ],
        "round_limits": {"backbones": 2, "sequences_per_backbone": 2, "maximum_rounds": 1},
        "gate_policy": "Missing/early-stopped predictions fail readiness. Confidence and structural recovery must be reported separately; neither establishes binding. No automatic score-driven iteration until this baseline is calibrated.",
        "remaining": [
            "GPU folding control",
            "ProteinMPNN adapter/output ingestion",
            "RF3 output ingestion and calibrated selection",
            "complete design/refold/cofold execution",
        ],
    }
    save(output / "plan.json", plan)
    (output / "fold_control.sh").write_text(
        '#!/usr/bin/env bash\nset -euo pipefail\ncd "$(dirname "$0")"\nrf3 fold inputs=fold_inputs.json out_dir=fold_control\n'
    )
    return plan
