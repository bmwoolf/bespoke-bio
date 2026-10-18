"""Local motif/topology compiler. Does not run or emulate a generative model."""

import json
import math
from itertools import combinations, pairwise, permutations
from pathlib import Path

import gemmi
import numpy as np

from .contracts import MotifPlanRequest
from .geometry import distances
from .pipeline import save
from .placement import protein_only
from .structures import atoms, protein_residues, ref, sha, write_chains

FOUNDRY_COMMIT = "d078b24d6d0e6f8607e9a0311128f4a897f83953"


def native_disulfides(chain):
    cyst = [
        r for r in protein_residues(chain) if r.name == "CYS" and r.find_atom("SG", "*") is not None
    ]
    return [
        {
            "first": ref(chain.name, a).model_dump(),
            "second": ref(chain.name, b).model_dump(),
            "evidence": "inferred from observed SG geometry, not imported bond annotations",
            "sg_distance_angstrom": a.find_atom("SG", "*").pos.dist(b.find_atom("SG", "*").pos),
        }
        for a, b in combinations(cyst, 2)
        if 1.8 <= a.find_atom("SG", "*").pos.dist(b.find_atom("SG", "*").pos) <= 2.3
    ]


def selected_blocks(chain, contact_refs, flank, preserve_disulfides=True):
    residues = protein_residues(chain)
    lookup = {(r.seqid.num, r.seqid.icode.strip()): i for i, r in enumerate(residues)}
    contacts = set()
    for reference in contact_refs:
        if reference["chain"] != chain.name:
            raise ValueError("Contact-chain mismatch")
        k = (reference["number"], reference["insertion"])
        if k not in lookup:
            raise ValueError(f"Unobserved contact residue {k}")
        contacts.add(lookup[k])
    ss = native_disulfides(chain)
    centers = set(contacts)
    if preserve_disulfides:
        for bond in ss:
            for end in ["first", "second"]:
                r = bond[end]
                centers.add(lookup[(r["number"], r["insertion"])])
    selected = {
        i for j in centers for i in range(max(0, j - flank), min(len(residues), j + flank + 1))
    }
    blocks = []
    for i in sorted(selected):
        # Never join across an unresolved peptide segment in a source structure.
        continuous = False
        if blocks and i == blocks[-1][-1] + 1:
            previous = residues[i - 1].find_atom("C", "*")
            current = residues[i].find_atom("N", "*")
            continuous = (
                previous is not None
                and current is not None
                and 1.1 <= previous.pos.dist(current.pos) <= 1.8
            )
        if continuous:
            blocks[-1].append(i)
        else:
            blocks.append([i])
    if not blocks or len(blocks) > 4:
        raise ValueError("This bounded prototype supports 1..4 coherent motif blocks per module")
    return blocks, ss, contacts


def topology_orders(blocks):
    if len(blocks) > 6:
        raise ValueError("Topology enumeration capped at six total blocks")
    ids = [b["id"] for b in blocks]
    for order in permutations(ids):
        # Preserve sequence order within each module; change only how modules interleave.
        if all(
            [i for i in order if i.startswith(chain)]
            == [b["id"] for b in blocks if b["chain"] == chain]
            for chain in {b["chain"] for b in blocks}
        ):
            yield order


def plan_topology(order, blocks, request, target_tokens, matched_length=None):
    by_id = {b["id"]: b for b in blocks}
    bridges = []
    reasons = []
    for left, right in pairwise(order):
        gap = float(np.linalg.norm(np.array(by_id[left]["end_c"]) - by_id[right]["start_n"]))
        necessary = max(1, math.ceil(gap / 3.8) - 1)
        chosen = max(request.bridge_min_residues, necessary + request.bridge_slack_residues)
        if chosen > request.bridge_max_residues:
            reasons.append(f"{left}->{right}:bridge_budget")
        bridges.append(
            {
                "from": left,
                "to": right,
                "anchor_distance_angstrom": gap,
                "necessary_contour_minimum_residues": necessary,
                "proposed_residues": chosen,
                "status": "unmodeled",
            }
        )
    length = (
        sum(b["count"] for b in blocks)
        + 2 * request.terminal_residues
        + sum(b["proposed_residues"] for b in bridges)
    )
    if matched_length is not None:
        if matched_length < length:
            reasons.append("target_length_below_necessary_budget")
        extra = matched_length - length
        while extra > 0:
            available = [b for b in bridges if b["proposed_residues"] < request.bridge_max_residues]
            if not available:
                reasons.append("matched_length_bridge_budget")
                break
            min(available, key=lambda b: b["proposed_residues"])["proposed_residues"] += 1
            extra -= 1
            length += 1
    if length > request.maximum_scaffold_length:
        reasons.append("scaffold_length_budget")
    tokens = [str(request.terminal_residues)]
    spans = []
    position = request.terminal_residues + 1
    for i, bid in enumerate(order):
        block = by_id[bid]
        tokens.append(block["token"])
        spans.append(
            {
                "block_id": bid,
                "output_chain": "A",
                "output_first": position,
                "output_last": position + block["count"] - 1,
                "input_chain": block["chain"],
                "input_first": block["input_first"],
                "input_last": block["input_last"],
            }
        )
        position += block["count"]
        if i < len(bridges):
            n = bridges[i]["proposed_residues"]
            tokens.append(str(n))
            position += n
    tokens.append(str(request.terminal_residues))
    for token in target_tokens:
        tokens.extend(["/0", token])
    runs = 1 + sum(a[0] != b[0] for a, b in pairwise(order))
    return {
        "id": "topology_" + "_".join(order),
        "order": list(order),
        "provenance_order": "interleaved" if runs > 2 else "tandem_control",
        "status": "rejected" if reasons else "locally_preflighted",
        "failed_checks": reasons,
        "bridges": bridges,
        "scaffold_length": length,
        "contig": ",".join(tokens),
        "output_motif_spans": spans,
        "note": "Sequence provenance interleaves modules; spatial interdigitation and a shared folded core are not established",
    }


def prepare_plan(assembly, output, request=None):
    request = request or MotifPlanRequest()
    assembly, out = Path(assembly), Path(output)
    if out.exists() and any(out.iterdir()):
        raise ValueError("Use an empty output directory")
    result = json.loads((assembly / "placement_search.json").read_text())
    if not result["top_layout_ids"]:
        raise ValueError("No geometry-screened layout to prepare")
    top = next(r for r in result["layouts"] if r["id"] == result["top_layout_ids"][0])
    st = gemmi.read_structure(str(assembly / top["artifact"]))
    if len(st) != 1:
        raise ValueError("Explicit single model required")
    chains = {c.name: protein_only(c) for c in st[0]}
    for name in ["K", "L", "T", "U"]:
        if name not in chains:
            raise ValueError("Expected benchmark chains K,L,T,U")
    out.mkdir(parents=True, exist_ok=True)
    mapping = []
    blocks = []
    disulfides = []
    motif_chains = []
    contact_refs = []
    normalized = {}
    for name in ["K", "L", "T", "U"]:
        chain = chains[name]
        residues = protein_residues(chain)
        norm = chain.clone()
        for i, r in enumerate(norm, 1):
            mapping.append(
                {
                    "input_chain": name,
                    "input_number": i,
                    "source": ref(name, residues[i - 1]).model_dump(),
                    "residue_name": r.name,
                }
            )
            r.seqid = gemmi.SeqId(i, " ")
            r.label_seq = i
            r.subchain = name
        normalized[name] = norm
    for name, target_name in [("K", "T"), ("L", "U")]:
        chain = chains[name]
        partner = chains[target_name]
        target_xyz = atoms(partner)[0]
        refs = []
        for residue in protein_residues(chain):
            one = gemmi.Chain(name)
            one.add_residue(residue.clone())
            if distances(atoms(one)[0], target_xyz).min() < 5.0:
                refs.append(ref(name, residue).model_dump())
        contact_refs.extend(refs)
        raw_blocks, ss, contacts = selected_blocks(
            chain, refs, request.flank_residues, request.preserve_native_disulfides
        )
        residues = protein_residues(chain)
        for n, indices in enumerate(raw_blocks, 1):
            block_chain = gemmi.Chain(name + str(n))
            for i in indices:
                block_chain.add_residue(residues[i].clone())
            motif_chains.append(block_chain)
            start = residues[indices[0]].find_atom("N", "*")
            end = residues[indices[-1]].find_atom("C", "*")
            if start is None or end is None:
                raise ValueError("Motif boundary backbone atoms missing")
            blocks.append(
                {
                    "id": name + str(n),
                    "chain": name,
                    "token": f"{name}{indices[0] + 1}-{indices[-1] + 1}",
                    "input_first": indices[0] + 1,
                    "input_last": indices[-1] + 1,
                    "count": len(indices),
                    "source_residues": [ref(name, residues[i]).model_dump() for i in indices],
                    "contact_source_positions": [i + 1 for i in indices if i in contacts],
                    "start_n": start.pos.tolist(),
                    "end_c": end.pos.tolist(),
                    "fixed_atoms": "all observed atoms; no coordinates fabricated",
                    "support_sufficiency": "unknown; flanks and disulfides do not establish an independently folding motif",
                }
            )
        disulfides.extend(ss)
    write_chains(list(normalized.values()), out / "model_input.cif")
    write_chains(motif_chains, out / "motifs.cif")
    save(out / "residue_mapping.json", mapping)
    target_tokens = [f"{c}1-{len(normalized[c])}" for c in ["T", "U"]]
    topologies = [plan_topology(o, blocks, request, target_tokens) for o in topology_orders(blocks)]
    viable = sorted(
        [t for t in topologies if t["status"] == "locally_preflighted"],
        key=lambda t: (t["scaffold_length"], t["id"]),
    )
    control = next((t for t in viable if t["provenance_order"] == "tandem_control"), None)
    selected = [control] if control is not None else []
    matched_length = request.target_scaffold_length or max(
        (t["scaffold_length"] for t in selected), default=None
    )
    selected = [
        plan_topology(t["order"], blocks, request, target_tokens, matched_length) for t in selected
    ]
    selected_by_id = {t["id"]: t for t in selected}
    topologies = [selected_by_id.get(t["id"], t) for t in topologies]
    selected = [t for t in selected if t["status"] == "locally_preflighted"]
    specs = {
        t["id"]: {
            "input": "model_input.cif",
            "contig": t["contig"],
            "select_fixed_atoms": True,
            "dialect": 2,
        }
        for t in selected
    }
    save(out / "rfd3_inputs.draft.json", specs)
    report = {
        "schema_version": "0.1.0",
        "request": request.model_dump(),
        "source_layout_id": top["id"],
        "source_sha256": sha((assembly / top["artifact"]).read_bytes()),
        "foundry_commit": FOUNDRY_COMMIT,
        "model_validation": "draft: own geometry/mapping checks only; native RFD3 parsing and inference not executed",
        "evidence": "geometry_only",
        "contact_residues": contact_refs,
        "blocks": blocks,
        "native_disulfides": disulfides,
        "disulfide_policy": "endpoints and flanks preserved in fixed motifs"
        if request.preserve_native_disulfides
        else "explicitly not preserved",
        "topologies": topologies,
        "locally_viable_topologies": len(viable),
        "matched_comparison_scaffold_length": matched_length,
        "selected_topology_ids": [t["id"] for t in selected],
        "original_two_key_observed_length": sum(len(normalized[c]) for c in ["K", "L"]),
        "fixed_motif_residues": sum(b["count"] for b in blocks),
        "target_chain_order": ["T", "U"],
        "predicted_output_chain_order": [
            "A:designed_scaffold",
            "B:target_from_T",
            "C:target_from_U",
        ],
        "unknowns": [
            "Folding core and framework support not established",
            "Generated bridges and terminal regions not modeled",
            "Disulfide pairing/folding after redesign unvalidated",
            "Motif order does not establish spatial interdigitation",
            "Relative target placement is hypothetical",
            "No affinity, selectivity or function computed",
        ],
        "gpu_acceptance": [
            "Verify native input mapping and exact protected atom identities before generating a batch",
            "Check motif geometry and each separate interface; reject missing or scrambled motifs",
            "Redesign generated regions with inverse folding, retaining protected identities",
            "Independently predict the discrete scaffold alone and the two-target complex",
            "Check C-N bonds across every motif/bridge boundary and disulfide geometry",
            "Reassess compactness, forbidden clashes and each chain-pair metric after reprediction",
        ],
    }
    save(out / "motif_plan.json", report)
    save(out / "request.json", request)
    save(out / "request.schema.json", MotifPlanRequest.model_json_schema())
    # Concrete run command for the future native GPU environment; not executed locally.
    (out / "gpu_smoke_test.sh").write_text("""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
# Activate an isolated environment from the pinned Foundry commit first.
# Check native parsing and fixed-residue identity before trusting any model output.
rfd3 design inputs=rfd3_inputs.draft.json out_dir=gpu_smoke prevalidate_inputs=true diffusion_batch_size=1 n_batches=1 inference_mode=legacy
""")
    lines = [
        "# Motif scaffold preflight",
        "",
        "Local geometry and topology preparation only. No backbone, sequence or binding prediction has run.",
        "",
        f"Protected blocks: {len(blocks)}; retained motif residues: {report['fixed_motif_residues']}; local topology options: {len(topologies)}; within specified budgets: {len(viable)}.",
        "",
        "| Topology | Provenance order | Scaffold residues | Local status |",
        "|---|---|---:|---|",
    ]
    lines += [
        f"| {t['id']} | {t['provenance_order']} | {t['scaffold_length']} | {t['status']} |"
        for t in topologies
    ]
    lines += [
        "",
        "Interleaved generation is retired; only a tandem topology is selected. Motifs are preserved as contiguous indexed blocks. Full native parser validation remains pending.",
        "",
        "A shorter planned chain is not a demonstrated smaller or functional protein. Framework removal remains an explicit scientific risk.",
    ]
    (out / "REVIEW.md").write_text("\n".join(lines) + "\n")
    save(
        out / "manifest.json",
        {
            "files": {p.name: sha(p.read_bytes()) for p in out.iterdir() if p.is_file()},
            "foundry_commit": FOUNDRY_COMMIT,
            "evidence": "geometry_only",
        },
    )
    return report


def check_motif_geometry(plan_dir, candidate_path, topology_id, rmsd_tolerance=0.5):
    """CPU post-check for a future model output; missing coordinates never pass.

    Uses ONE joint target-frame fit. Does not independently fit motifs, which
    would hide movement relative to a required partner. No binding score.
    """
    from .geometry import apply, fit

    plan_dir = Path(plan_dir)
    plan = json.loads((plan_dir / "motif_plan.json").read_text())
    topology = next(t for t in plan["topologies"] if t["id"] == topology_id)
    reference = gemmi.read_structure(str(plan_dir / "model_input.cif"))
    candidate = gemmi.read_structure(str(candidate_path))
    if len(candidate) != 1:
        raise ValueError("Explicit single model required")
    model = candidate[0]
    if not {"A", "B", "C"} <= {c.name for c in model}:
        return {
            "status": "missing",
            "reason": "Expected scaffold A and target chains B,C",
            "evidence": "geometry_only",
        }
    expected_targets = []
    observed_targets = []
    for src, dst in [("T", "B"), ("U", "C")]:
        source_res = protein_residues(reference[0][src])
        output_res = protein_residues(model[dst])
        if len(source_res) != len(output_res) or any(
            a.name != b.name for a, b in zip(source_res, output_res)
        ):
            return {
                "status": "fail",
                "reason": "Target identity/length mismatch",
                "evidence": "geometry_only",
            }
        for a, b in zip(source_res, output_res):
            aa, bb = a.find_atom("CA", "*"), b.find_atom("CA", "*")
            if aa is None or bb is None:
                return {
                    "status": "missing",
                    "reason": "Target Cα atom missing",
                    "evidence": "geometry_only",
                }
            expected_targets.append(aa.pos.tolist())
            observed_targets.append(bb.pos.tolist())
    tr, target_rmsd = fit(observed_targets, expected_targets)
    output_by_number = {(r.seqid.num, r.seqid.icode.strip()): r for r in model["A"]}
    findings = []
    for span in topology["output_motif_spans"]:
        expected = []
        observed = []
        missing = []
        mismatches = []
        for source_number, output_number in zip(
            range(span["input_first"], span["input_last"] + 1),
            range(span["output_first"], span["output_last"] + 1),
            strict=True,
        ):
            source = reference[0][span["input_chain"]][source_number - 1]
            residue = output_by_number.get((output_number, ""))
            if residue is None:
                missing.append(f"A{output_number}")
                continue
            if residue.name != source.name:
                mismatches.append(f"A{output_number}:{residue.name}!={source.name}")
            for atom in source:
                if atom.element.is_hydrogen:
                    continue
                other = residue.find_atom(atom.name, "*")
                if other is None:
                    missing.append(f"A{output_number}:{atom.name}")
                    continue
                expected.append(atom.pos.tolist())
                observed.append(other.pos.tolist())
        rmsd = (
            None
            if missing or mismatches or not expected
            else float(np.sqrt(np.mean(np.sum((apply(observed, tr) - expected) ** 2, axis=1))))
        )
        maximum = (
            None
            if missing or mismatches or not expected
            else float(np.linalg.norm(apply(observed, tr) - expected, axis=1).max())
        )
        findings.append(
            {
                "block_id": span["block_id"],
                "source_interface": span["input_chain"]
                + ":"
                + ("T" if span["input_chain"] == "K" else "U"),
                "status": "missing"
                if missing
                else (
                    "fail"
                    if mismatches
                    or rmsd is None
                    or rmsd > rmsd_tolerance
                    or maximum > rmsd_tolerance
                    else "pass"
                ),
                "fixed_atom_rmsd_angstrom": rmsd,
                "maximum_fixed_atom_displacement_angstrom": maximum,
                "missing": missing,
                "identity_mismatches": mismatches,
            }
        )
    failed = target_rmsd > rmsd_tolerance or any(f["status"] != "pass" for f in findings)
    return {
        "status": "fail" if failed else "pass",
        "evidence": "geometry_only",
        "joint_target_frame_ca_rmsd_angstrom": target_rmsd,
        "fixed_atom_rmsd_tolerance_angstrom": rmsd_tolerance,
        "blocks": findings,
        "scope": "Motif identity and target-relative geometry only; not folded scaffold acceptance",
        "remaining": [
            "all generated residues and peptide bonds",
            "full-atom inventory",
            "disulfide pairing",
            "binding/function",
            "independent prediction",
        ],
    }
