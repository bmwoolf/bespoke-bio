"""CPU-only compact rigid assembly search with explicit geometry and evidence."""

import json
import platform
import time
from pathlib import Path

import gemmi
import numpy as np

from .contracts import DesignResponse, PlacementSearchRequest, Transform
from .geometry import apply, distances, pair_metrics
from .pipeline import renamed, save
from .retrieval import move
from .structures import atoms, ports, protein_residues, sequence, sha, write_chains


def radii(chain):
    return np.array([a.element.vdw_r for r in chain for a in r if not a.element.is_hydrogen])


def cross_screen(a, b, ra, rb, tolerance):
    """Element-aware overlaps. Negative gap = overlap of tabulated vdW spheres."""
    overlap_count = 0
    contact_count = 0
    minimum_gap = float("inf")
    minimum_distance = float("inf")
    for offset in range(0, len(a), 128):
        d = distances(a[offset : offset + 128], b)
        gap = d - ra[offset : offset + 128, None] - rb[None, :]
        overlap_count += int((gap < -tolerance).sum())
        contact_count += int((d < 5.0).sum())
        minimum_gap = min(minimum_gap, float(gap.min()))
        minimum_distance = min(minimum_distance, float(d.min()))
    return {
        "vdw_overlap_pairs": overlap_count,
        "minimum_vdw_gap_angstrom": minimum_gap,
        "minimum_heavy_atom_distance_angstrom": minimum_distance,
        "contact_atom_pairs_5A": contact_count,
    }


def protein_only(chain):
    result = gemmi.Chain(chain.name)
    for r in protein_residues(chain):
        result.add_residue(r.clone())
    if not len(result):
        raise ValueError("No protein residues")
    return result


def layout_transforms(request, center):
    yield "overlap_control", Transform()
    for separation in request.offsets_angstrom:
        for axis in range(3):
            for sign in [1, -1]:
                delta = np.zeros(3)
                delta[axis] = sign * separation
                for angle in request.z_rotations_degrees:
                    radians = np.deg2rad(angle)
                    c, s = np.cos(radians), np.sin(radians)
                    rotation = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
                    translation = center - rotation @ center + delta
                    yield (
                        "sampled",
                        Transform(rotation=rotation.tolist(), translation=translation.tolist()),
                    )


def evaluate_layout(key, target, transform, request):
    k, _ = atoms(key)
    t, _ = atoms(target)
    k2, t2 = apply(k, transform), apply(t, transform)
    kr, tr = radii(key), radii(target)
    forbidden = {
        name: cross_screen(a, b, ra, rb, request.vdw_overlap_tolerance_angstrom)
        for name, a, b, ra, rb in [
            ("key1:key2", k, k2, kr, kr),
            ("target1:target2", t, t2, tr, tr),
            ("key1:target2", k, t2, kr, tr),
            ("key2:target1", k2, t, kr, tr),
        ]
    }
    # Own partners are always co-moved: preserve the full observed binding pose.
    ca = atoms(key, ca=True)[0]
    all_ca = np.concatenate([ca, apply(ca, transform)])
    rg = float(np.sqrt(np.mean(np.sum((all_ca - all_ca.mean(0)) ** 2, axis=1))))
    res = protein_residues(key)
    anchor_c = np.array(res[-1].find_atom("C", "*").pos.tolist())
    anchor_n = np.array(res[0].find_atom("N", "*").pos.tolist())
    if request.attachment == "C":
        end1, end2 = anchor_c, apply([anchor_n], transform)[0]
    else:
        end1, end2 = apply([anchor_c], transform)[0], anchor_n
    gap = float(np.linalg.norm(end2 - end1))
    contour = 3.8 * (len(request.linker_sequence) + 1)
    failed = [name for name, m in forbidden.items() if m["vdw_overlap_pairs"]]
    if gap > contour:
        failed.append("linker_contour_upper_bound")
    return {
        "status": "rejected" if failed else "passes_geometry_screen",
        "failed_checks": failed,
        "transform": transform.model_dump(),
        "forbidden_pairs": forbidden,
        "compactness": {
            "key_ca_radius_of_gyration_angstrom": rg,
            "key_ca_bounding_box_extents_angstrom": np.ptp(all_ca, axis=0).tolist(),
        },
        "connection": {
            "anchor_distance_angstrom": gap,
            "contour_upper_bound_angstrom": contour,
            "extension_fraction": gap / contour,
            "status": "unmodeled",
            "linker_sequence": request.linker_sequence,
        },
        "required_interfaces": {
            "key1:target1": {"status": "pose_preserved", "metrics": pair_metrics(k, t)},
            "key2:target2": {"status": "pose_preserved", "metrics": pair_metrics(k2, t2)},
        },
        "confidence": None,
        "function": None,
    }


def run_search(bundle, output, request=None):
    request = request or PlacementSearchRequest()
    bundle, out = Path(bundle), Path(output)
    if out.exists() and any(out.iterdir()):
        raise ValueError("Use an empty output directory")
    response = DesignResponse.model_validate_json((bundle / "response.json").read_text())
    source_component = response.components[0]
    state = response.target_states[0]
    original_key = gemmi.read_structure(str(bundle / source_component.structure))[0][0].clone()
    key = protein_only(original_key)
    target = gemmi.read_structure(str(bundle / state.prepared))[0][0].clone()
    source_ports = ports(key)
    if any(p.frame is None for p in source_ports):
        raise ValueError("Resolved terminal backbone anchors required")
    out.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    write_chains([key], out / "key.cif")
    write_chains([target], out / "target.cif")
    cid = "component_" + sha((sequence(key) + sha((out / "key.cif").read_bytes())).encode())[:16]
    center = atoms(key)[0].mean(0)
    rows = []
    for label, transform in layout_transforms(request, center):
        record = evaluate_layout(key, target, transform, request)
        record["id"] = (
            "layout_"
            + sha(
                json.dumps(
                    {
                        "component": cid,
                        "state": state.id,
                        "transform": transform.model_dump(),
                        "request": request.model_dump(),
                    },
                    sort_keys=True,
                ).encode()
            )[:16]
        )
        record["kind"] = label
        rows.append(record)
    eligible = sorted(
        [r for r in rows if r["status"] == "passes_geometry_screen"],
        key=lambda r: (
            r["compactness"]["key_ca_radius_of_gyration_angstrom"],
            r["connection"]["anchor_distance_angstrom"],
            r["id"],
        ),
    )
    for rank, row in enumerate(eligible[: request.top_k], 1):
        tr = Transform(**row["transform"])
        file = f"compact_{rank}.cif"
        write_chains(
            [
                renamed(target, "T"),
                renamed(key, "K"),
                renamed(move(target, tr), "U"),
                renamed(move(key, tr), "L"),
            ],
            out / file,
        )
        row["artifact"] = file
        row["instances"] = [
            {
                "id": "key1",
                "component_id": cid,
                "chain": "K",
                "transform": Transform().model_dump(),
            },
            {"id": "key2", "component_id": cid, "chain": "L", "transform": tr.model_dump()},
        ]
        row["target_placements"] = [
            {
                "id": "target1",
                "target_state_id": state.id,
                "chain": "T",
                "transform": Transform().model_dump(),
            },
            {
                "id": "target2",
                "target_state_id": state.id,
                "chain": "U",
                "transform": tr.model_dump(),
            },
        ]
    # Include the spread-out reference used for an honest same-atoms compactness comparison.
    reference_transform = Transform(translation=[0.0, 100.0, 0.0])
    reference = evaluate_layout(key, target, reference_transform, request)
    reference["kind"] = "spread_out_reference"
    reference["artifact"] = "spread_out_reference.cif"
    write_chains(
        [
            renamed(target, "T"),
            renamed(key, "K"),
            renamed(move(target, reference_transform), "U"),
            renamed(move(key, reference_transform), "L"),
        ],
        out / reference["artifact"],
    )
    result = {
        "schema_version": "0.1.0",
        "request": request.model_dump(),
        "component_id": cid,
        "source_component_id": source_component.id,
        "component_revision_note": "Protein-only derivative: crystallization nonprotein atoms omitted; sequence and observed polypeptide coordinates preserved",
        "target_state_id": state.id,
        "sampled_layouts": len(rows) - 1,
        "controls": 1,
        "eligible_layouts": len(eligible),
        "ranking": "Filter forbidden vdW overlaps and necessary linker reach; ascending key Cα radius of gyration, then anchor distance, then ID",
        "overlap_rule": "distance < sum of Gemmi element vdW radii minus overlap tolerance; not a force-field energy",
        "context": "Selected full target author chain including its nonprotein context; no whole-biological-assembly or membrane expansion",
        "relative_target_pose": "hypothesized: each target co-moves with its key; not compatible-site fitting on a prepositioned target assembly",
        "evidence": "geometry_only",
        "backbone_or_sequence_redesign": False,
        "conformational_sampling": False,
        "interdigitation": "not established: intact domains remain rigid; shared-scaffold integration requires redesign",
        "unmodeled": [
            "linker",
            "fusion folding",
            "independent complex prediction",
            "alternate-state assembly",
            "binding and function",
        ],
        "top_layout_ids": [r["id"] for r in eligible[: request.top_k]],
        "layouts": rows,
        "spread_out_reference": reference,
        "elapsed_seconds": time.monotonic() - start,
        "device": "CPU",
        "python": platform.python_version(),
    }
    save(out / "request.json", request)
    save(out / "placement_search.json", result)
    save(out / "request.schema.json", PlacementSearchRequest.model_json_schema())
    protected = []
    polymer_target_xyz = atoms(protein_only(target))[0]
    from .structures import ref

    for residue in protein_residues(key):
        residue_chain = gemmi.Chain(key.name)
        residue_chain.add_residue(residue.clone())
        if distances(atoms(residue_chain)[0], polymer_target_xyz).min() < 5.0:
            protected.append(ref(key.name, residue).model_dump())
    save(
        out / "redesign_brief.json",
        {
            "status": "proposal_not_executed",
            "architecture": "compact_two_recognition_modules",
            "component_id": cid,
            "candidate_layout_ids": result["top_layout_ids"],
            "protected_binding_residues_per_instance": protected,
            "selection_definition": "observed component residues within 5 Å heavy-atom distance of target polypeptide",
            "geometry_inputs": [r["artifact"] for r in eligible[: request.top_k]],
            "instances": [
                {"id": "key1", "output_chain": "K", "source_chain": key.name},
                {"id": "key2", "output_chain": "L", "source_chain": key.name},
            ],
            "redesign_goal": "Explore a shared compact scaffold supporting both protected binding motifs",
            "editable_regions": None,
            "unknowns": [
                "Required folding core and disulfides beyond the binding surface have not been assigned",
                "Non-interface residues are not automatically safe to delete",
                "Precise generative model conditioning must be validated before execution",
            ],
            "future_routes": {
                "motif_scaffolding": "preferred assembly-redesign experiment; executable model input not compiled",
                "Raygun_or_SCISOR": "optional single-module shortening experiment; not run or installed",
            },
            "acceptance_requirements": [
                "preserve each required interface after independent full-atom reprediction",
                "evaluate the joint complex and forbidden cross-contacts",
                "model sequence connectivity and linkers",
                "compare full-construct compactness at equal evidence level",
            ],
            "conformational_design": "deferred by user",
        },
    )
    seq = sequence(key)
    (out / "proposed_fusion.fasta").write_text(
        ">tandem_observed_only_unmodeled\n" + seq + request.linker_sequence + seq + "\n"
    )
    write_review(out, result)
    save(
        out / "manifest.json",
        {
            "source_hashes": {
                n: sha((bundle / n).read_bytes())
                for n in ["response.json", source_component.structure, state.prepared]
            },
            "files": {p.name: sha(p.read_bytes()) for p in out.iterdir() if p.is_file()},
            "evidence": "geometry_only",
            "gemmi_version": gemmi.__version__,
            "numpy_version": np.__version__,
        },
    )
    return result


def write_review(out, result):
    import html

    selected = [r for r in result["layouts"] if r["id"] in result["top_layout_ids"]]
    lines = [
        '<!doctype html><meta charset="utf-8"><title>Compact assembly search</title>',
        "<style>body{font:15px system-ui;margin:32px;max-width:1200px}td,th{padding:7px;border:1px solid #ddd}table{border-collapse:collapse}svg{background:#f3f6f8}</style>",
        "<h1>Compact rigid assembly search</h1>",
        f"<p>{result['sampled_layouts']} sampled placements; {result['eligible_layouts']} pass the geometry screen. No sequence or backbone design, linker folding, or conformational sampling.</p>",
        "<p>Each target co-moves with its key. Relative target placement is hypothetical. Compact packing is a proposal for subsequent redesign, not an interdigitated shared scaffold.</p>",
        '<p><a href="placement_search.json">Full results, transforms and rejected placements</a></p>',
    ]
    reference = result["spread_out_reference"]
    for row in [reference] + selected:
        label = row.get("id", "Spread-out reference")
        rg = row["compactness"]["key_ca_radius_of_gyration_angstrom"]
        lines.append(
            f'<h2>{html.escape(label)}</h2><p>Cα radius of gyration: {rg:.2f} Å; linker anchor distance: {row["connection"]["anchor_distance_angstrom"]:.2f} Å; screen: {row["status"]}. <a href="{row["artifact"]}">mmCIF</a></p>'
        )
        if row is not reference:
            change = 100 * (1 - rg / reference["compactness"]["key_ca_radius_of_gyration_angstrom"])
            lines.append(
                f"<p>{change:.1f}% smaller key Cα radius of gyration than spread-out reference, with the same polypeptide atoms. This is spatial packing, not a reduction in sequence length.</p>"
            )
        st = gemmi.read_structure(str(out / row["artifact"]))
        chain_data = [(c.name, atoms(c, ca=True)[0]) for c in st[0]]
        # Shared projection scale makes compactness visible across layouts.
        for axes, name in [([0, 1], "XY"), ([0, 2], "XZ"), ([1, 2], "YZ")]:
            xyz = np.concatenate([x for _, x in chain_data])
            center = xyz.mean(0)
            lines.append(
                f'<svg width="360" height="360"><text x="10" y="20">{name}: 2 pixels / Å</text>'
            )
            for i, (chain, coords) in enumerate(chain_data):
                color = ["#1685aa", "#df731e", "#824cc5", "#46994e"][i]
                for x, y in (coords - center)[:, axes] * 2 + 180:
                    lines.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="2" fill="{color}"/>')
            lines.append("</svg>")
    lines.append(
        "<h2>All sampled layouts and overlap control</h2><table><tr><th>ID</th><th>Status</th><th>Key Cα Rg (Å)</th><th>Anchor gap (Å)</th><th>Failed checks</th></tr>"
    )
    for row in result["layouts"]:
        values = [
            row["id"],
            row["status"],
            round(row["compactness"]["key_ca_radius_of_gyration_angstrom"], 2),
            round(row["connection"]["anchor_distance_angstrom"], 2),
            ", ".join(row["failed_checks"]),
        ]
        lines.append("<tr>" + "".join(f"<td>{html.escape(str(v))}</td>" for v in values) + "</tr>")
    lines.append("</table>")
    (out / "review.html").write_text("\n".join(lines))
