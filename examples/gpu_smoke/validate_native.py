"""Run under the pinned native Foundry environment; loads no model checkpoint."""

import argparse
import importlib.metadata
import json
from pathlib import Path

import numpy as np
from rfd3.inference.input_parsing import DesignInputSpecification


def validate(plan_dir, output):
    root = Path(plan_dir).resolve()
    plan = json.loads((root / "motif_plan.json").read_text())
    inputs = json.loads((root / "rfd3_inputs.draft.json").read_text())
    reports = []
    for name, data in inputs.items():
        spec = DesignInputSpecification(**{**data, "input": str(root / data["input"])})
        source = spec.atom_array_input.copy()
        aa, metadata = spec.build(return_metadata=True)
        topology = next(t for t in plan["topologies"] if t["id"] == name)
        expected = {}
        for span in topology["output_motif_spans"]:
            for src, dst in zip(
                range(span["input_first"], span["input_last"] + 1),
                range(span["output_first"], span["output_last"] + 1),
                strict=True,
            ):
                expected[(span["input_chain"] + str(src))] = ("A", dst)
        for original, assembled in [("T", "B"), ("U", "C")]:
            for res_id in np.unique(source.res_id[source.chain_id == original]):
                expected[original + str(res_id)] = (assembled, int(res_id))
        counts = {
            str(c): len(np.unique(aa.res_id[aa.chain_id == c])) for c in np.unique(aa.chain_id)
        }
        assert counts == {"A": topology["scaffold_length"], "B": 115, "C": 115}, (
            f"Unexpected chain lengths {counts}"
        )
        source_map = {
            (str(c) + str(r), str(a)): (str(n), coord)
            for c, r, a, n, coord in zip(
                source.chain_id,
                source.res_id,
                source.atom_name,
                source.res_name,
                source.coord,
                strict=True,
            )
        }
        points_src = []
        points_built = []
        seen = set()
        generated = set()
        for i, component in enumerate(aa.src_component):
            component = str(component)
            if component not in expected:
                assert component.endswith("P"), f"Unrecognized source token {component}"
                assert (
                    not aa.is_motif_atom_with_fixed_coord[i]
                    and not aa.is_motif_atom_with_fixed_seq[i]
                ), "Generated residue accidentally fixed"
                generated.add((str(aa.chain_id[i]), int(aa.res_id[i])))
                continue
            chain, number = expected[component]
            assert str(aa.chain_id[i]) == chain and int(aa.res_id[i]) == number, (
                "Native output residue map differs from compiler"
            )
            assert bool(aa.is_motif_atom_with_fixed_coord[i]) and bool(
                aa.is_motif_atom_with_fixed_seq[i]
            ), "Protected atom lost fixed coordinate/sequence status"
            token = (component, str(aa.atom_name[i]))
            identity, coord = source_map[token]
            assert str(aa.res_name[i]) == identity, "Protected residue identity changed"
            points_src.append(coord)
            points_built.append(aa.coord[i])
            seen.add(token)
        expected_atoms = {k for k in source_map if k[0] in expected}
        assert seen == expected_atoms, "Protected input atoms were lost or duplicated"
        assert len(generated) == topology["scaffold_length"] - plan["fixed_motif_residues"], (
            "Unexpected generated-region size"
        )
        a, b = np.array(points_built, dtype=float), np.array(points_src, dtype=float)
        u, _, vt = np.linalg.svd((a - a.mean(0)).T @ (b - b.mean(0)))
        d = np.eye(3)
        d[-1, -1] = np.linalg.det(vt.T @ u.T)
        r = vt.T @ d @ u.T
        fitted = np.einsum("nj,ij->ni", a - a.mean(0), r) + b.mean(0)
        maximum = float(np.linalg.norm(fitted - b, axis=1).max())
        assert maximum < 0.001, "Native parser changed protected relative geometry"
        reports.append(
            {
                "id": name,
                "status": "passed",
                "chain_lengths": counts,
                "residue_count": sum(counts.values()),
                "atom_count": len(aa),
                "protected_atoms_checked": len(seen),
                "protected_residues": len(expected),
                "generated_scaffold_residues": len(generated),
                "maximum_fixed_atom_displacement_after_joint_alignment_angstrom": maximum,
                "sampled_contig": metadata["extra"]["sampled_contig"],
                "scaffold_motif_residue_map": [
                    {"source": k, "chain": v[0], "residue": v[1]}
                    for k, v in expected.items()
                    if v[0] == "A"
                ],
            }
        )
    result = {
        "status": "passed" if reports else "failed",
        "foundry_commit": plan["foundry_commit"],
        "method": "Unmodified native DesignInputSpecification constructor and build(return_metadata=True)",
        "device": "CPU",
        "checkpoint_loaded": False,
        "cases": reports,
        "packages": {
            p: importlib.metadata.version(p)
            for p in ["rc-foundry", "atomworks", "biotite", "torch", "numpy", "pydantic"]
        },
        "scope": "Native input construction and fixed-atom mapping; no generation/folding validation",
    }
    Path(output).write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps({"status": result["status"], "cases": len(reports), "checkpoint_loaded": False})
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("plan")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    validate(args.plan, args.output)
