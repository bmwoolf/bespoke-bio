"""Pinned-source Mosaic adapter. GPU execution remains unverified until a real run."""

import importlib.util
import subprocess
import time
from pathlib import Path

import numpy as np

from .geometry import distances
from .pipeline import save
from .structures import atoms, ports, select_indices, sequence

MOSAIC_COMMIT = "1eb4af64655870b45a715a5a7e61570e8b098701"


def compile_request(request, target):
    if request.target_state_ids[0] != target.id:
        raise ValueError("Primary target state mismatch")
    chain = request.site[0].chain
    idx = select_indices(target, request.site, chain)
    return {
        "engine": "mosaic",
        "commit": MOSAIC_COMMIT,
        "model": "Boltz2",
        "target_sequence": target.chains[chain],
        "target_state_id": target.id,
        "target_chain": chain,
        "binder_length": request.binder_length,
        "epitope_indices_target_relative_zero_based": idx,
        "seed": request.seed,
        "steps": request.steps,
        "objective_weights": {
            "site_contact": 2.0,
            "within_binder_contact": 1.0,
            "mpnn_sequence_recovery": 5.0,
        },
        "chain_order": ["binder", "target"],
        "full_atom_reprediction_required": True,
        "request": request.model_dump(),
        "constraint_enforcement": {
            "site": "during_generation: BinderTargetContact(epitope_idx)",
            "length": "during_generation",
            "attachment_clearance": "after_generation: nearest target atom to selected terminal anchor",
            "selected_target_states": "postgeneration evaluation required; only primary state modeled by this adapter",
            "exclusions": "unsupported by this single-target adapter",
            "preserved_second_module": "assembly stage: immutable component revision",
            "connected_construct": "unsupported: linker folding and whole-fusion reprediction pending",
            "affinity_or_function": "unsupported",
        },
        "execution_status": "source_verified_not_gpu_executed",
    }


def run(plan, output):
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    if any(out.iterdir()):
        raise ValueError("Design output directory must be empty")
    save(out / "plan.json", plan)
    if not importlib.util.find_spec("mosaic") or not importlib.util.find_spec("jax"):
        save(
            out / "status.json",
            {
                "status": "blocked",
                "reason": "Pinned Mosaic, GPU JAX runtime and model weights are not installed",
            },
        )
        raise RuntimeError(
            "Mosaic GPU environment unavailable; see status.json and docs/COMPUTE.md"
        )
    import jax
    import jax.numpy as jnp
    import mosaic
    import mosaic.losses.structure_prediction as sp
    from mosaic.common import TOKENS
    from mosaic.losses.protein_mpnn import InverseFoldingSequenceRecovery, inverse_fold
    from mosaic.models.boltz2 import Boltz2
    from mosaic.optimizers import gradient_MCMC, simplex_APGM
    from mosaic.proteinmpnn.mpnn import load_mpnn_sol
    from mosaic.structure_prediction import TargetChain

    repo = Path(mosaic.__file__).resolve().parents[2]
    revision = subprocess.check_output(
        ["git", "-C", str(repo), "rev-parse", "HEAD"], text=True
    ).strip()
    if revision != MOSAIC_COMMIT:
        raise RuntimeError(f"Mosaic commit mismatch: {revision}")
    if not any(d.platform in ("gpu", "tpu") for d in jax.devices()):
        raise RuntimeError("GPU/TPU required for this bounded live runner")
    start = time.monotonic()
    model = Boltz2()
    mpnn = load_mpnn_sol(0.05)
    n = plan["binder_length"]
    target = TargetChain(plan["target_sequence"], use_msa=False)
    features, writer = model.binder_features(binder_length=n, chains=[target])
    terms = [
        sp.BinderTargetContact(epitope_idx=plan["epitope_indices_target_relative_zero_based"]),
        sp.WithinBinderContact(),
        InverseFoldingSequenceRecovery(mpnn, temp=jnp.array(0.01)),
    ]
    weights = list(plan["objective_weights"].values())
    loss = model.build_loss(
        loss=sum((w * t for w, t in zip(weights[1:], terms[1:])), weights[0] * terms[0]),
        features=features,
        recycling_steps=2,
    )
    key = jax.random.key(plan["seed"])
    keys = jax.random.split(key, 8)
    x = jax.nn.softmax(0.5 * jax.random.gumbel(keys[0], (n, 20)), axis=-1)
    _, best = simplex_APGM(
        loss_function=loss, x=x, n_steps=plan["steps"], stepsize=0.1, momentum=0.0, key=keys[1]
    )
    _, best = simplex_APGM(
        loss_function=loss, x=best, n_steps=10, stepsize=0.1, scale=1.5, key=keys[2]
    )
    design = model.predict(
        PSSM=best, features=features, writer=writer, recycling_steps=3, key=keys[3]
    )
    design.st.make_mmcif_document().write_file(
        str(out / "design_backbone_not_for_atomistic_checks.cif")
    )
    ids = inverse_fold(
        mpnn=mpnn, binder_length=n, output=design.model_output, temp=0.1, key=keys[4]
    )
    ids = gradient_MCMC(loss=loss, sequence=ids, steps=10, key=keys[5])
    seq = "".join(TOKENS[int(i)] for i in ids)
    # target_only_features rebuilds amino-acid-specific atom channels for both discrete chains.
    ff, fw = model.target_only_features(
        [TargetChain(seq, use_msa=False), TargetChain(plan["target_sequence"], use_msa=False)]
    )
    final = model.predict(features=ff, writer=fw, recycling_steps=3, key=keys[6])
    if [sequence(c) for c in final.st[0]] != [seq, plan["target_sequence"]]:
        raise RuntimeError("Final predicted chain identities/order differ from binder-first plan")
    final.st.make_mmcif_document().write_file(str(out / "final_complex.cif"))
    mf, mw = model.target_only_features([TargetChain(seq, use_msa=False)])
    monomer = model.predict(features=mf, writer=mw, recycling_steps=3, key=keys[7])
    monomer.st.make_mmcif_document().write_file(str(out / "final_monomer.cif"))
    final_output = final.model_output
    per_pair = {f"{a}:{b}": v for (a, b), v in final_output.chain_pair_iptm().items()}
    hard = jax.nn.one_hot(ids, 20)
    objective = []
    for name, w, term in zip(plan["objective_weights"], weights, terms):
        v, aux = term(hard, final_output, key=keys[0])
        objective.append(
            {
                "name": name,
                "weight": w,
                "value": float(v),
                "scope": "final full-atom prediction; raw objective and weight recorded separately",
                "aux": {k: np.asarray(v).tolist() for k, v in aux.items()},
            }
        )
    port = next(p for p in ports(final.st[0][0]) if p.terminus == plan["request"]["attachment"])
    clearance = (
        float(distances([port.frame.translation], atoms(final.st[0][1])[0]).min())
        if port.frame
        else None
    )
    (out / "candidate.fasta").write_text(">candidate\n" + seq + "\n")
    save(
        out / "result.json",
        {
            "sequence": seq,
            "chain_order": ["binder", "target"],
            "mosaic_commit": revision,
            "evidence": "predicted",
            "elapsed_seconds": time.monotonic() - start,
            "chain_pair_iptm": per_pair,
            "confidence_scale": "iPTM 0..1; not affinity",
            "terminal_clearance_angstrom": clearance,
            "terminal_clearance_pass": None
            if clearance is None
            else clearance >= plan["request"]["terminal_clearance_angstrom"],
            "objective": objective,
            "remaining_evaluations": plan["constraint_enforcement"],
            "status": "primary_state_predicted; multi-state matrix and fusion reprediction still required",
        },
    )
