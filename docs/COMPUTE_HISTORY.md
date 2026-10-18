> Historical development notes. Superseded by [current status](STATUS.md), [integration](INTEGRATION.md), and [compute](COMPUTE.md). Statements about the active engine, pending work and resource state describe earlier stages.

# Live compute handoff

Local execution used macOS ARM64, Python 3.12.9, Gemmi 0.7.3, NumPy 2.2.6 and Pydantic 2.12.3. Six pytest checks pass. Foldseek downloaded from the official universal macOS archive reports commit `463739e0014a1549a527de589102cde98f802f37`; archive and binary hashes are recorded in `research/runtime.json`. No CUDA GPU, Mosaic installation or model checkpoints were configured in this package environment. No remote compute was launched.

The live adapter is source-checked at Mosaic commit `1eb4af64655870b45a715a5a7e61570e8b098701`, not execution-verified. Use an isolated Linux GPU checkout and the upstream dependency lock; avoid mixing it into this lightweight geometry environment. Exact GPU choice/VRAM and runtime must be measured on our selected target; no unmeasured runtime claim is made.

```sh
git clone https://github.com/escalante-bio/mosaic.git /workspace/mosaic
git -C /workspace/mosaic checkout 1eb4af64655870b45a715a5a7e61570e8b098701
cd /workspace/mosaic
uv sync --group jax-cuda
# Copy this BioCAD package to /workspace/bespoke_bio, then use the Mosaic environment:
uv pip install --python .venv/bin/python --no-deps -e /workspace/bespoke_bio
.venv/bin/python -c 'import jax; print(jax.devices())'
.venv/bin/biocad design /workspace/bespoke_bio/examples/mosaic_plan.json --output /workspace/design-001
```

Upstream package/model dependencies may require additional downloads or credentials; verify them on the selected machine before a campaign. The runner requires a Git source checkout at the exact commit and checks for a GPU/TPU before model construction. It attempts only one primary-state candidate. Its output is not a completed campaign: named alternate/exclusion predictions, acceptance/diversity filtering and fusion reprediction remain to be wired.

To compile a fresh plan locally:

```sh
.venv/bin/biocad compile examples/replay/request.json examples/replay/response.json --output examples/mosaic_plan.json
```

The plan declares each constraint's enforcement and unsupported features. The script runs site contact + within-binder contact + soluble ProteinMPNN recovery objectives, then discrete sequence refinement and full-atom complex/monomer prediction. Attachment access is assessed after generation. No existing scaffold is silently treated as editable or preserved by the live generation step.

A local preflight was executed; it writes `outputs/live_preflight/status.json` and exits with the specific missing-runtime blocker. No neural scores were substituted. Next live acceptance test: one bounded candidate, verify final atom identities and chain order, store weight hashes/device/software versions, predict both named states, and import those results into the same evaluation matrix before any campaign expansion.

## Current RFD3 preparation (2026-10-17)

Runpod API authentication is verified through a read-only account request. No pod has been created. The approved phase cap is $100. A live catalog read reports A100 SXM 80GB Secure Cloud at $1.79/hour; a three-hour initial GPU allowance is $5.37 before storage/other charges. Availability and the final all-in quote must be refreshed at launch. See `research/runpod_launch_plan.json`.

The active smoke-test bundle is `examples/gpu_smoke`, not the older Mosaic plan. It includes an isolated worker bootstrap, native parser/configuration gates, exact Foundry commit, two matched 170-residue scaffold cases, runtime/checkpoint provenance capture and local post-check commands. Native CPU input construction passed for both cases (400 residues including targets). No checkpoint or inference output is present locally. See its REVIEW.md for execution and cleanup.

Credentials live in ignored `.env.runpod.local` and are loaded by `scripts/runpod.py`; do not upload this file, account snapshots or the vendor directory in the design bundle. The CLI can read current infrastructure through the authenticated wrapper. No password is needed.

## Completed smoke run

The authorized two-case run completed and its pod was deleted after verified download. Use the 3.0.1 AtomWorks worker constraint; 3.0.3 passed parser-only checks but failed full model initialization. See `examples/gpu_smoke_results` for actual outputs/provenance and `docs/STATUS.md` for measured runtime and cleanup. No active compute remains from this run.

## Off-the-shelf BoltzGen run

The complete seven-stage unconstrained PD-L1 test finished on A100 80GB and was downloaded
before pod deletion. See `examples/boltzgen_results/review.html` and its status/provenance
records. Use `examples/boltzgen_smoke/bootstrap.sh` and `run.sh` for the tested runtime
setup. This path requires a CUDA-13-capable host for upstream's separate ESMFold2 scorer,
even though the main BoltzGen environment uses CUDA 12.8. All cloud resources created
for this run were removed. The source/runtime/model cache references are retained; large
re-downloadable checkpoint caches were not copied into the results archive.
