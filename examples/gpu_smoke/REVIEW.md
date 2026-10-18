# Two RFD3 GPU smoke tests

Prepared, not GPU-executed. The native RFD3 parser and installed CLI configuration composition passed locally without weights. `native_validation.json` records exact chain lengths, protected atom identity/mapping and joint-frame geometry.

Both cases preserve 82 residues in four blocks and both 115-residue PD-L1 target chains. Each proposed scaffold has 170 residues: 82 protected and 88 generated. One order interleaves the two modules (K1 L1 L2 K2); the matched control is tandem (K1 K2 L1 L2). Two outputs total, one per case, seed 7. The longer scaffold is an engineering allowance for rebuilding framework support, not a demonstrated minimum functional size.

The core audit finds nine omitted core-like packing residues per module: author 36W, 37F, 51L, 52L, 68F, 79V, 81L, 83M, 124V. These are geometric candidates for support that must be replaced, not measured essential residues. Twenty-two retained buried hydrophobic/Cys positions receive native `select_buried` conditioning. Preserve recognition geometry and disulfide compatibility while rebuilding packing and a coherent backbone. Native acceptance does not establish design feasibility.

## Worker preparation

Use the launch profile in `research/runpod_launch_plan.json`. Upload this directory to `/workspace/gpu_smoke`. The official PyTorch template is only the base image. Install uv, then run `bash bootstrap.sh` inside the worker. It checks out the pinned Foundry source, creates an isolated Python 3.12 environment, installs CUDA 12.8 Torch 2.11 and constrained native dependencies, downloads only the RFD3 checkpoint, records its hash/device/runtime, reruns native validation and configuration composition, then runs the two cases sequentially with diffusion batch size one. Local parser testing used macOS Torch 2.14.1; the GPU runtime must pass its own validation before inference.

Retain the complete `gpu_smoke/` model output directory, worker provenance, frozen requirements, resolved configuration, native worker validation and console log. Download these before terminating the task-created pod. A process timeout does not terminate the pod or end billing. The $100 phase cap includes setup/storage; do not leave the worker running between experiments.

## Local review after download

For each final, nontrajectory CIF, use the actual matching topology ID:

```sh
.venv/bin/python scripts/check_rfd3_output.py /absolute/path/to/output.cif.gz topology_K1_L1_L2_K2 --output outputs/interleaved-check.json
.venv/bin/python scripts/check_rfd3_output.py /absolute/path/to/control.cif.gz topology_K1_K2_L1_L2 --output outputs/control-check.json
```

These check protected atom identities/geometry against one joint target fit, scaffold length and peptide continuity. Passing is coordinate triage only. Sequence design, independent refolding, clashes/packing/disulfide assessment, and functional testing remain downstream. A smoke test is successful when the model loads, both jobs produce inspectable complete artifacts, and discrepancies are reported; it need not produce a viable binder.

## GPU compatibility finding

First GPU startup exposed an AtomWorks API change in 3.0.2/3.0.3 (`AnnotateCovalentModifications`). The worker constraint now pins 3.0.1, which retains the constructor required by this exact Foundry commit. The original local parser report used 3.0.3; the worker reruns native parsing with 3.0.1. Preserve the failed first-attempt log alongside the retry provenance.
