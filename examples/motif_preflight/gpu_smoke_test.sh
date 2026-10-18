#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
# Activate an isolated environment from the pinned Foundry commit first.
# Check native parsing and fixed-residue identity before trusting any model output.
rfd3 design inputs=rfd3_inputs.draft.json out_dir=gpu_smoke prevalidate_inputs=true diffusion_batch_size=1 n_batches=1 inference_mode=legacy
