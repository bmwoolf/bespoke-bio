#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python check_config.py
python validate_native.py . --output native_validation_on_worker.json
rfd3 design inputs=rfd3_inputs.draft.json out_dir=gpu_smoke prevalidate_inputs=true diffusion_batch_size=1 n_batches=1 inference_mode=legacy seed=7
