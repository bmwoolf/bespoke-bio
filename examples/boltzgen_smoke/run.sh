#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export HF_HOME=/workspace/cache/huggingface
export UV_CACHE_DIR=/tmp/boltzgen-uv-cache
export UV_LINK_MODE=hardlink
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=4
boltzgen check pdl1.yaml --output checked_unconstrained --cache /workspace/cache/boltzgen
esm_python=$(cat esm_python.txt)
boltzgen run pdl1.yaml --output results_unconstrained --protocol protein-anything --reuse --num_designs 2 --budget 2 --diffusion_batch_size 1 --inverse_fold_num_sequences 1 --num_workers 0 --devices 1 --esmfold2_python "$esm_python" --use_kernels false --cache /workspace/cache/boltzgen
