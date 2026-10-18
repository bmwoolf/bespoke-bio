#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
rf3 fold inputs=fold_inputs.json out_dir=fold_control
