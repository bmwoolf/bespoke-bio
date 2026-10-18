#!/usr/bin/env bash
# Run inside a Linux CUDA 12.8 GPU pod, with this directory already uploaded.
set -euo pipefail
cd "$(dirname "$0")"
command -v uv >/dev/null || python3 -m pip install 'uv==0.8.22'
commit=d078b24d6d0e6f8607e9a0311128f4a897f83953
if [ ! -d /workspace/foundry/.git ]; then
  git clone https://github.com/RosettaCommons/foundry.git /workspace/foundry
fi
git -C /workspace/foundry checkout "$commit"
uv venv --python 3.12 /workspace/rfd3-env
uv pip install --python /workspace/rfd3-env/bin/python 'torch==2.11.0' --index-url https://download.pytorch.org/whl/cu128
uv pip install --python /workspace/rfd3-env/bin/python '/workspace/foundry[rfd3]' --constraint constraints.txt
source /workspace/rfd3-env/bin/activate
python -c 'import torch; assert torch.cuda.is_available(); print(torch.__version__, torch.cuda.get_device_name(0))'
export FOUNDRY_CHECKPOINT_DIRS=/workspace/checkpoints
foundry install rfd3 --checkpoint-dir /workspace/checkpoints
uv pip freeze --python /workspace/rfd3-env/bin/python > worker_requirements.txt
python - <<'PY'
import hashlib, json, pathlib, subprocess, torch
weights={str(p):hashlib.file_digest(p.open('rb'), 'sha256').hexdigest() for p in pathlib.Path('/workspace/checkpoints').rglob('*.ckpt')}
pathlib.Path('worker_provenance.json').write_text(json.dumps({'foundry_commit':subprocess.check_output(['git','-C','/workspace/foundry','rev-parse','HEAD'],text=True).strip(),'torch':torch.__version__,'cuda':torch.version.cuda,'gpu':torch.cuda.get_device_name(0),'weights':weights},indent=2))
PY
bash gpu_smoke_test.sh
