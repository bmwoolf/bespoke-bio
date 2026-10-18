#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
# Model files live on the network volume; interpreters execute on the pod's local disk.
export HF_HOME=/workspace/cache/huggingface
export UV_CACHE_DIR=/tmp/boltzgen-uv-cache
export UV_LINK_MODE=hardlink
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=4
export PYTHONUNBUFFERED=1
command -v uv >/dev/null || python3 -m pip install uv==0.12.24
revision=a07fea1153111ef453aaea7a15697ee05de9291f
snapshot=/workspace/cache/runtimes/boltzgen-${revision}-torch280-py312
mkdir -p "$snapshot" /workspace/cache/boltzgen
# Snapshots are trusted account-local artifacts, restored to their original absolute paths.
if [[ -f "$snapshot/runtime.tar" && -f "$snapshot/runtime.tar.sha256" ]]; then
  (cd "$snapshot" && sha256sum -c runtime.tar.sha256)
  tar -xf "$snapshot/runtime.tar" -C /
  cp "$snapshot/esm_python.txt" esm_python.txt
  echo 'Restored versioned local runtime snapshot'
else
  uv venv --system-site-packages --python 3.12 /opt/boltzgen-env
  curl -fsSL "https://codeload.github.com/HannesStark/boltzgen/tar.gz/$revision" -o /tmp/boltzgen-source.tar.gz
  mkdir -p /tmp/boltzgen-source
  tar -xzf /tmp/boltzgen-source.tar.gz -C /tmp/boltzgen-source --strip-components=1
  echo 'torch==2.8.0' > constraints.txt
  uv pip install --python /opt/boltzgen-env/bin/python /tmp/boltzgen-source --constraint constraints.txt
fi
source /opt/boltzgen-env/bin/activate
python -c 'import torch; assert torch.__version__.startswith("2.8.0"); assert torch.cuda.is_available(); torch.ones(1,device="cuda"); print(torch.__version__, torch.cuda.get_device_name(0))'
python - <<'PY'
from pathlib import Path
from boltzgen.task.esmfold2.runtime import resolve_python
p = Path('esm_python.txt')
candidate = p.read_text().strip() if p.exists() else None
resolved = resolve_python(candidate if candidate and Path(candidate).is_file() else None, require_cuda=True)
p.write_text(resolved + '\n')
PY
uv pip freeze --python /opt/boltzgen-env/bin/python > worker_requirements.txt
uv pip freeze --python "$(cat esm_python.txt)" > esmfold2_requirements.txt
# Commit the snapshot only after both environments pass real CUDA allocation checks.
if [[ ! -f "$snapshot/runtime.tar" || ! -f "$snapshot/runtime.tar.sha256" ]]; then
  esm_environment=$(dirname "$(dirname "$(cat esm_python.txt)")")
  managed_python=$("$(cat esm_python.txt)" -c 'import sys; print(sys.base_prefix)')
  snapshot_paths=(opt/boltzgen-env "${esm_environment#/}")
  if [[ "$managed_python" == /root/.local/share/uv/python/* ]]; then
    snapshot_paths+=("${managed_python#/}")
  fi
  tar -cf "$snapshot/runtime.tar.partial" -C / "${snapshot_paths[@]}"
  mv "$snapshot/runtime.tar.partial" "$snapshot/runtime.tar"
  cp esm_python.txt "$snapshot/esm_python.txt"
  (cd "$snapshot" && sha256sum runtime.tar > runtime.tar.sha256.partial && mv runtime.tar.sha256.partial runtime.tar.sha256)
fi
printf '%s\n' 'Runtime checked; model-cache readiness still requires a completed pipeline.'
if [[ "${1:-run}" == warmup ]]; then
  boltzgen check pdl1.yaml --output checked_unconstrained --cache /workspace/cache/boltzgen
  python - <<'PYWARM'
import json
from pathlib import Path
from types import SimpleNamespace
from huggingface_hub import hf_hub_download, snapshot_download
from boltzgen.cli.boltzgen import ARTIFACTS, get_artifact_path
from boltzgen.task.esmfold2.contract import MODEL_REPO, MODEL_REVISION, ESMC_REPO, ESMC_REVISION
args = SimpleNamespace(cache='/workspace/cache/boltzgen', force_download=False, models_token=None)
files = []
for name in ['design-diverse', 'design-adherence', 'solublempnn', 'folding', 'moldir']:
    artifact, kind = ARTIFACTS[name]
    files.append(get_artifact_path(args, artifact, repo_type=kind))
files.append(Path(hf_hub_download(MODEL_REPO, 'ccd.pkl', revision=MODEL_REVISION)))
for repo, revision in [(MODEL_REPO, MODEL_REVISION), (ESMC_REPO, ESMC_REVISION)]:
    directory = Path(snapshot_download(repo, revision=revision, allow_patterns=['*.json', '*.safetensors']))
    files.extend(p for p in directory.rglob('*') if p.is_file())
manifest = [{'path': str(p), 'bytes': p.stat().st_size} for p in files]
assert all(p['bytes'] > 0 for p in manifest)
Path('/workspace/cache/warm_manifest.json').write_text(json.dumps({'status': 'downloaded; runtime CUDA checked; inference not tested', 'files': manifest}, indent=2))
print('Model files downloaded and runtime CUDA checked; inference readiness still requires an actual job.')
PYWARM
else
  bash run.sh
fi
