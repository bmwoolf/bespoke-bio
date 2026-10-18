"""Compose the actual installed CLI configuration without constructing a model."""

import json
from pathlib import Path

import rfd3.cli
from hydra import compose, initialize_config_dir
from omegaconf import OmegaConf

module = Path(rfd3.cli.__file__)
directory = module.parent.parent.parent / "configs"
if not directory.exists():
    directory = module.parent / "configs"
with initialize_config_dir(config_dir=str(directory), version_base="1.3"):
    config = compose(
        config_name="inference",
        overrides=[
            "inference_engine=rfdiffusion3",
            "inputs=rfd3_inputs.draft.json",
            "out_dir=gpu_smoke",
            "prevalidate_inputs=true",
            "diffusion_batch_size=1",
            "n_batches=1",
            "inference_mode=legacy",
            "seed=7",
        ],
    )
assert config.inference_mode == "legacy" and config.diffusion_batch_size == 1
Path("resolved_config.yaml").write_text(OmegaConf.to_yaml(config))
print(json.dumps({"status": "passed", "scope": "CLI configuration composition; no model loaded"}))
