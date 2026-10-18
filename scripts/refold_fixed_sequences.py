#!/usr/bin/env python3
"""Refold fixed protein sequences using the installed BoltzGen folding checkpoint.

Input JSON: {"tasks": [{"id":"gene_reference_b1", "target_sequence":"ACD...",
"binder_sequence":"EFG..."}]}. Omit binder_sequence for target-alone ensembles.
Run inside BoltzGen's environment, passing --checkpoint and --moldir explicitly.
--validate-only performs native parsing/featurization without model/GPU execution.
Every seed retains one structure per task; no best-of-N structural selection.
The seed controls Lightning/Torch RNG; native feature construction also creates its
own numpy.default_rng(None), so bitwise reproducibility is not claimed.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


class FixedSequenceDataset:
    """Label scoring groups after featurization, as native FromGenerated does."""

    def __init__(self, dataset):
        self.dataset = dataset

    def __len__(self):
        return len(self.dataset)

    def __getitem__(self, index):
        features = self.dataset[index]
        # Two-chain inputs: A target, B binder. Monomer: score the whole chain.
        asym = features["asym_id"]
        mask = asym == asym.max()
        features["design_mask"] = mask.bool()
        features["chain_design_mask"] = mask.bool()
        return features


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--moldir", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--seeds", nargs="+", type=int, default=[11, 22, 33])
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    import hydra
    import numpy as np
    import pytorch_lightning as pl
    import yaml
    from omegaconf import OmegaConf

    tasks = json.loads(args.manifest.read_text())["tasks"]
    ids = [t["id"] for t in tasks]
    if not ids or len(ids) != len(set(ids)):
        raise ValueError("Task IDs must be nonempty and unique")
    inputs = args.output.resolve() / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    paths = []
    for task in tasks:
        if not re.fullmatch(r"[A-Za-z0-9_-]+", task["id"]):
            raise ValueError("Unsafe task id")
        sequences = [task["target_sequence"]]
        if task.get("binder_sequence"):
            sequences.append(task["binder_sequence"])
        if any(not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", s) for s in sequences):
            raise ValueError("Only fully specified canonical protein sequences are accepted")
        path = inputs / (task["id"] + ".yaml")
        path.write_text(yaml.safe_dump({"entities": [
            {"protein": {"id": chr(65 + i), "sequence": seq, "msa": "empty"}}
            for i, seq in enumerate(sequences)]}))
        paths.append(str(path))
    receipt = {"mode": "validate_only" if args.validate_only else "fixed_sequence_refolding",
               "checkpoint": args.checkpoint, "seeds": args.seeds, "tasks": tasks,
               "sequence_redesign": False, "target_templates": False, "results": []}
    for seed in args.seeds:
        out = args.output.resolve() / f"seed_{seed}"
        cfg = {
            "_target_": "boltzgen.task.predict.predict.Predict",
            "data": {
                "_target_": "boltzgen.task.predict.data_from_yaml.FromYamlDataModule",
                "cfg": {
                    "_target_": "boltzgen.task.predict.data_from_yaml.DataConfig",
                    "moldir": args.moldir, "multiplicity": 1, "yaml_path": paths,
                    "tokenizer": {"_target_": "boltzgen.data.tokenize.tokenizer.Tokenizer",
                                  "atomize_modified_residues": False},
                    "featurizer": {"_target_": "boltzgen.data.feature.featurizer.Featurizer"},
                    "design": True, "atom14": True, "disulfide_on": True,
                    "disulfide_prob": 1.0, "diffusion_samples": 1,
                }, "batch_size": 1, "num_workers": 0, "pin_memory": True,
            },
            "writer": {"_target_": "boltzgen.task.predict.writer.FoldingWriter",
                       "design_dir": str(out)},
            "checkpoint": args.checkpoint, "output": str(out), "name": None,
            "recycling_steps": 3, "sampling_steps": 200, "diffusion_samples": 1,
            "keys_dict_out": ["iptm", "ptm", "protein_iptm", "design_ptm", "target_ptm",
                              "design_to_target_iptm", "min_design_to_target_pae"],
            "trainer": {"accelerator": "gpu", "devices": 1, "logger": False,
                        "precision": "bf16-mixed"},
            "override": {"validators": None, "use_kernels": False},
        }
        pl.seed_everything(seed, workers=True)
        native = hydra.utils.instantiate(OmegaConf.create(cfg))
        native.data.predict_set = FixedSequenceDataset(native.data.predict_set)
        out.mkdir(parents=True, exist_ok=True)
        (out / "config.yaml").write_text(yaml.safe_dump(cfg))
        if args.validate_only:
            for i, task in enumerate(tasks):
                feature = native.data.predict_set[i]
                expected = len(task["target_sequence"]) + len(task.get("binder_sequence", ""))
                assert len(feature["asym_id"]) == expected
                assert feature["design_mask"].sum().item() == len(
                    task.get("binder_sequence") or task["target_sequence"])
                assert not feature["template_mask"].any()
            receipt["validated_tasks"] = len(tasks)
            break
        native.run()
        archives = sorted(out.rglob("*.npz"))
        if {p.stem for p in archives} != set(ids):
            raise RuntimeError("Native prediction did not produce one archive per task")
        for archive in archives:
            metrics = {}
            with np.load(archive, allow_pickle=False) as data:
                for key in cfg["keys_dict_out"]:
                    if key in data and data[key].size == 1:
                        value = float(data[key].reshape(-1)[0])
                        metrics[key] = value if np.isfinite(value) else None
            receipt["results"].append({"id": archive.stem, "seed": seed,
                                       "archive": str(archive), "metrics": metrics})
        (args.output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    (args.output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"mode": receipt["mode"], "tasks": len(tasks),
                      "scored_outputs": len(receipt["results"])}))


if __name__ == "__main__":
    main()
