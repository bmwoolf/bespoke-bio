import json
from pathlib import Path

import gemmi

from biocad.baseline import prepare

ROOT = Path(__file__).resolve().parents[1]


def test_baseline_separates_sequence_only_folding_from_motif_free_binder_design(tmp_path):
    plan = prepare(ROOT / "examples", tmp_path)
    assert plan["interleaving"] is False
    fold = json.loads((tmp_path / "fold_inputs.json").read_text())[0]
    assert set(fold["components"][0]) == {"seq", "chain_id"}
    assert len(fold["components"][0]["seq"]) == 127
    design = json.loads((tmp_path / "binder_inputs.json").read_text())["pdl1_baseline"]
    assert design["contig"] == "75,/0,A1-115"
    target = gemmi.read_structure(str(tmp_path / "target.cif"))[0]["A"]
    assert [r.seqid.num for r in target] == list(range(1, 116))
    assert all(r.het_flag == "A" for r in target)
    assert plan["status"] == "prepared_not_executed"
