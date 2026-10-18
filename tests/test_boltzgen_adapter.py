import json
from pathlib import Path

import gemmi

from biocad.boltzgen_adapter import prepare

ROOT = Path(__file__).resolve().parents[1]


def test_target_entity_sequence_and_normalized_indices_survive_adapter(tmp_path):
    plan = prepare(ROOT / "examples", tmp_path)
    cif = gemmi.cif.read_file(str(tmp_path / "target.cif")).sole_block()
    assert len(cif.find_values("_entity_poly_seq.mon_id")) == 115
    structure = gemmi.read_structure(str(tmp_path / "target.cif"))
    assert [r.label_seq for r in structure[0]["A"]] == list(range(1, 116))
    spec = json.loads((tmp_path / "pdl1.yaml").read_text())
    assert spec["entities"][0]["protein"]["sequence"] == 75
    assert plan["num_designs"] == 2 and plan["interleaving"] is False
    assert plan["site_selection"]["method"] == "unconstrained"
    assert "binding_types" not in spec["entities"][1]["file"]
    assert not (tmp_path / "binder_inputs.json").exists()
    assert "--skip_inverse_folding" not in plan["command"]
    assert plan["status"] == "prepared_not_executed"


def test_summary_keeps_failed_candidates_and_requires_complete_filter_output(tmp_path):
    import pytest

    from biocad.boltzgen_adapter import summarize

    results = tmp_path / "results"
    with pytest.raises(FileNotFoundError):
        summarize(results, tmp_path / "summary.json")
    folder = results / "final_ranked_designs"
    folder.mkdir(parents=True)
    (folder / "all_designs_metrics.csv").write_text(
        "id,file_name,pass_filters,pass_filter_rmsd_filter,filter_rmsd\nok,ok.cif,True,True,1.0\nbad,bad.cif,False,False,8.0\n"
    )
    report = summarize(results, tmp_path / "summary.json")
    assert report["total"] == 2 and report["passed"] == 1
    assert report["candidates"][1]["failed_filters"] == ["pass_filter_rmsd_filter"]
