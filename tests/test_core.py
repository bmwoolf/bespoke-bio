from pathlib import Path

from biocad.core import audit_core

ROOT = Path(__file__).resolve().parents[1]


def test_core_audit_rejects_assumption_that_motif_mask_preserves_support(tmp_path):
    report = audit_core(
        ROOT / "examples/compact_assembly", ROOT / "examples/gpu_smoke", tmp_path / "audit.json"
    )
    assert report["evidence"].startswith("geometry_derived")
    assert len(report["modules"]) == 2
    for module in report["modules"]:
        assert module["omitted_packing_support"]
        assert any(r["residue"]["number"] == 36 for r in module["omitted_packing_support"])
