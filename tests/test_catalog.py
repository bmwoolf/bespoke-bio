from pathlib import Path

import pytest

from biocad.catalog import CatalogEntry, search_catalog

ROOT = Path(__file__).resolve().parents[1]
FOLDSEEK = ROOT / "vendor/foldseek/bin/foldseek"


def entry():
    return CatalogEntry(
        target_structure=str(ROOT / "examples/5JDS.cif"),
        target_chain="A",
        binder_structure=str(ROOT / "examples/5JDS.cif"),
        binder_chain="B",
        reference="https://www.rcsb.org/structure/5JDS",
        evidence={"type": "experimental_complex", "construct": "observed only"},
        conditions={"fixture": "5JDS"},
    )


@pytest.mark.skipif(not FOLDSEEK.exists(), reason="Bounded native Foldseek binary unavailable")
def test_live_catalog_analogy_and_exact_match(tmp_path):
    analogous = search_catalog(
        ROOT / "examples/5J89.cif",
        "A",
        [entry()],
        tmp_path / "analogy",
        conditions={"fixture": "5J89"},
        executable=str(FOLDSEEK),
    )
    exact = search_catalog(
        ROOT / "examples/5JDS.cif",
        "A",
        [entry()],
        tmp_path / "exact",
        conditions={"fixture": "5JDS"},
        executable=str(FOLDSEEK),
    )
    result = analogous["results"][0]
    exact_result = exact["results"][0]
    assert result["template_id"] == exact_result["template_id"]
    assert result["generated_candidate"] is False
    assert result["conditions_match"] is False
    assert exact_result["conditions_match"] is True
    assert exact_result["exact_observed_target_sequence_match"] is True
    assert exact_result["automatic_candidate_reuse"] is False
    assert result["pose_mapping"]["ca_rmsd_angstrom"] < 3
    assert result["pose_mapping"]["evidence"] == "transferred_analogy"
    assert (tmp_path / "analogy" / result["pose"]).exists()
    assert (tmp_path / "exact/catalog_results.json").exists()
    with pytest.raises(ValueError, match="empty output"):
        search_catalog(ROOT / "examples/5JDS.cif", "A", [entry()], tmp_path / "exact")


def test_catalog_bounded(tmp_path):
    with pytest.raises(ValueError, match="1..100"):
        search_catalog("unused", "A", [], tmp_path)
