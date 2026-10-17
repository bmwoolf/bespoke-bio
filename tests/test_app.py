from pathlib import Path

import pytest
from fastapi.testclient import TestClient

import bespoke_bio.api as service
from bespoke_bio.api import (
    TranslateRequest,
    Variant,
    VariantError,
    app,
    apply_variants,
    build_bundle,
    translate_cds,
)

client = TestClient(app)


def variant(position: int, ref: str, alt: str) -> Variant:
    return Variant(variant_id=f"v{position}", cds_position=position, ref=ref, alt=alt)


def payload() -> dict:
    transcript = {
        "gene_id": "gene-1",
        "gene_symbol": "GENE1",
        "transcript_id": "tx-1.1",
        "haplotype": "paternal",
        "coding_sequence": "ATGGAACTGAAATAA",
        "expression": {"tpm": 10.5},
        "variants": [{"variant_id": "v1", "cds_position": 4, "ref": "G", "alt": "C"}],
    }
    return {
        "case_id": "case-1",
        "reference": {"genome_build": "GRCh38", "annotation": "GENCODE-demo"},
        "transcripts": [transcript, {**transcript, "transcript_id": "tx-2.1"}],
    }


def test_translate_endpoint() -> None:
    response = client.post("/v1/translate", json=payload())

    assert response.status_code == 200
    body = response.json()
    assert body["target_count"] == 2
    assert body["unique_sequence_count"] == 1
    assert body["targets"][0]["amino_acid_sequence"] == "MQLK"
    assert body["targets"][0]["first_changed_residue"] == 2


def test_reference_mismatch_returns_422() -> None:
    request = payload()
    request["transcripts"][0]["variants"][0]["ref"] = "A"

    response = client.post("/v1/translate", json=request)

    assert response.status_code == 422
    assert "expected A" in response.json()["detail"]


def test_indels_keep_reference_coordinates() -> None:
    personalized = apply_variants(
        "ATGAAACCCGGGTAA",
        [variant(4, "A", "AT"), variant(10, "G", "A")],
    )

    assert personalized == "ATGATAACCCAGGTAA"


def test_overlapping_variants_are_rejected() -> None:
    with pytest.raises(VariantError, match="overlaps"):
        apply_variants("ATGAAATAA", [variant(4, "AA", "C"), variant(5, "A", "G")])


def test_translation_stops_at_internal_stop() -> None:
    protein, warnings = translate_cds("ATGGAATAAAAATAA")

    assert protein == "ME"
    assert warnings == ["internal_stop"]


def test_unphased_variants_generate_possible_proteins() -> None:
    assert service.genotype_state("./.") == "uncertain"
    request = payload()
    request["transcripts"] = [
        {
            **request["transcripts"][0],
            "haplotype": "unknown",
            "variants": [
                {"variant_id": "v1", "cds_position": 4, "ref": "G", "alt": "C", "genotype": "0/1"},
                {"variant_id": "v2", "cds_position": 7, "ref": "C", "alt": "A", "genotype": "0/1"},
            ],
        }
    ]

    response = client.post("/v1/translate", json=request)

    assert response.status_code == 200
    body = response.json()
    assert body["target_count"] == 4
    assert body["unique_sequence_count"] == 4
    assert {target["amino_acid_sequence"] for target in body["targets"]} == {
        "MELK",
        "MEMK",
        "MQLK",
        "MQMK",
    }
    assert all("unphased_variants" in target["warnings"] for target in body["targets"])


def test_bundle_contains_step_3_handoff_files(tmp_path) -> None:
    response = build_bundle(TranslateRequest.model_validate(payload()), tmp_path)

    assert response.target_count == 2
    assert (tmp_path / "proteins.fasta").read_text().count(">") == 1
    assert len((tmp_path / "targets.jsonl").read_text().splitlines()) == 2
    assert (tmp_path / "provenance.json").is_file()
    assert all(target.structure.status == "unresolved" for target in response.targets)


def test_alphafold_match_is_downloaded(monkeypatch, tmp_path) -> None:
    request = payload()
    request["transcripts"] = [{**request["transcripts"][0], "uniprot_id": "P12345", "variants": []}]
    protein = service.translate_case(TranslateRequest.model_validate(request)).targets[0]

    monkeypatch.setattr(
        service,
        "fetch_json",
        lambda url: [
            {
                "sequence": protein.amino_acid_sequence,
                "cifUrl": "https://example.test/model.cif",
                "paeDocUrl": "https://example.test/pae.json",
                "entryId": "AF-P12345-F1",
                "globalMetricValue": 91.2,
                "latestVersion": 6,
            }
        ],
    )
    monkeypatch.setattr(service, "download", lambda url, path: path.write_text(url))
    (tmp_path / "targets").mkdir()

    structure = service.find_alphafold(protein, tmp_path)

    assert structure is not None
    assert structure.status == "complete"
    assert structure.source == "AlphaFold DB"
    assert (tmp_path / structure.file).is_file()


def test_colabfold_result_is_collected(monkeypatch, tmp_path) -> None:
    protein = service.translate_case(TranslateRequest.model_validate(payload())).targets[0]

    def fake_run(command, check):
        assert check is True
        predictions = Path(command[2])
        predictions.mkdir()
        job = Path(command[1]).read_text().splitlines()[0].removeprefix(">")
        (predictions / f"{job}_unrelaxed_rank_001_model.pdb").write_text("MODEL\nEND\n")
        (predictions / f"{job}_scores_rank_001_model.json").write_text('{"plddt":[80.0,90.0]}')

    monkeypatch.setattr(service.shutil, "which", lambda name: "/fake/colabfold_batch")
    monkeypatch.setattr(service.subprocess, "run", fake_run)
    (tmp_path / "targets").mkdir()

    structure = service.run_colabfold([protein], tmp_path)[protein.sequence_hash]

    assert structure.status == "complete"
    assert structure.source == "ColabFold"
    assert structure.mean_plddt == 85.0
    assert (tmp_path / structure.file).is_file()
