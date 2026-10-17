import pytest
from fastapi.testclient import TestClient

from bespoke_bio.api import Variant, VariantError, app, apply_variants, translate_cds

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
