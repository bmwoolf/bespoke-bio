from fastapi.testclient import TestClient

from bespoke_bio.api import app

client = TestClient(app)


def payload() -> dict:
    return {
        "case_id": "case-1",
        "reference": {"genome_build": "GRCh38", "annotation": "GENCODE-demo"},
        "transcripts": [
            {
                "gene_id": "gene-1",
                "gene_symbol": "GENE1",
                "transcript_id": "tx-1.1",
                "haplotype": "paternal",
                "coding_sequence": "ATGGAACTGAAATAA",
                "expression": {"tpm": 10.5, "raw_count": 100},
                "variants": [
                    {
                        "variant_id": "v1",
                        "cds_position": 4,
                        "ref": "G",
                        "alt": "C",
                    }
                ],
            },
            {
                "gene_id": "gene-2",
                "gene_symbol": "GENE2",
                "transcript_id": "tx-2.1",
                "haplotype": "maternal",
                "coding_sequence": "ATGGAACTGAAATAA",
                "expression": {"tpm": 2.0},
                "variants": [
                    {
                        "variant_id": "v2",
                        "cds_position": 4,
                        "ref": "G",
                        "alt": "C",
                    }
                ],
            },
        ],
    }


def test_translate_endpoint_returns_targets_and_deduplication_count() -> None:
    response = client.post("/v1/translate", json=payload())

    assert response.status_code == 200
    body = response.json()
    assert body["target_count"] == 2
    assert body["unique_sequence_count"] == 1
    assert body["reference"]["genome_build"] == "GRCh38"
    assert body["targets"][0]["amino_acid_sequence"] == "MQLK"
    assert body["targets"][0]["first_changed_residue"] == 2
    assert "structure" not in body["targets"][0]
    assert "keys" not in body["targets"][0]


def test_translate_endpoint_rejects_reference_mismatch() -> None:
    request = payload()
    request["transcripts"][0]["variants"][0]["ref"] = "A"

    response = client.post("/v1/translate", json=request)

    assert response.status_code == 422
    assert "expected A" in response.json()["detail"]
