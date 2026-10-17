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


def test_build_endpoint_uses_full_pipeline(monkeypatch, tmp_path) -> None:
    monkeypatch.setenv("BESPOKE_BIO_OUTPUT_DIR", str(tmp_path))

    def fake_build(request, output, folding_backend, folding_options):
        assert folding_backend == "esmfold"
        assert folding_options == []
        output.mkdir(parents=True)
        structure_file = Path("targets") / "model.pdb"
        (output / structure_file).parent.mkdir()
        (output / structure_file).write_text("MODEL\nEND\n")
        translated = service.translate_case(request)
        for target in translated.targets:
            target.structure = service.Structure(
                status="complete",
                source="ESMFold",
                file=structure_file.as_posix(),
            )
        return translated

    monkeypatch.setattr(service, "build_bundle", fake_build)
    response = client.post("/v1/build?folding_backend=esmfold", json=payload())

    assert response.status_code == 200
    body = response.json()
    assert body["bundle_id"]
    assert body["target_count"] == 2
    assert body["unique_sequence_count"] == 1
    assert body["targets"][0]["amino_acid_sequence"] == "MQLK"
    assert body["targets"][0]["first_changed_residue"] == 2
    assert body["targets"][0]["structure"]["source"] == "ESMFold"
    file_response = client.get(
        f"/v1/bundles/{body['bundle_id']}/{body['targets'][0]['structure']['file']}"
    )
    assert file_response.status_code == 200
    assert file_response.text == "MODEL\nEND\n"


def test_reference_mismatch_returns_422(monkeypatch, tmp_path) -> None:
    monkeypatch.setenv("BESPOKE_BIO_OUTPUT_DIR", str(tmp_path))
    request = payload()
    request["transcripts"][0]["variants"][0]["ref"] = "A"

    response = client.post("/v1/build?folding_backend=esmfold", json=request)

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


def test_unknown_phase_must_be_explicit(monkeypatch, tmp_path) -> None:
    monkeypatch.setenv("BESPOKE_BIO_OUTPUT_DIR", str(tmp_path))
    request = payload()
    request["transcripts"] = [
        {
            **request["transcripts"][0],
            "haplotype": "unknown",
        }
    ]

    response = client.post("/v1/build?folding_backend=esmfold", json=request)

    assert response.status_code == 422
    assert "Step 1 must supply an explicit phase_scenario" in response.json()["detail"]


def test_bundle_folds_missing_structures(monkeypatch, tmp_path) -> None:
    monkeypatch.setattr(
        service,
        "resolve_structure",
        lambda protein, output: service.Structure(status="unresolved"),
    )

    def fake_fold(proteins, output, options):
        assert options == ["--demo"]
        structures = {}
        for protein in proteins:
            path = Path("targets") / f"{service.structure_key(protein)}.pdb"
            (output / path).write_text("MODEL\nEND\n")
            structures[protein.sequence_hash] = service.Structure(
                status="complete",
                source="test-folder",
                file=path.as_posix(),
            )
        return structures

    monkeypatch.setitem(service.FOLDING_BACKENDS, "esmfold", fake_fold)
    response = build_bundle(
        TranslateRequest.model_validate(payload()),
        tmp_path,
        folding_backend="esmfold",
        folding_options=["--demo"],
    )

    assert response.target_count == 2
    assert (tmp_path / "proteins.fasta").read_text().count(">") == 1
    assert len((tmp_path / "targets.jsonl").read_text().splitlines()) == 2
    assert (tmp_path / "provenance.json").is_file()
    assert '"folding_backend": "esmfold"' in (tmp_path / "provenance.json").read_text()
    assert all(target.structure.status == "complete" for target in response.targets)


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


def pdb_with_plddt(value: float) -> str:
    return (
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00"
        f"{value:6.2f}           C\n"
    )


def test_esmfold_result_is_collected(monkeypatch, tmp_path) -> None:
    protein = service.translate_case(TranslateRequest.model_validate(payload())).targets[0]

    def fake_run(command, check):
        assert check is True
        assert command[0] == "/fake/esm-fold"
        predictions = Path(command[command.index("-o") + 1])
        predictions.mkdir()
        fasta = Path(command[command.index("-i") + 1])
        job = fasta.read_text().splitlines()[0].removeprefix(">")
        (predictions / f"{job}.pdb").write_text(pdb_with_plddt(85.0))

    monkeypatch.setattr(service.shutil, "which", lambda name: f"/fake/{name}")
    monkeypatch.setattr(service.subprocess, "run", fake_run)
    (tmp_path / "targets").mkdir()

    structure = service.run_esmfold([protein], tmp_path, ["--chunk-size", "64"])[
        protein.sequence_hash
    ]

    assert structure.status == "complete"
    assert structure.source == "ESMFold"
    assert structure.mean_plddt == 85.0
    assert (tmp_path / structure.file).is_file()


def test_boltz2_result_is_collected(monkeypatch, tmp_path) -> None:
    protein = service.translate_case(TranslateRequest.model_validate(payload())).targets[0]

    def fake_run(command, check):
        assert check is True
        assert command[:2] == ["/fake/boltz", "predict"]
        inputs = Path(command[2])
        predictions = Path(command[command.index("--out_dir") + 1])
        job = next(inputs.glob("*.yaml")).stem
        result_dir = predictions / "predictions" / job
        result_dir.mkdir(parents=True)
        (result_dir / f"{job}_model_0.cif").write_text("data_model\n")
        (result_dir / f"confidence_{job}_model_0.json").write_text(
            '{"complex_plddt":0.91}'
        )

    monkeypatch.setattr(service.shutil, "which", lambda name: f"/fake/{name}")
    monkeypatch.setattr(service.subprocess, "run", fake_run)
    (tmp_path / "targets").mkdir()

    structure = service.run_boltz2([protein], tmp_path, ["--use_msa_server"])[
        protein.sequence_hash
    ]

    assert structure.status == "complete"
    assert structure.source == "Boltz-2"
    assert structure.mean_plddt == 91.0
    assert (tmp_path / structure.file).is_file()
    assert (tmp_path / structure.confidence_file).is_file()


def test_openfold_result_is_collected(monkeypatch, tmp_path) -> None:
    protein = service.translate_case(TranslateRequest.model_validate(payload())).targets[0]
    script = tmp_path / "run_pretrained_openfold.py"
    templates = tmp_path / "templates"
    script.write_text("")
    templates.mkdir()
    monkeypatch.setenv("BESPOKE_BIO_OPENFOLD_SCRIPT", str(script))
    monkeypatch.setenv("BESPOKE_BIO_OPENFOLD_TEMPLATE_MMCIF_DIR", str(templates))

    def fake_run(command, check, cwd):
        assert check is True
        assert cwd == script.parent
        assert command[1] == str(script)
        inputs = Path(command[2])
        predictions = Path(command[command.index("--output_dir") + 1])
        job = next(inputs.glob("*.fasta")).stem
        result_dir = predictions / "predictions"
        result_dir.mkdir(parents=True)
        (result_dir / f"{job}_model_1_unrelaxed.pdb").write_text(pdb_with_plddt(88.0))

    monkeypatch.setattr(service.subprocess, "run", fake_run)
    (tmp_path / "targets").mkdir()

    structure = service.run_openfold([protein], tmp_path, ["--model_device", "cuda:0"])[
        protein.sequence_hash
    ]

    assert structure.status == "complete"
    assert structure.source == "OpenFold"
    assert structure.mean_plddt == 88.0
    assert (tmp_path / structure.file).is_file()
