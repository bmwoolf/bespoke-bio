import json
import shutil
from pathlib import Path

import gemmi

from bespoke_bio.api import CODON_TABLE, Structure, TranslateRequest, build_bundle
from biocad.structures import sha
from biocad.upstream import import_targets

ROOT = Path(__file__).resolve().parents[1]


def make_upstream(output):
    structure = gemmi.read_structure(str(ROOT / "examples/5JDS.cif"))
    structure.setup_entities()
    seq = gemmi.one_letter_code(next(e for e in structure.entities if e.name == "2").full_sequence)
    reverse = {aa: codon for codon, aa in CODON_TABLE.items()}
    cds = "".join(reverse[aa] for aa in seq) + "TAA"
    request = TranslateRequest.model_validate(
        {
            "case_id": "fixture-import",
            "reference": {
                "genome_build": "synthetic-coding-fixture",
                "annotation": "5JDS deposited protein construct",
            },
            "transcripts": [
                {
                    "gene_id": "CD274",
                    "gene_symbol": "CD274",
                    "transcript_id": tid,
                    "coding_sequence": cds,
                    "haplotype": "reference",
                    "expression": {"tpm": value},
                    "variants": [],
                }
                for tid, value in [("tx1", 10.0), ("tx2", 20.0)]
            ],
        }
    )
    response = build_bundle(request, output)
    shutil.copy2(ROOT / "examples/5JDS.cif", output / "targets/5JDS.cif")
    for r in response.targets:
        r.structure = Structure(
            status="complete", source="RCSB PDB", source_id="5JDS", file="targets/5JDS.cif"
        )
    (output / "targets.jsonl").write_text(
        "".join(r.model_dump_json() + "\n" for r in response.targets)
    )
    return response


def test_real_upstream_contract_to_design_state(tmp_path):
    response = make_upstream(tmp_path / "upstream")
    result = import_targets(tmp_path / "upstream", tmp_path / "imported")
    assert len(result.records) == 2 and all(r["status"] == "imported" for r in result.records)
    assert len(result.target_states) == 1
    state = result.target_states[0]
    assert state.context["exact_protein_sequence"] == response.targets[0].amino_acid_sequence
    assert len(state.chains["A"]) < len(state.context["exact_protein_sequence"])
    assert {m["value"] for m in state.measurements} == {10.0, 20.0}
    assert all(
        m["kind"] == "measured_RNA" and m["protein_abundance"] is None for m in state.measurements
    )
    assert len(state.context["upstream_records"]) == 2
    assert all(m["protein_sequence_position"] is not None for m in state.residue_mapping)


def test_import_retains_unresolved_and_rejects_wrong_sequence_and_escape(tmp_path):
    make_upstream(tmp_path / "upstream")
    path = tmp_path / "upstream/targets.jsonl"
    records = [json.loads(x) for x in path.read_text().splitlines()]
    unresolved = dict(records[0], target_id="unresolved", structure={"status": "unresolved"})
    wrong = dict(
        records[0],
        target_id="wrong",
        amino_acid_sequence="M" + records[0]["amino_acid_sequence"][1:].replace("A", "G", 1),
    )
    wrong["sequence_hash"] = "sha256:" + sha(wrong["amino_acid_sequence"].encode())
    escape = dict(
        records[0],
        target_id="escape",
        structure={**records[0]["structure"], "file": "../secret.cif"},
    )
    path.write_text("\n".join(json.dumps(x) for x in [unresolved, wrong, escape]) + "\n")
    result = import_targets(tmp_path / "upstream", tmp_path / "imported")
    assert [r["status"] for r in result.records] == ["unresolved", "failed", "failed"]
    assert result.target_states == []
    assert "sequence" in result.records[1]["reason"]
    assert "escapes" in result.records[2]["reason"]


def test_pdb_upstream_model_can_be_prepared(tmp_path):
    make_upstream(tmp_path / "upstream")
    original = gemmi.read_structure(str(ROOT / "examples/replay/target.cif"))
    st = gemmi.Structure()
    st.add_model(original[0].clone())
    st.write_pdb(str(tmp_path / "upstream/targets/target.pdb"))
    seq = "".join(
        gemmi.find_tabulated_residue(r.name).one_letter_code
        for r in st[0]["A"]
        if r.het_flag == "A"
    )
    path = tmp_path / "upstream/targets.jsonl"
    record = json.loads(path.read_text().splitlines()[0])
    record.update(
        amino_acid_sequence=seq,
        sequence_hash="sha256:" + sha(seq.encode()),
        amino_acid_length=len(seq),
    )
    record["structure"] = {
        "status": "complete",
        "source": "ColabFold",
        "file": "targets/target.pdb",
    }
    path.write_text(json.dumps(record) + "\n")
    imported = import_targets(tmp_path / "upstream", tmp_path / "imported")
    assert imported.records[0]["status"] == "imported"
    assert imported.target_states[0].evidence == "predicted"
