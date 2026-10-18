import pytest

from bespoke_bio.stage1.references import protein_sequences


def test_gencode_transcript_name_may_repeat_unversioned_identifier(tmp_path):
    path = tmp_path / "proteins.fa"
    path.write_text(
        ">ENSP00000465763.1|ENST00000484532.6|ENSG00000285629.1|OTTHUMG1|OTTHUMT1|ENST00000484532|ENSG00000285629|3\nMEK\n"
    )
    assert protein_sequences(path) == {"ENST00000484532.6": "MEK"}


def test_duplicate_versioned_protein_records_are_rejected(tmp_path):
    path = tmp_path / "proteins.fa"
    entry = ">ENSP1.1|ENST00000484532.6|ENSG1.1\nMEK\n"
    path.write_text(entry + entry)
    with pytest.raises(ValueError, match="ambiguous"):
        protein_sequences(path)
