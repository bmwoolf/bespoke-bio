"""Exercise sequence-only publication and confidence-boundary rejection with native VCFs."""

from pathlib import Path

import pysam
import pytest

from bespoke_bio.stage1.io import artifact, decode_json, read_json, write_jsonl
from bespoke_bio.stage1.precomputed import PrecomputedSource, qualify_source
from bespoke_bio.stage1.proteins import export_proteins, validate_protein_package
from bespoke_bio.stage1.references import CodingTranscript, Segment


def fixture(tmp_path, *, boundary=False, spanning=False, annotation_id="ENST.1"):
    paths = []
    for annotated in (False, True):
        path = tmp_path / ("annotated.vcf.gz" if annotated else "original.vcf.gz")
        header = pysam.VariantHeader()
        header.contigs.add("chr1", length=100)
        header.formats.add("GT", 1, "String", "Genotype")
        header.add_sample("HG002")
        if annotated:
            header.info.add(
                "CSQ", ".", "String", "Format: Allele|Consequence|Feature_type|Feature|ALLELE_NUM"
            )
        with pysam.VariantFile(str(path), "wz", header=header) as stream:
            row = stream.new_record(
                contig="chr1", start=13, alleles=("G", "C", "*") if spanning else ("G", "C")
            )
            row.samples["HG002"]["GT"] = (1, 2) if spanning else (0, 1)
            row.samples["HG002"].phased = True
            if annotated:
                row.info["CSQ"] = (f"C|missense_variant|Transcript|{annotation_id}|1",)
            stream.write(row)
            if boundary:
                row = stream.new_record(contig="chr1", start=18, alleles=("AAAA", "A"))
                row.samples["HG002"]["GT"] = (0, 1)
                stream.write(row)
        pysam.tabix_index(str(path), preset="vcf")
        paths.append(path)
    bed = tmp_path / "confidence.bed"
    bed.write_text("chr1\t10\t19\n")
    readme = tmp_path / "README"
    readme.write_text("Synthetic fixture, not a benchmark.\n")
    source = PrecomputedSource(
        donor="HG002",
        sample="HG002",
        release="fixture",
        variants=artifact(paths[0], "variants", "VCF"),
        variant_index=artifact(Path(str(paths[0]) + ".tbi"), "index", "TBI"),
        benchmark_regions=artifact(bed, "regions", "BED"),
        source_readme=artifact(readme, "readme", "text"),
    )
    tx = CodingTranscript(
        "ENST.1",
        "ENSG.1",
        "G",
        "protein_coding",
        "chr1",
        "+",
        segments=[Segment(10, 19, "exon")],
        sequence="ATGGAATAA",
    )
    records = tmp_path / "references.jsonl"
    write_jsonl(records, [tx.record()])
    return source, records, paths[0], paths[1], bed


def test_actual_sequences_without_expression_or_variant_payload(tmp_path):
    source, records, native, annotated, bed = fixture(tmp_path)
    result = export_proteins(source, records, native, annotated, bed, tmp_path / "output")
    assert result["counts"]["target_records"] == 2
    assert result["unique_sequences"] == 2
    targets = [
        decode_json(line)
        for line in (tmp_path / "output/protein-targets.jsonl").read_text().splitlines()
    ]
    assert [r["amino_acid_sequence"] for r in targets] == ["ME", "MQ"]
    assert all(r["expression"] is None and "variants" not in r for r in targets)
    assert result["existing_step2_api_compatible"] is False
    assert read_json(tmp_path / "output/manifest.json")["status"] == "sequence_export_complete"
    assert validate_protein_package(tmp_path / "output")["target_records"] == 2
    with (tmp_path / "output/proteins.fasta").open("a") as stream:
        stream.write("corruption\n")
    with pytest.raises(ValueError, match="artifact integrity"):
        validate_protein_package(tmp_path / "output")


def test_crossing_event_removes_confidence_instead_of_inventing_reference(tmp_path):
    source, _, _, _, _ = fixture(tmp_path, boundary=True)
    result = qualify_source(source, tmp_path / "qualified")
    assert result == {"qualified": 1, "outside_confidence": 1, "filtered": 0, "missing_genotype": 0}
    assert (tmp_path / "qualified/effective-confidence.bed").read_text() == "chr1\t10\t18\n"


def test_empty_release_does_not_publish_complete_manifest(tmp_path):
    source, records, native, annotated, _bed = fixture(tmp_path)
    bed2 = tmp_path / "narrow.bed"
    bed2.write_text("chr1\t10\t18\n")
    with pytest.raises(ValueError, match="no complete"):
        export_proteins(source, records, native, annotated, bed2, tmp_path / "output")
    assert not (tmp_path / "output/manifest.json").exists()


def test_annotation_loss_stops_sequence_publication(tmp_path):
    source, records, native, annotated, bed = fixture(tmp_path)
    with pysam.VariantFile(str(annotated)) as stream:
        header = stream.header.copy()
    empty = tmp_path / "empty.vcf"
    with pysam.VariantFile(str(empty), "w", header=header):
        pass
    with pytest.raises(ValueError, match="dropped"):
        export_proteins(source, records, native, empty, bed, tmp_path / "output")
    assert not (tmp_path / "output/manifest.json").exists()


def test_spanning_deletion_is_deferred_without_losing_supported_partner(tmp_path):
    source, records, native, annotated, bed = fixture(tmp_path, spanning=True)
    result = export_proteins(source, records, native, annotated, bed, tmp_path / "output")
    assert result["counts"]["target_records"] == 1
    assert result["counts"]["partially_resolved"] == 1
    assert result["reason_counts"]["unsupported_allele"] == 1
    assert result["allele_reason_counts"]["unsupported_allele"] == 1
    disposition = decode_json((tmp_path / "output/evidence/dispositions.jsonl").read_text())
    assert disposition["reason"] == "unsupported_allele"
    assert validate_protein_package(tmp_path / "output")["target_records"] == 1
    target = decode_json((tmp_path / "output/protein-targets.jsonl").read_text())
    assert target["amino_acid_sequence"] == "MQ"


def test_versioned_annotation_mismatch_is_an_integrity_failure(tmp_path):
    source, records, native, annotated, bed = fixture(tmp_path, annotation_id="ENST.2")
    with pytest.raises(ValueError, match="missing versioned coding annotation"):
        export_proteins(source, records, native, annotated, bed, tmp_path / "output")
    assert not (tmp_path / "output/manifest.json").exists()
