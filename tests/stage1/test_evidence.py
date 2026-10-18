import pytest

from bespoke_bio.stage1.contracts import Policy
from bespoke_bio.stage1.evidence import (
    callability,
    crosscheck_identity,
    depth_mask,
    intersect,
    read_expression,
    subtract,
)


def test_salmon_preserves_fractional_counts_and_rejects_missing_columns(tmp_path):
    path = tmp_path / "quant.sf"
    path.write_text("Name\tLength\tEffectiveLength\tTPM\tNumReads\nENST.1\t30\t20.1\t2\t1.75\n")
    assert read_expression(path)["ENST.1"].num_reads == 1.75
    path.write_text("Name\tLength\tTPM\tNumReads\n")
    with pytest.raises(ValueError):
        read_expression(path)


@pytest.mark.parametrize("value", ["NaN", "inf", "-1"])
def test_bad_expression_fails(tmp_path, value):
    path = tmp_path / "quant.sf"
    path.write_text(f"Name\tLength\tEffectiveLength\tTPM\tNumReads\nENST.1\t30\t20\t{value}\t1\n")
    with pytest.raises(ValueError):
        read_expression(path)


def test_depth_requires_zeros_and_complete_ordered_bed(tmp_path):
    path = tmp_path / "depth.tsv"
    path.write_text("chr1\t11\t10\nchr1\t12\t0\nchr1\t13\t20\n")
    assert depth_mask(path, {"chr1": [(10, 13)]}, 10) == {"chr1": [(10, 11), (12, 13)]}
    path.write_text("chr1\t11\t10\nchr1\t13\t20\n")
    with pytest.raises(ValueError, match="missing"):
        depth_mask(path, {"chr1": [(10, 13)]}, 10)


def test_interval_algebra_keeps_holes():
    assert intersect([(0, 10), (20, 30)], [(5, 25)]) == [(5, 10), (20, 25)]
    assert subtract([(0, 10), (20, 30)], [(5, 25)]) == [(0, 5), (25, 30)]


HEADER = """##fileformat=VCFv4.3
##contig=<ID=chr1,length=100>
##FILTER=<ID=PASS,Description="Passed">
##FILTER=<ID=LowQual,Description="Low quality">
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Quality">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDNA
"""


def test_original_deletion_span_supplies_evidence_without_base_depth(tmp_path):
    original, gvcf = tmp_path / "raw.vcf", tmp_path / "raw.g.vcf"
    original.write_text(HEADER + "chr1\t5\t.\tAAA\tA\t60\tPASS\t.\tGT:DP:GQ\t0/1:30:60\n")
    gvcf.write_text(
        HEADER
        + "chr1\t1\t.\tA\t<*>\t.\t.\tEND=4\tGT:GQ:MIN_DP\t0/0:60:30\n"
        + "chr1\t8\t.\tA\t<*>\t.\t.\tEND=10\tGT:GQ:MIN_DP\t0/0:60:30\n"
    )
    mask = callability(original, gvcf, "DNA", Policy(), {"chr1": [(0, 4), (7, 10)]})
    assert mask["chr1"].covers(0, 10)
    # A failed call cannot inherit confidence from a reference block crossing its span.
    original.write_text(original.read_text().replace("60\tPASS", "60\tLowQual"))
    gvcf.write_text(HEADER + "chr1\t1\t.\tA\t<*>\t.\t.\tEND=10\tGT:GQ:MIN_DP\t0/0:60:30\n")
    mask = callability(original, gvcf, "DNA", Policy(), {"chr1": [(0, 10)]})
    assert not mask["chr1"].covers(4, 7)


def test_identity_needs_actual_off_diagonal_lod(tmp_path):
    path = tmp_path / "fingerprints.tsv"
    path.write_text("LEFT_GROUP_VALUE\tRIGHT_GROUP_VALUE\tLOD_SCORE\nDNA\tDNA\t100\n")
    with pytest.raises(ValueError, match="off-diagonal"):
        crosscheck_identity(path, "DNA", "RNA", 5)
    path.write_text(path.read_text() + "DNA\tRNA\t4.99\n")
    with pytest.raises(ValueError, match="insufficient"):
        crosscheck_identity(path, "DNA", "RNA", 5)
