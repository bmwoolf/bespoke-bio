from copy import deepcopy

import pytest

from bespoke_bio.api import Reference
from bespoke_bio.stage1.biology import IntervalMask, Site
from bespoke_bio.stage1.contracts import ExpressionRow, Policy
from bespoke_bio.stage1.evidence import SiteIndex
from bespoke_bio.stage1.handoff import evaluate, preflight
from bespoke_bio.stage1.references import CodingTranscript, Segment


def inputs():
    tx = CodingTranscript(
        "ENST.1",
        "ENSG.1",
        "G",
        "protein_coding",
        "chr1",
        "+",
        segments=[Segment(10, 19, "e1")],
        sequence="ATGGAATAA",
    )
    expression = ExpressionRow(
        transcript_id=tx.transcript_id, length=9, effective_length=5, tpm=10, num_reads=3.5
    )
    return tx, expression


def test_identical_proteins_keep_two_biological_records():
    tx, expression = inputs()
    result = evaluate(
        "case",
        Reference(genome_build="GRCh38", annotation="test"),
        {tx.transcript_id: tx},
        {tx.transcript_id: expression},
        {"chr1": IntervalMask([(10, 19)])},
        [("chr1", SiteIndex([]))],
        Policy(),
    )
    response = preflight(result.request, result.alleles)
    assert response.target_count == 2
    assert response.unique_sequence_count == 1
    assert len({t.target_id for t in response.targets}) == 2
    assert result.summary["released_expressed"] == 1


def test_missing_coverage_never_becomes_reference():
    tx, expression = inputs()
    result = evaluate(
        "case",
        Reference(genome_build="GRCh38", annotation="test"),
        {tx.transcript_id: tx},
        {tx.transcript_id: expression},
        {"chr1": IntervalMask([(10, 18)])},
        [("chr1", SiteIndex([]))],
        Policy(),
    )
    assert result.request is None
    assert result.dispositions[0].reasons == ["coverage"]


def test_full_output_retains_every_supported_transcript():
    tx, row = inputs()
    transcripts, expression = {}, {}
    for i in range(12):
        value = deepcopy(tx)
        value.transcript_id = f"ENST{i:02}.1"
        transcripts[value.transcript_id] = value
        expression[value.transcript_id] = row.model_copy(
            update={"transcript_id": value.transcript_id}
        )
    result = evaluate(
        "case",
        Reference(genome_build="GRCh38", annotation="test"),
        transcripts,
        expression,
        {"chr1": IntervalMask([(10, 19)])},
        [("chr1", SiteIndex([]))],
        Policy(),
    )
    assert result.summary["released_expressed"] == 12
    assert result.summary["transcript_counts"] == {"selected": 12}
    assert len(result.request.transcripts) == 24


def test_missing_expression_is_not_zero():
    tx, _ = inputs()
    with pytest.raises(ValueError, match="expression join"):
        evaluate(
            "case",
            Reference(genome_build="GRCh38", annotation="test"),
            {tx.transcript_id: tx},
            {},
            {"chr1": IntervalMask([(10, 19)])},
            [("chr1", SiteIndex([]))],
            Policy(),
        )


def test_zero_expression_is_evidence_not_an_output_filter():
    tx, row = inputs()
    result = evaluate(
        "case",
        Reference(genome_build="GRCh38", annotation="test"),
        {tx.transcript_id: tx},
        {tx.transcript_id: row.model_copy(update={"tpm": 0, "num_reads": 0})},
        {"chr1": IntervalMask([(10, 19)])},
        [("chr1", SiteIndex([]))],
        Policy(),
    )
    assert len(result.request.transcripts) == 2
    assert all(record.expression.tpm == 0 for record in result.request.transcripts)
    assert result.summary["transcript_counts"] == {"selected": 1}


def test_one_unsupported_allele_does_not_discard_resolved_partner():
    tx, row = inputs()
    site = Site(
        "chr1",
        13,
        "G",
        ("T",),
        (0, 1),
        False,
        None,
        "sample",
        ("PASS",),
        20,
        30,
        "record",
        annotated_transcripts=frozenset({tx.transcript_id}),
    )
    result = evaluate(
        "case",
        Reference(genome_build="GRCh38", annotation="test"),
        {tx.transcript_id: tx},
        {tx.transcript_id: row},
        {"chr1": IntervalMask([(10, 19)])},
        [("chr1", SiteIndex([site]))],
        Policy(),
    )
    assert len(result.request.transcripts) == 1
    assert result.dispositions[0].status == "partially_resolved"
    assert result.dispositions[0].reasons == ["internal_stop"]
    assert result.summary["allele_counts"] == {"resolved": 1, "deferred": 1}
    assert preflight(result.request, result.alleles).target_count == 1
