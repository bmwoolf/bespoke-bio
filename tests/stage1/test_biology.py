"""Independent sequence expectations for phase, strand, indels and complete allele pairs."""

from dataclasses import replace

import pytest

from bespoke_bio.stage1.biology import (
    Deferred,
    IntervalMask,
    Site,
    complete_protein,
    construct_pair,
    diploid_slots,
    edited_sequence,
    project,
)
from bespoke_bio.stage1.contracts import ExpressionRow
from bespoke_bio.stage1.references import CodingTranscript, Segment, coding_segments


def transcript(strand="+"):
    return CodingTranscript(
        "ENST_TEST.1",
        "ENSG_TEST.1",
        "TEST",
        "protein_coding",
        "chr1",
        strand,
        segments=[Segment(100, 115, "exon1")],
        sequence="ATGGAATTCCCCTAA",
    )


def site(start=103, ref="G", alts=("C",), genotype=(0, 1), **kwargs):
    return Site(
        "chr1",
        start,
        ref,
        alts,
        genotype,
        False,
        None,
        "DNA",
        ("PASS",),
        30,
        60,
        "source:1",
        **kwargs,
    )


EXPRESSION = ExpressionRow(
    transcript_id="ENST_TEST.1", length=15, effective_length=10, tpm=4.2, num_reads=8.25
)


def test_unphased_pair_preserves_reference_and_fractional_expression():
    records, evidence = construct_pair(transcript(), [site()], EXPRESSION, "solution")
    assert [e["protein"] for e in evidence] == ["MEFP", "MQFP"]
    assert all(r.coding_sequence == "ATGGAATTCCCCTAA" for r in records)
    assert records[1].variants[0].cds_position == 4
    assert [r.expression.raw_count for r in records] == [8.25, 8.25]
    assert records[0].phase_scenario != records[1].phase_scenario


@pytest.mark.parametrize(
    "genotype,expected",
    [
        ((0, 0), ["MEFP", "MEFP"]),
        ((1, 1), ["MQFP", "MQFP"]),
        ((1, 2), ["MQFP", "MKFP"]),
        ((2, 1), ["MQFP", "MKFP"]),
    ],
)
def test_actual_diploid_alleles(genotype, expected):
    _, evidence = construct_pair(
        transcript(), [site(alts=("C", "A"), genotype=genotype)], EXPRESSION, "solution"
    )
    assert [e["protein"] for e in evidence] == expected


def test_negative_strand_reverses_both_full_anchored_alleles():
    tx = transcript("-")
    variant = project(tx, site(start=109, ref="TTC", alts=("G",)), 1)
    assert (variant.cds_position, variant.ref, variant.alt) == (4, "GAA", "C")


def test_full_ref_span_must_fit_one_exon():
    tx = transcript()
    tx.segments = [Segment(100, 106, "one"), Segment(120, 129, "two")]
    with pytest.raises(Deferred, match="boundary_edit"):
        project(tx, site(start=105, ref="AT", alts=("A",)), 1)


def test_ref_mismatch_is_integrity_failure():
    with pytest.raises(ValueError, match="REF mismatch"):
        project(transcript(), site(ref="A"), 1)


def test_multiple_hets_require_same_supported_block():
    first = site()
    second = site(start=106, ref="T", alts=("C",))
    with pytest.raises(Deferred, match="unresolved_phase"):
        diploid_slots([first, second])
    first = replace(first, phased=True, phase_set="100")
    second = replace(second, phased=True, phase_set="101")
    with pytest.raises(Deferred, match="unresolved_phase"):
        diploid_slots([first, second])
    assert diploid_slots([first, replace(second, phase_set="100")]) == [(0, 1), (0, 1)]


@pytest.mark.parametrize("gt", [(None, 1), (1,), (0, 1, 1), (0, 2)])
def test_missing_or_invalid_gt_is_never_reference(gt):
    with pytest.raises(Deferred):
        diploid_slots([site(genotype=gt)])


def test_joint_edits_can_rescue_single_stop():
    first = replace(site(alts=("T",)), phased=True, phase_set="100")
    second = replace(site(start=104, ref="A", alts=("C",)), phased=True, phase_set="100")
    _, evidence = construct_pair(transcript(), [first, second], EXPRESSION, "solution")
    assert [e["protein"] for e in evidence] == ["MEFP", "MSFP"]


def test_resolved_partner_survives_one_deferred_allele():
    records, evidence = construct_pair(transcript(), [site(alts=("T",))], EXPRESSION, "solution")
    assert len(records) == 1
    assert records[0].phase_scenario.endswith(":allele-1")
    assert [e["status"] for e in evidence] == ["resolved", "deferred"]
    assert evidence[1]["reason"] == "internal_stop"


def test_indels_keep_subsequent_original_coordinates():
    tx = transcript()
    deletion = project(tx, site(start=102, ref="GGAA", alts=("G",)), 1)
    later = project(tx, site(start=109, ref="C", alts=("A",)), 1)
    assert edited_sequence(tx.sequence, [later, deletion]) == "ATGTTCACCTAA"
    assert complete_protein("ATGTTCACCTAA") == "MFT"


def test_overlapping_edits_defer_instead_of_ordering_arbitrarily():
    edit = project(transcript(), site(), 1)
    with pytest.raises(Deferred, match="overlapping"):
        edited_sequence(transcript().sequence, [edit, edit])


@pytest.mark.parametrize(
    "sequence,reason",
    [
        ("ATGGAA", "missing_terminal_stop"),
        ("ATGGAATA", "partial_codon"),
        ("CTGGAATAA", "missing_start"),
        ("ATGNNNTAA", "ambiguous_allele"),
    ],
)
def test_complete_orf_required(sequence, reason):
    with pytest.raises(Deferred, match=reason):
        complete_protein(sequence)


def test_terminal_stop_added_exactly_once():
    tx = transcript()
    tx.cds = [Segment(100, 115, "one", 0)]
    tx.stops = [Segment(112, 115, "one", 0)]
    assert coding_segments(tx) == [Segment(100, 115, "one")]


def test_mask_requires_every_base_and_rejects_unsorted_intervals():
    mask = IntervalMask([(10, 20), (21, 30)])
    assert mask.covers(10, 20)
    assert not mask.covers(19, 22)
    with pytest.raises(ValueError):
        IntervalMask([(10, 20), (15, 30)])
