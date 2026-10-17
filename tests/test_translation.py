import pytest

from bespoke_bio.models import CdsVariant
from bespoke_bio.translation import VariantApplicationError, apply_variants, translate_cds


def variant(position: int, ref: str, alt: str) -> CdsVariant:
    return CdsVariant(
        variant_id=f"v{position}",
        cds_position=position,
        ref=ref,
        alt=alt,
    )


def test_translates_a_complete_cds_and_omits_terminal_stop() -> None:
    result = translate_cds("ATGGAACTGAAATAA")

    assert result.protein == "MELK"
    assert result.warnings == ()


def test_applies_substitution_using_one_based_coordinates() -> None:
    personalized = apply_variants("ATGGAACTGAAATAA", [variant(7, "C", "T")])

    assert personalized == "ATGGAATTGAAATAA"
    assert translate_cds(personalized).protein == "MELK"


def test_applies_indel_without_shifting_later_reference_coordinates() -> None:
    personalized = apply_variants(
        "ATGAAACCCGGGTAA",
        [variant(4, "A", "AT"), variant(10, "G", "A")],
    )

    assert personalized == "ATGATAACCCAGGTAA"


def test_rejects_a_reference_allele_mismatch() -> None:
    with pytest.raises(VariantApplicationError, match="expected C"):
        apply_variants("ATGAAATAA", [variant(4, "C", "T")])


def test_stops_translation_at_an_internal_stop() -> None:
    result = translate_cds("ATGGAATAAAAATAA")

    assert result.protein == "ME"
    assert [warning.code for warning in result.warnings] == ["internal_stop"]
