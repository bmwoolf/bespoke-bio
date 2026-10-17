from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256

from bespoke_bio.models import CdsVariant, TranslationWarning

CODON_TABLE = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


class VariantApplicationError(ValueError):
    """Raised when variants cannot be applied unambiguously to a reference CDS."""


@dataclass(frozen=True)
class TranslationResult:
    protein: str
    warnings: tuple[TranslationWarning, ...]


def apply_variants(reference_cds: str, variants: list[CdsVariant]) -> str:
    """Apply non-overlapping variants using coordinates from the original CDS."""
    ordered = sorted(variants, key=lambda item: item.cds_position)
    result: list[str] = []
    cursor = 0

    for variant in ordered:
        start = variant.cds_position - 1
        end = start + len(variant.ref)
        if start < cursor:
            raise VariantApplicationError(
                f"variant {variant.variant_id} overlaps an earlier variant"
            )
        if end > len(reference_cds):
            raise VariantApplicationError(
                f"variant {variant.variant_id} extends beyond the reference CDS"
            )
        observed = reference_cds[start:end]
        if observed != variant.ref:
            raise VariantApplicationError(
                f"variant {variant.variant_id} expected {variant.ref} at CDS position "
                f"{variant.cds_position}, found {observed}"
            )
        result.extend((reference_cds[cursor:start], variant.alt))
        cursor = end

    result.append(reference_cds[cursor:])
    return "".join(result)


def translate_cds(cds: str) -> TranslationResult:
    warnings: list[TranslationWarning] = []
    if not cds.startswith("ATG"):
        warnings.append(
            TranslationWarning(code="missing_start", message="CDS does not begin with ATG")
        )
    if len(cds) % 3:
        warnings.append(
            TranslationWarning(
                code="partial_codon",
                message=f"CDS has {len(cds) % 3} trailing nucleotide(s)",
            )
        )

    amino_acids: list[str] = []
    complete_codon_count = len(cds) // 3
    for index in range(complete_codon_count):
        codon = cds[index * 3 : index * 3 + 3]
        amino_acid = CODON_TABLE.get(codon, "X")
        if amino_acid == "X":
            warnings.append(
                TranslationWarning(
                    code="ambiguous_codon",
                    message=f"Codon {index + 1} contains an ambiguous nucleotide",
                )
            )
        if amino_acid == "*":
            if index != complete_codon_count - 1:
                warnings.append(
                    TranslationWarning(
                        code="internal_stop",
                        message=f"Translation stopped at amino-acid position {index + 1}",
                    )
                )
            break
        amino_acids.append(amino_acid)

    return TranslationResult(protein="".join(amino_acids), warnings=tuple(warnings))


def protein_hash(sequence: str) -> str:
    return f"sha256:{sha256(sequence.encode()).hexdigest()}"


def first_changed_residue(reference: str, alternate: str) -> int | None:
    for index, (left, right) in enumerate(zip(reference, alternate), start=1):
        if left != right:
            return index
    if len(reference) != len(alternate):
        return min(len(reference), len(alternate)) + 1
    return None
