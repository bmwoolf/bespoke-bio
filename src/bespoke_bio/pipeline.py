from __future__ import annotations

from hashlib import sha256

from bespoke_bio.models import (
    TargetProtein,
    TranslateRequest,
    TranslateResponse,
    TranslationWarning,
)
from bespoke_bio.translation import (
    apply_variants,
    first_changed_residue,
    protein_hash,
    translate_cds,
)


def translate_case(request: TranslateRequest) -> TranslateResponse:
    targets: list[TargetProtein] = []

    for transcript in request.transcripts:
        reference_translation = translate_cds(transcript.coding_sequence)
        personalized_cds = apply_variants(transcript.coding_sequence, transcript.variants)
        personalized_translation = translate_cds(personalized_cds)

        warnings = list(personalized_translation.warnings)
        if transcript.haplotype == "unknown" and transcript.variants:
            warnings.append(
                TranslationWarning(
                    code="unphased_variants",
                    message="Variants were translated together without known phasing",
                )
            )

        target_seed = "|".join(
            (
                request.case_id,
                transcript.transcript_id,
                transcript.haplotype,
                protein_hash(personalized_translation.protein),
            )
        )
        target_id = f"target_{sha256(target_seed.encode()).hexdigest()[:12]}"
        targets.append(
            TargetProtein(
                target_id=target_id,
                gene_id=transcript.gene_id,
                gene_symbol=transcript.gene_symbol,
                transcript_id=transcript.transcript_id,
                haplotype=transcript.haplotype,
                source_variant_ids=[variant.variant_id for variant in transcript.variants],
                expression=transcript.expression,
                reference_amino_acid_sequence=reference_translation.protein,
                amino_acid_sequence=personalized_translation.protein,
                amino_acid_length=len(personalized_translation.protein),
                sequence_hash=protein_hash(personalized_translation.protein),
                first_changed_residue=first_changed_residue(
                    reference_translation.protein, personalized_translation.protein
                ),
                translation_status=("complete_with_warnings" if warnings else "complete"),
                warnings=warnings,
            )
        )

    return TranslateResponse(
        case_id=request.case_id,
        targets=targets,
        target_count=len(targets),
        unique_sequence_count=len({target.sequence_hash for target in targets}),
    )
