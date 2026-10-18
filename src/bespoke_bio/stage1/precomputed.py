"""Import already called variants without claiming to have processed their source reads."""

from pathlib import Path

import pysam

from .biology import IntervalMask
from .contracts import Artifact, Identifier, Record
from .evidence import read_bed, require_sample, subtract
from .io import canonical, verify, write_new
from .references import AUTOSOMES, attributes, merged_intervals, open_text
from .variants import validate_cardinalities


class PrecomputedSource(Record):
    """A benchmark-derived demonstration input, kept separate from caller evaluation.

    The matching BED defines the source's small-variant confidence scope. Its evidence
    is assembly based; missing read depth or genotype quality must remain missing.
    """

    donor: Identifier
    sample: Identifier
    release: Identifier
    variants: Artifact
    variant_index: Artifact
    benchmark_regions: Artifact
    source_readme: Artifact


def phase_identity(record, sample: str, source_sha256: str) -> str | None:
    """Preserve explicit phase, including VCF's documented default phase-set semantics.

    A phased GT without a PS subfield belongs to the source's default phase set on that
    chromosome. An unphased GT never acquires phase from this rule. Namespace the set
    by source and sample, and retain the original FORMAT fields in the native VCF.
    """
    call = record.samples[sample]
    if not call.phased:
        return None
    if "PS" not in record.format:
        return f"{source_sha256}:{sample}:{record.contig}:vcf-default"
    if call.get("PS") is None:
        return None
    return f"{source_sha256}:{sample}:{record.contig}:PS:{call['PS']}"


def qualify_source(source: PrecomputedSource, output: Path) -> dict:
    """Separate supported small-variant evidence before reference-dependent normalization.

    The source VCF also contains calls outside its matching benchmark BED, including
    structural events. An event crossing a confidence boundary cannot establish a
    reference allele inside that boundary. Remove its full span from the usable mask.
    The complete original callset remains immutable and separately registered.
    """
    regions = read_bed(Path(source.benchmark_regions.location))
    masks = {chrom: IntervalMask(spans) for chrom, spans in regions.items()}
    holes, counts = (
        {},
        {"qualified": 0, "outside_confidence": 0, "filtered": 0, "missing_genotype": 0},
    )
    output.mkdir(parents=True, exist_ok=False)
    with pysam.VariantFile(source.variants.location) as original:
        require_sample(original, source.sample)
        with pysam.VariantFile(
            str(output / "qualified.bcf"), "wb", header=original.header
        ) as selected:
            for row in original:
                gt = row.samples[source.sample].get("GT") or ()
                reason = None
                if not masks.get(row.contig, IntervalMask([])).covers(row.start, row.stop):
                    reason = "outside_confidence"
                elif set(row.filter) - {"PASS"}:
                    reason = "filtered"
                elif len(gt) != 2 or any(a is None or a < 0 or a >= len(row.alleles) for a in gt):
                    reason = "missing_genotype"
                if reason:
                    counts[reason] += 1
                    holes.setdefault(row.contig, []).append((row.start, row.stop))
                else:
                    selected.write(row)
                    counts["qualified"] += 1
    effective = {
        chrom: subtract(spans, merged_intervals(holes.get(chrom, [])))
        for chrom, spans in regions.items()
    }
    bed = "".join(f"{chrom}\t{a}\t{b}\n" for chrom, spans in effective.items() for a, b in spans)
    write_new(output / "effective-confidence.bed", bed.encode())
    write_new(output / "qualification.json", canonical(counts) + b"\n")
    return counts


def coding_regions(gtf: Path, fasta: Path, flank: int) -> dict[str, list[tuple[int, int]]]:
    """Select annotated coding intervals and explicit local splice context for annotation.

    This is an annotation work reduction, not a definition of all protein-altering
    biology. The complete native callset remains available; novel splice reconstruction,
    regulatory effects and structural variants require their own qualified analyses.
    """
    if flank < 0:
        raise ValueError("negative annotation flank")
    spans = {}
    with pysam.FastaFile(str(fasta)) as genome, open_text(gtf) as stream:
        lengths = dict(zip(genome.references, genome.lengths, strict=True))
        for line in stream:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError("expected nine GTF fields")
            chrom, _, feature, start, end, _, _, _, raw = fields
            if feature not in {"CDS", "stop_codon"} or chrom not in AUTOSOMES:
                continue
            if attributes(raw).get("transcript_type") != ["protein_coding"]:
                continue
            left, right = int(start) - 1, int(end)
            if chrom not in lengths or not 0 <= left < right <= lengths[chrom]:
                raise ValueError("coding interval outside the reference")
            spans.setdefault(chrom, []).append(
                (max(0, left - flank), min(lengths[chrom], right + flank))
            )
    return {chrom: merged_intervals(rows) for chrom, rows in spans.items()}


def register_precomputed(
    source: PrecomputedSource, gtf: Path, fasta: Path, output: Path, flank: int = 50
) -> dict:
    """Verify source identity and genotype encodings, then prepare native annotation regions.

    No ready handoff is published here. Reference confidence, annotation and complete
    allele reconstruction still need validation. This source cannot validate itself as truth.
    """
    for item in (
        source.variants,
        source.variant_index,
        source.benchmark_regions,
        source.source_readme,
    ):
        verify(item)
    output.mkdir(parents=True, exist_ok=False)
    regions = coding_regions(gtf, fasta, flank)
    bed = "".join(
        f"{chrom}\t{start}\t{end}\n" for chrom in sorted(regions) for start, end in regions[chrom]
    )
    write_new(output / "annotation-regions.bed", bed.encode())
    counts = {"records": 0, "phased_records": 0, "records_without_dp": 0, "records_without_gq": 0}
    with pysam.VariantFile(source.variants.location) as variants:
        require_sample(variants, source.sample)
        for record in variants:
            validate_cardinalities(record, source.sample)
            call = record.samples[source.sample]
            counts["records"] += 1
            counts["phased_records"] += call.phased
            counts["records_without_dp"] += call.get("DP") is None
            counts["records_without_gq"] += call.get("GQ") is None
    result = {
        "input_mode": "precomputed_benchmark_demo",
        "source": source.model_dump(mode="json"),
        "annotation_flank_bases": flank,
        "annotation_region_bases": sum(b - a for rows in regions.values() for a, b in rows),
        "counts": counts,
        "handoff_ready": False,
        "caller_accuracy_evaluated": False,
        "expression_status": "not_measured",
    }
    write_new(output / "intake.json", canonical(result) + b"\n")
    return result
