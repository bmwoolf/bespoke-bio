"""Deterministic diploid allele construction. No reads, processes or cloud state here."""

from bisect import bisect_right
from dataclasses import dataclass
from hashlib import sha256

from Bio.Seq import Seq

from bespoke_bio.api import Expression, Transcript, Variant

from .contracts import ExpressionRow
from .references import CodingTranscript


class Deferred(ValueError):
    """Supported input that cannot establish a complete transcript under the pilot policy."""


@dataclass(frozen=True)
class Site:
    """One unsplit normalized VCF record; start is zero-based, genotype indexes alleles."""

    chrom: str
    start: int
    ref: str
    alts: tuple[str, ...]
    genotype: tuple[int | None, ...]
    phased: bool
    phase_set: str | None
    sample: str
    filters: tuple[str, ...]
    dp: int | None
    gq: int | None
    source_record: str
    splice_transcripts: frozenset[str] = frozenset()
    annotated_transcripts: frozenset[str] = frozenset()

    @property
    def end(self):
        return self.start + len(self.ref)

    def allele_key(self, index: int) -> str:
        return f"GRCh38:{self.chrom}:{self.start + 1}:{self.ref}:{self.alts[index - 1]}"


class IntervalMask:
    """Sorted disjoint callable intervals, supporting logarithmic full-interval checks."""

    def __init__(self, intervals: list[tuple[int, int]]):
        self.intervals = intervals
        self.starts = [left for left, _ in intervals]
        previous = -1
        for start, end in intervals:
            if start < 0 or start <= previous or end <= start:
                raise ValueError("mask must contain sorted disjoint nonadjacent intervals")
            previous = end

    def covers(self, start: int, end: int) -> bool:
        index = bisect_right(self.starts, start) - 1
        return index >= 0 and self.intervals[index][1] >= end


def identifier(key: str) -> str:
    return key if len(key) <= 200 else "allele_" + sha256(key.encode()).hexdigest()


def project(tx: CodingTranscript, site: Site, index: int) -> Variant:
    """Project a full anchored REF span inside one coding exon, reversing both minus alleles."""
    alt = site.alts[index - 1]
    if not site.ref or not alt or set(site.ref + alt) - set("ACGT"):
        raise Deferred("unsupported_allele")
    offset = 0
    for segment in tx.segments:
        if segment.start <= site.start and site.end <= segment.end:
            # Insertions anchored on either exon/CDS boundary need a splice-aware model.
            if len(alt) > len(site.ref) and (
                site.start == segment.start or site.end == segment.end
            ):
                raise Deferred("boundary_edit")
            if tx.strand == "+":
                position = offset + site.start - segment.start
                ref = site.ref
            else:
                position = offset + segment.end - site.end
                ref = str(Seq(site.ref).reverse_complement())
                alt = str(Seq(alt).reverse_complement())
            if tx.sequence[position : position + len(ref)] != ref:
                raise ValueError(
                    f"projected REF mismatch: {tx.transcript_id}/{site.allele_key(index)}"
                )
            separator = "|" if site.phased else "/"
            return Variant(
                variant_id=identifier(site.allele_key(index)),
                cds_position=position + 1,
                ref=ref,
                alt=alt,
                genotype=separator.join("." if x is None else str(x) for x in site.genotype),
                phase_set=site.phase_set,
                quality=";".join(site.filters) or ".",
            )
        offset += segment.end - segment.start
    raise Deferred("boundary_edit")


def diploid_slots(sites: list[Site]) -> list[tuple[int, int]]:
    """Use explicit genotype slots; multiple heterozygous sites require one phase block."""
    heterozygous = []
    for site in sites:
        gt = site.genotype
        if len(gt) != 2 or any(x is None or x < 0 or x > len(site.alts) for x in gt):
            raise Deferred("missing_or_nondiploid_genotype")
        if gt[0] != gt[1]:
            heterozygous.append(site)
    if len(heterozygous) > 1:
        blocks = {(s.sample, s.chrom, s.phase_set) for s in heterozygous}
        if any(not s.phased or s.phase_set is None for s in heterozygous) or len(blocks) != 1:
            raise Deferred("unresolved_phase")
    return [tuple(s.genotype) if s.phased else tuple(sorted(s.genotype)) for s in sites]


def edited_sequence(reference: str, edits: list[Variant]) -> str:
    """Apply edits independently of Step 2, preserving original coordinates after indels."""
    cursor = 0
    chunks = []
    for edit in sorted(edits, key=lambda e: e.cds_position):
        start = edit.cds_position - 1
        end = start + len(edit.ref)
        if start < cursor:
            raise Deferred("overlapping_edits")
        if start < 0 or reference[start:end] != edit.ref:
            raise ValueError("edit REF does not match reference CDS")
        chunks.extend((reference[cursor:start], edit.alt))
        cursor = end
    chunks.append(reference[cursor:])
    return "".join(chunks)


def complete_protein(sequence: str) -> str:
    """Require a complete standard-code ORF, then translate with Biopython."""
    if set(sequence) - set("ACGT"):
        raise Deferred("ambiguous_allele")
    if not sequence.startswith("ATG"):
        raise Deferred("missing_start")
    if len(sequence) % 3:
        raise Deferred("partial_codon")
    if sequence[-3:] not in {"TAA", "TAG", "TGA"}:
        raise Deferred("missing_terminal_stop")
    protein = str(Seq(sequence).translate())
    if "*" in protein[:-1]:
        raise Deferred("internal_stop")
    return protein[:-1]


def construct_alleles(
    tx: CodingTranscript, sites: list[Site], solution_id: str
) -> list[dict]:
    """Resolve each allele independently after shared phase checks; retain failed-allele evidence."""
    if any(tx.transcript_id in s.splice_transcripts for s in sites):
        raise Deferred("splice_consequence")
    coding = [s for s in sites if any(s.start < e.end and s.end > e.start for e in tx.segments)]
    slots = diploid_slots(coding)
    evidence = []
    for allele in range(2):
        phase_scenario = f"{solution_id}:allele-{allele + 1}"
        item = {
            "transcript_id": tx.transcript_id,
            "solution_id": solution_id,
            "phase_scenario": phase_scenario,
            "allele": allele + 1,
            "reference_cds_sha256": sha256(tx.sequence.encode()).hexdigest(),
            "source_records": [s.source_record for s in coding],
            "status": "deferred",
            "reason": None,
            "edits": [],
            "cds_sha256": None,
            "protein_sha256": None,
            "protein": None,
        }
        evidence.append(item)
        try:
            edits = [
                project(tx, site, gt[allele])
                for site, gt in zip(coding, slots, strict=True)
                if gt[allele] != 0
            ]
            item["edits"] = [v.model_dump(mode="json") for v in edits]
            sequence = edited_sequence(tx.sequence, edits)
            item["cds_sha256"] = sha256(sequence.encode()).hexdigest()
            protein = complete_protein(sequence)
        except Deferred as error:
            item["reason"] = str(error)
            continue
        item.update(
            status="resolved", protein=protein, protein_sha256=sha256(protein.encode()).hexdigest()
        )
    return evidence


def construct_pair(
    tx: CodingTranscript, sites: list[Site], expression: ExpressionRow, solution_id: str
) -> tuple[list[Transcript], list[dict]]:
    """Adapt resolved alleles to the current CDS-based Step 2 contract."""
    evidence = construct_alleles(tx, sites, solution_id)
    records = []
    for item in evidence:
        item["expression_transcript_id"] = expression.transcript_id
        if item["status"] != "resolved":
            continue
        records.append(
            Transcript(
                gene_id=tx.gene_id,
                gene_symbol=tx.gene_symbol,
                transcript_id=tx.transcript_id,
                haplotype="unknown",
                phase_scenario=item["phase_scenario"],
                coding_sequence=tx.sequence,
                expression=Expression(tpm=expression.tpm, raw_count=expression.num_reads),
                variants=[Variant.model_validate(edit) for edit in item["edits"]],
            )
        )
    return records, evidence
