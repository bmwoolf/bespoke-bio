"""Read native evidence with explicit missing-data semantics and batched interval operations."""

import csv
import math
from bisect import bisect_left
from pathlib import Path

import pysam

from .biology import IntervalMask, Site
from .contracts import ExpressionRow, Policy
from .references import merged_intervals


def read_expression(path: Path) -> dict[str, ExpressionRow]:
    result = {}
    with path.open() as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if reader.fieldnames != ["Name", "Length", "EffectiveLength", "TPM", "NumReads"]:
            raise ValueError("unexpected Salmon quant.sf columns")
        for row in reader:
            item = ExpressionRow(
                transcript_id=row["Name"],
                length=row["Length"],
                effective_length=row["EffectiveLength"],
                tpm=row["TPM"],
                num_reads=row["NumReads"],
            )
            if item.transcript_id in result:
                raise ValueError(f"duplicate expression row: {item.transcript_id}")
            result[item.transcript_id] = item
    if not result:
        raise ValueError("empty quantification")
    return result


def intersect(left: list[tuple[int, int]], right: list[tuple[int, int]]) -> list[tuple[int, int]]:
    result = []
    i = j = 0
    while i < len(left) and j < len(right):
        start, end = max(left[i][0], right[j][0]), min(left[i][1], right[j][1])
        if start < end:
            result.append((start, end))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return merged_intervals(result)


def subtract(
    intervals: list[tuple[int, int]], holes: list[tuple[int, int]]
) -> list[tuple[int, int]]:
    result = []
    index = 0
    for start, end in intervals:
        while index < len(holes) and holes[index][1] <= start:
            index += 1
        cursor = start
        current = index
        while current < len(holes) and holes[current][0] < end:
            left, right = holes[current]
            if left > cursor:
                result.append((cursor, min(left, end)))
            cursor = max(cursor, right)
            current += 1
        if cursor < end:
            result.append((cursor, end))
    return result


def read_bed(path: Path) -> dict[str, list[tuple[int, int]]]:
    result = {}
    with path.open() as stream:
        for line in stream:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            chrom, start, end, *_ = line.split()
            result.setdefault(chrom, []).append((int(start), int(end)))
    return {chrom: merged_intervals(spans) for chrom, spans in result.items()}


def depth_mask(path: Path, bed: dict[str, list[tuple[int, int]]], minimum: int):
    """Consume one native depth pass. Every requested base, including zero depth, must exist."""
    result = {chrom: [] for chrom in bed}
    cursors = {chrom: (0, spans[0][0]) for chrom, spans in bed.items()}
    with path.open() as stream:
        for line in stream:
            chrom, position, depth = line.rstrip().split("\t")
            position, depth = int(position) - 1, int(depth)
            if chrom not in bed or depth < 0:
                raise ValueError("unexpected depth contig or negative value")
            index, expected = cursors[chrom]
            if index >= len(bed[chrom]) or position != expected:
                raise ValueError(
                    f"missing, duplicate or out-of-order depth base: {chrom}:{position + 1}"
                )
            if depth >= minimum:
                spans = result[chrom]
                if spans and spans[-1][1] == position:
                    spans[-1] = (spans[-1][0], position + 1)
                else:
                    spans.append((position, position + 1))
            if position + 1 == bed[chrom][index][1]:
                index += 1
                expected = bed[chrom][index][0] if index < len(bed[chrom]) else -1
            else:
                expected += 1
            cursors[chrom] = (index, expected)
    if any(index != len(bed[chrom]) for chrom, (index, _) in cursors.items()):
        raise ValueError("depth output omits requested bases")
    return result


def good_call(record, sample: str, policy: Policy) -> bool:
    call = record.samples[sample]
    genotype = call.get("GT") or ()
    return (
        list(record.filter) == ["PASS"]
        and len(genotype) == 2
        and all(a is not None and 0 <= a < len(record.alleles) for a in genotype)
        and call.get("DP") is not None
        and call["DP"] >= policy.min_depth
        and call.get("GQ") is not None
        and call["GQ"] >= policy.min_genotype_quality
    )


def require_sample(stream, sample: str) -> None:
    if list(stream.header.samples) != [sample]:
        raise ValueError(f"expected one VCF sample {sample}; found {list(stream.header.samples)}")


def callability(
    original: Path, gvcf: Path, sample: str, policy: Policy, depth: dict[str, list[tuple[int, int]]]
) -> dict[str, IntervalMask]:
    """Overlay original event spans on reference blocks; failed calls cannot be rescued."""
    good, bad, reference = {}, {}, {}
    seen = set()
    with pysam.VariantFile(str(original)) as stream:
        require_sample(stream, sample)
        for record in stream:
            if record.contig not in depth:
                continue
            key = (record.contig, record.start, record.ref, record.alts)
            if key in seen:
                raise ValueError("duplicate original variant record")
            seen.add(key)
            target = good if good_call(record, sample, policy) else bad
            target.setdefault(record.contig, []).append(
                (record.start, record.start + len(record.ref))
            )
    with pysam.VariantFile(str(gvcf)) as stream:
        require_sample(stream, sample)
        for record in stream:
            if record.contig not in depth:
                continue
            call = record.samples[sample]
            is_block = all(a in {"<*>", "<NON_REF>", "."} for a in (record.alts or (".",)))
            if not is_block or call.get("GT") != (0, 0):
                continue
            if (
                list(record.filter) not in ([], ["PASS"])
                or call.get("GQ") is None
                or call["GQ"] < policy.min_genotype_quality
                or call.get("MIN_DP") is None
                or call["MIN_DP"] < policy.min_depth
            ):
                continue
            reference.setdefault(record.contig, []).append((record.start, record.stop))
    result = {}
    for chrom, sufficient_depth in depth.items():
        supported_reference = intersect(
            merged_intervals(reference.get(chrom, [])), sufficient_depth
        )
        supported = merged_intervals(supported_reference + good.get(chrom, []))
        result[chrom] = IntervalMask(subtract(supported, merged_intervals(bad.get(chrom, []))))
    return result


def crosscheck_identity(path: Path, dna_file: str, rna_file: str, threshold: float) -> float:
    """Require the actual off-diagonal normal-normal LOD, not self-comparisons or exit zero."""
    with path.open() as stream:
        rows = [line for line in stream if line.strip() and not line.startswith("#")]
    scores = []
    for row in csv.DictReader(rows, delimiter="\t"):
        groups = {row["LEFT_GROUP_VALUE"], row["RIGHT_GROUP_VALUE"]}
        if groups == {dna_file, rna_file} and dna_file != rna_file:
            score = float(row["LOD_SCORE"])
            if not math.isfinite(score) or score <= threshold:
                raise ValueError("DNA/RNA identity mismatch or insufficient fingerprint evidence")
            scores.append(score)
    if not scores:
        raise ValueError("missing off-diagonal DNA/RNA fingerprint comparison")
    return min(scores)


def csq_fields(header) -> list[str]:
    if "CSQ" not in header.info:
        raise ValueError("VEP CSQ header missing")
    description = header.info["CSQ"].description
    fields = description.split("Format: ", 1)[-1].split("|")
    if not {"ALLELE_NUM", "Feature", "Feature_type", "Consequence"} <= set(fields):
        raise ValueError("VEP consequence schema lacks required fields")
    return fields


def consequences(
    record, fields: list[str], *, allow_spanning_deletions: bool = False
) -> list[dict]:
    """Read CSQ by header field names and numeric ALT index, never by trimmed allele strings."""
    result = []
    for value in record.info.get("CSQ", ()):
        parts = value.split("|")
        if len(parts) != len(fields):
            raise ValueError("VEP consequence width differs from header")
        row = dict(zip(fields, parts, strict=True))
        index = int(row["ALLELE_NUM"])
        if not 1 <= index <= len(record.alts or ()):
            raise ValueError("VEP ALLELE_NUM outside native ALT array")
        result.append(row)
    missing = set(range(1, len(record.alts or ()) + 1)) - {int(r["ALLELE_NUM"]) for r in result}
    # VEP does not annotate VCF's '*' placeholder for an upstream deletion. Keep
    # its genotype slot intact; sequence construction must explicitly defer it.
    allowed = {
        i for i, alt in enumerate(record.alts or (), 1) if allow_spanning_deletions and alt == "*"
    }
    if missing - allowed:
        raise ValueError("VEP dropped an ALT allele")
    return result


class SiteIndex:
    """Index one chromosome's normalized calls for batched transcript interval queries."""

    def __init__(self, sites: list[Site]):
        self.sites = sorted(sites, key=lambda s: s.start)
        self.starts = [s.start for s in self.sites]
        self.max_span = max((len(s.ref) for s in sites), default=1)
        self.splice_sites: dict[str, list[Site]] = {}
        for site in sites:
            for transcript_id in site.splice_transcripts:
                self.splice_sites.setdefault(transcript_id, []).append(site)

    def overlapping(self, start: int, end: int) -> list[Site]:
        first = bisect_left(self.starts, max(0, start - self.max_span + 1))
        last = bisect_left(self.starts, end)
        return [s for s in self.sites[first:last] if s.end > start]
