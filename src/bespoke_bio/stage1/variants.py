"""Join original, normalized, annotated and phased calls without losing allele identity."""

import math
from itertools import zip_longest
from pathlib import Path

import pysam

from .biology import Site
from .evidence import SiteIndex, consequences, csq_fields, good_call, require_sample
from .io import canonical


def key(record) -> tuple:
    return record.contig, record.pos, record.ref, tuple(record.alts or ())


def original_key(record) -> tuple:
    """ORIG links only changed rows; unchanged records retain their native identity."""
    if "ORIG" not in record.info:
        return key(record)
    value = record.info["ORIG"]
    if isinstance(value, tuple):
        if len(value) != 1:
            raise ValueError("normalization merged multiple original records")
        value = value[0]
    parts = value.split("|")
    if len(parts) not in (4, 5):
        raise ValueError("unexpected BCFtools ORIG encoding")
    return parts[0], int(parts[1]), parts[2], tuple(parts[3].split(","))


def validate_cardinalities(record, sample: str) -> None:
    """Check allele/genotype arrays before trusting any normalization or phase adapter."""
    count = len(record.alleles)
    ploidy = len(record.samples[sample].get("GT") or ())
    for collection, declarations in [
        (record.info, record.header.info),
        (record.samples[sample], record.header.formats),
    ]:
        for field, value in collection.items():
            number = declarations[field].number
            expected = {
                "A": count - 1,
                "R": count,
                "G": math.comb(count + ploidy - 1, ploidy) if ploidy else 0,
            }.get(number)
            if expected is None or value is None or value == (None,):
                continue
            values = value if isinstance(value, tuple) else (value,)
            if len(values) != expected:
                raise ValueError(
                    f"invalid Number={number} cardinality for {field} at {key(record)}"
                )


def qualify_phase_input(source: Path, destination: Path, sample: str, policy) -> dict:
    """Keep all quality-qualified genomic calls, including noncoding phase bridges."""
    kept = rejected = 0
    with pysam.VariantFile(str(source)) as reader:
        require_sample(reader, sample)
        with pysam.VariantFile(str(destination), "wz", header=reader.header) as writer:
            for record in reader:
                validate_cardinalities(record, sample)
                if good_call(record, sample, policy):
                    writer.write(record)
                    kept += 1
                else:
                    rejected += 1
    pysam.tabix_index(str(destination), preset="vcf", force=False)
    return {"qualified": kept, "rejected": rejected}


def chromosome_sites(
    original: Path,
    normalized: Path,
    annotated: Path,
    phased: Path,
    sample: str,
    policy,
    transcript_ids: set[str],
    source_sha256: str,
    evidence_path: Path,
):
    """Yield one chromosome at a time; native callsets retain every consequence and FORMAT field."""
    with (
        pysam.VariantFile(str(original)) as raw,
        pysam.VariantFile(str(normalized)) as norm,
        pysam.VariantFile(str(annotated)) as annotations,
        pysam.VariantFile(str(phased)) as phases,
        evidence_path.open("xb") as evidence,
    ):
        for stream in (raw, norm, annotations, phases):
            require_sample(stream, sample)
        fields = csq_fields(annotations.header)
        ordinal = 0
        total = 0
        for chrom in raw.header.contigs:
            original_rows = {}
            for record in raw.fetch(chrom):
                ordinal += 1
                identity = key(record)
                if identity in original_rows:
                    raise ValueError(f"duplicate original variant: {identity}")
                validate_cardinalities(record, sample)
                original_rows[identity] = (ordinal, record)
            if not original_rows:
                continue
            phased_rows = {}
            for record in phases.fetch(chrom):
                identity = key(record)
                if identity in phased_rows:
                    raise ValueError(f"duplicate phased variant: {identity}")
                phased_rows[identity] = record
            sites = []
            seen_sources, seen_normalized, seen_phase = set(), set(), set()
            for normalized_record, annotation in zip_longest(
                norm.fetch(chrom), annotations.fetch(chrom)
            ):
                if (
                    normalized_record is None
                    or annotation is None
                    or key(normalized_record) != key(annotation)
                ):
                    raise ValueError("VEP changed, dropped or reordered native variant records")
                identity = key(normalized_record)
                if identity in seen_normalized:
                    raise ValueError("duplicate/conflicting normalized variant")
                seen_normalized.add(identity)
                origin = original_key(normalized_record)
                if origin not in original_rows or origin in seen_sources:
                    raise ValueError("normalization lost the one-to-one original-record link")
                seen_sources.add(origin)
                row_number, raw_record = original_rows[origin]
                if len(normalized_record.alts or ()) != len(raw_record.alts or ()):
                    raise ValueError("normalization split or dropped multiallelic alleles")
                native = normalized_record.samples[sample]
                for field in ("GT", "DP", "GQ", "AD", "PL"):
                    if native.get(field) != raw_record.samples[sample].get(field):
                        raise ValueError(f"normalization changed genotype evidence {field}")
                rows = consequences(annotation, fields)
                called = {a for a in (native.get("GT") or ()) if a is not None and a > 0}
                relevant = [
                    r
                    for r in rows
                    if r["Feature_type"] == "Transcript" and r["Feature"] in transcript_ids
                ]
                splice = frozenset(
                    r["Feature"]
                    for r in relevant
                    if int(r["ALLELE_NUM"]) in called
                    and any(term.startswith("splice_") for term in r["Consequence"].split("&"))
                )
                selected = normalized_record
                if good_call(normalized_record, sample, policy):
                    if identity not in phased_rows:
                        raise ValueError("qualified call missing from phase output")
                    selected = phased_rows[identity]
                    seen_phase.add(identity)
                    if sorted(selected.samples[sample]["GT"]) != sorted(native["GT"]):
                        raise ValueError("phasing changed the genotype")
                call = selected.samples[sample]
                source_id = f"{source_sha256}:{row_number}"
                sites.append(
                    Site(
                        chrom,
                        selected.start,
                        selected.ref,
                        tuple(selected.alts or ()),
                        tuple(call.get("GT") or ()),
                        call.phased,
                        str(call["PS"]) if call.get("PS") is not None else None,
                        sample,
                        tuple(selected.filter),
                        call.get("DP"),
                        call.get("GQ"),
                        source_id,
                        splice,
                        frozenset(r["Feature"] for r in relevant),
                    )
                )
                if relevant:
                    evidence.write(
                        canonical(
                            {
                                "original_record": source_id,
                                "original_key": origin,
                                "normalized_key": identity,
                                "qual": raw_record.qual,
                                "filter": list(raw_record.filter) or ["."],
                                "sample": dict(raw_record.samples[sample]),
                                "phased": call.phased,
                                "phase_set": call.get("PS"),
                                "consequences": relevant,
                            }
                        )
                        + b"\n"
                    )
                total += 1
            if seen_sources != set(original_rows) or seen_phase != set(phased_rows):
                raise ValueError("variant/phase record reconciliation failed")
            yield chrom, SiteIndex(sites)
        if total != ordinal:
            raise ValueError("original/normalized/annotated record totals differ")
