"""Extract coding-orientation sequences from one exact primary-assembly GTF/FASTA pair."""

import gzip
import re
from dataclasses import asdict, dataclass, field
from hashlib import sha256
from pathlib import Path

import pysam
from Bio.Seq import Seq

from .io import decode_json, write_jsonl, write_new

AUTOSOMES = frozenset(f"chr{i}" for i in range(1, 23))


@dataclass(frozen=True, order=True)
class Segment:
    """A zero-based half-open genomic CDS or stop interval in one annotated exon."""

    start: int
    end: int
    exon: str
    phase: int | None = None


@dataclass
class CodingTranscript:
    """Versioned GTF identity with intervals ordered in coding orientation."""

    transcript_id: str
    gene_id: str
    gene_symbol: str
    biotype: str
    chrom: str
    strand: str
    tags: set[str] = field(default_factory=set)
    cds: list[Segment] = field(default_factory=list)
    stops: list[Segment] = field(default_factory=list)
    segments: list[Segment] = field(default_factory=list)
    sequence: str = ""
    reference_reasons: list[str] = field(default_factory=list)

    @property
    def scope_reason(self) -> str | None:
        if self.biotype != "protein_coding":
            return "biotype"
        if self.chrom not in AUTOSOMES:
            return "contig"
        if "selenoprotein" in self.tags or "Selenocysteine" in self.tags:
            return "unsupported_genetic_code"
        return None

    @property
    def mane(self) -> bool:
        return "MANE_Select" in self.tags

    def record(self) -> dict:
        result = asdict(self)
        result["tags"] = sorted(self.tags)
        result["sequence_sha256"] = sha256(self.sequence.encode()).hexdigest()
        return result


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def attributes(text: str) -> dict[str, list[str]]:
    result: dict[str, list[str]] = {}
    for key, value in re.findall(r'(\S+)\s+"([^"]*)"\s*;', text):
        result.setdefault(key, []).append(value)
    if not result:
        raise ValueError("malformed GTF attributes")
    return result


def read_gtf(path: Path) -> dict[str, CodingTranscript]:
    """Read every transcript; retain exclusions instead of shrinking the input universe."""
    result: dict[str, CodingTranscript] = {}
    with open_text(path) as stream:
        for number, line in enumerate(stream, 1):
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"GTF line {number}: expected nine fields")
            chrom, _, feature, start, end, _, strand, phase, raw = fields
            if feature not in {"transcript", "CDS", "stop_codon", "Selenocysteine"}:
                continue
            attrs = attributes(raw)
            txid = attrs["transcript_id"][0]
            tx = result.get(txid)
            identity = (
                attrs["gene_id"][0],
                attrs["gene_name"][0],
                attrs["transcript_type"][0],
                chrom,
                strand,
            )
            if tx is None:
                if strand not in {"+", "-"}:
                    raise ValueError(f"invalid strand: {txid}")
                tx = CodingTranscript(txid, *identity)
                result[txid] = tx
            if identity != (tx.gene_id, tx.gene_symbol, tx.biotype, tx.chrom, tx.strand):
                raise ValueError(f"conflicting GTF identity: {txid}")
            tx.tags.update(attrs.get("tag", []))
            if feature == "Selenocysteine":
                tx.tags.add(feature)
            if feature in {"CDS", "stop_codon"}:
                left, right = int(start) - 1, int(end)
                if left < 0 or right <= left:
                    raise ValueError(f"invalid GTF interval: {txid}")
                exon = attrs.get("exon_id", attrs.get("exon_number"))
                if not exon:
                    raise ValueError(f"missing exon identity: {txid}")
                segment = Segment(left, right, exon[0], None if phase == "." else int(phase))
                (tx.cds if feature == "CDS" else tx.stops).append(segment)
    return result


def coding_segments(tx: CodingTranscript) -> list[Segment]:
    """Union CDS and stop bases within each exon, preserving exon boundaries."""
    result: list[Segment] = []
    for part in sorted(tx.cds + tx.stops):
        if result and part.start < result[-1].end and part.exon != result[-1].exon:
            raise ValueError(f"overlapping exons: {tx.transcript_id}")
        if result and part.exon == result[-1].exon and part.start <= result[-1].end:
            previous = result.pop()
            result.append(Segment(previous.start, max(previous.end, part.end), part.exon))
        else:
            result.append(Segment(part.start, part.end, part.exon))
    return sorted(result, reverse=tx.strand == "-")


def protein_sequences(path: Path) -> dict[str, str]:
    """Read GENCODE protein FASTA, keyed by its versioned transcript identifier."""
    result = {}
    with pysam.FastxFile(str(path)) as stream:
        for entry in stream:
            fields = entry.name.split("|")
            # GENCODE fixes the versioned transcript ID in column two. Later name
            # fields may repeat an unversioned ENST identifier.
            ids = fields[1:2]
            if not ids or not re.fullmatch(r"ENST[0-9]+\.[0-9]+(?:_PAR_Y)?", ids[0]):
                raise ValueError(f"invalid protein transcript identity: {entry.name}")
            if len(ids) != 1 or ids[0] in result:
                raise ValueError(f"ambiguous protein/transcript identity: {entry.name}")
            result[ids[0]] = entry.sequence.rstrip("*")
    return result


def reference_sequence(tx: CodingTranscript, genome: pysam.FastaFile) -> str:
    parts = []
    for segment in tx.segments:
        value = genome.fetch(tx.chrom, segment.start, segment.end).upper()
        if len(value) != segment.end - segment.start:
            raise ValueError(f"reference interval outside contig: {tx.transcript_id}")
        parts.append(str(Seq(value).reverse_complement()) if tx.strand == "-" else value)
    return "".join(parts)


def validate_reference(tx: CodingTranscript, expected: str | None) -> list[str]:
    """Classify unsupported reference ORFs; unexpected translation disagreement is an error."""
    reasons = []
    if tx.tags & {"cds_start_NF", "cds_end_NF"} or not tx.cds:
        return ["incomplete_reference"]
    count = 0
    for segment in sorted(tx.cds, reverse=tx.strand == "-"):
        if segment.phase != (3 - count % 3) % 3:
            raise ValueError(f"inconsistent CDS phase: {tx.transcript_id}")
        count += segment.end - segment.start
    sequence = tx.sequence
    if set(sequence) - set("ACGT"):
        reasons.append("ambiguous_reference")
    if not sequence.startswith("ATG"):
        reasons.append("unsupported_start")
    if len(sequence) % 3 or sequence[-3:] not in {"TAA", "TAG", "TGA"}:
        reasons.append("incomplete_reference")
    if reasons:
        return reasons
    actual = str(Seq(sequence).translate(to_stop=False))
    if "*" in actual[:-1]:
        return ["reference_internal_stop"]
    if expected is None or actual[:-1] != expected:
        raise ValueError(f"GENCODE protein disagreement: {tx.transcript_id}")
    return []


def prepare_references(gtf: Path, fasta: Path, proteins: Path) -> dict[str, CodingTranscript]:
    transcripts = read_gtf(gtf)
    expected = protein_sequences(proteins)
    with pysam.FastaFile(str(fasta)) as genome:
        for tx in transcripts.values():
            if tx.chrom not in genome.references:
                raise ValueError(f"annotation contig absent from reference: {tx.chrom}")
            if tx.scope_reason:
                continue
            tx.segments = coding_segments(tx)
            tx.sequence = reference_sequence(tx, genome)
            tx.reference_reasons = validate_reference(tx, expected.get(tx.transcript_id))
    return transcripts


def merged_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Return a sorted union; adjacent intervals form one continuous covered interval."""
    result = []
    for start, end in sorted(intervals):
        if start < 0 or end <= start:
            raise ValueError("invalid half-open interval")
        if result and start <= result[-1][1]:
            result[-1] = (result[-1][0], max(result[-1][1], end))
        else:
            result.append((start, end))
    return result


def write_reference_records(transcripts: dict[str, CodingTranscript], output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    write_jsonl(
        output / "transcripts.jsonl",
        (tx.record() for tx in sorted(transcripts.values(), key=lambda tx: tx.transcript_id)),
    )
    by_chrom: dict[str, list[tuple[int, int]]] = {}
    for tx in transcripts.values():
        if not tx.scope_reason:
            by_chrom.setdefault(tx.chrom, []).extend((s.start, s.end) for s in tx.segments)
    lines = [
        f"{chrom}\t{start}\t{end}\n"
        for chrom in sorted(by_chrom)
        for start, end in merged_intervals(by_chrom[chrom])
    ]
    write_new(output / "merged_cds.bed", "".join(lines).encode())


def read_reference_records(path: Path) -> dict[str, CodingTranscript]:
    """Restore registered reference records, preserving exact sequence and interval identities."""
    result = {}
    required = set(CodingTranscript.__dataclass_fields__) | {"sequence_sha256"}
    with path.open() as stream:
        for line in stream:
            item = decode_json(line)
            if set(item) != required:
                raise ValueError("reference record has missing or unknown fields")
            expected = item.pop("sequence_sha256")
            if sha256(item["sequence"].encode()).hexdigest() != expected:
                raise ValueError("reference sequence hash mismatch")
            item["tags"] = set(item["tags"])
            for field in ("cds", "stops", "segments"):
                item[field] = [Segment(**part) for part in item[field]]
            tx = CodingTranscript(**item)
            if tx.transcript_id in result:
                raise ValueError("duplicate reference transcript")
            result[tx.transcript_id] = tx
    if not result:
        raise ValueError("empty transcript reference")
    return result


def compare_gffread(transcripts: dict[str, CodingTranscript], cds: Path) -> dict:
    """Compare every supported complete reference CDS with independently extracted gffread CDS."""
    expected = {
        name for name, tx in transcripts.items() if not tx.scope_reason and not tx.reference_reasons
    }
    seen = set()
    with pysam.FastxFile(str(cds)) as records:
        for record in records:
            if record.name not in expected:
                continue
            if record.name in seen:
                raise ValueError("duplicate native CDS record")
            if record.sequence.upper() != transcripts[record.name].sequence:
                raise ValueError(f"native gffread/reference CDS disagreement: {record.name}")
            seen.add(record.name)
    if seen != expected:
        raise ValueError(f"native CDS extraction omitted {len(expected - seen)} supported records")
    return {"complete_reference_cds": len(expected), "matching_native_cds": len(seen)}
