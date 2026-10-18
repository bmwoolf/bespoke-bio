"""Export complete protein sequences from confidence-scoped, precomputed human variants.

This export has no expression requirement and does not call a folding service. The
current CDS-based Step 2 API remains a separate contract. Native variants and allele
edits stay in the evidence directory, outside the sequence payload.
"""

import os
from collections import Counter
from contextlib import ExitStack
from hashlib import sha256
from itertools import zip_longest
from pathlib import Path

import pysam
from pydantic import model_validator

from bespoke_bio.api import Variant, apply_variants, translate_cds

from .biology import Deferred, IntervalMask, Site, construct_alleles
from .contracts import Artifact, Digest, Identifier, Record
from .evidence import SiteIndex, consequences, csq_fields, read_bed, require_sample
from .io import (
    artifact,
    canonical,
    decode_json,
    digest,
    file_hash,
    read_json,
    source_hash,
    sync_directory,
    verify,
    write_new,
)
from .precomputed import PrecomputedSource, phase_identity
from .references import read_reference_records
from .variants import key


class ProteinTarget(Record):
    """One biological allele with its actual, complete amino-acid sequence.

    Expression is unknown in this input mode. Sequence hashes deduplicate computation;
    target records retain the donor, transcript and allele scenario independently.
    """

    target_id: Identifier
    donor: Identifier
    gene_id: str
    gene_symbol: str
    transcript_id: str
    phase_scenario: str
    amino_acid_sequence: str
    sequence_sha256: Digest
    reference_sequence_sha256: Digest
    expression: None = None

    @model_validator(mode="after")
    def complete_sequence(self):
        if not self.amino_acid_sequence or set(self.amino_acid_sequence) - set(
            "ACDEFGHIKLMNPQRSTVWY"
        ):
            raise ValueError("expected a nonempty complete standard-code protein")
        if sha256(self.amino_acid_sequence.encode()).hexdigest() != self.sequence_sha256:
            raise ValueError("protein sequence hash mismatch")
        return self


def annotated_sites(normalized: Path, annotated: Path, source: PrecomputedSource, evidence):
    """Require every selected native record and genotype to survive VEP unchanged."""
    by_chrom = {}
    with pysam.VariantFile(str(normalized)) as native, pysam.VariantFile(str(annotated)) as vep:
        require_sample(native, source.sample)
        require_sample(vep, source.sample)
        fields = csq_fields(vep.header)
        seen = set()
        for row, annotation in zip_longest(native, vep):
            if row is None or annotation is None or key(row) != key(annotation):
                raise ValueError("VEP changed, dropped or reordered selected records")
            identity = key(row)
            if identity in seen:
                raise ValueError("duplicate normalized variant")
            seen.add(identity)
            call = row.samples[source.sample]
            other = annotation.samples[source.sample]
            if dict(call) != dict(other) or call.phased != other.phased:
                raise ValueError("VEP changed genotype evidence")
            rows = consequences(annotation, fields, allow_spanning_deletions=True)
            if any(r.get("CHECK_REF") for r in rows):
                raise ValueError("VEP reference check failed")
            called = {a for a in (call.get("GT") or ()) if a is not None and a > 0}
            relevant = [r for r in rows if r["Feature_type"] == "Transcript"]
            splice = frozenset(
                r["Feature"]
                for r in relevant
                if int(r["ALLELE_NUM"]) in called
                and any(t.startswith("splice_") for t in r["Consequence"].split("&"))
            )
            source_record = source.variants.sha256 + ":" + digest(identity)
            site = Site(
                row.contig,
                row.start,
                row.ref,
                tuple(row.alts or ()),
                tuple(call.get("GT") or ()),
                call.phased,
                phase_identity(row, source.sample, source.variants.sha256),
                source.sample,
                tuple(row.filter),
                call.get("DP"),
                call.get("GQ"),
                source_record,
                splice,
                frozenset(r["Feature"] for r in relevant),
            )
            by_chrom.setdefault(row.contig, []).append(site)
            evidence.write(
                canonical(
                    {
                        "source_record": source_record,
                        "normalized_key": identity,
                        "original_key": row.info.get("ORIG") if "ORIG" in row.info else None,
                        "sample": dict(call),
                        "phased": call.phased,
                        "phase_set": site.phase_set,
                        "filter": list(row.filter),
                        "annotated_transcript_count": len(site.annotated_transcripts),
                        "splice_transcripts": sorted(site.splice_transcripts),
                        "spanning_deletion_alleles": [
                            i for i, alt in enumerate(row.alts or (), 1) if alt == "*"
                        ],
                    }
                )
                + b"\n"
            )
    return {chrom: SiteIndex(sites) for chrom, sites in by_chrom.items()}


def transcript_sites(tx, index):
    """Query coding exons and annotated splice effects, without scanning a whole chromosome."""
    sites = set(index.splice_sites.get(tx.transcript_id, []))
    for segment in tx.segments:
        sites.update(index.overlapping(segment.start, segment.end))
    return sorted(sites, key=lambda s: (s.start, s.ref, s.alts))


def eligibility(tx, sites, mask):
    """Keep source confidence distinct from read-based DP/GQ qualification."""
    if not tx.segments or not tx.sequence:
        raise ValueError("missing validated reference CDS")
    if not all(mask.covers(s.start, s.end) for s in tx.segments):
        return "outside_source_confidence"
    for site in sites:
        if set(site.filters) - {"PASS"}:
            return "source_filtered_variant"
        if not mask.covers(site.start, site.end):
            return "variant_outside_source_confidence"
        if (
            any(
                a is not None and 0 < a <= len(site.alts) and set(site.alts[a - 1]) <= set("ACGT")
                for a in site.genotype
            )
            and any(site.start < s.end and site.end > s.start for s in tx.segments)
            and tx.transcript_id not in site.annotated_transcripts
        ):
            raise ValueError(
                f"missing versioned coding annotation: {tx.transcript_id} at "
                f"{site.chrom}:{site.start + 1}"
            )
    return None


def protein_target(tx, allele, donor):
    """Check the independent Step 2 edit/translation implementation before export."""
    edits = [Variant.model_validate(edit) for edit in allele["edits"]]
    cds = apply_variants(tx.sequence, edits)
    protein, warnings = translate_cds(cds)
    if (
        warnings
        or protein != allele["protein"]
        or sha256(cds.encode()).hexdigest() != allele["cds_sha256"]
    ):
        raise ValueError("independent protein reconstruction disagreement")
    reference, warnings = translate_cds(tx.sequence)
    if warnings:
        raise ValueError("invalid reference translation")
    return ProteinTarget(
        target_id="protein_" + digest([donor, allele["phase_scenario"]])[:32],
        donor=donor,
        gene_id=tx.gene_id,
        gene_symbol=tx.gene_symbol,
        transcript_id=tx.transcript_id,
        phase_scenario=allele["phase_scenario"],
        amino_acid_sequence=protein,
        sequence_sha256=allele["protein_sha256"],
        reference_sequence_sha256=sha256(reference.encode()).hexdigest(),
    )


def export_proteins(
    source: PrecomputedSource,
    records: Path,
    normalized: Path,
    annotated: Path,
    confidence: Path,
    output: Path,
) -> dict:
    """Export every supported allele and account for every reference transcript.

    A final manifest appears only after all checks and files are complete. Interrupted
    attempts have no manifest and are never reused as finished sequence packages.
    """
    for item in (
        source.variants,
        source.variant_index,
        source.benchmark_regions,
        source.source_readme,
    ):
        verify(item)
    inputs = [
        artifact(records, "validated_reference_records", "JSONL"),
        artifact(normalized, "coding_context", "BCF"),
        artifact(annotated, "consequences", "VCF"),
        artifact(confidence, "qualified_confidence_regions", "BED"),
    ]
    transcripts = read_reference_records(records)
    original = read_bed(Path(source.benchmark_regions.location))
    regions = read_bed(confidence)
    for chrom, spans in regions.items():
        declared = IntervalMask(original.get(chrom, []))
        if not all(declared.covers(start, end) for start, end in spans):
            raise ValueError("qualified confidence exceeds source confidence")
    masks = {chrom: IntervalMask(spans) for chrom, spans in regions.items()}
    output.mkdir(parents=True, exist_ok=False)
    evidence_dir = output / "evidence"
    evidence_dir.mkdir()
    counts, reasons, allele_reasons, sequences = Counter(), Counter(), Counter(), {}
    with ExitStack() as stack:
        variants = stack.enter_context((evidence_dir / "variants.jsonl").open("xb"))
        indices = annotated_sites(normalized, annotated, source, variants)
        targets = stack.enter_context((output / "protein-targets.jsonl").open("xb"))
        alleles = stack.enter_context((evidence_dir / "alleles.jsonl").open("xb"))
        dispositions = stack.enter_context((evidence_dir / "dispositions.jsonl").open("xb"))
        empty = SiteIndex([])
        for tx in transcripts.values():
            status, reason, released = "deferred", None, 0
            transcript_reasons = []
            if tx.scope_reason:
                status, reason = "out_of_scope", tx.scope_reason
            elif tx.reference_reasons:
                reason = tx.reference_reasons[0]
                transcript_reasons.extend(tx.reference_reasons)
            else:
                sites = transcript_sites(tx, indices.get(tx.chrom, empty))
                reason = eligibility(tx, sites, masks.get(tx.chrom, IntervalMask([])))
                if reason is None:
                    solution = (
                        "solution_"
                        + digest(
                            [
                                source.variants.sha256,
                                tx.transcript_id,
                                sha256(tx.sequence.encode()).hexdigest(),
                            ]
                        )[:24]
                    )
                    try:
                        pair = construct_alleles(tx, sites, solution)
                    except Deferred as error:
                        reason = str(error)
                    else:
                        for allele in pair:
                            alleles.write(
                                canonical({k: v for k, v in allele.items() if k != "protein"})
                                + b"\n"
                            )
                            if allele["status"] != "resolved":
                                allele_reasons[allele["reason"]] += 1
                                transcript_reasons.append(allele["reason"])
                                continue
                            target = protein_target(tx, allele, source.donor)
                            targets.write(canonical(target.model_dump(mode="json")) + b"\n")
                            sequences[target.sequence_sha256] = target.amino_acid_sequence
                            released += 1
                        status = (
                            "selected"
                            if released == 2
                            else "partially_resolved"
                            if released
                            else "deferred"
                        )
            counts[status] += 1
            counts["target_records"] += released
            if reason:
                transcript_reasons.append(reason)
            transcript_reasons = sorted(set(transcript_reasons))
            reasons.update(transcript_reasons)
            dispositions.write(
                canonical(
                    {
                        "transcript_id": tx.transcript_id,
                        "status": status,
                        "reason": reason or (transcript_reasons[0] if transcript_reasons else None),
                        "reasons": transcript_reasons,
                        "released_alleles": released,
                    }
                )
                + b"\n"
            )
        for stream in (variants, targets, alleles, dispositions):
            stream.flush()
            os.fsync(stream.fileno())
    if sum(
        counts[s] for s in ("out_of_scope", "selected", "partially_resolved", "deferred")
    ) != len(transcripts):
        raise ValueError("transcript disposition accounting failed")
    if not sequences:
        raise ValueError("no complete supported proteins")
    fasta = "".join(
        f">sha256:{identity}\n" + "\n".join(seq[i : i + 80] for i in range(0, len(seq), 80)) + "\n"
        for identity, seq in sorted(sequences.items())
    )
    write_new(output / "proteins.fasta", fasta.encode())
    artifacts = [
        artifact(
            path, str(path.relative_to(output)), "FASTA" if path.suffix == ".fasta" else "JSONL"
        ).model_dump(mode="json")
        for path in sorted(output.rglob("*"))
        if path.is_file()
    ]
    result = {
        "schema_version": "stage1-proteins-1",
        "status": "sequence_export_complete",
        "input_mode": "precomputed_benchmark_demo",
        "source": source.model_dump(mode="json"),
        "inputs": [item.model_dump(mode="json") for item in inputs],
        "adapter_sha256": source_hash(Path(__file__).parent),
        "step2_source_sha256": file_hash(Path(__file__).parents[1] / "api.py"),
        "counts": dict(counts),
        "reason_counts": dict(reasons),
        "reason_count_scope": "transcripts; multiple reasons can apply",
        "allele_reason_counts": dict(allele_reasons),
        "unique_sequences": len(sequences),
        "reference_transcripts": len(transcripts),
        "expression_status": "not_measured",
        "caller_accuracy_evaluated": False,
        "existing_step2_api_compatible": False,
        "artifacts": artifacts,
        "scope": "Autosomal annotated standard-code complete CDS; source small-variant confidence regions. Splice consequences, unresolved phase, structural and incomplete ORFs are deferred.",
    }
    staging = output / ".manifest.partial"
    write_new(staging, canonical(result) + b"\n")
    os.link(staging, output / "manifest.json")
    staging.unlink()
    sync_directory(output)
    return result


def validate_protein_package(root: Path) -> dict:
    """Recheck a completed or restored package without requiring its original machine paths."""
    root = root.resolve()
    manifest = read_json(root / "manifest.json")
    if (
        manifest["schema_version"] != "stage1-proteins-1"
        or manifest["status"] != "sequence_export_complete"
    ):
        raise ValueError("not a completed protein package")
    for raw in manifest["artifacts"]:
        item = Artifact.model_validate(raw)
        path = (root / item.role).resolve()
        if not path.is_relative_to(root):
            raise ValueError("artifact path escapes package")
        verify(item.model_copy(update={"location": str(path)}))
    sequences = set()
    with pysam.FastxFile(str(root / "proteins.fasta")) as fasta:
        for record in fasta:
            identity = sha256(record.sequence.encode()).hexdigest()
            if record.name != "sha256:" + identity or identity in sequences:
                raise ValueError("FASTA identity mismatch or duplicate")
            sequences.add(identity)
    targets, seen_ids, seen_sequences = {}, set(), set()
    with (root / "protein-targets.jsonl").open() as stream:
        for line in stream:
            target = ProteinTarget.model_validate(decode_json(line))
            if target.target_id in seen_ids or target.phase_scenario in targets:
                raise ValueError("duplicate biological target")
            if target.sequence_sha256 not in sequences:
                raise ValueError("target missing from FASTA")
            seen_ids.add(target.target_id)
            seen_sequences.add(target.sequence_sha256)
            targets[target.phase_scenario] = target.sequence_sha256
    if len(targets) != manifest["counts"]["target_records"] or seen_sequences != sequences:
        raise ValueError("target/FASTA accounting mismatch")
    if len(sequences) != manifest["unique_sequences"]:
        raise ValueError("unique sequence count mismatch")
    resolved = {}
    with (root / "evidence/alleles.jsonl").open() as stream:
        for line in stream:
            item = decode_json(line)
            if item["status"] == "resolved":
                scenario = item["phase_scenario"]
                if scenario in resolved:
                    raise ValueError("duplicate allele evidence")
                resolved[scenario] = item["protein_sha256"]
    if resolved != targets:
        raise ValueError("allele/target reconciliation failed")
    transcripts, released = set(), 0
    with (root / "evidence/dispositions.jsonl").open() as stream:
        for line in stream:
            item = decode_json(line)
            if item["transcript_id"] in transcripts:
                raise ValueError("duplicate transcript disposition")
            transcripts.add(item["transcript_id"])
            released += item["released_alleles"]
            if item["status"] != "selected" and not item["reason"]:
                raise ValueError("unresolved transcript lacks a disposition reason")
    if len(transcripts) != manifest["reference_transcripts"] or released != len(targets):
        raise ValueError("transcript disposition accounting mismatch")
    return {
        "target_records": len(targets),
        "unique_sequences": len(sequences),
        "status": "verified",
    }
