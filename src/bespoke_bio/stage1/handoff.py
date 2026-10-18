"""Gate transcript pairs, preserve dispositions and publish the existing Step 2 contract."""

from collections import Counter
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path

from bespoke_bio.api import Reference, TranslateRequest, translate_case

from .biology import Deferred, IntervalMask, construct_pair
from .contracts import Artifact, Check, Disposition, ExpressionRow, Manifest, Policy
from .evidence import SiteIndex
from .io import artifact, canonical, digest, publish_ready, verify, write_jsonl, write_new
from .references import CodingTranscript

PHASE_REASONS = {"unresolved_phase", "missing_or_nondiploid_genotype"}
PROJECTION_REASONS = {"boundary_edit", "unsupported_allele", "overlapping_edits"}


@dataclass
class Evaluation:
    """Released alleles plus complete disposition accounting for the fixed GTF universe."""

    request: TranslateRequest | None
    dispositions: list[Disposition]
    alleles: list[dict]
    summary: dict


def reason_priority(reason: str) -> int:
    if reason == "coverage":
        return 1
    if reason == "splice_consequence":
        return 2
    if reason in PHASE_REASONS:
        return 3
    if reason in PROJECTION_REASONS:
        return 4
    return 5


def evaluate(
    case_id: str,
    reference: Reference,
    transcripts: dict[str, CodingTranscript],
    expression: dict[str, ExpressionRow],
    masks: dict[str, IntervalMask],
    chromosomes,
    policy: Policy,
) -> Evaluation:
    """Evaluate every in-scope transcript; release all supported results regardless of TPM."""
    dispositions = {}
    candidates = []
    all_alleles = []
    eligible_expressed = expressed = 0
    by_chrom: dict[str, list[CodingTranscript]] = {}
    for tx in transcripts.values():
        if tx.scope_reason:
            dispositions[tx.transcript_id] = Disposition(
                transcript_id=tx.transcript_id,
                status="out_of_scope",
                reasons=[tx.scope_reason],
                primary_reason=tx.scope_reason,
            )
        else:
            by_chrom.setdefault(tx.chrom, []).append(tx)
    seen_chromosomes = set()
    for chrom, index in chromosomes:
        if chrom in seen_chromosomes:
            raise ValueError("chromosome evaluated twice")
        seen_chromosomes.add(chrom)
        for tx in by_chrom.get(chrom, []):
            if tx.transcript_id not in expression:
                raise ValueError(f"missing versioned expression join: {tx.transcript_id}")
            row = expression[tx.transcript_id]
            expressed += row.tpm > 0
            if tx.reference_reasons:
                dispositions[tx.transcript_id] = Disposition(
                    transcript_id=tx.transcript_id,
                    status="deferred",
                    reasons=tx.reference_reasons,
                    primary_reason=tx.reference_reasons[0],
                    not_evaluated=["coverage", "splice", "phase", "projection", "final_orf"],
                )
                continue
            if not tx.segments or not tx.sequence:
                raise ValueError("reference record lacks validated coding sequence")
            left = min(s.start for s in tx.segments)
            right = max(s.end for s in tx.segments)
            sites = set(index.overlapping(left, right))
            sites.update(index.splice_sites.get(tx.transcript_id, []))
            sites = sorted(sites, key=lambda s: (s.start, s.ref, s.alts))
            for site in sites:
                if (
                    any(a not in (0, None) for a in site.genotype)
                    and any(site.start < s.end and site.end > s.start for s in tx.segments)
                    and tx.transcript_id not in site.annotated_transcripts
                ):
                    raise ValueError(f"unexplained VEP transcript join: {tx.transcript_id}")
            reasons = []
            not_evaluated = []
            mask = masks.get(chrom, IntervalMask([]))
            if not all(mask.covers(s.start, s.end) for s in tx.segments):
                reasons.append("coverage")
            if any(tx.transcript_id in s.splice_transcripts for s in sites):
                reasons.append("splice_consequence")
            solution = (
                "solution_"
                + digest(
                    {
                        "case": case_id,
                        "transcript": tx.transcript_id,
                        "reference": reference.model_dump(),
                        "sites": [
                            (s.source_record, s.genotype, s.phased, s.phase_set) for s in sites
                        ],
                    }
                )[:24]
            )
            records = allele_evidence = None
            try:
                records, allele_evidence = construct_pair(tx, sites, row, solution)
            except Deferred as error:
                reason = str(error)
                if reason not in reasons:
                    reasons.append(reason)
                if reason in PHASE_REASONS:
                    not_evaluated = ["projection", "final_orf"]
                elif reason in PROJECTION_REASONS or reason == "splice_consequence":
                    not_evaluated = ["final_orf"]
            if allele_evidence is not None:
                if reasons:
                    # Shared evidence gates constrain both alleles, even when translation works.
                    for item in allele_evidence:
                        item["status"] = "deferred"
                        item["reason"] = min(reasons, key=reason_priority)
                all_alleles.extend(allele_evidence)
            if not reasons and records is not None and not records:
                reasons = sorted({item["reason"] for item in allele_evidence}, key=reason_priority)
            if reasons:
                reasons.sort(key=reason_priority)
                dispositions[tx.transcript_id] = Disposition(
                    transcript_id=tx.transcript_id,
                    status="deferred",
                    reasons=reasons,
                    primary_reason=reasons[0],
                    not_evaluated=not_evaluated,
                )
                continue
            if records is None or allele_evidence is None:
                raise ValueError("eligible transcript has no allele evidence")
            eligible_expressed += row.tpm > 0
            rank = (not any(r.variants for r in records), not tx.mane, -row.tpm, tx.transcript_id)
            candidates.append((rank, records, allele_evidence))
    missing_chromosomes = set(by_chrom) - seen_chromosomes
    if missing_chromosomes:
        # Absence of variants is valid, but callers must explicitly supply an empty index.
        raise ValueError(f"missing chromosome evidence: {sorted(missing_chromosomes)}")
    selected_records = []
    for _, records, evidence in sorted(candidates, key=lambda candidate: candidate[0]):
        txid = records[0].transcript_id
        reasons = sorted({item["reason"] for item in evidence if item["status"] != "resolved"})
        dispositions[txid] = Disposition(
            transcript_id=txid,
            status="selected" if len(records) == 2 else "partially_resolved",
            reasons=reasons,
            primary_reason=reasons[0] if reasons else None,
        )
        selected_records.extend(records)
    if set(dispositions) != set(transcripts):
        raise ValueError("dispositions do not partition the primary-GTF universe")
    request = (
        TranslateRequest(case_id=case_id, reference=reference, transcripts=selected_records)
        if selected_records
        else None
    )
    counts = dict(Counter(d.status for d in dispositions.values()))
    reasons = dict(Counter(reason for d in dispositions.values() for reason in d.reasons))
    summary = {
        "transcript_counts": counts,
        "allele_counts": dict(Counter(item["status"] for item in all_alleles)),
        "reason_counts": reasons,
        "in_scope_expressed": expressed,
        "released_expressed": eligible_expressed,
        "released_expressed_fraction": eligible_expressed / expressed if expressed else None,
    }
    return Evaluation(request, list(dispositions.values()), all_alleles, summary)


def empty_chromosomes(chromosomes, expected: set[str]):
    """Account explicitly for chromosomes with no calls rather than inventing reference evidence."""
    seen = set()
    for chrom, index in chromosomes:
        seen.add(chrom)
        yield chrom, index
    for chrom in sorted(expected - seen):
        yield chrom, SiteIndex([])


def preflight(request: TranslateRequest, alleles: list[dict]):
    """Compare Step 2 translation against independently constructed complete protein sequences."""
    response = translate_case(request)
    alleles = [item for item in alleles if item["status"] == "resolved"]
    if len(response.targets) != len(alleles):
        raise ValueError("allele/target count mismatch")
    if len({t.target_id for t in response.targets}) != len(response.targets):
        raise ValueError("duplicate Step 2 target identity")
    expected = {e["phase_scenario"]: e for e in alleles}
    if len(expected) != len(alleles):
        raise ValueError("duplicate allele scenario identity")
    for target in response.targets:
        allele = expected[target.phase_scenario]
        if (
            target.warnings
            or target.amino_acid_sequence != allele["protein"]
            or target.sequence_hash != "sha256:" + allele["protein_sha256"]
        ):
            raise ValueError("independent sequence/Step 2 translation disagreement")
    return response


def publish(
    evaluation: Evaluation,
    attempt: Path,
    ready: Path,
    identities: dict,
    evidence: list[Artifact],
    checks: list[Check],
    durable_roots: tuple[Path, ...] = (),
) -> Manifest:
    """Flush immutable evidence and publish readiness last; no_targets never produces ready."""
    attempt.mkdir(parents=True, exist_ok=False)
    for item in evidence:
        verify(item)
        roots = durable_roots or (attempt.parent,)
        if not any(Path(item.location).resolve().is_relative_to(root.resolve()) for root in roots):
            raise ValueError("published evidence must reside in the declared durable storage roots")
    dispositions = attempt / "dispositions.jsonl"
    write_jsonl(dispositions, (d.model_dump(mode="json") for d in evaluation.dispositions))
    allele_path = attempt / "alleles.jsonl"
    write_jsonl(allele_path, evaluation.alleles)
    summary = attempt / "summary.json"
    write_new(summary, canonical(evaluation.summary) + b"\n")
    evidence = evidence + [
        artifact(dispositions, "dispositions", "jsonl"),
        artifact(allele_path, "allele_targets", "jsonl"),
        artifact(summary, "summary", "json"),
    ]
    case_sha = None
    target_count = unique_count = 0
    if evaluation.request:
        response = preflight(evaluation.request, evaluation.alleles)
        payload = canonical(evaluation.request.model_dump(mode="json")) + b"\n"
        path = attempt / "case.json"
        write_new(path, payload)
        case_sha = sha256(payload).hexdigest()
        evidence.append(artifact(path, "case", "json"))
        target_count, unique_count = response.target_count, response.unique_sequence_count
    manifest = Manifest(
        **identities,
        status="ready" if target_count else "no_targets",
        case_sha256=case_sha,
        artifacts=evidence,
        checks=checks,
        transcript_counts=evaluation.summary["transcript_counts"],
        target_count=target_count,
        unique_sequence_count=unique_count,
    )
    if target_count:
        publish_ready(attempt, ready, manifest.model_dump(mode="json"))
    else:
        write_new(attempt / "manifest.json", canonical(manifest.model_dump(mode="json")) + b"\n")
    return manifest
