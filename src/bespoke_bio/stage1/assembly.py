"""Join registered native evidence and publish every supported Stage 1 allele."""

from pathlib import Path

import yaml

from bespoke_bio.api import Reference, TranslateRequest

from .biology import IntervalMask
from .contracts import Artifact, Check, Intake, Policy, Record
from .evidence import callability, depth_mask, intersect, read_bed, read_expression
from .handoff import empty_chromosomes, evaluate, preflight, publish
from .io import (
    artifact,
    canonical,
    digest,
    read_json,
    run_identity,
    run_lock,
    source_hash,
    write_new,
)
from .references import read_reference_records
from .variants import chromosome_sites


class AssemblyInputs(Record):
    """Registered immutable files, identified by role rather than guessed filename patterns."""

    intake: Path
    policy: Path
    tools: Path
    reference_records: Path
    cds_bed: Path
    depth: Path
    original_vcf: Path
    gvcf: Path
    normalized_vcf: Path
    annotated_vcf: Path
    phased_vcf: Path
    expression: Path
    sample: str
    reference: Reference
    reference_sha256: str
    evidence_registry: Path
    release_checks: Path
    durable_roots: list[Path]


def assemble(inputs: AssemblyInputs, repository: Path, output: Path, attempt_id: str):
    """Validate joins, retain dispositions and release only after all named checks pass."""
    policy = Policy.model_validate(yaml.safe_load(inputs.policy.read_text()))
    intake = Intake.model_validate(read_json(inputs.intake))
    if {r.specimen for r in intake.reads if r.assay == "germline_wgs"} != {inputs.sample}:
        raise ValueError("variant sample does not match the DNA intake specimen")
    if Path(attempt_id).name != attempt_id or attempt_id in {".", ".."}:
        raise ValueError("unsafe assembly attempt identifier")
    tools = read_json(inputs.tools)
    evidence = [Artifact.model_validate(item) for item in read_json(inputs.evidence_registry)]
    by_path = {}
    for item in evidence:
        path = Path(item.location).resolve()
        if not path.is_file() or path.stat().st_size != item.bytes:
            raise ValueError(f"native evidence absent or truncated: {item.role}")
        # Publication verifies each content hash once, after the native joins succeed.
        if path in by_path:
            raise ValueError("native evidence registered more than once")
        by_path[path] = item
    required = (
        inputs.reference_records,
        inputs.cds_bed,
        inputs.depth,
        inputs.original_vcf,
        inputs.gvcf,
        inputs.normalized_vcf,
        inputs.annotated_vcf,
        inputs.phased_vcf,
        inputs.expression,
        inputs.intake,
        inputs.policy,
        inputs.tools,
    )
    if any(path.resolve() not in by_path for path in required):
        raise ValueError("assembly input missing from immutable evidence registry")
    checks = [Check.model_validate(item) for item in read_json(inputs.release_checks)]
    adapter = source_hash(repository / "src/bespoke_bio/stage1")
    step2 = source_hash(repository / "src/bespoke_bio")
    contract = digest(TranslateRequest.model_json_schema())
    config = digest(policy.model_dump(mode="json"))
    run_id = run_identity(
        intake,
        reference=inputs.reference_sha256,
        config=config,
        tools=digest(tools),
        adapter=adapter,
        step2_source=step2,
        contract=contract,
    )
    output.mkdir(parents=True, exist_ok=True)
    attempt = output / attempt_id
    with run_lock(output / "assembly.lock"):
        if (output / "ready.json").exists():
            raise FileExistsError("handoff already published; validate or create a new run")
        working = output / (attempt_id + "-evidence")
        working.mkdir(exist_ok=False)
        transcripts = read_reference_records(inputs.reference_records)
        expression = read_expression(inputs.expression)
        bed = read_bed(inputs.cds_bed)
        depth = depth_mask(inputs.depth, bed, policy.min_depth)
        masks = callability(inputs.original_vcf, inputs.gvcf, inputs.sample, policy, depth)
        masks = {
            chrom: IntervalMask(intersect(mask.intervals, bed.get(chrom, [])))
            for chrom, mask in masks.items()
        }
        mask_path = working / "callable-cds.bed"
        write_new(
            mask_path,
            "".join(
                f"{chrom}\t{start}\t{end}\n"
                for chrom in sorted(masks)
                for start, end in masks[chrom].intervals
            ).encode(),
        )
        links = working / "variant-transcript-links.jsonl"
        chromosomes = chromosome_sites(
            inputs.original_vcf,
            inputs.normalized_vcf,
            inputs.annotated_vcf,
            inputs.phased_vcf,
            inputs.sample,
            policy,
            set(transcripts),
            by_path[inputs.original_vcf.resolve()].sha256,
            links,
        )
        evaluation = evaluate(
            intake.case_id,
            inputs.reference,
            transcripts,
            expression,
            masks,
            empty_chromosomes(
                chromosomes, {tx.chrom for tx in transcripts.values() if not tx.scope_reason}
            ),
            policy,
        )
        response = preflight(evaluation.request, evaluation.alleles) if evaluation.request else None
        report = working / "assembly-checks.json"
        write_new(
            report,
            canonical(
                {
                    "transcripts": len(transcripts),
                    "expression_rows": len(expression),
                    "dispositions": len(evaluation.dispositions),
                    "targets": response.target_count if response else 0,
                    "unique_proteins": response.unique_sequence_count if response else 0,
                    "summary": evaluation.summary,
                }
            )
            + b"\n",
        )
        evidence.extend(
            [
                artifact(mask_path, "callable_cds", "bed"),
                artifact(links, "variant_transcript_links", "jsonl"),
                artifact(report, "assembly_checks", "json"),
            ]
        )
        for name in (
            "callability",
            "variant_reconciliation",
            "handoff_translation",
            "disposition_accounting",
            "artifact_integrity",
        ):
            if any(check.name == name for check in checks):
                raise ValueError("external release checks cannot replace an assembly-owned check")
            checks.append(
                Check(
                    name=name,
                    status="passed",
                    command=["bespoke-stage1", "assemble"],
                    evidence_roles=["assembly_checks", "callable_cds", "variant_transcript_links"],
                )
            )
        identities = {
            "run_id": run_id,
            "attempt_id": attempt_id,
            "case_id": intake.case_id,
            "adapter_sha256": adapter,
            "contract_sha256": contract,
            "reference_sha256": inputs.reference_sha256,
            "config_sha256": config,
            "intake_sha256": digest(intake.model_dump(mode="json")),
        }
        return publish(
            evaluation,
            attempt,
            output / "ready.json",
            identities,
            evidence,
            checks,
            tuple(inputs.durable_roots) + (output,),
        )
