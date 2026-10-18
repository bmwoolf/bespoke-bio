"""Explicit Stage 1 operations; all paths and machine settings are supplied by the operator."""

import argparse
from pathlib import Path

from .assembly import AssemblyInputs, assemble
from .contracts import Intake
from .intake import fetch_locked, pilot_intake
from .io import canonical, read_json, write_new
from .precomputed import PrecomputedSource, register_precomputed
from .precomputed_run import PrecomputedRun, process_precomputed
from .proteins import export_proteins, validate_protein_package
from .references import (
    compare_gffread,
    prepare_references,
    read_reference_records,
    write_reference_records,
)
from .runtime import Runtime
from .workflows import ReferenceAssets, inspect, launch, prepare


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    fetch = commands.add_parser("fetch", help="Acquire and register pinned public sources")
    fetch.add_argument("--lock", required=True, type=Path)
    fetch.add_argument("--root", required=True, type=Path)
    fetch.add_argument("--workers", type=int, default=4)
    refs = commands.add_parser("references", help="Validate every in-scope reference CDS")
    refs.add_argument("--gtf", required=True, type=Path)
    refs.add_argument("--fasta", required=True, type=Path)
    refs.add_argument("--proteins", required=True, type=Path)
    refs.add_argument("--output", required=True, type=Path)
    workflows = commands.add_parser("prepare-workflows", help="Write full-input workflow settings")
    workflows.add_argument("--intake", required=True, type=Path)
    workflows.add_argument("--references", required=True, type=Path)
    workflows.add_argument("--runtime", required=True, type=Path)
    workflows.add_argument("--output", required=True, type=Path)
    workflows.add_argument("--run-id", required=True)
    inspection = commands.add_parser(
        "inspect-workflow", help="Resolve reachable native tool images"
    )
    inspection.add_argument("--runtime", required=True, type=Path)
    inspection.add_argument("--repository", required=True, type=Path)
    inspection.add_argument("--root", required=True, type=Path)
    inspection.add_argument("--branch", choices=["dna", "rna"], required=True)
    run = commands.add_parser("run-workflows", help="Run complete DNA and RNA workflows")
    run.add_argument("--runtime", required=True, type=Path)
    run.add_argument("--repository", required=True, type=Path)
    run.add_argument("--root", required=True, type=Path)
    run.add_argument("--run-id", required=True)
    run.add_argument("--attempt-id", required=True)
    run.add_argument("--resume", action="store_true")
    native_refs = commands.add_parser(
        "validate-reference-native", help="Compare all complete CDS with gffread"
    )
    native_refs.add_argument("--records", required=True, type=Path)
    native_refs.add_argument("--gffread-cds", required=True, type=Path)
    native_refs.add_argument("--output", required=True, type=Path)
    handoff = commands.add_parser(
        "assemble", help="Publish all validated allele records and evidence"
    )
    handoff.add_argument("--inputs", required=True, type=Path)
    handoff.add_argument("--repository", required=True, type=Path)
    handoff.add_argument("--output", required=True, type=Path)
    handoff.add_argument("--attempt-id", required=True)
    imported = commands.add_parser(
        "import-precomputed", help="Register precomputed benchmark evidence for a downstream demo"
    )
    imported.add_argument("--source", required=True, type=Path)
    imported.add_argument("--gtf", required=True, type=Path)
    imported.add_argument("--fasta", required=True, type=Path)
    imported.add_argument("--output", required=True, type=Path)
    imported.add_argument("--flank", type=int, default=50)
    proteins = commands.add_parser(
        "export-proteins", help="Export complete precomputed-evidence protein sequences"
    )
    proteins.add_argument("--source", required=True, type=Path)
    proteins.add_argument("--records", required=True, type=Path)
    proteins.add_argument("--normalized", required=True, type=Path)
    proteins.add_argument("--annotated", required=True, type=Path)
    proteins.add_argument("--confidence", required=True, type=Path)
    proteins.add_argument("--output", required=True, type=Path)
    prepared = commands.add_parser(
        "process-precomputed", help="Run pinned native annotation and protein export"
    )
    prepared.add_argument("--settings", required=True, type=Path)
    prepared.add_argument("--output", required=True, type=Path)
    verify_proteins = commands.add_parser(
        "validate-proteins", help="Verify a completed or restored protein package"
    )
    verify_proteins.add_argument("--root", required=True, type=Path)
    arguments = parser.parse_args()
    if arguments.command == "fetch":
        lock = read_json(arguments.lock)
        artifacts = fetch_locked(lock, arguments.root, arguments.workers)
        intake = pilot_intake(lock, artifacts)
        write_new(arguments.root / "intake.json", canonical(intake.model_dump(mode="json")) + b"\n")
        write_new(
            arguments.root / "sources.json",
            canonical([a.model_dump(mode="json") for a in artifacts]) + b"\n",
        )
        print(f"Registered {len(artifacts)} source files.")
    elif arguments.command == "references":
        transcripts = prepare_references(arguments.gtf, arguments.fasta, arguments.proteins)
        write_reference_records(transcripts, arguments.output)
        print(f"Resolved {len(transcripts)} annotated transcripts.")
    elif arguments.command == "prepare-workflows":
        prepare(
            Intake.model_validate(read_json(arguments.intake)),
            ReferenceAssets.model_validate(read_json(arguments.references)),
            Runtime.model_validate(read_json(arguments.runtime)),
            arguments.output,
            arguments.run_id,
        )
    elif arguments.command == "inspect-workflow":
        inspect(
            Runtime.model_validate(read_json(arguments.runtime)),
            arguments.repository,
            arguments.root,
            arguments.branch,
        )
    elif arguments.command == "run-workflows":
        launch(
            Runtime.model_validate(read_json(arguments.runtime)),
            arguments.repository,
            arguments.root,
            arguments.run_id,
            arguments.attempt_id,
            arguments.resume,
        )
    elif arguments.command == "validate-reference-native":
        result = compare_gffread(read_reference_records(arguments.records), arguments.gffread_cds)
        write_new(arguments.output, canonical(result) + b"\n")
        print(result)
    elif arguments.command == "assemble":
        manifest = assemble(
            AssemblyInputs.model_validate(read_json(arguments.inputs)),
            arguments.repository,
            arguments.output,
            arguments.attempt_id,
        )
        print(
            f"{manifest.status}: {manifest.target_count} biological targets; {manifest.unique_sequence_count} unique proteins"
        )
    elif arguments.command == "import-precomputed":
        result = register_precomputed(
            PrecomputedSource.model_validate(read_json(arguments.source)),
            arguments.gtf,
            arguments.fasta,
            arguments.output,
            arguments.flank,
        )
        print(result["counts"])
    elif arguments.command == "export-proteins":
        result = export_proteins(
            PrecomputedSource.model_validate(read_json(arguments.source)),
            arguments.records,
            arguments.normalized,
            arguments.annotated,
            arguments.confidence,
            arguments.output,
        )
        print(result["counts"])
    elif arguments.command == "process-precomputed":
        result = process_precomputed(
            PrecomputedRun.model_validate(read_json(arguments.settings)),
            arguments.output,
        )
        print(result["counts"])
    elif arguments.command == "validate-proteins":
        print(validate_protein_package(arguments.root))


if __name__ == "__main__":
    main()
