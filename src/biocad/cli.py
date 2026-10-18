import argparse
import json
from pathlib import Path

from .contracts import DesignRequest, DesignResponse
from .pipeline import run_example, save


def main():
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="command", required=True)
    demo = sub.add_parser("demo")
    demo.add_argument("--fixtures", default="examples")
    demo.add_argument("--output", required=True)
    demo.add_argument("--request")
    demo.add_argument("--foldseek", default="vendor/foldseek/bin/foldseek")
    sub.add_parser("schema").add_argument("--output", default="examples/response.schema.json")
    compile_ = sub.add_parser("compile")
    compile_.add_argument("request")
    compile_.add_argument("response")
    compile_.add_argument("--output", required=True)
    live = sub.add_parser("design")
    live.add_argument("plan")
    live.add_argument("--output", required=True)
    placement = sub.add_parser("placements")
    placement.add_argument("--bundle", default="examples/replay")
    placement.add_argument("--request")
    placement.add_argument("--output", required=True)
    motifs = sub.add_parser("motifs")
    motifs.add_argument("--assembly", default="examples/compact_assembly")
    motifs.add_argument("--request")
    motifs.add_argument("--output", required=True)
    importer = sub.add_parser("import-targets")
    importer.add_argument("bundle")
    importer.add_argument("--output", required=True)
    importer.add_argument("--chains", help="JSON mapping upstream target IDs to author-chain lists")
    baseline = sub.add_parser("baseline")
    baseline.add_argument("--fixtures", default="examples")
    baseline.add_argument("--output", required=True)
    boltzgen = sub.add_parser("boltzgen-prepare")
    boltzgen.add_argument("--fixtures", default="examples")
    boltzgen.add_argument("--output", required=True)
    boltzgen.add_argument("--source", help="Single FASTA, PDB or mmCIF target")
    boltzgen.add_argument("--chain")
    boltzgen.add_argument("--name", default="design")
    boltzgen.add_argument("--binder-length", type=int, default=75)
    boltzgen.add_argument("--num-designs", type=int, default=2)
    export = sub.add_parser("export-binders")
    export.add_argument("results")
    export.add_argument("target_state")
    export.add_argument("--output", required=True)
    catalog = sub.add_parser("catalog-search")
    catalog.add_argument("query")
    catalog.add_argument("--chain", required=True)
    catalog.add_argument("--entries", required=True, help="JSON list of declared local templates")
    catalog.add_argument("--output", required=True)
    catalog.add_argument("--foldseek", default="vendor/foldseek/bin/foldseek")
    submit = sub.add_parser("submit-bundle", help="Send an upstream package to the design API")
    submit.add_argument("bundle")
    submit.add_argument("--url", default="http://127.0.0.1:8000")
    submit.add_argument("--output", required=True, help="Local receipt directory")
    submit.add_argument("--chains", help="JSON mapping target IDs to author-chain lists")
    submit.add_argument("--binder-length", type=int, default=75)
    submit.add_argument("--num-designs", type=int, default=2)
    submit.add_argument(
        "--execute",
        action="store_true",
        help="Explicitly start each prepared child on the configured GPU worker",
    )
    submit.add_argument(
        "--wait", action="store_true", help="Wait for submitted jobs; requires --execute"
    )
    submit.add_argument("--timeout", type=int, default=3600)
    a = p.parse_args()
    if a.command == "submit-bundle":
        from .client import DesignClient

        if a.wait and not a.execute:
            p.error("--wait requires --execute")
        destination = Path(a.output)
        if destination.exists() and any(destination.iterdir()):
            p.error("Use an empty receipt directory")
        destination.mkdir(parents=True, exist_ok=True)
        client = DesignClient(a.url)
        receipt = client.submit_bundle(
            a.bundle,
            binder_length=a.binder_length,
            num_designs=a.num_designs,
            chain_selections=json.loads(Path(a.chains).read_text()) if a.chains else None,
        )
        save(destination / "intake.json", receipt)
        outcomes = []
        for record in receipt.get("records", []):
            child = record.get("child_run_id")
            if not child:
                outcomes.append(record)
                continue
            try:
                outcome = client.start(child) if a.execute else client.status(child)
                if a.wait:
                    outcome = client.wait(child, timeout=a.timeout)
                    if outcome["status"] == "completed":
                        client.download(child, destination / (child + ".zip"))
                        save(destination / (child + "_results.json"), client.results(child))
                outcomes.append(outcome)
            except (RuntimeError, OSError, TimeoutError) as exc:
                outcomes.append({"run_id": child, "client_error": str(exc)})
            save(destination / "jobs.json", outcomes)
        print(
            json.dumps(
                {"run_id": receipt["run_id"], "jobs": outcomes, "receipts": str(destination)},
                indent=2,
            )
        )
        if not receipt.get("prepared_targets") or any(
            "client_error" in o or o.get("status") in {"failed", "failed_preparation"}
            for o in outcomes
        ):
            raise SystemExit(1)
        return
    if a.command == "boltzgen-prepare":
        from .boltzgen_adapter import prepare, prepare_target

        result = (
            prepare_target(
                a.source,
                a.output,
                chain=a.chain,
                name=a.name,
                binder_length=a.binder_length,
                num_designs=a.num_designs,
            )
            if a.source
            else prepare(a.fixtures, a.output)
        )
        print(json.dumps(result, indent=2))
    elif a.command == "export-binders":
        from .design_bundle import export_bundle

        print(json.dumps(export_bundle(a.results, a.target_state, a.output), indent=2))
    elif a.command == "catalog-search":
        from .catalog import search_catalog

        print(
            json.dumps(
                search_catalog(
                    a.query,
                    a.chain,
                    json.loads(Path(a.entries).read_text()),
                    a.output,
                    executable=a.foldseek,
                ),
                indent=2,
            )
        )
    elif a.command == "baseline":
        from .baseline import prepare

        print(json.dumps(prepare(a.fixtures, a.output), indent=2))
    elif a.command == "demo":
        req = DesignRequest.model_validate_json(Path(a.request).read_text()) if a.request else None
        r = run_example(a.fixtures, a.output, req, a.foldseek)
        print(
            json.dumps(
                {"response": str(Path(a.output) / "response.json"), "stages": r.stages}, indent=2
            )
        )
    elif a.command == "schema":
        save(a.output, DesignResponse.model_json_schema())
    elif a.command == "compile":
        from .mosaic_adapter import compile_request

        r = DesignRequest.model_validate_json(Path(a.request).read_text())
        response = DesignResponse.model_validate_json(Path(a.response).read_text())
        save(a.output, compile_request(r, response.target_states[0]))
    elif a.command == "placements":
        from .contracts import PlacementSearchRequest
        from .placement import run_search

        req = (
            PlacementSearchRequest.model_validate_json(Path(a.request).read_text())
            if a.request
            else None
        )
        result = run_search(a.bundle, a.output, req)
        print(
            json.dumps(
                {
                    k: result[k]
                    for k in [
                        "sampled_layouts",
                        "eligible_layouts",
                        "top_layout_ids",
                        "elapsed_seconds",
                    ]
                },
                indent=2,
            )
        )
    elif a.command == "motifs":
        from .contracts import MotifPlanRequest
        from .motifs import prepare_plan

        request = (
            MotifPlanRequest.model_validate_json(Path(a.request).read_text()) if a.request else None
        )
        result = prepare_plan(a.assembly, a.output, request)
        print(
            json.dumps(
                {
                    k: result[k]
                    for k in [
                        "fixed_motif_residues",
                        "locally_viable_topologies",
                        "selected_topology_ids",
                        "model_validation",
                    ]
                },
                indent=2,
            )
        )
    elif a.command == "import-targets":
        from .upstream import import_targets

        selections = json.loads(Path(a.chains).read_text()) if a.chains else None
        result = import_targets(a.bundle, a.output, selections)
        print(
            json.dumps({"states": len(result.target_states), "records": result.records}, indent=2)
        )
    elif a.command == "design":
        from .mosaic_adapter import run

        run(json.loads(Path(a.plan).read_text()), a.output)


if __name__ == "__main__":
    main()
