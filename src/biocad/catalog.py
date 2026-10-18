"""Optional bounded retrieval of declared target/binder templates, not binder generation."""

import json
from pathlib import Path

import gemmi
from pydantic import Field

from .contracts import Record
from .retrieval import search, transfer
from .structures import sequence, sha, write_chains


class CatalogEntry(Record):
    target_structure: str
    target_chain: str
    binder_structure: str
    binder_chain: str
    reference: str
    evidence: dict
    conditions: dict = Field(default_factory=dict)


def search_catalog(
    query, query_chain, entries, output, *, conditions=None, executable="foldseek", top_k=5
):
    """Return known templates and mapped poses; never certify affinity or a cache hit.

    Entries must describe explicit local structures. Exact observed sequence and condition
    matches are reported separately from structural analogy; even those do not establish
    identical full constructs, modifications, or experimental validation of truncations.
    """
    if not 1 <= len(entries) <= 100:
        raise ValueError("Declare 1..100 catalog entries")
    output = Path(output)
    if output.exists() and any(output.iterdir()):
        raise ValueError("Use an empty output directory")
    library = output / "library"
    library.mkdir(parents=True)
    query_structure = gemmi.read_structure(str(query))
    query_structure.remove_alternative_conformations()
    query_c = query_structure[0][query_chain]
    write_chains([query_c], output / "query.cif")
    declared = {}
    for index, raw in enumerate(entries):
        entry = CatalogEntry.model_validate(raw)
        target_path, binder_path = Path(entry.target_structure), Path(entry.binder_structure)
        target_st = gemmi.read_structure(str(target_path))
        binder_st = gemmi.read_structure(str(binder_path))
        target_st.remove_alternative_conformations()
        binder_st.remove_alternative_conformations()
        target, binder = target_st[0][entry.target_chain], binder_st[0][entry.binder_chain]
        provenance = {
            "entry": entry.model_dump(),
            "target_source_sha256": sha(target_path.read_bytes()),
            "binder_source_sha256": sha(binder_path.read_bytes()),
        }
        identity = {
            "target_sha256": provenance["target_source_sha256"],
            "binder_sha256": provenance["binder_source_sha256"],
            "target_chain": entry.target_chain,
            "binder_chain": entry.binder_chain,
            "reference": entry.reference,
            "conditions": entry.conditions,
        }
        entry_id = "template_" + sha(json.dumps(identity, sort_keys=True).encode())[:20]
        name = f"entry{index:03d}"
        write_chains([target], library / f"{name}.cif")
        declared[f"{name}_{target.name}"] = (entry_id, entry, target, binder, provenance)
    raw_search = search(output / "query.cif", library, output / "foldseek", executable, top_k)
    results = []
    for rank, hit in enumerate(raw_search["hits"]):
        entry_id, entry, target, binder, provenance = declared[hit["target"]]
        mapped, mapping = transfer(hit, query_c, target, binder)
        # Use distinct output chains even when the source partner shares query chain name.
        target_copy = query_c.clone()
        target_copy.name = "T"
        mapped.name = "B"
        pose_path = f"retrieved_{rank:02d}.cif"
        write_chains([target_copy, mapped], output / pose_path)
        results.append(
            {
                "template_id": entry_id,
                "kind": "retrieved_known_binder_template",
                "generated_candidate": False,
                "reference": entry.reference,
                "provenance": provenance,
                "binder_sequence": sequence(binder),
                "exact_observed_target_sequence_match": sequence(query_c) == sequence(target),
                "conditions_match": conditions is not None and conditions == entry.conditions,
                "automatic_candidate_reuse": False,
                "structural_match": hit,
                "pose_mapping": mapping,
                "pose": pose_path,
                "binding_validation": "not established for query by structural similarity",
                "construct_scope": "observed coordinates only; full construct identity not inferred",
            }
        )
    report = {
        "schema_version": "0.1.0",
        "route": "optional_local_template_retrieval",
        "query_sha256": sha(Path(query).read_bytes()),
        "query_chain": query_chain,
        "conditions": conditions,
        "catalog_size": len(entries),
        "foldseek": raw_search,
        "results": results,
        "next_step": "generate and independently evaluate binders when no suitable validated component exists",
    }
    (output / "catalog_results.json").write_text(json.dumps(report, indent=2) + "\n")
    return report
