# Optional local template retrieval

The hackathon's primary route is a supplied target structure → BoltzGen design → sequence design → independent structure predictions → candidate filtering and downloadable evidence. A sequence-only input can use the native BoltzGen sequence-conditioned route; preparation is tested locally and a fresh GPU run is needed to establish its behavior here. Structural Foldseek retrieval requires coordinates. The website/catalog configurator can be deferred: the original implementation prompt explicitly assigns the website and infrastructure to other engineers and asks for portable scientific artifacts.

`biocad.catalog.search_catalog` adds a small optional branch before design. It searches **1–100 explicitly declared local mmCIF target templates** with the existing native Foldseek wrapper. It returns source hashes, a stable template identifier, raw structural match scores, residue correspondence and a rigidly transferred known-binder pose. It never treats a structural hit as evidence of binding to the query. It does not download a database or substitute a known experimental binder for a newly generated design.

```python
from biocad.catalog import CatalogEntry, search_catalog

report = search_catalog(
    "examples/5J89.cif", "A",
    [CatalogEntry(
        target_structure="examples/5JDS.cif", target_chain="A",
        binder_structure="examples/5JDS.cif", binder_chain="B",
        reference="https://www.rcsb.org/structure/5JDS",
        evidence={"type": "experimental_complex", "scope": "observed coordinates"},
        conditions={"fixture": "5JDS"},
    )],
    "outputs/template_lookup",
    conditions={"fixture": "5J89"},
    executable="vendor/foldseek/bin/foldseek",
)
```

Use an empty output directory. `catalog_results.json` contains `kind="retrieved_known_binder_template"` and `generated_candidate=false`; the emitted mmCIFs contain target chain T and transferred binder chain B. Coordinates are in angstroms. `pose_mapping` records the transform and aligned author residues.

An exact **observed** target-sequence match and exact declared conditions are reported separately. Neither implies full-construct identity: unresolved residues, modifications, assay conditions and sequence truncation still matter. Automatic candidate reuse is deliberately disabled. The 5JDS binder is an observed-only KN035 construct, not a validated newly designed candidate. This keeps the demonstration honest even when a real experimental template is available.

## Smallest scope that satisfies the demo now

- Keep: exact target preparation, upstream package import, reproducible BoltzGen generation and score ingestion, candidate/target evidence, structure/sequence exports and run provenance.
- Optional now: this bounded known-template branch; it reconnects the reusable key/lock vocabulary to actual coordinate retrieval without requiring a new search platform.
- Defer: production catalog UI, public structural database indexing, motif interleaving, compact co-key redesign, flexible assembly, functional switching and multistate selectivity optimization. Existing rigid assembly geometry remains a separate hypothesis rather than a validated connected protein.
- Scientific gaps remain explicit: model confidence does not measure affinity; strong binder claims need independent assessment and experimental validation. Two completed designs demonstrate executable plumbing, not a successful discovery campaign.

Local verification: `pytest tests/test_catalog.py` executes native Foldseek on the public 5JDS/5J89 fixtures when the binary is present; it checks real mapped poses, provenance, exact-match versus condition mismatch and the no-automatic-reuse policy. If the binary is absent, that native integration test is clearly skipped.
