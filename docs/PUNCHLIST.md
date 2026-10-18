# Hackathon delivery punchlist

The deliverable is an inspectable scientific pipeline with a teammate-facing HTTP interface. Production catalog/site work is outside this milestone.

## Critical path

- [x] Execute the complete off-the-shelf BoltzGen protocol on a public PD-L1 target without specifying an epitope.
- [x] Download every candidate, native filter outcome, confidence metric, log and model provenance; delete the original GPU pod.
- [x] Prepare arbitrary single-chain PDB/mmCIF and single-sequence FASTA inputs with residue mappings and input hashes.
- [x] Export actual predictions as typed target-linked binder components, attachment frames, FASTA, structures and evaluations.
- [x] Implement bounded Foldseek retrieval of declared target/binder templates, preserving experimental context and distinguishing analogy from validated reuse.
- [x] Finish and test the HTTP intake for direct targets and upstream target packages; expose status and artifact downloads.
- [x] Publish runnable teammate request examples and replay the real cached results through the API.
- [x] Create a 100 GB persistent Runpod model/runtime cache volume in US-MO-1.
- [x] Populate the cache on a GPU, verify runtime/weight readiness, and prove restoration on another disposable pod.
- [x] Connect prepared HTTP jobs to a static worker or authenticated managed service with persistent queue/budget, automatic GPU startup/deletion, cancellation and status.
- [x] Verify that HTTP live dispatch completes a real fresh GPU job and returns its bundle. Pod lifecycle remains operator-owned with a cost guard.
- [x] Refresh the portable handoff archive with the API and validated examples.

## After the integration works

- [ ] Run one novel incoming teammate target end to end; the public PD-L1 sequence-only path has now executed; a novel teammate target remains to be supplied.
- [ ] Increase candidate count and retain diversity; the two-candidate smoke test establishes execution, not reliable hit discovery.
- [ ] Add independent complex validation or a BindCraft comparison if time permits. Keep BoltzGen as the working primary engine.
- [ ] Add explicit target-state/off-target panels and selection rules.

## Deferred

Interdigitated motif redesign, a folded connected multi-module assembly, full negative-design optimization, large external structure databases, production website/database/scheduling, and experimental binding validation.

## Current scientific result

Two 75-residue candidates were generated against 115 observed PD-L1 residues. One passes all nine native default filters; the other fails complex backbone recovery (2.61 Å versus a 2.5 Å threshold). Both have weak interface confidence (ipTM about 0.30). Neither is an experimentally validated binder. Filter pass and evidence for binding remain separate in the exported bundle.

Storage remains billable after GPU deletion: the retained 100 GB network volume is approximately $7/month until deleted. GPU setup and testing remain subject to the authorized $100 cap.

## Latest integration proof

A fresh sequence-only HTTP job completed all seven stages and returned one native-filter-passing candidate. Complex recovery RMSD is 0.83 Å, binder-alone recovery 1.33 Å, ipTM 0.41 and ESMFold2 minimum ipSAE 0.172. The predicted target recovers the 115 observed experimental PD-L1 Cα positions at 1.36 Å RMSD after proper rigid alignment; those coordinates were not supplied to generation. This is still computational evidence. See `examples/boltzgen_api_results/`.

A second disposable pod restored the persistent cache in 76 seconds and passed both CUDA runtime checks. Both GPU pods were deleted and deletion verified. This phase's estimated GPU cost was about $1.14 plus storage. The code and proposed feed contract are pushed on `codex/binder-pipeline-demo` and open for team review in [PR #1](https://github.com/bmwoolf/bespoke-bio/pull/1). The portable ZIP is also available.
