> Historical development notes. Superseded by [current status](STATUS.md), [integration](INTEGRATION.md), and [compute](COMPUTE.md). Statements about the active engine, pending work and resource state describe earlier stages.

# Execution status

Current milestone: [hackathon punchlist](PUNCHLIST.md). The active engine is BoltzGen;
older Mosaic and interdigitated assembly plans below describe research history.

**Confirmed now:** full seven-stage no-epitope PD-L1 run completed, two candidates downloaded,
one native-filter pass, weak interface confidence for both, original GPU pod deleted.
Generic sequence/structure preparation, typed binder export, bounded native Foldseek
catalog retrieval, upstream HTTP intake, explicit existing-worker dispatch and ZIP
downloads pass local tests. The API is running locally and real cached replay/download
has been verified over HTTP. Fresh sequence-only HTTP dispatch also completed all seven native stages, exported one
filter-passing candidate and downloaded its ZIP. A fresh GPU restored the cache in
76 seconds. Both disposable GPUs were deleted; the persistent volume remains.
See `examples/boltzgen_api_results/summary.json` for the new result.


| Component | Status | Evidence or blocker |
|---|---|---|
| Original prompt, example JSON and handwritten diagram | Read | Scientific ownership remains between upstream structures and downstream package |
| Typed contracts and package CLI | Executed | Installed local package; strict JSON Schema emitted |
| Structure preparation and residue maps | Executed | Two public real structures; author/label IDs separated; missing standard heavy atoms recorded |
| Candidate/state matrix | Executed for geometry | Same observed-only component across two conditions; confidence explicitly missing |
| Site annotation and mapping | Executed | Observed contact site and alignment-mapped alternate site |
| Bounded Foldseek search | Executed | One-template declared library, full native alignment retained |
| Partner pose transfer | Executed | Cα Kabsch refit; alternate-state transfer rejects severe overlaps |
| Reusable key and terminal ports | Executed | Observed-only computational construct; native unresolved terminal residues omitted |
| Two-module tandem placement | Executed | Two identical modules, artificial relative target placement |
| Required interfaces and clashes | Executed | Each intended edge and unintended component/target cross-pair assessed separately |
| Connection sequence and reach bound | Executed | Fusion FASTA explicit; linker remains structurally unmodeled |
| Offline scientific review | Exported | HTML with XY/XZ/YZ Cα projections plus per-interface evaluations; mmCIF for all atoms |
| Real cached de novo output | Executed | Two BoltzGen PD-L1 candidates; full logs, structures and explicit filter flags retained |
| Mosaic live adapter | Source-checked; preflight executed | No installed GPU runtime/weights; no model inference run |
| Multistate negative-design optimization | Unsupported in live runner | Matrix contract exists; per-state model execution/selection not wired |
| Whole connected-construct/joint prediction | Pending | No folded linker or full assembly confidence |
| Key-side retrieval / assembly-template retrieval | Pending | Generic bounded search helper reusable, but no route-C example demonstrated |
| Interface-search/Folddisco | Deferred | Installed Foldseek interface command inspected; motif adapter not implemented |
| Full biological-assembly/context expansion | Pending | Source retained; selected author chains prepared; no inferred membrane/glycan context |
| Measurements/clinical effect | Not computed | No affinity, degradation, selectivity or therapeutic efficacy claim |

Verification covers known-transform pose transfer, reflections/collinearity, cropped/insertion-coded residue mappings, a full-context clash control, actual fixture identity and sequence export, candidate/state linkage, missing prediction semantics, separate required interfaces, reuse of identical component IDs, overlapping assembly rejection, and real configuration-dependent clearance/reach failures. Six focused tests passed. GPU behavior is not covered by those tests.

Next two priorities:

1. Identify authorized GPU compute and budget; execute the pinned Mosaic primary-state smoke test, add exact weight/device provenance, then re-predict the same discrete sequence across the supplied state/exclusion panel and ingest those artifacts.
2. Supply the intended distinct second module and architecture; build a small template library, search attachment-compatible layouts, and model the actual linker/whole construct. Add ensemble/state-selective objectives only after this loop is reproducible.

## Assembly-focused update, 2026-10-17

The user's priority is now compact assembly first; conformational modeling is deferred.
A CPU-only `biocad placements` command executed a 48-layout rigid search plus an overlap
control, with element-aware vdW overlap checks and a linker reach rejection. Five layouts
passed these necessary checks. It exports three compact proposals, exact transforms, an
offline report and an unexecuted motif-redesign brief. Shared-scaffold interdigitation,
sequence shortening and folded linker modeling have not run. The combined repository
suite now has 18 passing tests, including the upstream translation API tests.

Immediate next priority is scaffold/connection definition and generative redesign around
protected interfaces; model-based validation still needs an identified runtime. Details
and the Raygun/SCISOR assessment are in [ASSEMBLY.md](ASSEMBLY.md).

## Motif-redesign preparation update

`biocad motifs` now executes a local, weight-free motif/topology preflight. It preserves
82 residues in four indexed blocks, including endpoints of four disulfides inferred from
observed coordinates. Six orders pass the current coarse budgets; one interleaved order
and one tandem control are exported at a matched planned length of 127 residues. Native
RFD3 parser validation and model inference have not run. CPU controls test exact atom
mapping, support retention, impossible bridge budgets, and detection of moved/missing
motifs relative to a joint target frame. See [MOTIF_REDESIGN.md](MOTIF_REDESIGN.md).

## Local readiness update (2026-10-17)

This update supersedes the earlier pending-native-parser and compute-discovery statements above. Native Foundry/RFD3 input construction has now passed on CPU for both the original 127-residue pair and the expanded matched 170-residue smoke pair. Installed CLI Hydra configuration composition also passed without model loading. The current bundle is `examples/gpu_smoke`.

The new core audit identifies omitted hydrophobic support under retained motifs, so this remains scaffold replacement, not validated framework deletion. Upstream `targets.jsonl` packages now import through `biocad import-targets` with exact sequence/hash matching, author-chain selection, provenance, transcript deduplication, missing-structure status and explicit RNA/protein-abundance separation. Both mmCIF and PDB paths are tested.

Runpod authentication is confirmed; the $100 phase budget is approved. Two single-output GPU jobs, their isolated runtime bootstrap and result checks are prepared. No GPU instance, model inference, redesigned sequence, independent refold or biological validation has run. Conformation remains deferred.

## GPU smoke results (2026-10-17)

Both prepared RFD3 cases now ran on a single Runpod A100 SXM 80GB with Torch 2.11.0+cu128. First model initialization exposed an AtomWorks constructor change in 3.0.2/3.0.3; pinning 3.0.1 resolved it without patching Foundry/model code. Worker-native parsing passed again. Inference completed in 24.29 seconds (interleaved) and 23.81 seconds (tandem control).

Both actual outputs have 170 scaffold residues, pass motif identities/geometry (maximum block RMSD about 0.0504 Å in one joint target frame), and pass the 1.1–1.8 Å peptide C–N continuity screen. These are generation and geometry checks, not independent folding or binder validation. Outputs, successful/failed logs, checkpoint SHA256 and frozen runtime are in `examples/gpu_smoke_results`; open its review.html.

The results archive was downloaded and SHA256-verified before pod `02l2bw7tyfpw7n` was deleted. A subsequent API read returned pod-not-found; the temporary account SSH key was removed. GPU elapsed-time cost is estimated at $0.19 plus small disk charges, well below the $100 cap. The billing endpoint had not posted a record immediately after deletion. No network volume was created or retained.

## Baseline pivot

Interleaving is retired from active motif selection; new motif plans select only tandem. Archived comparisons remain available. The active work is now ordinary folding followed by ordinary single-target binder design, not shared-scaffold redesign. `biocad baseline` prepares a sequence-only RF3 KN035 control and motif-free RFD3 PD-L1 binder inputs with explicit normalized residue mapping. Both native input constructors passed locally. Twenty-six repository tests pass. See `examples/baseline_loop/README.md` for staged execution, model choices and remaining adapters. No new GPU prediction has run for this baseline.

## Full BoltzGen pipeline completed (2026-10-17)

The active baseline now uses the off-the-shelf BoltzGen 0.3.2 `protein-anything` protocol,
pinned at `a07fea1153111ef453aaea7a15697ee05de9291f`. The user removed epitope constraints;
the completed test uses no hotspot/binding-site mask, a 115-residue observed PD-L1 target
from 5JDS and two 75-residue designed binders. All seven stages completed: design,
SolubleMPNN inverse folding, Boltz-2 complex folding, binder-alone folding, ESMFold2
scoring, analysis and filtering. This supersedes the unexecuted custom RF3 baseline.

One candidate (`pdl1_1`) passes the upstream default geometry/composition filters.
`pdl1_0` fails complex-recovery RMSD (2.6117 Å, threshold 2.5 Å). Both have low interface
confidence: Boltz-2 ipTM about 0.30/0.29 and ESMFold2 minimum ipSAE about 0.054/0.040.
Neither is a confirmed binder or a strong computational binding hit. The current native
protocol ranks by interface-related metrics but does not hard-filter for high interface
confidence. Keep this distinction when consuming `pass_filters`.

Results, candidate sequences/coordinates, native CSV/PDF, all intermediate artifacts,
configuration, failed-attempt logs, frozen main/ESMFold2 runtimes and model references
are retained under `examples/boltzgen_results`; open review.html. The complete downloaded
archive was SHA256-verified before deleting pod `g4zwiy9n5grhlq`; subsequent API read
returned not-found and the temporary account SSH key was removed. Estimated GPU cost
was $0.846 plus disk charges; no billing record was posted immediately after deletion.

Runtime lessons incorporated into the reusable scripts: use local container disk for
Python environments/uv cache on network-backed pods, cap OMP/MKL threads at four, and
use standard Torch kernels (`--use_kernels false`) with the image's Torch 2.8.0. Optional
cuEquivariance 0.12 called a newer Torch API and failed; model code/weights and scientific
filters were not patched. ESMFold2 runs in upstream's separate pinned Torch 2.11/CUDA13
runtime. The adapter and retained rejection semantics are covered by 28 local tests.
