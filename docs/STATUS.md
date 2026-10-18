# Current implementation status

The shipped path is **target intake → BoltzGen → predicted candidates and evidence**.
The feed format is proposed for team review. Production CAD/catalog UI, interdigitation,
conformational design and folded multi-module assembly are deferred.

| Component | Verified behavior | Limits |
|---|---|---|
| Upstream translation/package | CDS/variant translation, sequence hashing, RNA provenance, optional structure lookup/folding adapters | Raw sequencing processing belongs upstream; optional external folding is not automatically run |
| Target intake | Sequence, PDB/mmCIF chain or upstream package; exact hash/sequence checks and residue mapping | One polypeptide target per child job; no automatic biological-assembly expansion |
| Feed integration | Python client, `submit-bundle` CLI, independent child jobs, status/progress, typed results and ZIP downloads | Proposed team contract |
| Binder generation | Full pinned BoltzGen seven-stage protocol, no specified epitope | No bespoke target prioritization or experimental hit guarantee |
| Candidate export | Exact sequences, binder/complex mmCIF, native filter metrics/flags, confidence, target linkage and provenance | Filter pass is separate from binding evidence |
| Foldseek | Tested bounded search of 1–100 declared local templates and transferred partner poses | Optional library/CLI branch; not automatically invoked by the design API and not a complementarity detector |
| GPU cache | Fresh-pod restore, checksum verification, both CUDA checks and model-cache reuse in 76 seconds | Setup timing excludes inference; volume persists after pod deletion |
| Verification | Automated tests; real structure-input run and fresh sequence-only HTTP run | Synthetic feed tests and public PD-L1 benchmarks, not novel teammate-target or assay validation |

## Retained scientific results

The structure-input PD-L1 run produced two 75-residue candidates: one passes native
filters and one fails complex recovery RMSD (2.61 Å versus a 2.5 Å threshold). Both have
low interface confidence (ipTM about 0.30). See [the report](../examples/boltzgen_results/review.html).

The fresh sequence-only HTTP run produced one native-filter-passing candidate. Complex
recovery RMSD is 0.83 Å, binder-alone recovery 1.33 Å, ipTM 0.41 and ESMFold2 minimum ipSAE
0.172. A post-run alignment gives target Cα RMSD 1.36 Å against the 115 observed PD-L1
residues; those experimental coordinates were not input to generation. See
[the result](../examples/boltzgen_api_results/summary.json). Neither run establishes measured binding.

## Compute state and shared service

The original disposable benchmark GPU pods were deleted, with deletion verified. The 100 GB persistent
cache volume remains in US-MO-1 at approximately $7/month. The setup, fresh API run and
restore-test phase cost approximately $1.14 in estimated GPU charges, plus storage; these
are estimates rather than finalized billing. The shared service now uses an authenticated CPU gateway, separate persistent job storage,
a durable serial queue and managed disposable GPUs. It reserves $10 per executed job
against the phase cap, with $13 reserved for prior/service costs. Fresh jobs can be
submitted explicitly without manually configuring a worker. Queued cancellations release
their reservation; running jobs retain their conservative reservation. The gateway expires
after 24 hours. See [shared operation](SHARED_SERVICE.md).

## Remaining work

Agree the team feed/deployment contract; run a novel incoming target; define any required
selectivity/state panels; and validate promising candidates independently and experimentally.
See [the punchlist](PUNCHLIST.md). Earlier CAD and model-selection notes are retained as
[history](LEGACY.md), not current instructions.

## Shared service proof

Managed startup, runtime restore, full native inference, artifact download and verified GPU
deletion also completed through the shared HTTPS API. Queued cancellation and persistent
status/budget after restart were checked against the deployed service. See
[the sanitized proof](../examples/shared_service_smoke/summary.json). A CPU gateway remains
active for the 24-hour demo lifetime; the managed test GPU is deleted.
