# Research history and deferred CAD experiments

These files preserve earlier experiments and decision history. They are not the active
API contract, current model choice or next required implementation steps.

| Historical document | Scope |
|---|---|
| [Earlier integration/CAD notes](CAD_HISTORY.md) | Mosaic objectives, experimental replay and geometry contracts |
| [Assembly experiments](ASSEMBLY.md) | Rigid packing and proposed compact/interdigitated assembly |
| [Motif redesign](MOTIF_REDESIGN.md) | Protected motifs, RFD3 preparation and scaffold hypotheses |
| [Research survey](RESEARCH.md) | Earlier model choices and literature leads; not a current recommendation list |
| [Compute development log](COMPUTE_HISTORY.md) | Earlier local, Mosaic and RFD3 environment state |
| [Chronological status notes](STATUS_HISTORY.md) | Development snapshots, including superseded pending/next-step claims |

Historical artifacts remain in `examples/replay`, `compact_assembly`, `motif_preflight`,
`gpu_smoke`, `gpu_smoke_results`, and `baseline_loop`. Their CLI routes remain available
for reproducibility. The `biocad` namespace is still used by the active binder pipeline;
its name does not imply a shipped CAD frontend or validated molecular assembly.

For current work, use [the README](../README.md), [status](STATUS.md),
[integration](INTEGRATION.md), and [feed proposal](FEED_CONTRACT.md).
