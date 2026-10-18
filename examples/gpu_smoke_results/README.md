# Actual RFD3 GPU smoke outputs

Two single-sample generations completed on 2026-10-17 (America/Los_Angeles). Open review.html for projections and coordinate links. These are generated scaffold backbones around protected motifs, not sequence-designed/refolded or experimentally validated binders.

Both 170-residue outputs pass exact retained identity, joint target-frame motif geometry and peptide continuity checks. Do not interpret these checks as independent folding confidence. The shared PD-L1 target arrangement is a computational fixture.

The first startup failed with AtomWorks 3.0.3; the successful runtime pins 3.0.1. Both logs are retained. The checkpoint SHA256, exact Foundry commit, Torch/CUDA/device and frozen packages are recorded. The Runpod pod was deleted after a checksum-verified download, and its temporary SSH account key was removed. Estimated GPU charge is approximately $0.19; finalized billing was not yet posted.
