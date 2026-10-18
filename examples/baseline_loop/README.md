# Ordinary folding and single-target binder baseline

Interleaving and assembly redesign are retired from the active path. Historical smoke outputs remain unchanged. `biocad baseline --output examples/baseline_loop` prepares these inputs.

1. RF3 folds the 127-residue observed KN035 construct from sequence alone, with no structure/template supplied. The experimental coordinate reference is for evaluation only, not model input. This is an observed-only computational truncation and a plumbing control, not a leakage-free scientific benchmark. Start single-sequence; record early stopping/missing output as a failure requiring diagnosis.
2. RFD3 generates ordinary 75-residue binders against one observed PD-L1 domain, with three public contact-site hotspots. It preserves no binder motif and rearranges no modules. Target author IDs are mapped explicitly to normalized model IDs.
3. ProteinMPNN assigns sequences to the generated binder backbone with the target held fixed. Prefer the original implementation for this baseline: the pinned Foundry MPNN documentation cautions that its reimplementation is still being validated.
4. RF3 predicts each discrete binder sequence alone and together with PD-L1, from sequence inputs rather than the designed coordinates. Compare scaffold recovery, genuine model confidence, interface location and clashes. Keep monomer and complex predictions distinct.
5. Review a bounded batch (two backbones, two sequences each). Reuse a promising backbone for sequence redesign, or generate another batch if it fails. Do not convert missing model evidence into a passing score or treat confidence as measured affinity. Negative-panel/multistate design comes later.

Native RF3 sequence-input construction and native RFD3 binder-input construction passed locally. Neither model has run on these baseline inputs yet. Sequence-design integration, prediction output ingestion and calibrated selection remain unbuilt; plan.json records those gaps. Earlier GPU runs demonstrated motif scaffold generation only.

Primary implementation references: [Foundry](https://github.com/RosettaCommons/foundry), [original ProteinMPNN](https://github.com/dauparas/ProteinMPNN). Model inputs follow the locally pinned Foundry source rather than unpinned current examples. The existing $100 phase cap remains the ceiling for a future GPU execution; no new pod is active.
