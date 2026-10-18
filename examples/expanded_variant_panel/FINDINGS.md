# Expanded binder and variant study — completed

Across both campaigns: **22 proteins, 28 sequence states, 132 generated binder candidates and 243 repeated comparison predictions**. The expansion contributed 18 proteins, 21 states, 104 candidates and 168 comparisons. All 21 expansion jobs completed. IL10 and EPO were prepared but not admitted to the timed queue.

**15/132 candidates passed native BoltzGen filters** (12/104 in the expansion). This is a structural-screening count, not 15 validated binders. Native filters can pass candidates with low interface confidence. No binding, affinity, variant selectivity or therapeutic activity was measured.

## Findings worth retaining

| Target / candidate | Native Boltz interface | Native ESM ipSAE | Mean repeated interface score | Interpretation |
| --- | ---: | ---: | ---: | --- |
| PF4 / binder_ec321e7a89fa50d1 | 0.670 | 0.744 | 0.638 | Repeatable monomer-model interface, but severe clash with the experimental tetramer. |
| CXCL12 / binder_249d6d7357137754 | 0.691 | 0.834 | 0.639 | Repeatable monomer-model interface, but severe overlap with the experimental dimer. |
| HMGB1 / binder_08f2d30ba3c804fb | 0.621 | 0.454 | 0.401 | Moderate computational support; biological context and binding untested. |
| FGB R478K / binder_6cbcee43955fdbd5 | 0.541 | 0.252 | 0.088 reference / 0.092 R478K | Native pass does not persist as a confident interface in the untemplated comparison protocol. |
| MSTN / binder_46c8fff417a05c9a | 0.479 | 0.345 | 0.115 | Weak repeated interface support; native dimer context not assessed. |

Native evaluation and repeated comparisons are different protocols: native evaluation uses a generated target template, while repeated comparisons use fixed sequences without target templates or MSAs. Their absolute score differences are **not variant effects**. The table shows model robustness across protocols, not comparable affinity estimates.

## What the variant experiment shows

- SERPINA1 reference, E288V and E366K each produced 0/8 native passes. Seven selected exploratory binders were compared against all three target states, three times per pair. Mean interface scores remained low (approximately 0.08–0.15).
- FGB reference produced 1/8 native passes and R478K produced 2/8. Four selected binders, including mutation-contact failures, were compared against both states. Mean scores were low for both states (approximately 0.09–0.12).
- The earlier APOE, F7 and GDF15 comparisons likewise did not establish a convincing variant-selective binder.
- These are small exploratory samples. Different design cohorts contain different binder sequences; their pass-rate differences are not fixed-binder selectivity. Matched seed labels do not imply identical random draws, since native feature randomness remains uncontrolled.
- SERPINA1 E366K is a controlled single substitution on the canonical sequence background, not a reconstruction of the complete natural Z haplotype.

**Conclusion: the pipeline can generate and compare candidates, but this experiment does not establish variant-specific binding.** Low target-fold confidence limits interpretation; a small variant effect cannot be distinguished reliably from model limitations here.

## Folding calibration

The repeated target-only predictions were also compared with the exact experimental fragments used for ten design inputs. FGF2 reproduced its fragment at 2.07–2.28 Å Cα RMSD; other targets were substantially less reliable (for example IL1RN: 13.75–15.55 Å). SERPINA1 target-only mean local confidence was approximately 0.38 and FGB approximately 0.53 on the native normalized 0–1 scale. This argues for validating the folding/comparison protocol before making biological claims about structural variant effects. Replicate predictions are model samples, not an equilibrium conformational ensemble.

## Completion and next step

All 2,726 expansion files matched their remote SHA256 hashes; all 168 repeated outputs passed exact sequence, finite-coordinate and confidence-range checks. Including the earlier campaign, 3,561 files were verified. Both campaign GPUs were deleted and deletion confirmed. Combined listed-rate estimate: **$20.87 excluding storage**, not a provider invoice, within the $100 authorization. The team API and persistent model cache were retained.

Before more design volume, prioritize a calibrated target-folding protocol and assembly-aware target inputs. Then test a small set of fixed binders on matched reference/variant structures with an independent validation model. Wet-lab binding remains necessary.
