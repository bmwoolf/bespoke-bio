# Folding validation scientific handoff

Completed computational evidence; scientific interpretation remains required.

See science.json for scores and structure links. Native filters are preserved separately from standard fixed-sequence validation. All native bundle candidates are retained, including failures; counts are not restricted to successful examples.

- Computational predictions only; no experimental binding, selectivity, safety, or benefit demonstrated.
- Native design filters and independent fixed-sequence validation are separate evidence.
- When provided, posthoc target reference geometry is evaluated after prediction; those reference coordinates are not conditioning inputs to the predictor.
- Template-assisted recovery is not independent fold recovery; historical structures may be training data.
- Three seeds measure limited prediction variability, not equilibrium ensembles or statistical significance.
- SERPINA1 primary comparison uses shared experimental core positions, excluding variant sites; all-residue comparisons can be dominated by the unresolved amino-terminal tail.
- SERPINA1 variants are single substitutions on the canonical background, not full natural haplotypes.
- The experimental positive binder and shuffled decoy calibrate this workflow only; the decoy is not a validated nonbinder and may differ in folding.
- Oligomer context and competing biological partners require separate interpretation; a confidence score does not establish binding.
