# Folding validation and focused design

Status: the 64-design focused run and all 18 standard fixed-sequence predictions are complete.
Three native-filter passes remain computational candidates; none has experimental binding evidence.

## What changed

The validation path now uses standard Boltz 2.2.1 with its supported
`boltz2_conf.ckpt`, frozen target MSAs, exact sequence checks, three seeds, three
recycles and 200 sampling steps. It evaluates geometry as well as confidence.
The previous empty-MSA comparisons and BoltzGen's native refolding results remain
separate historical evidence, not interchangeable confidence measurements.

For each natural target, a full MSA × template control grid was run. The following
are three-seed medians; pLDDT uses the 0–100 scale and RMSD is mapped Cα distance
to a withheld experimental structure in Å.

| Target | Sequence only pLDDT / RMSD | MSA only | Template only | MSA + template |
|---|---:|---:|---:|---:|
| FGF2, 124 aa | 92.3 / 2.10 | 96.2 / 0.33 | 96.3 / 0.31 | 96.5 / 0.29 |
| IL1RN, 145 aa | 32.0 / 14.66 | 94.3 / 0.54 | 95.0 / 0.56 | 95.8 / 0.55 |
| SERPINA1, 394 aa | 83.1 / 9.30 | 86.7 / 1.65 | 87.8 / 1.76 | 85.8 / 1.51 |

MSAs substantially improve these natural targets. SERPINA1 also demonstrates why
high confidence alone is insufficient: sequence-only prediction is confident but
geometrically inaccurate. Withheld evaluation coordinates are not necessarily
absent from model training. These results do not establish reliable folding of
arbitrary de novo targets without homologs.

## Existing candidates and positive control

Three previously selected designed binders retain strong binder-fold confidence.
PF4 and CXCL12 monomer complexes look promising computationally, but their
tetramer/dimer tests do not preserve the tested assemblies reliably. These are
not evidence of binding intact native oligomers. HMGB1 preserves its individual
HMG boxes while changing their relative orientation; whole-protein RMSD should
not be interpreted as unfolding of those boxes. Construct and redox caveats remain.

The experimental IGF1–IGFBP4 control with official paired and unpaired MSAs for
both natural partners recovers 86–93% of experimental contacts across three seeds.
Its reference-contact Cα joint-fit RMSD is 0.39–0.40 Å, and binder pLDDT is 96.4.
Whole-target-aligned binder RMSD remains 4–5 Å. The empty-binder-MSA condition
recovers only 5–16% of contacts. This is a different MSA treatment from a de novo
binder and cannot establish experimental binding of our new sequences. A shuffled
sequence is a putative negative control, not a proven nonbinder.

## Variant result

SERPINA1 S and Z were modeled as single substitutions on the canonical background,
with unchanged reference homolog alignments and changed query sequences. These
are not full phased haplotypes. On 367 shared experimentally observed core
positions, reference-to-variant median RMSDs are 0.705 Å and 0.807 Å. The
within-reference median is 0.693 Å; between-allele and within-seed ranges overlap.
No convincing variant-dependent structural difference is established. Flexible
terminal residues are reported separately, and seeds are not physical ensembles.

## Completed focused design and validation

The run produced 32 designs per target, 75 residues each, without a specified epitope.
Two FGF2 and one IL1RN candidate passed native filters. All three received standard
Boltz complex and binder-alone predictions over three seeds: 18/18 completed.
Complexes used the frozen target MSA, an empty binder MSA and **no coordinate templates**.

| Candidate | Binder-alone pLDDT range | Both directional interface scores | Held-out target RMSD, Å | Pose RMSD to design, Å |
|---|---:|---:|---:|---:|
| FGF2 `446e69445d74fbf2` | 94.47–94.71 | 0.868–0.910 | 0.283–0.311 | 51.05–51.42 |
| FGF2 `8d7f64798b00c4dc` | 92.22–92.89 | 0.433–0.617 | 0.311–0.493 | 5.45–43.57 |
| IL1RN `70c5ec0eca4403d4` | 97.08–97.14 | 0.261–0.515 | 0.563–0.605 | 45.32–47.55 |

FGF2 `446e` is the strongest computational follow-up: its alternate docking pose
repeats across seeds (0.93–1.50 Å binder RMSD after target alignment), with no severe
heavy-atom overlap. However, it recovers zero original design contacts and its native
ESMFold2 minimum ipSAE is only 0.0138. Standard Boltz confidence does not resolve this
cross-method disagreement. This is an alternative predicted pose, not a validated binder.

The other FGF2 candidate and IL1RN retain credible standalone folds but variable
complex poses (seed-pair RMSDs 19.62–44.08 Å and 9.85–23.30 Å). Their ESMFold2 minimum
ipSAE scores are zero. Native pass status remains historical and separate from these
findings. All three predicted targets remain close to their held-out experimental
structures; these comparisons were performed after prediction.

See the [public scientific findings](../examples/folding_validation/FINAL_FINDINGS.md)
and [metrics figures](../examples/folding_validation/summary_figures/README.md).
FGF2 receptor/GAG context and IL1RN receptor competition remain biological questions.
No affinity, activity, variant selectivity or therapeutic utility is established.

Preparation and analysis are reproducible with `prepare_standard_validation.py`,
`analyze_standard_validation.py`, `prepare_candidate_validation.py`,
`analyze_candidate_validation.py`, and `prepare_focused_validation.py` in `scripts/`.
Live logs, full local structures and checksum receipts are under
`outputs/folding_validation/`; that private execution directory is not shipped
by Git. The allowlisted public package is in `examples/folding_validation/`; operational logs,
credentials and runtime caches are excluded.

Final cleanup: all scientific results downloaded and SHA-verified; campaign GPU `3oih9qanxyyj4k` deleted and API not_found verified. Estimated GPU compute $9.78, excluding storage. Team API and network volumes retained. Full local scientific archive: `outputs/folding_validation/scientific_campaign.tar.gz` (~177 MB). Public scientific handoff: `examples/folding_validation/`.
