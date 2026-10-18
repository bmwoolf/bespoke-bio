# Final scientific findings

The workflow completed 64 new designs and all 18 planned fixed-sequence validation predictions. Three designs passed native filters: two of 32 FGF2 candidates and one of 32 IL1RN candidates. Native filters assess structural recovery and composition; passing is not an interface-confidence or experimental-binding threshold.

## Focused candidates

All three are 75 residues. Standard Boltz 2.2.1 used a frozen target MSA, empty binder MSA, no coordinate templates, and seeds 11/22/33 for both complex and binder-alone prediction. Experimental target structures were compared only after prediction. Values below are ranges across seeds; chain-pair scores include both raw directions.

| Candidate | Binder-alone pLDDT | Directional interface score | Target vs held-out Cα RMSD, Å | Binder pose vs design RMSD, Å | Interpretation |
|---|---:|---:|---:|---:|---|
| FGF2 `binder_446e69445d74fbf2` | 94.47–94.71 | 0.868–0.910 | 0.283–0.311 | 51.05–51.42 | Strongest computational follow-up; repeatable alternate pose, original design pose not recovered |
| FGF2 `binder_8d7f64798b00c4dc` | 92.22–92.89 | 0.433–0.617 | 0.311–0.493 | 5.45–43.57 | Binder fold credible; docking varies substantially |
| IL1RN `binder_70c5ec0eca4403d4` | 97.08–97.14 | 0.261–0.515 | 0.563–0.605 | 45.32–47.55 | Strong standalone fold; interface and pose remain weak/inconsistent |

FGF2 `446e` predicts a similar alternate docking pose in all three seeds: binder Cα RMSD between predictions is 0.93–1.50 Å after fitting the target. It has 26–27 residue-contact pairs, no heavy-atom overlaps below 2 Å, and zero recovery of original design contacts. The target itself remains close to experimental structure. This supports a reproducible **alternative computational pose**, not successful recovery of the designed interface or measured binding.

Method disagreement remains material. Native Boltz ipTM / ESMFold2 minimum ipSAE are 0.157 / 0.0138 for `446e`, 0.249 / 0 for the other FGF2 candidate, and 0.213 / 0 for IL1RN. These are different pipelines and metrics; they are not interchangeable with standard Boltz directional scores. Strong standard Boltz confidence does not resolve the ESM disagreement. Prioritize `446e` for additional independent evaluation, without calling it a validated binder.

The second FGF2 candidate's seed-pair pose RMSDs are 19.62–44.08 Å; original-contact recall is 9–45%. One prediction has one heavy-atom pair below 2 Å. IL1RN seed-pair pose RMSDs are 9.85–23.30 Å and original-contact recall is zero. These two provide good evidence for foldable predicted scaffolds, but limited support for a reproducible binding pose. Seed-pair RMSDs use one proper rigid Cα fit on the target followed by binder comparison; the three pairs are not independent observations.

## What the controls established

All 36 target-only MSA × template control predictions completed. MSA-only median held-out Cα RMSDs are 0.33 Å for FGF2, 0.54 Å for IL1RN and 1.65 Å for SERPINA1, versus 2.10, 14.66 and 9.30 Å with neither MSA nor template. SERPINA1 illustrates why high local confidence alone is inadequate. This supports using homolog evidence for these natural targets; it does not validate arbitrary de novo target folding. Historical experimental structures are withheld from these inputs, not necessarily from model training.

The experimental IGF1–IGFBP4 control with MSAs for both natural partners recovers 86–93% of reference contacts and a 0.39–0.40 Å reference-interface joint Cα fit. Empty-binder-MSA predictions recover only 5–16% of contacts. This differing MSA treatment limits transfer of calibration to de novo binders. The shuffled control is a putative negative, not an experimentally established nonbinder.

SERPINA1 S/Z single-substitution models show reference-to-variant median core RMSDs of 0.705/0.807 Å versus within-reference 0.693 Å, with overlapping ranges. No convincing variant-dependent structural difference or selectivity is established. The primary comparison uses 367 shared experimentally observed core positions excluding variant sites; flexible termini are separate. These are canonical-background substitutions, not complete natural haplotypes, and three seeds are not physical ensembles.

Earlier PF4/CXCL12 monomer results do not establish compatibility with intact native oligomers. HMGB1 domain orientation and redox/construct context remain unresolved. FGF2 receptor/glycosaminoglycan context and IL1RN receptor competition still require biological interpretation. No new sequence has measured affinity, activity, selectivity, safety, or therapeutic benefit.
