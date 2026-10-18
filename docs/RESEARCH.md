> Historical research / deferred experiment. This is not the active design pipeline. Decisions and "next experiment" statements below reflect earlier work. The shipped path is [BoltzGen generation and API integration](INTEGRATION.md); see [current status](STATUS.md).

# Research and first implementation choices

Reviewed 2026-10-17. This is a targeted primary-source scan, not a comprehensive benchmark or a claim that one method is universally best. Recent preprints and project-page success rates are author reports; they do not establish performance on our targets. Software inspected does not imply GPU execution.

## First pipeline

Use Mosaic with Boltz2 and soluble ProteinMPNN: site-conditioned contact loss, internal-contact and inverse-folding recovery objectives; soft optimization, sharpening, inverse folding, discrete refinement, then fresh full-atom complex and monomer predictions. Start with one sequence/primary target and a terminal-clearance postfilter. This follows the requested custom-objective direction while keeping scope bounded. Mosaic is an optimization framework, not a complete screening pipeline. Its official example itself recommends BindCraft for a simpler minibinder workflow. [Pinned Mosaic source](https://github.com/escalante-bio/mosaic/tree/1eb4af64655870b45a715a5a7e61570e8b098701).

Source inspection verified binder-first features, target-relative epitope indices, and `chain_pair_iptm()`. The generic iPTM loss pools all target chains, so it cannot stand in for separate required interfaces. Design features do not provide final amino-acid-specific side chains; `target_only_features` must rebuild both discrete chains before atomistic checks. The implementation preserves objective weights and raw terms. [Structure losses](https://github.com/escalante-bio/mosaic/blob/1eb4af64655870b45a715a5a7e61570e8b098701/src/mosaic/losses/structure_prediction.py), [Boltz2 wrapper](https://github.com/escalante-bio/mosaic/blob/1eb4af64655870b45a715a5a7e61570e8b098701/src/mosaic/models/boltz2.py).

## Other viable approaches

| Approach | Why it matters | Decision |
|---|---|---|
| BindCraft2 | Current documentation exposes multistate targeting/detargeting, termini accessibility and several binder formats | Most relevant simpler comparison; documentation reviewed, not installed or validated |
| BoltzGen | Integrated generative design workflow | Viable alternative, deferred to avoid another engine |
| RFdiffusion3 → MPNN → RF3 | All-atom constrained backbone generation with documented pipeline | Strong route for explicit atom/geometry conditioning; more setup than the first custom-objective loop |
| Proteina-Complexa | Joint sequence/structure generation plus inference-time search; Steinegger collaboration | Highest-priority follow-up generative comparison; compare on the same state panel and compute budget |
| Biohub ESMFold2 / ESMC | New open prediction/design ecosystem released May 2026 | Candidate later predictor or objective; requires target-specific validation |

Primary sources: [BindCraft2](https://github.com/PacesaLab/BindCraft2), [BoltzGen](https://github.com/HannesStark/boltzgen), [official RFD3 tutorial](https://github.com/RosettaCommons/foundry/blob/production/models/rfd3/docs/tutorials/binder_design_tutorial.md), [Proteina-Complexa code](https://github.com/NVIDIA-BioNeMo/Proteina-Complexa), [NVIDIA project and validation](https://research.nvidia.com/labs/genair/proteina-complexa/), [Biohub release](https://biohub.org/news/world-model-of-protein-biology/).

## “Fluid” modeling: separate three layers

1. **Retrieve and rearrange parts.** Foldseek supplies domain/complex analogy; align the target then transfer the associated partner using the same fitted transform. Preserve the original interaction provenance. Its multimer/interface commands concern input complexes, not discovery of a biological interface from an isolated monomer. The downloaded build exposes `easy-interfacesearch`. [Foldseek](https://github.com/steineggerlab/foldseek).
2. **Find discontinuous local chemistry.** Folddisco indexes spatial motifs without inventing a continuous chain from discontinuous residues. Keep it optional until a concrete motif-retrieval question justifies another index. Whole-proteome indexes can be large; begin with a declared small library. [Folddisco](https://github.com/steineggerlab/folddisco).
3. **Model states and flexibility.** A layout graph stores alternative rigid arrangements, not thermodynamic populations. Start with named experimental/predicted conformations and evaluate every candidate against each. Later, add ensembles, linker modeling and local relaxation. Foldseek alone provides neither dynamics nor binding complementarity.

AlloGen (June 2026) specifically treats conformation-selective design through differential state scoring and reports calmodulin experimental validation. This is closely aligned with positive/negative state panels; availability and generalization still need code-level investigation before adoption. [AlloGen preprint](https://arxiv.org/abs/2606.05474).

Correia and collaborators' AmyBind work (June 2026) reports polymorph-specific α-synuclein fibril binders, with binding/labeling behavior separated from downstream fibril effects. This supports explicit state identity and separate functional evidence in our contracts. It does not establish a general molecular machine. [AmyBind preprint](https://www.biorxiv.org/content/10.64898/2026.06.13.732066v1).

Two recent leads are research watch items, not dependencies: ComplexDesign (June 2026) addresses flexible relative target placement for ternary complexes, but its abstract says code will be released upon publication; T-REX (September 2026) studies budget-aware campaign orchestration. Neither is necessary for the first bounded loop. [ComplexDesign](https://www.biorxiv.org/content/10.64898/2026.06.21.733655v1), [T-REX primary metadata/abstract](https://pubmed.ncbi.nlm.nih.gov/42818266/).

An October 13 Correia-associated metamorphic-landscape preprint surfaced in search, but the original bioRxiv page could not be retrieved. Its full methods and implementation were not verified, so it does not support an implementation choice here. [Original preprint link](https://www.biorxiv.org/content/10.64898/2026.10.13.758820v1).

## Fixture and biological interpretation

The initial target is human PD-L1 with KN035 in 5JDS. 5J89 is an alternate inhibitor-bound crystal context, not a patient variant. The known complex provides real geometry without pretending a model has generated a binder. The example upstream TP53 transcript/variant records remain illustrative; phase, complete sequence, exact construct and structure are needed before using them as actual design inputs. RNA TPM is not protein abundance. [5JDS](https://www.rcsb.org/structure/5JDS), [5J89](https://www.rcsb.org/structure/5J89).

No new degradation mechanism is assumed. A recognition module plus an existing recruiter will need compartment, recruited machinery, connectivity and joint structure checks; affinity, structural confidence and downstream activity remain distinct evidence categories.
