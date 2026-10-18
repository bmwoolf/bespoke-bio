> Historical research / deferred experiment. This is not the active design pipeline. Decisions and "next experiment" statements below reflect earlier work. The shipped path is [BoltzGen generation and API integration](INTEGRATION.md); see [current status](STATUS.md).

# Motif-preserving scaffold redesign: local/GPU boundary

The next experiment replaces two complete recognition domains with one connected scaffold carrying both sets of protected recognition geometry. This is a motif-grafting hypothesis, not a validated way to remove the nanobody framework. Conformational landscape design remains deferred.

## Local preparation already executed

```sh
uv run --extra design biocad motifs --assembly examples/compact_assembly --output outputs/motif-plan
```

The local fixture yields:

- Four coherent blocks: source K20–35, K94–118, L20–35 and L94–118, retaining 82 residues in total.
- Both contact patches per module, two flanking residues around selected positions, and the endpoints of the two inferred native disulfides per module. SG distances are approximately 2.03 and 2.08 Å. These bonds are inferred from deposited coordinates, not imported bond records.
- Six block orders preserving the original order within each module; four interleave module provenance and two are tandem controls.
- One interleaved order and one tandem control, matched to a planned 127-residue scaffold, including generated segments. This is a design length specification, not a generated 127-residue protein.
- Source-to-normalized residue maps, fixed-coordinate atom selections and predicted motif spans for each proposed output chain.
- Draft RFdiffusion3 indexed-contig inputs pinned to Foundry commit `d078b24d6d0e6f8607e9a0311128f4a897f83953`.

The existing packed layout supplies relative motif/target coordinates. The two targets still occupy a hypothetical relative pose. The new motif input explicitly contains protein atoms only; full original context remains in the source bundle and must be rechecked after prediction.

`model_input.cif` preserves the original atoms and coordinates under explicit new residue numbering. `motifs.cif` displays the four protected segments as separate chains. No discontinuous patch is silently concatenated. `motif_plan.json` lists unsupported assumptions, selected bridge lengths, rejected topologies and acceptance checks. A one-residue bridge cap rejects every tested order.

## What to isolate and test

| Component | Local CPU tests now | Model/GPU work |
|---|---|---|
| Motif definition | Exact residue/atom identity; contiguous blocks; flanks; inferred disulfide endpoints | Whether those blocks retain enough scaffold support to bind |
| Topology/compiler | Source/output mapping; matched chain-length comparison; explicit peptide gaps and chain breaks; rejection of impossible contour budgets | Generate a backbone that connects those motifs without broken chains or strained geometry |
| Geometry validator | Known global rotation/translation; per-motif target-relative deviations; missing residues/atoms; contact identity | Apply unchanged validator to real generated and independently predicted outputs |
| Sequence assignment | Fixed/editable index masks and output contract can be checked locally | ProteinMPNN on generated regions; CPU inference is possible with weights, but GPU is preferable for iterative use |
| Structural validation | Parse real outputs, check bond lengths/disulfides/clashes, calculate compactness and interface geometry | Independently predict final discrete sequence alone and in the full two-target complex |
| Biological interpretation | Keep confidence, affinity and function separate | Models do not establish affinity or biological function; those require experiments |

The CPU motif checker uses **one joint target-frame alignment**, then assesses every protected block. It does not independently align each motif, which could conceal loss of its binding pose. It checks RMSD and maximum atom displacement and rejects missing results. Its passing result is explicitly limited to motif geometry; it does not pass a connected construct with absent bridges.

The synthetic tests contain only target and fixed-motif geometry, with deliberate movements or omissions. They are not model outputs and supply no learned scores.

## Model choice and native validation

RFdiffusion3's documented indexed contigs can connect fixed sequence/coordinate fragments with generated residues and retain separate target chains. This is the proposed generative route for assembly redesign; Mosaic remains the existing binder/objective adapter. We are adding one specific assembly experiment, not an interchangeable-engine framework. [Pinned input documentation](https://github.com/RosettaCommons/foundry/blob/d078b24d6d0e6f8607e9a0311128f4a897f83953/models/rfd3/docs/input.md).

For the first test, enumerate explicit block orders instead of relying on unindexed motifs to discover their order. Maintainers state that adjacency in unindexed motifs is conditioning rather than a hard constraint, so native output mappings must be inspected. [Maintainer discussion](https://github.com/RosettaCommons/foundry/issues/395).

The draft JSON has passed our own geometry/mapping tests only. **Native parser validation has not run.** That validation can run on CPU in an appropriate Foundry environment without model weights; it is not intrinsically a GPU task. The current lightweight environment does not include Foundry/AtomWorks. The supplied `gpu_smoke_test.sh` uses `prevalidate_inputs=true` before its inference run; that flag is not a standalone CPU-only validation command.

Full heavy-atom motif preservation can overconstrain generation. Reported examples show the importance of checking that generated residues actually connect to constrained atoms. If the first smoke test fails, inspect constraints and chain geometry before increasing a campaign budget. [Maintainer analysis](https://github.com/RosettaCommons/foundry/issues/372).

## Smallest useful GPU experiment

1. In an isolated environment at the pinned commit, validate the two native inputs and confirm source selections, output chain numbering and fixed sequence/coordinate semantics. Pin and hash the actual model checkpoint.
2. Generate **one backbone for each of the two matched-length topologies**. These are two structural smoke tests, not a success-rate estimate. They model the 127-residue proposed scaffold plus two 115-residue target contexts.
3. Use `check_motif_geometry` with the exported mapping and verify generated-region coverage and every new peptide boundary. Reject missing, scrambled or disconnected motifs; do not interpret `select_fixed_atoms` as proof that the output is chemically connected.
4. If the geometry passes, inverse-fold generated regions while retaining protected identities. Independently predict the discrete construct alone and with both targets. Evaluate the two interfaces separately, the joint assembly, disulfide pairing, clashes and compactness after reprediction.
5. Expand to a few independent samples only after the source mapping, connectivity and validation checks work on actual outputs.

No GPU runtime/VRAM estimate is claimed before profiling this exact input. A GPU is the practical route for repeated backbone generation and complex prediction; it is not required for the local preparation and checking code.

Framework/core support, off-targets, folding, affinity and simultaneous occupancy remain unknown. All four block orders labeled “interleaved” describe sequence provenance only. Actual structural interdigitation and a shared folded core must be demonstrated in generated, repredicted and eventually experimental structures.
