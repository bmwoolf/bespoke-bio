> Historical research / deferred experiment. This is not the active design pipeline. Decisions and "next experiment" statements below reflect earlier work. The shipped path is [BoltzGen generation and API integration](INTEGRATION.md); see [current status](STATUS.md).

# Compact assembly first

User direction, 2026-10-17: prioritize compact, interdigitated assembly; defer conformational modeling. Keep the public PD-L1 benchmark. GPU provisioning and budget remain undecided.

## What runs without a GPU

```sh
uv sync --extra design --extra dev --locked
uv run --extra design biocad placements --bundle examples/replay --output outputs/compact-assembly
uv run --extra design --extra dev pytest -q
```

The placement command reads local fixtures only. It does not need model weights, remote compute, internet, or Foldseek. The broader `biocad demo` command additionally uses Foldseek for the bounded search and partner-pose transfer.

The current search places two copies of the complete observed nanobody while preserving each nanobody's relative pose to its own target. It samples 48 combinations of translation and rotation, plus an explicit overlapping control. The second target co-moves with its module. Relative target geometry is therefore hypothetical, rather than a fit to two sites on an already observed assembly.

The source nanobody contains crystallization nonprotein atoms. This search exports a new protein-only component revision with its own ID; it retains the exact observed polypeptide sequence and coordinates and records the parent component ID. Target author-chain context remains present. Full biological assemblies, membranes and glycans are not added automatically.

## Screening and results

Forbidden pairs are key–key, target–target, and each key against the other target. A pair fails if a heavy-atom distance is below the sum of Gemmi's element-specific van der Waals radii minus 0.4 Å. This is stricter than the prototype's 1.5 Å coincident-atom screen, but is still geometric triage rather than an energy model.

A connection also fails if its terminal C-to-N anchor distance exceeds the linker's necessary contour upper bound, `3.8 × (linker residues + 1)` Å. Being below that bound does not establish a viable linker; the report preserves extension fraction and marks the connection unmodeled.

Among surviving layouts, rank by the two keys' Cα radius of gyration, then terminal-anchor distance, then stable ID. Raw clash counts, gaps, compactness, transforms and required-interface metrics remain visible. Contact counts are not treated as binding energies. Targets are excluded from the compactness objective.

The default local run produced **5 layouts passing the screen out of 48 samples**. The best has a key Cα radius of gyration of **23.96 Å**, compared with **51.76 Å** for the same intact keys separated by 100 Å. Its linker anchors are **46.66 Å** apart, about **77%** of the 15-residue contour estimate. This is a reduction in spatial spread, not protein sequence length. The overlap control fails; replacing the linker with one glycine leaves no surviving layout in the tested grid.

Outputs include `placement_search.json` with all accepted/rejected placements, exact transforms, `compact_*.cif` for the top three, a spread-out reference, an offline `review.html`, proposed fusion FASTA, request schema and a hashed manifest. `redesign_brief.json` identifies observed binding-contact residues for each instance and carries the selected geometry forward as an explicitly unexecuted redesign proposal.

## Shrinkers versus integrated scaffold design

**Raygun** is a plausible interpretation of “protein shrinker.” It takes a template sequence and generates length-controlled variants with insertions, deletions and substitutions. The authors report miniaturization experiments, including fluorescent proteins. This makes it worth testing on an individual oversized module. Its documented task does not establish joint preservation of two binding poses, connectivity, or ternary-complex compatibility. Sequence shortening alone also does not imply a smaller spatial envelope. [Raygun primary repository](https://github.com/rohitsinghlab/raygun), [original preprint](https://www.biorxiv.org/content/10.1101/2024.08.13.607858v1).

**SCISOR** proposes sequence deletions using a diffusion model and evaluates deletion effects and motif retention. It is another single-module shortening candidate, not a joint geometry or interface assembly validator. Neither shrinker is installed or executed here. [SCISOR primary repository](https://github.com/baronet2/SCISOR), [paper](https://arxiv.org/abs/2511.07390).

For **interdigitated shared-scaffold design**, my first proposed experiment is motif scaffolding: keep both required recognition geometries, then rebuild the scaffold/connections around them. RFdiffusion demonstrates structural motif scaffolding, and RFdiffusion3 exposes atom-level fixed-coordinate and sequence-conditioning selections. This is a better conceptual match to joint integration. It is not guaranteed to work for two complete recognition surfaces, and an executable, version-pinned model configuration is still required. [RFdiffusion primary paper](https://pmc.ncbi.nlm.nih.gov/articles/10468394/), [RFD3 input parser](https://github.com/RosettaCommons/foundry/blob/production/models/rfd3/src/rfd3/inference/input_parsing.py).

The exported binding-contact mask is only a starting constraint. It does not establish that all other residues can be removed: folding cores, disulfides and framework support need explicit preservation/redesign decisions. No loop is cut, no residue is deleted, and no arbitrary insertion is claimed to be connected by the local search.

## Next experiment

Start from a compact layout that passes the geometric screen; retain its exact module/target transforms and protected binding regions. Define the essential cores and permitted connection topology, then compare a motif-preserving shared-scaffold redesign against a simpler terminal fusion. Use separate full-atom predictions for the isolated construct and joint complex, assess both binding interfaces, and recheck compactness after prediction. Evaluate conformation/state ensembles after this assembly path is working.
