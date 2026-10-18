import hashlib
from pathlib import Path

import gemmi
import numpy as np

from .atom_inventory import missing_atoms
from .contracts import Port, ResidueRef, TargetState
from .geometry import frame


def sha(data):
    return hashlib.sha256(data).hexdigest()


def protein_residues(chain):
    return [
        r
        for r in chain
        if r.het_flag == "A" and gemmi.find_tabulated_residue(r.name).is_amino_acid()
    ]


def ref(chain, residue):
    return ResidueRef(chain=chain, number=residue.seqid.num, insertion=residue.seqid.icode.strip())


def key(reference):
    return reference.chain, reference.number, reference.insertion


def sequence(chain):
    return "".join(
        gemmi.find_tabulated_residue(r.name).one_letter_code.upper()
        for r in protein_residues(chain)
    )


def atoms(chain, ca=False):
    xyz = []
    labels = []
    for r in chain:
        for a in r:
            if a.element.is_hydrogen or (ca and a.name != "CA"):
                continue
            # preparation selects one conformer; refuse accidental duplicates upstream.
            xyz.append([a.pos.x, a.pos.y, a.pos.z])
            labels.append((r.seqid.num, r.seqid.icode.strip(), a.name))
    return np.asarray(xyz), labels


def prepare(source, out, chain_names, evidence="experimental", context=None):
    original = gemmi.read_structure(str(source))
    original.setup_entities()
    if len(original) != 1:
        raise ValueError("Explicit model selection required for multi-model input")
    st = original.clone()
    model = st[0]
    for name in [c.name for c in model]:
        if name not in chain_names:
            del model[name]
    if {c.name for c in model} != set(chain_names):
        raise ValueError("Missing requested chain")
    # Preserve original on disk. Retain nonprotein partners, ligands and waters in selected chains.
    st.remove_alternative_conformations()
    mapping = []
    chains = {}
    unresolved = []
    for c in model:
        chains[c.name] = sequence(c)
        for i, r in enumerate(protein_residues(c)):
            mapping.append(
                {
                    "structure": ref(c.name, r).model_dump(),
                    "model_input_index": i,
                    "residue_name": r.name,
                    "label_seq_id": r.label_seq,
                    "missing_standard_heavy_atoms": missing_atoms(r),
                    "observed_sequence_position": i + 1,
                    "missing_backbone_atoms": [
                        n for n in ("N", "CA", "C", "O") if r.find_atom(n, "*") is None
                    ],
                }
            )
    # Deposited unobserved table retained as raw rows; do not infer UniProt numbering from auth IDs.
    is_cif = str(source).lower().endswith((".cif", ".mmcif", ".cif.gz", ".mmcif.gz"))
    block = gemmi.cif.read(str(source)).sole_block() if is_cif else None
    table = (
        block.find_mmcif_category("_pdbx_unobs_or_zero_occ_residues.")
        if block is not None
        else None
    )
    if table:
        tags = list(table.tags)
        for row in table:
            unresolved.append(dict(zip(tags, list(row))))
    out = Path(out)
    st.setup_entities()
    st.make_mmcif_document().write_file(str(out))
    return TargetState(
        id="state_" + sha((sha(Path(source).read_bytes()) + ",".join(chain_names)).encode())[:16],
        source=str(source),
        prepared=str(out),
        source_sha256=sha(Path(source).read_bytes()),
        chains=chains,
        residue_mapping=mapping,
        unresolved=unresolved,
        evidence=evidence,
        context={
            "selected_model": 0,
            "chain_identity": "author chain ID; label_seq_id kept separately",
            "sequence_scope": "observed construct; full deposited polymer sequence in original mmCIF",
            "deposited_entity_sequences": {e.name: e.full_sequence for e in original.entities},
            "altloc_policy": "gemmi.remove_alternative_conformations",
            **(context or {}),
        },
    )


def select_indices(state, references, chain):
    lookup = {
        key(ResidueRef(**m["structure"])): m["model_input_index"] for m in state.residue_mapping
    }
    if any(r.chain != chain for r in references):
        raise ValueError("Site must belong to selected target chain")
    try:
        return [lookup[key(r)] for r in references]
    except KeyError as e:
        raise ValueError(f"Unobserved or absent site residue {e}") from e


def write_chains(chains, path):
    st = gemmi.Structure()
    model = gemmi.Model("1")
    for c in chains:
        model.add_chain(c.clone())
    st.add_model(model)
    st.setup_entities()
    st.make_mmcif_document().write_file(str(path))
    return st


def ports(chain):
    residues = protein_residues(chain)
    result = []
    for end, r, origin, x, p in [
        ("N", residues[0], "N", "CA", "C"),
        ("C", residues[-1], "C", "CA", "N"),
    ]:
        aa = [r.find_atom(n, "*") for n in (origin, x, p)]
        f = (
            frame(*[[a.pos.x, a.pos.y, a.pos.z] for a in aa])
            if all(a is not None for a in aa)
            else None
        )
        result.append(
            Port(terminus=end, residue=ref(chain.name, r), anchors=[origin, x, p], frame=f)
        )
    return result
