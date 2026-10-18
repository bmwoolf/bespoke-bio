"""Standard amino-acid heavy atom inventory; excludes optional terminal OXT."""

SIDECHAINS = dict(
    zip(
        [
            "ALA",
            "ARG",
            "ASN",
            "ASP",
            "CYS",
            "GLN",
            "GLU",
            "GLY",
            "HIS",
            "ILE",
            "LEU",
            "LYS",
            "MET",
            "PHE",
            "PRO",
            "SER",
            "THR",
            "TRP",
            "TYR",
            "VAL",
        ],
        [
            "CB",
            "CB CG CD NE CZ NH1 NH2",
            "CB CG OD1 ND2",
            "CB CG OD1 OD2",
            "CB SG",
            "CB CG CD OE1 NE2",
            "CB CG CD OE1 OE2",
            "",
            "CB CG ND1 CD2 CE1 NE2",
            "CB CG1 CG2 CD1",
            "CB CG CD1 CD2",
            "CB CG CD CE NZ",
            "CB CG SD CE",
            "CB CG CD1 CD2 CE1 CE2 CZ",
            "CB CG CD",
            "CB OG",
            "CB OG1 CG2",
            "CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2",
            "CB CG CD1 CD2 CE1 CE2 CZ OH",
            "CB CG1 CG2",
        ],
    )
)


def missing_atoms(residue):
    if residue.name not in SIDECHAINS:
        return None
    expected = set(("N CA C O " + SIDECHAINS[residue.name]).split())
    return sorted(expected - {a.name for a in residue if not a.element.is_hydrogen})
