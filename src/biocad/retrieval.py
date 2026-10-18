import csv
import subprocess
from pathlib import Path

import gemmi

from .geometry import apply, fit
from .structures import protein_residues

FIELDS = (
    "query,target,qstart,qend,tstart,tend,qaln,taln,qcov,tcov,alntmscore,qtmscore,ttmscore,evalue"
)


def search(query, library, output, executable="foldseek", top_k=5):
    if not 1 <= top_k <= 25:
        raise ValueError("top_k must be 1..25")
    library = Path(library)
    files = sorted(library.glob("*.cif"))
    if not files or len(files) > 100:
        raise ValueError("Declare a collection of 1..100 mmCIF files")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    version = subprocess.check_output([executable, "version"], text=True).strip()
    commands = [
        [executable, "createdb", str(library), str(output / "db"), "--chain-name-mode", "1"],
        [
            executable,
            "easy-search",
            str(query),
            str(output / "db"),
            str(output / "hits.tsv"),
            str(output / "tmp"),
            "--format-output",
            FIELDS,
            "--alignment-type",
            "1",
            "--exhaustive-search",
            "1",
            "--max-seqs",
            str(top_k),
            "--threads",
            "1",
            "-e",
            "10",
        ],
    ]
    with (output / "foldseek.log").open("w") as log:
        for cmd in commands:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
    with (output / "hits.tsv").open() as f:
        hits = list(csv.DictReader(f, fieldnames=FIELDS.split(","), delimiter="\t"))
    hits = sorted(hits, key=lambda h: float(h["alntmscore"]), reverse=True)[:top_k]
    return {
        "version": version,
        "commands": commands,
        "mode": "TMalign; exhaustive bounded library",
        "score_semantics": "q/t TM scores normalized by respective chain lengths; evalue retained raw, not ranked",
        "hits": hits,
    }


def aligned_pairs(hit, query_chain, template_chain):
    qr, tr = protein_residues(query_chain), protein_residues(template_chain)
    qi, ti = int(hit["qstart"]) - 1, int(hit["tstart"]) - 1
    pairs = []
    for q, t in zip(hit["qaln"], hit["taln"], strict=True):
        if q != "-" and t != "-":
            qa, ta = qr[qi].find_atom("CA", "*"), tr[ti].find_atom("CA", "*")
            if qa is not None and ta is not None:
                pairs.append(
                    (qi, ti, [qa.pos.x, qa.pos.y, qa.pos.z], [ta.pos.x, ta.pos.y, ta.pos.z])
                )
        qi += q != "-"
        ti += t != "-"
    return pairs


def transfer(hit, query_chain, template_chain, partner_chain):
    pairs = aligned_pairs(hit, query_chain, template_chain)
    transform, rmsd = fit([p[3] for p in pairs], [p[2] for p in pairs])
    partner = move(partner_chain, transform)
    return partner, {
        "transform": transform.model_dump(),
        "ca_rmsd_angstrom": rmsd,
        "residue_pairs": [
            {
                "query_index": p[0],
                "template_index": p[1],
                "query_auth": str(protein_residues(query_chain)[p[0]].seqid),
                "template_auth": str(protein_residues(template_chain)[p[1]].seqid),
            }
            for p in pairs
        ],
        "evidence": "transferred_analogy",
        "direction": "template target -> query target; same transform applied to partner",
    }


def move(chain, transform):
    result = chain.clone()
    for r in result:
        for a in r:
            p = apply([[a.pos.x, a.pos.y, a.pos.z]], transform)[0]
            a.pos = gemmi.Position(*p)
    return result
