"""Explicit publication of exact-sequence target structures for an immutable catalog."""

import argparse
import json
import math
import re
import shutil
import sqlite3
import tempfile
from contextlib import closing
from hashlib import sha256
from pathlib import Path

from .catalog import Catalog

DIGEST = re.compile(r"^[a-f0-9]{64}$")
ARTIFACT = re.compile(r"^[a-f0-9]{64}\.(cif|json|npz)$")


def checksum(path):
    with path.open("rb") as stream:
        digest = sha256()
        while block := stream.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def verify_coordinates(path, sequence):
    # Publication uses the optional Stage 1 runtime; serving remains dependency-light.
    from Bio.PDB import MMCIFParser
    from Bio.SeqUtils import seq1

    structure = MMCIFParser(QUIET=True, auth_chains=False, auth_residues=False).get_structure(
        "target", str(path)
    )
    models = list(structure)
    chains = list(structure.get_chains())
    residues = list(structure.get_residues())
    if (
        len(models) != 1
        or len(chains) != 1
        or "".join(seq1(r.resname) for r in residues) != sequence
        or [r.id[1] for r in residues] != list(range(1, len(sequence) + 1))
        or any(not {"N", "CA", "C", "O"} <= set(r.child_dict) for r in residues)
        or any(not math.isfinite(float(v)) for a in structure.get_atoms() for v in a.coord)
    ):
        raise ValueError("Structure must cover the exact complete sequence in one chain")


def publish_structures(index, source, destination):
    """Publish a reviewed manifest and its files; never infer a sequence or overwrite a release."""
    catalog = Catalog(Path(index))
    source = Path(source).resolve()
    destination = Path(destination).resolve()
    if destination.exists():
        raise FileExistsError(destination)
    supplied = json.loads(source.read_text())
    catalog_hash = catalog.metadata()["manifest_sha256"]
    if supplied["catalog_sha256"] != catalog_hash or supplied["version"] != 1:
        raise ValueError("Structure input does not match this catalog")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=".structures-", dir=destination.parent))
    manifest = {"version": 1, "catalog_sha256": catalog_hash, "structures": {}, "files": {}}
    sequences = {}
    try:
        (temporary / "files").mkdir()

        def copy_file(relative, expected, extension):
            path = (source.parent / relative).resolve()
            if not path.is_relative_to(source.parent) or not DIGEST.fullmatch(expected):
                raise ValueError("Invalid structure source path or checksum")
            if checksum(path) != expected:
                raise ValueError("Structure source checksum mismatch")
            name = expected + extension
            shutil.copyfile(path, temporary / "files" / name)
            manifest["files"][name] = {"sha256": expected, "size": path.stat().st_size}
            return "/api/structures/" + name

        for entry in supplied["structures"]:
            sequence = entry["amino_acid_sequence"]
            digest = sha256(sequence.encode()).hexdigest()
            if digest != entry["sequence_sha256"] or digest in sequences:
                raise ValueError("Invalid or duplicate sequence hash")
            url = copy_file(entry["file"], entry["file_sha256"], ".cif")
            verify_coordinates(temporary / "files" / url.rsplit("/", 1)[1], sequence)
            confidence = entry["confidence"]
            extension = Path(confidence["file"]).suffix
            if extension not in {".npz", ".json"}:
                raise ValueError("Unsupported confidence artifact")
            confidence_url = copy_file(confidence["file"], confidence["sha256"], extension)
            mean = entry.get("mean_plddt")
            if mean is not None and (not math.isfinite(mean) or not 0 <= mean <= 100):
                raise ValueError("pLDDT must use the 0–100 scale")
            flags = entry["provenance"].get("peptide_geometry_screen_flags")
            if flags is not None and (type(flags) is not int or flags < 0):
                raise ValueError("Geometry screen flags must be a nonnegative count or null")
            provenance = json.dumps(entry["provenance"], allow_nan=False).encode()
            provenance_hash = sha256(provenance).hexdigest()
            name = provenance_hash + ".json"
            (temporary / "files" / name).write_bytes(provenance)
            manifest["files"][name] = {"sha256": provenance_hash, "size": len(provenance)}
            manifest["structures"][digest] = {
                "sequence_sha256": digest,
                "source": entry["source"],
                "residues": len(sequence),
                "mean_plddt": mean,
                "geometry_screen_flags": flags,
                "structure_url": url,
                "confidence_url": confidence_url,
                "provenance_url": "/api/structures/" + name,
            }
            sequences[digest] = sequence
        seen = set()
        with (
            closing(sqlite3.connect(temporary / "ready.sqlite3")) as db,
            closing(catalog.connect()) as original,
        ):
            db.execute("CREATE TABLE targets (target_id TEXT PRIMARY KEY, sequence_sha256 TEXT)")
            for row in original.execute("SELECT target_id,payload FROM targets"):
                payload = json.loads(row["payload"])
                digest = payload["sequence_sha256"]
                if digest in sequences:
                    if payload["amino_acid_sequence"] != sequences[digest]:
                        raise ValueError("Catalog sequence disagrees with structure identity")
                    db.execute("INSERT INTO targets VALUES (?,?)", (row["target_id"], digest))
                    seen.add(digest)
            count = db.execute("SELECT count(*) FROM targets").fetchone()[0]
            db.commit()
        if seen != set(sequences):
            raise ValueError("A published structure has no target in this catalog")
        manifest["structured_records"] = count
        manifest["index_sha256"] = checksum(temporary / "ready.sqlite3")
        (temporary / "manifest.json").write_text(json.dumps(manifest, allow_nan=False) + "\n")
        StructureLibrary(temporary, catalog_hash)
        temporary.rename(destination)
        return {"structured_records": count, "structured_sequences": len(sequences)}
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)


class StructureLibrary:
    def __init__(self, root, catalog_hash):
        self.root = Path(root).resolve()
        self.index = self.root / "ready.sqlite3"
        manifest = json.loads((self.root / "manifest.json").read_text())
        if manifest["version"] != 1 or manifest["catalog_sha256"] != catalog_hash:
            raise ValueError("Structure library does not match this catalog")
        if self.index.is_symlink() or checksum(self.index) != manifest["index_sha256"]:
            raise ValueError("Structure index checksum mismatch")
        self.records = manifest["structured_records"]
        self.structures = manifest["structures"]
        self.files = manifest["files"]
        with closing(sqlite3.connect(self.index.as_uri() + "?mode=ro", uri=True)) as db:
            count = db.execute("SELECT count(*) FROM targets").fetchone()[0]
            sequences = {
                row[0] for row in db.execute("SELECT DISTINCT sequence_sha256 FROM targets")
            }
        if count != self.records or sequences != set(self.structures):
            raise ValueError("Structure coverage disagrees with its index")
        for name, entry in self.files.items():
            if not ARTIFACT.fullmatch(name) or name.split(".")[0] != entry["sha256"]:
                raise ValueError("Invalid structure artifact name")
            path = self.root / "files" / name
            if path.is_symlink() or path.stat().st_size != entry["size"]:
                raise ValueError("Structure artifact size or location changed")
            if checksum(path) != entry["sha256"]:
                raise ValueError("Structure artifact checksum mismatch")
        for digest, entry in self.structures.items():
            if not DIGEST.fullmatch(digest) or entry["sequence_sha256"] != digest:
                raise ValueError("Invalid structure sequence identity")
            for field in ("structure_url", "confidence_url", "provenance_url"):
                url = entry[field]
                if (
                    not url.startswith("/api/structures/")
                    or url.removeprefix("/api/structures/") not in self.files
                ):
                    raise ValueError("Unlisted structure artifact")

    def file(self, name):
        if name not in self.files:
            raise KeyError(name)
        return self.root / "files" / name


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--index", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(publish_structures(args.index, args.source, args.output)))


if __name__ == "__main__":
    main()
