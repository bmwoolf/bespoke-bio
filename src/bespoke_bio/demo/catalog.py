"""Index a complete, verified Stage 1 package without changing biological records."""

import json
import os
import sqlite3
import tempfile
from pathlib import Path


def build_index(package: Path, destination: Path) -> dict:
    # Import the optional Stage 1 runtime only when building an index.
    from bespoke_bio.stage1.io import file_hash
    from bespoke_bio.stage1.proteins import validate_protein_package

    package = package.resolve()
    manifest = package / "manifest.json"
    manifest_hash = file_hash(manifest)
    source = package / "protein-targets.jsonl"
    source_hash = file_hash(source)
    verified = validate_protein_package(package)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=".catalog-", dir=destination.parent)
    os.close(fd)
    try:
        with sqlite3.connect(temporary) as db:
            db.executescript("""
                CREATE TABLE targets (
                    target_id TEXT PRIMARY KEY, gene TEXT NOT NULL,
                    transcript TEXT NOT NULL, donor TEXT NOT NULL,
                    changed INTEGER NOT NULL, length INTEGER NOT NULL, payload TEXT NOT NULL);
                CREATE TABLE metadata (payload TEXT NOT NULL);
            """)
            with source.open() as stream:
                for line in stream:
                    record = json.loads(line)
                    db.execute(
                        "INSERT INTO targets VALUES (?,?,?,?,?,?,?)",
                        (
                            record["target_id"],
                            record["gene_symbol"],
                            record["transcript_id"],
                            record["donor"],
                            record["sequence_sha256"] != record["reference_sequence_sha256"],
                            len(record["amino_acid_sequence"]),
                            line.strip(),
                        ),
                    )
            count, changed = db.execute("SELECT count(*),sum(changed) FROM targets").fetchone()
            if count != verified["target_records"]:
                raise ValueError("Index lost target records")
            if file_hash(source) != source_hash or file_hash(manifest) != manifest_hash:
                raise ValueError("Source package changed while indexing")
            metadata = {
                **verified,
                "nonreference_records": changed,
                "manifest_sha256": manifest_hash,
                "targets_sha256": source_hash,
                "schema_version": "protein-catalog-1",
            }
            db.execute("INSERT INTO metadata VALUES (?)", (json.dumps(metadata),))
            db.executescript("""
                CREATE INDEX gene_index ON targets(gene COLLATE NOCASE);
                CREATE INDEX transcript_index ON targets(transcript COLLATE NOCASE);
                CREATE INDEX browse_index ON targets(gene,transcript,target_id,donor,changed,length);
                CREATE INDEX changed_index ON targets(changed,gene,transcript,target_id,donor,length);
            """)
        # Publish only a complete index; never replace a prior one implicitly.
        os.link(temporary, destination)
        return metadata
    finally:
        Path(temporary).unlink(missing_ok=True)


class Catalog:
    def __init__(self, path: Path):
        self.path = path.resolve()

    def connect(self):
        db = sqlite3.connect(self.path.as_uri() + "?mode=ro", uri=True)
        db.row_factory = sqlite3.Row
        return db

    def metadata(self):
        with self.connect() as db:
            return json.loads(db.execute("SELECT payload FROM metadata").fetchone()[0])

    def search(self, query="", changed=False, limit=25, offset=0):
        # Literal prefix search: user text cannot become SQL or LIKE wildcards.
        pattern = query.replace("\\", "\\\\").replace("%", "\\%").replace("_", "\\_") + "%"
        where = "(gene LIKE ? ESCAPE '\\' OR transcript LIKE ? ESCAPE '\\' OR target_id = ?)"
        args = [pattern, pattern, query]
        if changed:
            where += " AND changed=1"
        with self.connect() as db:
            if query:
                total = db.execute("SELECT count(*) FROM targets WHERE " + where, args).fetchone()[0]
            else:
                meta = json.loads(db.execute("SELECT payload FROM metadata").fetchone()[0])
                total = meta["nonreference_records" if changed else "target_records"]
            rows = db.execute(
                "SELECT target_id,gene,transcript,donor,changed,length FROM targets WHERE "
                + where
                + " ORDER BY gene,transcript,target_id LIMIT ? OFFSET ?",
                [*args, limit, offset],
            )
            return {
                "total": total,
                "offset": offset,
                "limit": limit,
                "records": [dict(row) for row in rows],
            }

    def get(self, target_id):
        with self.connect() as db:
            row = db.execute(
                "SELECT payload FROM targets WHERE target_id=?", (target_id,)
            ).fetchone()
        if row is None:
            raise KeyError(target_id)
        return json.loads(row[0])
