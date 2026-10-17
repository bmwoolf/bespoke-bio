import argparse
import json
import shutil
import subprocess
import tempfile
from datetime import UTC, datetime
from hashlib import sha256
from pathlib import Path
from typing import Annotated, Any, Literal
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, ConfigDict, Field, field_validator

# API schema

Identifier = Annotated[str, Field(min_length=1, max_length=200)]


def normalize_dna(value: str) -> str:
    value = value.upper().replace("U", "T")
    invalid = sorted(set(value) - set("ACGTN"))
    if invalid:
        raise ValueError(f"unsupported nucleotide symbols: {', '.join(invalid)}")
    return value


class Model(BaseModel):
    model_config = ConfigDict(extra="forbid")


class Reference(Model):
    genome_build: Identifier
    annotation: Identifier


class Expression(Model):
    tpm: float = Field(ge=0)
    raw_count: float | None = Field(default=None, ge=0)


class Variant(Model):
    variant_id: Identifier
    cds_position: int = Field(ge=1, description="One-based position in the reference CDS")
    ref: str = Field(min_length=1)
    alt: str = ""
    genotype: str | None = None
    phase_set: str | None = None
    quality: str | None = None

    @field_validator("ref", "alt")
    @classmethod
    def validate_allele(cls, value: str) -> str:
        return normalize_dna(value)


class Transcript(Model):
    gene_id: Identifier
    gene_symbol: Identifier
    transcript_id: Identifier
    uniprot_id: Identifier | None = None
    haplotype: Identifier = "unknown"
    phase_scenario: Identifier | None = None
    coding_sequence: str = Field(min_length=3)
    expression: Expression
    variants: list[Variant] = Field(default_factory=list)

    @field_validator("coding_sequence")
    @classmethod
    def validate_sequence(cls, value: str) -> str:
        return normalize_dna(value)


class TranslateRequest(Model):
    case_id: Identifier
    reference: Reference
    transcripts: list[Transcript] = Field(min_length=1)


class Structure(Model):
    status: Literal["complete", "unresolved"]
    source: str | None = None
    source_id: str | None = None
    file: str | None = None
    mean_plddt: float | None = None
    confidence_file: str | None = None
    model_version: str | None = None
    error: str | None = None


class Protein(Model):
    target_id: Identifier
    gene_id: Identifier
    gene_symbol: Identifier
    transcript_id: Identifier
    uniprot_id: Identifier | None = None
    haplotype: Identifier
    phase_scenario: str | None = None
    source_variant_ids: list[str]
    expression: Expression
    reference_amino_acid_sequence: str
    amino_acid_sequence: str
    amino_acid_length: int
    sequence_hash: str
    first_changed_residue: int | None
    translation_status: Literal["complete", "complete_with_warnings"]
    warnings: list[str]
    structure: Structure | None = None


class TranslateResponse(Model):
    case_id: Identifier
    reference: Reference
    targets: list[Protein]
    target_count: int
    unique_sequence_count: int


# Translation

# Standard genetic code in T, C, A, G base order.
_CODONS = (first + second + third for first in "TCAG" for second in "TCAG" for third in "TCAG")
_AMINO_ACIDS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = dict(zip(_CODONS, _AMINO_ACIDS, strict=True))


class VariantError(ValueError):
    pass


def apply_variants(reference: str, variants: list[Variant]) -> str:
    """Apply non-overlapping variants using positions from the reference CDS."""
    pieces: list[str] = []
    cursor = 0

    for variant in sorted(variants, key=lambda item: item.cds_position):
        start = variant.cds_position - 1
        end = start + len(variant.ref)
        if start < cursor:
            raise VariantError(f"{variant.variant_id} overlaps an earlier variant")
        if end > len(reference):
            raise VariantError(f"{variant.variant_id} extends beyond the reference CDS")
        observed = reference[start:end]
        if observed != variant.ref:
            raise VariantError(
                f"{variant.variant_id} expected {variant.ref} at CDS position "
                f"{variant.cds_position}, found {observed}"
            )
        pieces.extend((reference[cursor:start], variant.alt))
        cursor = end

    pieces.append(reference[cursor:])
    return "".join(pieces)


def translate_cds(cds: str) -> tuple[str, list[str]]:
    warnings: list[str] = []
    if not cds.startswith("ATG"):
        warnings.append("missing_start")
    if len(cds) % 3:
        warnings.append("partial_codon")

    protein: list[str] = []
    codon_count = len(cds) // 3
    for index in range(codon_count):
        amino_acid = CODON_TABLE.get(cds[index * 3 : index * 3 + 3], "X")
        if amino_acid == "X":
            warnings.append("ambiguous_codon")
        if amino_acid == "*":
            if index < codon_count - 1:
                warnings.append("internal_stop")
            break
        protein.append(amino_acid)

    return "".join(protein), warnings


def sequence_hash(sequence: str) -> str:
    return f"sha256:{sha256(sequence.encode()).hexdigest()}"


def first_change(reference: str, alternate: str) -> int | None:
    for index, (left, right) in enumerate(zip(reference, alternate), start=1):
        if left != right:
            return index
    return min(len(reference), len(alternate)) + 1 if len(reference) != len(alternate) else None


def translate_case(request: TranslateRequest) -> TranslateResponse:
    targets: list[Protein] = []

    for transcript in request.transcripts:
        if (
            transcript.haplotype == "unknown"
            and transcript.variants
            and not transcript.phase_scenario
        ):
            raise VariantError(
                f"{transcript.transcript_id} has variants with unknown phase; "
                "Step 1 must supply an explicit phase_scenario"
            )

        reference_protein, _ = translate_cds(transcript.coding_sequence)
        variants = sorted(
            transcript.variants,
            key=lambda item: (item.cds_position, item.variant_id),
        )
        protein, warnings = translate_cds(apply_variants(transcript.coding_sequence, variants))
        protein_hash = sequence_hash(protein)
        variant_ids = [variant.variant_id for variant in variants]
        identity = "|".join(
            (
                request.case_id,
                transcript.transcript_id,
                transcript.phase_scenario or transcript.haplotype,
                ",".join(variant_ids),
                protein_hash,
            )
        )
        targets.append(
            Protein(
                target_id=f"target_{sha256(identity.encode()).hexdigest()[:12]}",
                gene_id=transcript.gene_id,
                gene_symbol=transcript.gene_symbol,
                transcript_id=transcript.transcript_id,
                uniprot_id=transcript.uniprot_id,
                haplotype=transcript.haplotype,
                phase_scenario=transcript.phase_scenario,
                source_variant_ids=variant_ids,
                expression=transcript.expression,
                reference_amino_acid_sequence=reference_protein,
                amino_acid_sequence=protein,
                amino_acid_length=len(protein),
                sequence_hash=protein_hash,
                first_changed_residue=first_change(reference_protein, protein),
                translation_status="complete_with_warnings" if warnings else "complete",
                warnings=warnings,
            )
        )

    return TranslateResponse(
        case_id=request.case_id,
        reference=request.reference,
        targets=targets,
        target_count=len(targets),
        unique_sequence_count=len({target.sequence_hash for target in targets}),
    )


# Structure resolution

USER_AGENT = "bespoke-bio/0.1"
RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_DATA_URL = "https://data.rcsb.org/rest/v1/core/polymer_entity"
ALPHAFOLD_URL = "https://alphafold.ebi.ac.uk/api/prediction"


def fetch_json(url: str, payload: dict[str, Any] | None = None) -> Any:
    request = Request(
        url,
        data=json.dumps(payload).encode() if payload else None,
        headers={"Content-Type": "application/json", "User-Agent": USER_AGENT},
    )
    with urlopen(request, timeout=30) as response:
        return json.load(response)


def download(url: str, destination: Path) -> None:
    request = Request(url, headers={"User-Agent": USER_AGENT})
    with urlopen(request, timeout=60) as response:
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(response.read())


def structure_key(protein: Protein) -> str:
    return protein.sequence_hash.removeprefix("sha256:")


def find_rcsb(protein: Protein, output_dir: Path) -> Structure | None:
    query = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": 1,
                "identity_cutoff": 1,
                "sequence_type": "protein",
                "value": protein.amino_acid_sequence,
            },
        },
        "request_options": {
            "scoring_strategy": "sequence",
            "paginate": {"start": 0, "rows": 10},
        },
        "return_type": "polymer_entity",
    }
    for result in fetch_json(RCSB_SEARCH_URL, query).get("result_set", []):
        identifier = result["identifier"]
        entry_id, entity_id = identifier.split("_", 1)
        if len(entry_id) != 4:  # Exclude computed structure models.
            continue
        entity = fetch_json(f"{RCSB_DATA_URL}/{entry_id}/{entity_id}")
        observed = "".join(entity["entity_poly"]["pdbx_seq_one_letter_code_can"].split())
        if observed != protein.amino_acid_sequence:
            continue

        relative_file = Path("targets") / f"{structure_key(protein)}.cif"
        download(f"https://files.rcsb.org/download/{entry_id}.cif", output_dir / relative_file)
        return Structure(
            status="complete",
            source="RCSB PDB",
            source_id=identifier,
            file=relative_file.as_posix(),
        )
    return None


def find_alphafold(protein: Protein, output_dir: Path) -> Structure | None:
    if not protein.uniprot_id:
        return None

    records = fetch_json(f"{ALPHAFOLD_URL}/{protein.uniprot_id}")
    record = next(
        (item for item in records if item.get("sequence") == protein.amino_acid_sequence),
        None,
    )
    if not record:
        return None

    relative_file = Path("targets") / f"{structure_key(protein)}.cif"
    download(record["cifUrl"], output_dir / relative_file)
    confidence_file = None
    if record.get("paeDocUrl"):
        relative_confidence = Path("targets") / f"{structure_key(protein)}.pae.json"
        download(record["paeDocUrl"], output_dir / relative_confidence)
        confidence_file = relative_confidence.as_posix()

    return Structure(
        status="complete",
        source="AlphaFold DB",
        source_id=record["entryId"],
        file=relative_file.as_posix(),
        mean_plddt=record.get("globalMetricValue"),
        confidence_file=confidence_file,
        model_version=str(record.get("latestVersion")) if record.get("latestVersion") else None,
    )


def resolve_structure(protein: Protein, output_dir: Path) -> Structure:
    errors: list[str] = []
    for source, resolver in (("RCSB PDB", find_rcsb), ("AlphaFold DB", find_alphafold)):
        try:
            structure = resolver(protein, output_dir)
            if structure:
                return structure
        except (HTTPError, URLError, TimeoutError, KeyError, ValueError) as error:
            errors.append(f"{source}: {error}")
    return Structure(status="unresolved", error="; ".join(errors) or None)


def run_colabfold(proteins: list[Protein], output_dir: Path) -> dict[str, Structure]:
    binary = shutil.which("colabfold_batch")
    if not binary:
        raise RuntimeError("structure folding requires colabfold_batch on PATH")

    results: dict[str, Structure] = {}
    with tempfile.TemporaryDirectory(prefix="bespoke-colabfold-", dir=output_dir) as temporary:
        work_dir = Path(temporary)
        fasta = work_dir / "input.fasta"
        jobs = {f"seq_{structure_key(protein)[:16]}": protein for protein in proteins}
        fasta.write_text(
            "".join(f">{job}\n{protein.amino_acid_sequence}\n" for job, protein in jobs.items())
        )
        predictions = work_dir / "predictions"
        subprocess.run([binary, str(fasta), str(predictions)], check=True)

        for job, protein in jobs.items():
            pdb_files = sorted(predictions.glob(f"{job}*rank_001*.pdb"))
            if not pdb_files:
                pdb_files = sorted(predictions.glob(f"{job}*.pdb"))
            if not pdb_files:
                results[protein.sequence_hash] = Structure(
                    status="unresolved", error="ColabFold produced no PDB file"
                )
                continue

            relative_file = Path("targets") / f"{structure_key(protein)}.pdb"
            shutil.copyfile(pdb_files[0], output_dir / relative_file)
            score_files = sorted(predictions.glob(f"{job}*scores*rank_001*.json"))
            confidence_file = None
            mean_plddt = None
            if score_files:
                scores = json.loads(score_files[0].read_text())
                values = scores.get("plddt", [])
                mean_plddt = sum(values) / len(values) if values else None
                relative_confidence = Path("targets") / f"{structure_key(protein)}.scores.json"
                shutil.copyfile(score_files[0], output_dir / relative_confidence)
                confidence_file = relative_confidence.as_posix()

            results[protein.sequence_hash] = Structure(
                status="complete",
                source="ColabFold",
                source_id=job,
                file=relative_file.as_posix(),
                mean_plddt=mean_plddt,
                confidence_file=confidence_file,
            )
    return results


# Step 3 handoff package


def build_bundle(
    request: TranslateRequest,
    output_dir: Path,
) -> TranslateResponse:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "targets").mkdir(exist_ok=True)
    response = translate_case(request)
    unique = {target.sequence_hash: target for target in response.targets}

    structures = {
        sequence_hash: resolve_structure(protein, output_dir)
        for sequence_hash, protein in unique.items()
    }
    missing = [protein for key, protein in unique.items() if structures[key].status == "unresolved"]
    if missing:
        structures.update(run_colabfold(missing, output_dir))
    unresolved = [key for key, structure in structures.items() if structure.status == "unresolved"]
    if unresolved:
        raise RuntimeError(f"folding failed for {len(unresolved)} unique sequence(s)")

    for target in response.targets:
        target.structure = structures[target.sequence_hash]

    target_ids: dict[str, list[str]] = {}
    for target in response.targets:
        target_ids.setdefault(target.sequence_hash, []).append(target.target_id)
    (output_dir / "proteins.fasta").write_text(
        "".join(
            f">{sequence_hash.removeprefix('sha256:')} targets={','.join(target_ids[sequence_hash])}\n"
            f"{unique[sequence_hash].amino_acid_sequence}\n"
            for sequence_hash in sorted(unique)
        )
    )
    (output_dir / "targets.jsonl").write_text(
        "".join(
            json.dumps(target.model_dump(mode="json", exclude_none=True), separators=(",", ":"))
            + "\n"
            for target in response.targets
        )
    )
    provenance = {
        "case_id": request.case_id,
        "reference": request.reference.model_dump(mode="json"),
        "generated_at": datetime.now(UTC).isoformat(),
        "lookup_order": ["RCSB PDB", "AlphaFold DB"],
        "folding_software": "colabfold_batch",
    }
    (output_dir / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    return response


# HTTP API and CLI

app = FastAPI(
    title="Bespoke Bio",
    summary="Translate expressed rare variants into target proteins",
)


@app.get("/healthz")
def healthcheck() -> dict[str, str]:
    return {"status": "ok"}


@app.post("/v1/translate", response_model=TranslateResponse, response_model_exclude_none=True)
def translate(request: TranslateRequest) -> TranslateResponse:
    try:
        return translate_case(request)
    except VariantError as error:
        raise HTTPException(status_code=422, detail=str(error)) from error


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a Bespoke Bio Step 3 handoff package")
    parser.add_argument("input", type=Path, help="Step 1 JSON input")
    parser.add_argument("output", type=Path, help="Output directory")
    arguments = parser.parse_args()
    request = TranslateRequest.model_validate_json(arguments.input.read_text())
    response = build_bundle(request, arguments.output)
    print(
        json.dumps(
            {
                "targets": response.target_count,
                "unique_sequences": response.unique_sequence_count,
                "output": str(arguments.output),
            }
        )
    )


if __name__ == "__main__":
    main()
