"""Exercise indexed native VCF joins and publication with a small explicit DNA/RNA fixture."""

from pathlib import Path

import pysam

from bespoke_bio.api import Reference, TranslateRequest
from bespoke_bio.stage1.assembly import AssemblyInputs, assemble
from bespoke_bio.stage1.contracts import RELEASE_CHECKS, Check
from bespoke_bio.stage1.intake import pilot_intake
from bespoke_bio.stage1.io import artifact, canonical, read_json, write_jsonl
from bespoke_bio.stage1.references import CodingTranscript, Segment


def vcf(path, *, annotation=False, reference_blocks=False):
    header = pysam.VariantHeader()
    header.contigs.add("chr1", length=100)
    for name, number, kind in [
        ("GT", 1, "String"),
        ("DP", 1, "Integer"),
        ("GQ", 1, "Integer"),
        ("MIN_DP", 1, "Integer"),
    ]:
        header.formats.add(name, number, kind, name)
    header.add_sample("DNA")
    if annotation:
        header.info.add(
            "CSQ", ".", "String", "Format: Allele|Consequence|Feature_type|Feature|ALLELE_NUM"
        )
    with pysam.VariantFile(str(path), "wz", header=header) as out:
        if reference_blocks:
            for start, stop in [(10, 13), (14, 19)]:
                row = out.new_record(
                    contig="chr1", start=start, stop=stop, alleles=("A", "<NON_REF>")
                )
                row.samples["DNA"].update({"GT": (0, 0), "DP": 20, "GQ": 30, "MIN_DP": 20})
                out.write(row)
        else:
            row = out.new_record(contig="chr1", start=13, alleles=("G", "C"))
            row.filter.add("PASS")
            row.samples["DNA"].update({"GT": (0, 1), "DP": 20, "GQ": 30})
            if annotation:
                row.info["CSQ"] = ("C|missense_variant|Transcript|ENST.1|1",)
            out.write(row)
    pysam.tabix_index(str(path), preset="vcf")


def test_native_evidence_produces_uncapped_zero_tpm_alleles(tmp_path):
    inputs = {}
    registered = []
    for role in ["original_vcf", "normalized_vcf", "annotated_vcf", "phased_vcf", "gvcf"]:
        path = tmp_path / (role + ".vcf.gz")
        vcf(path, annotation=role == "annotated_vcf", reference_blocks=role == "gvcf")
        inputs[role] = path
        registered.extend(
            [
                artifact(path, role, "vcf.gz"),
                artifact(Path(str(path) + ".tbi"), role + "_index", "tbi"),
            ]
        )
    tx = CodingTranscript(
        "ENST.1",
        "ENSG.1",
        "G",
        "protein_coding",
        "chr1",
        "+",
        segments=[Segment(10, 19, "exon")],
        sequence="ATGGAATAA",
    )
    path = tmp_path / "references.jsonl"
    write_jsonl(path, [tx.record()])
    inputs["reference_records"] = path
    contents = {
        "cds_bed": "chr1\t10\t19\n",
        "depth": "".join(f"chr1\t{pos}\t20\n" for pos in range(11, 20)),
        "expression": "Name\tLength\tEffectiveLength\tTPM\tNumReads\nENST.1\t9\t5\t0\t0\n",
        "tools": "{}\n",
        "policy": "{}\n",
    }
    for role, content in contents.items():
        path = tmp_path / role
        path.write_text(content)
        inputs[role] = path
    sources = []
    for role in ("dna_r1", "dna_r2", "rna_r1", "rna_r2"):
        path = tmp_path / role
        path.write_text("@read\nATG\n+\nIII\n")
        sources.append(artifact(path, role, "fastq"))
    lock = {
        "dataset": {
            "donor": "donor",
            "sex": "XY",
            "dna_specimen": "DNA",
            "rna_specimen": "RNA",
            "dna_library": "dna-lib",
            "rna_library": "rna-lib",
        }
    }
    path = tmp_path / "intake.json"
    path.write_bytes(canonical(pilot_intake(lock, sources).model_dump(mode="json")))
    inputs["intake"] = path
    for role, path in inputs.items():
        if role not in {a.role for a in registered}:
            registered.append(artifact(path, role, path.suffix))
    registry = tmp_path / "registry.json"
    registry.write_bytes(canonical([a.model_dump(mode="json") for a in registered]))
    owned = {
        "callability",
        "variant_reconciliation",
        "handoff_translation",
        "disposition_accounting",
        "artifact_integrity",
    }
    checks = [
        Check(
            name=name,
            status="passed",
            command=["synthetic-test-fixture"],
            evidence_roles=["intake"],
        )
        for name in sorted(RELEASE_CHECKS - owned)
    ]
    check_file = tmp_path / "checks.json"
    check_file.write_bytes(canonical([c.model_dump(mode="json") for c in checks]))
    settings = AssemblyInputs(
        **inputs,
        sample="DNA",
        reference=Reference(genome_build="GRCh38", annotation="fixture"),
        reference_sha256="a" * 64,
        evidence_registry=registry,
        release_checks=check_file,
        durable_roots=[tmp_path],
    )
    output = tmp_path / "handoff"
    manifest = assemble(settings, Path(__file__).resolve().parents[2], output, "attempt-1")
    assert manifest.status == "ready"
    assert manifest.target_count == 2
    assert manifest.unique_sequence_count == 2
    request = TranslateRequest.model_validate(read_json(output / "attempt-1/case.json"))
    assert all(record.expression.tpm == 0 for record in request.transcripts)
    assert [len(record.variants) for record in request.transcripts] == [0, 1]
    assert request.transcripts[1].variants[0].cds_position == 4
    assert read_json(output / "ready.json")["target_count"] == 2
