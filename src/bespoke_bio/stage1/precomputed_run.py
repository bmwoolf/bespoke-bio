"""Run the prepared-evidence route with pinned native containers and explicit resources."""

import os
import re
from pathlib import Path

from .contracts import Record
from .io import canonical, verify, write_new
from .precomputed import PrecomputedSource, coding_regions, qualify_source
from .proteins import export_proteins
from .runtime import Allocation, DiskLimit, Events, disk_status, execute


class PrecomputedRun(Record):
    """Operator-supplied paths and container digests; no provider-specific settings."""

    source: PrecomputedSource
    genome: Path
    annotation: Path
    reference_records: Path
    assembly: str
    bcftools_image: str
    vep_image: str
    allocation: Allocation
    disks: list[DiskLimit]

    def validate_assets(self):
        for image in (self.bcftools_image, self.vep_image):
            if not re.fullmatch(r"[^\s]+@sha256:[0-9a-f]{64}", image):
                raise ValueError("native containers require immutable image digests")
        for path in (self.genome, self.annotation, self.reference_records):
            if not path.is_file():
                raise ValueError(f"required input missing: {path}")
        if not Path(str(self.genome) + ".fai").is_file():
            raise ValueError("indexed reference missing")
        if not self.disks:
            raise ValueError("disk headroom monitoring is required")
        disk_status(self.disks)


def process_precomputed(settings: PrecomputedRun, output: Path) -> dict:
    """Normalize supported source calls, annotate coding context and export all resolved alleles.

    Each attempt uses a fresh directory. Existing source data and indexes are immutable;
    native logs and intermediate evidence remain available when an attempt fails.
    """
    settings.validate_assets()
    source = settings.source
    for item in (
        source.variants,
        source.variant_index,
        source.benchmark_regions,
        source.source_readme,
    ):
        verify(item)
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    write_new(output / "run.json", canonical(settings.model_dump(mode="json")) + b"\n")
    run_id = "prepared-" + os.urandom(12).hex()
    events = Events(output / "events.jsonl", run_id, "initial")
    events.write("source-qualification", "started")
    counts = qualify_source(source, output / "qualified")
    events.write("source-qualification", "completed", counts=counts)
    regions = coding_regions(settings.annotation, settings.genome, 50)
    write_new(
        output / "annotation-regions.bed",
        "".join(
            f"{chrom}\t{a}\t{b}\n" for chrom, spans in regions.items() for a, b in spans
        ).encode(),
    )
    common = [
        "docker",
        "run",
        "--rm",
        "--network",
        "none",
        "--label",
        f"bespoke.run={run_id}",
        "--user",
        f"{os.getuid()}:{os.getgid()}",
        "--cpus",
        str(settings.allocation.cpus),
        "--memory",
        f"{settings.allocation.memory_gib}g",
        "-v",
        f"{output}:/output",
        "-v",
        f"{settings.genome.resolve().parent}:/reference:ro",
        "-v",
        f"{settings.annotation.resolve()}:/annotations/source.gtf:ro",
    ]
    reference = "/reference/" + settings.genome.name
    threads = str(max(0, settings.allocation.cpus - 1))
    # Build the indexed annotation from the exact GTF used for CDS extraction.
    # VEP enables offline mode itself for local GTF+FASTA input; explicitly adding
    # --offline also enables its separate Ensembl cache, which is not this source.
    annotation_script = """set -euo pipefail
awk '!/^#/' /annotations/source.gtf | LC_ALL=C sort -T /output -S 4G --parallel=4 -k1,1 -k4,4n -k5,5n | bgzip -@ 3 -c > /output/annotation.gtf.gz
tabix -p gff /output/annotation.gtf.gz
"""
    execute(
        common + [settings.bcftools_image, "bash", "-c", annotation_script],
        cwd=output,
        log=output / "annotation-index.log",
        events=events,
        stage="annotation-index",
        disks=settings.disks,
        timeout_seconds=settings.allocation.max_hours * 3600,
        poll_seconds=5,
        docker_cleanup=True,
    )
    native = [
        (
            "normalize",
            [
                "norm",
                "--check-ref",
                "e",
                "--old-rec-tag",
                "ORIG",
                "-f",
                reference,
                "--threads",
                threads,
                "-Ob",
                "-o",
                "/output/normalized.bcf",
                "/output/qualified/qualified.bcf",
            ],
        ),
        ("index", ["index", "--threads", threads, "/output/normalized.bcf"]),
        (
            "select",
            [
                "view",
                "--regions-file",
                "/output/annotation-regions.bed",
                "--regions-overlap",
                "1",
                "--threads",
                threads,
                "-Ob",
                "-o",
                "/output/coding-context.bcf",
                "/output/normalized.bcf",
            ],
        ),
        (
            "vep-input",
            ["view", "-Ov", "-o", "/output/coding-context.vcf", "/output/coding-context.bcf"],
        ),
    ]
    for name, arguments in native:
        execute(
            common + [settings.bcftools_image, "bcftools"] + arguments,
            cwd=output,
            log=output / f"{name}.log",
            events=events,
            stage=name,
            disks=settings.disks,
            timeout_seconds=settings.allocation.max_hours * 3600,
            poll_seconds=5,
            docker_cleanup=True,
        )
    arguments = [
        "vep",
        "--gtf",
        "/output/annotation.gtf.gz",
        "--assembly",
        settings.assembly,
        "--species",
        "homo_sapiens",
        "--fasta",
        reference,
        "--input_file",
        "/output/coding-context.vcf",
        "--format",
        "vcf",
        "--output_file",
        "/output/annotated.vcf",
        "--vcf",
        "--fork",
        str(settings.allocation.cpus),
        "--buffer_size",
        "1000",
        "--transcript_version",
        "--allele_number",
        "--check_ref",
        "--dont_skip",
        "--symbol",
        "--biotype",
        "--max_sv_size",
        "-1",
        "--stats_file",
        "/output/vep-summary.html",
        "--warning_file",
        "/output/vep-warnings.txt",
        "--skipped_variants_file",
        "/output/vep-skipped.txt",
    ]
    execute(
        common + [settings.vep_image] + arguments,
        cwd=output,
        log=output / "vep.log",
        events=events,
        stage="consequences",
        disks=settings.disks,
        timeout_seconds=settings.allocation.max_hours * 3600,
        poll_seconds=5,
        docker_cleanup=True,
    )
    events.write("protein-export", "started")
    result = export_proteins(
        source,
        settings.reference_records,
        output / "coding-context.bcf",
        output / "annotated.vcf",
        output / "qualified/effective-confidence.bed",
        output / "proteins",
    )
    events.write("protein-export", "completed", counts=result["counts"])
    return result
