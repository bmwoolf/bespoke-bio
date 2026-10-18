import pytest

from bespoke_bio.stage1.workflows import ReferenceAssets, groovy_string, parameters


def test_baseline_preserves_native_and_rna_alignment_evidence(tmp_path):
    assets = ReferenceAssets(
        fasta=tmp_path / "genome.fa",
        dictionary=tmp_path / "genome.dict",
        gtf=tmp_path / "annotation.gtf",
        transcript_fasta=tmp_path / "transcripts.fa",
        bwa_index=tmp_path / "bwa",
        star_index=tmp_path / "star",
        vep_cache=tmp_path / "vep",
        bundle_sha256="a" * 64,
    )
    dna = parameters("dna", tmp_path / "dna.csv", tmp_path / "dna", assets)
    rna = parameters("rna", tmp_path / "rna.csv", tmp_path / "rna", assets)
    assert dna["split_fastq"] == 0
    assert dna["skip_tools"] == "baserecalibrator"
    assert dna["save_output_as_bam"] and dna["save_mapped"]
    assert "--filter_common" not in dna["vep_custom_args"]
    assert "--per_gene" not in dna["vep_custom_args"]
    assert "--allele_number" in dna["vep_custom_args"]
    assert rna["save_align_intermeds"]
    assert rna["extra_salmon_quant_args"] == ""
    assert not rna["skip_trimming"]


def test_generated_config_cannot_inject_groovy_expressions():
    assert groovy_string("${value}") == "'${value}'"
    assert groovy_string("a'b") == "'a\\'b'"
    with pytest.raises(ValueError):
        groovy_string("x\nprocess.container='bad'")


def test_launch_rejects_incomplete_assets_and_wrong_run_identity(tmp_path):
    from bespoke_bio.stage1.io import canonical
    from bespoke_bio.stage1.runtime import Runtime
    from bespoke_bio.stage1.workflows import validate_prepared

    runtime = Runtime(
        root=tmp_path,
        scratch=tmp_path / "work",
        nextflow=tmp_path / "nextflow",
        nxf_home=tmp_path / "nxf",
        sarek=tmp_path / "sarek",
        rnaseq=tmp_path / "rnaseq",
        dna={"cpus": 2, "memory_gib": 4, "max_hours": 1},
        rna={"cpus": 2, "memory_gib": 4, "max_hours": 1},
        host_cpus=4,
        host_memory_gib=8,
        disks=[{"path": tmp_path, "min_free_bytes": 0, "min_free_inodes": 0}],
    )
    assets = ReferenceAssets(
        fasta=tmp_path / "genome.fa",
        dictionary=tmp_path / "genome.dict",
        gtf=tmp_path / "annotation.gtf",
        transcript_fasta=tmp_path / "transcripts.fa",
        bwa_index=tmp_path / "bwa",
        star_index=tmp_path / "star",
        vep_cache=tmp_path / "vep",
        bundle_sha256="a" * 64,
    )
    (tmp_path / "runtime.json").write_bytes(canonical(runtime.model_dump(mode="json")))
    (tmp_path / "reference-assets.json").write_bytes(canonical(assets.model_dump(mode="json")))
    (tmp_path / "run.json").write_bytes(canonical({"run_id": "prepared"}))
    with pytest.raises(ValueError, match="container labels"):
        validate_prepared(tmp_path, runtime, "different")
    with pytest.raises(ValueError, match="BWA-MEM2 index"):
        validate_prepared(tmp_path, runtime, "prepared")
    changed = runtime.model_copy(update={"host_cpus": 8})
    with pytest.raises(ValueError, match="budgets"):
        validate_prepared(tmp_path, changed, "prepared")
