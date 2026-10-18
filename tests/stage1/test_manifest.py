import pytest

from bespoke_bio.stage1.contracts import RELEASE_CHECKS, Artifact, Check, Manifest


def manifest():
    sha = "a" * 64
    return {
        "run_id": sha,
        "attempt_id": "attempt",
        "case_id": "case",
        "status": "ready",
        "case_sha256": sha,
        "adapter_sha256": sha,
        "contract_sha256": sha,
        "reference_sha256": sha,
        "config_sha256": sha,
        "intake_sha256": sha,
        "artifacts": [
            Artifact(role="case", location="/data/case.json", bytes=1, sha256=sha, format="json")
        ],
        "checks": [
            Check(name=name, status="passed", command=["validator"], evidence_roles=["case"])
            for name in sorted(RELEASE_CHECKS)
        ],
        "transcript_counts": {"selected": 1},
        "target_count": 2,
        "unique_sequence_count": 1,
    }


def test_ready_cannot_omit_an_entire_release_gate():
    value = manifest()
    Manifest(**value)
    value["checks"] = value["checks"][:-1]
    with pytest.raises(ValueError, match="required"):
        Manifest(**value)


def test_successful_check_requires_existing_evidence():
    value = manifest()
    value["checks"][0].evidence_roles = ["absent"]
    with pytest.raises(ValueError, match="missing evidence"):
        Manifest(**value)


def test_unique_protein_count_cannot_exceed_allele_count():
    value = manifest()
    value["unique_sequence_count"] = 3
    with pytest.raises(ValueError, match="exceed"):
        Manifest(**value)
