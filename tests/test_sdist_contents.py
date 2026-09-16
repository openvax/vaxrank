"""The source distribution contains one complete, intentional payload."""

from pathlib import Path, PurePosixPath
import subprocess
import sys
import tarfile


ROOT = Path(__file__).resolve().parents[1]
RAW_GENOMIC_SUFFIXES = (
    ".bam",
    ".bai",
    ".cram",
    ".fastq",
    ".fastq.gz",
    ".fq",
    ".fq.gz",
    ".sam",
    ".sam.gz",
)


def _build_sdist(output_directory):
    result = subprocess.run(
        [
            sys.executable,
            "setup.py",
            "--quiet",
            "sdist",
            "--dist-dir",
            str(output_directory),
        ],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    archives = list(output_directory.glob("vaxrank-*.tar.gz"))
    assert len(archives) == 1, archives
    return archives[0]


def test_sdist_omits_repository_only_tests_and_raw_genomic_data(tmp_path):
    """Inspect a built archive, not only setup metadata or the checkout.

    The suite cannot run without its repository-only helpers and fixtures.
    Publishing a subset makes the sdist look testable while imports and data
    access fail. Including the complete tree would publish raw genomic data,
    which is outside the package payload and requires separate approval.
    """
    archive = _build_sdist(tmp_path)
    with tarfile.open(archive, "r:gz") as sdist:
        members = [PurePosixPath(member.name) for member in sdist.getmembers()
                   if member.isfile()]

    # Every member is under the archive's versioned root directory.
    assert members and all(len(path.parts) > 1 for path in members)
    relative_paths = [PurePosixPath(*path.parts[1:]) for path in members]
    assert not any(path.parts[0] == "tests" for path in relative_paths)

    raw_genomic_data = sorted(
        str(path) for path in relative_paths
        if str(path).lower().endswith(RAW_GENOMIC_SUFFIXES)
    )
    assert raw_genomic_data == []

    # Guard against a vacuous archive that satisfies the exclusions by
    # accidentally dropping the actual installable project too.
    included = {str(path) for path in relative_paths}
    assert {
        "LICENSE",
        "MANIFEST.in",
        "README.md",
        "requirements.txt",
        "setup.py",
        "vaxrank/__init__.py",
        "vaxrank/version.py",
    } <= included
