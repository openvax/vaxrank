import io
import json
import subprocess
from pathlib import Path
from urllib.error import HTTPError

import pytest
import release_upload

from release_upload import (
    ReleaseUploadError,
    expected_release_filenames,
    main,
    pypi_release_filenames,
    publish_release,
    upload_distribution,
    wait_for_release_file,
)


EXPECTED = expected_release_filenames("vaxrank", "3.1.28")
PATHS = tuple(Path("dist") / filename for filename in EXPECTED)


def test_expected_release_filenames_are_exact_wheel_and_sdist():
    assert EXPECTED == {
        "vaxrank-3.1.28-py3-none-any.whl",
        "vaxrank-3.1.28.tar.gz",
    }


def test_pypi_release_filenames_reads_exact_release_metadata(monkeypatch):
    observed = {}

    def open_url(url, *, timeout):
        observed.update(url=url, timeout=timeout)
        return io.BytesIO(json.dumps({
            "urls": [{"filename": "vaxrank-3.1.28.tar.gz"}],
        }).encode("utf-8"))

    monkeypatch.setattr(release_upload, "urlopen", open_url)

    assert pypi_release_filenames(
        "vaxrank",
        "3.1.28",
        json_base_url="https://packages.example/pypi/",
        request_timeout_seconds=7,
    ) == {"vaxrank-3.1.28.tar.gz"}
    assert observed == {
        "url": "https://packages.example/pypi/vaxrank/3.1.28/json",
        "timeout": 7,
    }


def test_pypi_release_filenames_treats_missing_release_as_empty(monkeypatch):
    def open_url(url, *, timeout):
        raise HTTPError(url, 404, "Not Found", {}, None)

    monkeypatch.setattr(release_upload, "urlopen", open_url)
    assert pypi_release_filenames("vaxrank", "3.1.28") == frozenset()


def test_upload_distribution_is_bounded_and_disables_progress(monkeypatch):
    observed = {}

    def run(command, *, check, timeout):
        observed.update(command=command, check=check, timeout=timeout)

    monkeypatch.setattr(subprocess, "run", run)
    upload_distribution(
        "dist/vaxrank-3.1.28.tar.gz",
        python_executable="release-python",
        timeout_seconds=17,
    )

    assert observed == {
        "command": [
            "release-python",
            "-m",
            "twine",
            "upload",
            "--disable-progress-bar",
            "dist/vaxrank-3.1.28.tar.gz",
        ],
        "check": True,
        "timeout": 17,
    }


def test_wait_for_release_file_observes_eventual_metadata():
    responses = iter([set(), {"vaxrank-3.1.28.tar.gz"}])

    assert wait_for_release_file(
        "vaxrank-3.1.28.tar.gz",
        lambda: next(responses),
        timeout_seconds=1,
        poll_seconds=0,
    )


def test_timeout_after_server_acceptance_is_reconciled():
    published = set()

    def fetch_release_filenames():
        return published

    def upload_file(path):
        published.add(path.name)
        raise subprocess.TimeoutExpired(["twine", "upload"], timeout=60)

    result = publish_release(
        PATHS,
        expected_filenames=EXPECTED,
        fetch_release_filenames=fetch_release_filenames,
        upload_file=upload_file,
        verify_timeout_seconds=0,
    )

    assert EXPECTED <= result


def test_partial_release_uploads_only_missing_file_and_retry_is_safe():
    wheel = "vaxrank-3.1.28-py3-none-any.whl"
    sdist = "vaxrank-3.1.28.tar.gz"
    published = {wheel}
    uploaded = []

    def fetch_release_filenames():
        return published

    def upload_file(path):
        uploaded.append(path.name)
        published.add(path.name)

    publish_release(
        PATHS,
        expected_filenames=EXPECTED,
        fetch_release_filenames=fetch_release_filenames,
        upload_file=upload_file,
        verify_timeout_seconds=0,
    )
    publish_release(
        PATHS,
        expected_filenames=EXPECTED,
        fetch_release_filenames=fetch_release_filenames,
        upload_file=upload_file,
        verify_timeout_seconds=0,
    )

    assert uploaded == [sdist]


def test_failed_upload_without_server_artifact_fails_hard():
    def upload_file(path):
        raise subprocess.TimeoutExpired(["twine", "upload", path], timeout=60)

    with pytest.raises(ReleaseUploadError, match="did not publish"):
        publish_release(
            PATHS,
            expected_filenames=EXPECTED,
            fetch_release_filenames=frozenset,
            upload_file=upload_file,
            verify_timeout_seconds=0,
        )


def test_unexpected_distribution_set_is_rejected_before_upload():
    with pytest.raises(ReleaseUploadError, match="do not match"):
        publish_release(
            ["dist/vaxrank-3.1.28.tar.gz"],
            expected_filenames=EXPECTED,
            fetch_release_filenames=frozenset,
            upload_file=lambda path: None,
            verify_timeout_seconds=0,
        )


def test_main_passes_exact_artifact_contract_to_publisher(monkeypatch, capsys):
    observed = {}

    def publish_release_call(distribution_paths, **kwargs):
        observed.update(
            distribution_paths=tuple(distribution_paths),
            expected_filenames=kwargs["expected_filenames"],
        )
        return EXPECTED

    monkeypatch.setattr(release_upload, "publish_release", publish_release_call)

    assert main([
        "--project", "vaxrank",
        "--version", "3.1.28",
        *[str(path) for path in PATHS],
    ]) == 0
    assert observed == {
        "distribution_paths": tuple(str(path) for path in PATHS),
        "expected_filenames": EXPECTED,
    }
    assert "Verified PyPI release files" in capsys.readouterr().out
