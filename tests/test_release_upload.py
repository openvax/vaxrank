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
    file_sha256,
    main,
    publish_release,
    pypi_release_artifacts,
    upload_distribution,
    wait_for_complete_release,
    wait_for_release_file,
)


EXPECTED = expected_release_filenames("vaxrank", "3.1.28")
PATHS = tuple(Path("dist") / filename for filename in EXPECTED)


@pytest.fixture
def release_paths(tmp_path):
    paths = []
    for filename in sorted(EXPECTED):
        path = tmp_path / filename
        path.write_bytes((filename + "\n").encode("utf-8"))
        paths.append(path)
    return tuple(paths)


def _local_artifacts(paths):
    return {path.name: file_sha256(path) for path in paths}


def test_expected_release_filenames_are_exact_wheel_and_sdist():
    assert EXPECTED == {
        "vaxrank-3.1.28-py3-none-any.whl",
        "vaxrank-3.1.28.tar.gz",
    }


def test_pypi_release_artifacts_reads_project_release_map(monkeypatch):
    observed = {}

    def open_url(url, *, timeout):
        observed.update(url=url, timeout=timeout)
        return io.BytesIO(json.dumps({
            "releases": {
                "3.1.28": [{
                    "filename": "vaxrank-3.1.28.tar.gz",
                    "digests": {"sha256": "abc123"},
                }],
            },
        }).encode("utf-8"))

    monkeypatch.setattr(release_upload, "urlopen", open_url)
    monkeypatch.setattr(release_upload, "token_hex", lambda length: "fresh-token")

    assert pypi_release_artifacts(
        "vaxrank",
        "3.1.28",
        json_base_url="https://packages.example/pypi/",
        request_timeout_seconds=7,
    ) == {"vaxrank-3.1.28.tar.gz": "abc123"}
    assert observed == {
        "url": (
            "https://packages.example/pypi/vaxrank/json"
            "?cache_bust=fresh-token"
        ),
        "timeout": 7,
    }


def test_pypi_release_artifacts_uses_a_fresh_url_for_each_poll(monkeypatch):
    urls = []
    tokens = iter(["first-token", "second-token"])

    def open_url(url, *, timeout):
        urls.append(url)
        return io.BytesIO(b'{"releases": {}}')

    monkeypatch.setattr(release_upload, "urlopen", open_url)
    monkeypatch.setattr(release_upload, "token_hex", lambda length: next(tokens))

    pypi_release_artifacts("vaxrank", "3.1.28")
    pypi_release_artifacts("vaxrank", "3.1.28")

    assert urls == [
        "https://pypi.org/pypi/vaxrank/json?cache_bust=first-token",
        "https://pypi.org/pypi/vaxrank/json?cache_bust=second-token",
    ]


def test_pypi_release_artifacts_treats_missing_release_as_empty(monkeypatch):
    def open_url(url, *, timeout):
        raise HTTPError(url, 404, "Not Found", {}, None)

    monkeypatch.setattr(release_upload, "urlopen", open_url)
    assert pypi_release_artifacts("vaxrank", "3.1.28") == {}


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


def test_wait_for_release_file_observes_eventual_matching_metadata():
    filename = "vaxrank-3.1.28.tar.gz"
    responses = iter([{}, {filename: "abc123"}])

    assert wait_for_release_file(
        filename,
        "abc123",
        lambda: next(responses),
        timeout_seconds=1,
        poll_seconds=0,
    )


def test_wait_for_release_file_rejects_different_bytes_immediately():
    filename = "vaxrank-3.1.28.tar.gz"
    with pytest.raises(ReleaseUploadError, match="SHA-256"):
        wait_for_release_file(
            filename,
            "local-digest",
            lambda: {filename: "remote-digest"},
            timeout_seconds=1,
            poll_seconds=0,
        )


def test_complete_release_verification_retries_stale_metadata():
    expected = {"wheel": "wheel-digest", "sdist": "sdist-digest"}
    responses = iter([
        {"wheel": "wheel-digest"},
        expected,
    ])

    assert wait_for_complete_release(
        expected,
        lambda: next(responses),
        timeout_seconds=1,
        poll_seconds=0,
    ) == expected


def test_timeout_after_server_acceptance_is_reconciled(release_paths):
    published = {}

    def fetch_release_artifacts():
        return published

    def upload_file(path):
        published[path.name] = file_sha256(path)
        raise subprocess.TimeoutExpired(["twine", "upload"], timeout=60)

    result = publish_release(
        release_paths,
        expected_filenames=EXPECTED,
        fetch_release_artifacts=fetch_release_artifacts,
        upload_file=upload_file,
        verify_timeout_seconds=0,
    )

    assert EXPECTED <= set(result)


def test_partial_release_uploads_only_missing_file_and_retry_is_safe(
        release_paths):
    wheel = "vaxrank-3.1.28-py3-none-any.whl"
    sdist = "vaxrank-3.1.28.tar.gz"
    local = _local_artifacts(release_paths)
    published = {wheel: local[wheel]}
    uploaded = []

    def fetch_release_artifacts():
        return published

    def upload_file(path):
        uploaded.append(path.name)
        published[path.name] = file_sha256(path)

    for _ in range(2):
        publish_release(
            release_paths,
            expected_filenames=EXPECTED,
            fetch_release_artifacts=fetch_release_artifacts,
            upload_file=upload_file,
            verify_timeout_seconds=0,
        )

    assert uploaded == [sdist]


def test_existing_filename_with_different_bytes_fails_before_upload(
        release_paths):
    wheel = "vaxrank-3.1.28-py3-none-any.whl"
    uploaded = []

    with pytest.raises(ReleaseUploadError, match="SHA-256"):
        publish_release(
            release_paths,
            expected_filenames=EXPECTED,
            fetch_release_artifacts=lambda: {wheel: "0" * 64},
            upload_file=lambda path: uploaded.append(path.name),
            verify_timeout_seconds=0,
        )

    assert uploaded == []


def test_failed_upload_without_server_artifact_fails_hard(release_paths):
    def upload_file(path):
        raise subprocess.TimeoutExpired(["twine", "upload", path], timeout=60)

    with pytest.raises(ReleaseUploadError, match="did not publish"):
        publish_release(
            release_paths,
            expected_filenames=EXPECTED,
            fetch_release_artifacts=dict,
            upload_file=upload_file,
            verify_timeout_seconds=0,
        )


def test_unexpected_distribution_set_is_rejected_before_upload(tmp_path):
    path = tmp_path / "vaxrank-3.1.28.tar.gz"
    path.write_bytes(b"sdist")
    with pytest.raises(ReleaseUploadError, match="do not match"):
        publish_release(
            [path],
            expected_filenames=EXPECTED,
            fetch_release_artifacts=dict,
            upload_file=lambda artifact: None,
            verify_timeout_seconds=0,
        )


def test_main_passes_exact_artifact_contract_to_publisher(monkeypatch, capsys):
    observed = {}

    def publish_release_call(distribution_paths, **kwargs):
        observed.update(
            distribution_paths=tuple(distribution_paths),
            expected_filenames=kwargs["expected_filenames"],
        )
        return {filename: "digest" for filename in EXPECTED}

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
