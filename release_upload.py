"""Bounded, idempotent PyPI uploads for vaxrank releases."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from secrets import token_hex
import subprocess
import sys
import time
from pathlib import Path
from typing import Callable, Iterable, Mapping
from urllib.error import HTTPError
from urllib.parse import quote
from urllib.request import urlopen


PYPI_JSON_BASE_URL = "https://pypi.org/pypi"
DEFAULT_UPLOAD_TIMEOUT_SECONDS = 60.0
DEFAULT_VERIFY_TIMEOUT_SECONDS = 30.0
DEFAULT_VERIFY_POLL_SECONDS = 1.0


class ReleaseUploadError(RuntimeError):
    """Raised when the complete release artifact set cannot be verified."""


def expected_release_filenames(project: str, version: str) -> frozenset[str]:
    """Return the wheel and source-distribution names built by vaxrank."""
    wheel_project = project.replace("-", "_")
    return frozenset({
        f"{wheel_project}-{version}-py3-none-any.whl",
        f"{project}-{version}.tar.gz",
    })


def pypi_release_artifacts(
    project: str,
    version: str,
    *,
    json_base_url: str = PYPI_JSON_BASE_URL,
    request_timeout_seconds: float = 10.0,
) -> dict[str, str]:
    """Return one release's filename-to-SHA-256 map from PyPI metadata."""
    url = "%s/%s/json?cache_bust=%s" % (
        json_base_url.rstrip("/"),
        quote(project, safe=""),
        token_hex(8),
    )
    try:
        with urlopen(url, timeout=request_timeout_seconds) as response:
            payload = json.load(response)
    except HTTPError as error:
        if error.code == 404:
            return {}
        raise
    releases = payload.get("releases", {})
    return {
        item["filename"]: item.get("digests", {}).get("sha256", "")
        for item in releases.get(version, ())
    }


def file_sha256(path: str | Path) -> str:
    """Return the SHA-256 digest of an artifact's exact bytes."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as artifact:
        for chunk in iter(lambda: artifact.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _has_matching_artifact(
    filename: str,
    expected_sha256: str,
    published: Mapping[str, str],
) -> bool:
    """Return true for byte-identical artifacts and false when absent."""
    if filename not in published:
        return False
    observed_sha256 = published[filename]
    if observed_sha256 != expected_sha256:
        raise ReleaseUploadError(
            f"PyPI artifact {filename} has SHA-256 "
            f"{observed_sha256 or '<missing>'}, but the local artifact "
            f"has {expected_sha256}"
        )
    return True


def upload_distribution(
    distribution_path: str | Path,
    *,
    python_executable: str = sys.executable,
    timeout_seconds: float = DEFAULT_UPLOAD_TIMEOUT_SECONDS,
    repository_url: str | None = None,
) -> None:
    """Upload one distribution with a finite timeout and quiet progress."""
    command = [
        python_executable,
        "-m",
        "twine",
        "upload",
        "--disable-progress-bar",
    ]
    if repository_url:
        command.extend(["--repository-url", repository_url])
    command.append(str(distribution_path))
    subprocess.run(command, check=True, timeout=timeout_seconds)


def wait_for_release_file(
    filename: str,
    expected_sha256: str,
    fetch_release_artifacts: Callable[[], Mapping[str, str]],
    *,
    timeout_seconds: float = DEFAULT_VERIFY_TIMEOUT_SECONDS,
    poll_seconds: float = DEFAULT_VERIFY_POLL_SECONDS,
) -> bool:
    """Poll until ``filename`` is visible with the expected exact bytes."""
    deadline = time.monotonic() + timeout_seconds
    while True:
        published = dict(fetch_release_artifacts())
        if _has_matching_artifact(filename, expected_sha256, published):
            return True
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            return False
        time.sleep(min(poll_seconds, remaining))


def wait_for_complete_release(
    expected_artifacts: Mapping[str, str],
    fetch_release_artifacts: Callable[[], Mapping[str, str]],
    *,
    timeout_seconds: float = DEFAULT_VERIFY_TIMEOUT_SECONDS,
    poll_seconds: float = DEFAULT_VERIFY_POLL_SECONDS,
) -> dict[str, str]:
    """Return the latest metadata after waiting for a complete release."""
    deadline = time.monotonic() + timeout_seconds
    while True:
        published = dict(fetch_release_artifacts())
        complete = True
        for filename, expected_sha256 in expected_artifacts.items():
            if not _has_matching_artifact(
                filename,
                expected_sha256,
                published,
            ):
                complete = False
        if complete:
            return published
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            return published
        time.sleep(min(poll_seconds, remaining))


def publish_release(
    distribution_paths: Iterable[str | Path],
    *,
    expected_filenames: Iterable[str],
    fetch_release_artifacts: Callable[[], Mapping[str, str]],
    upload_file: Callable[[Path], None],
    verify_timeout_seconds: float = DEFAULT_VERIFY_TIMEOUT_SECONDS,
    verify_poll_seconds: float = DEFAULT_VERIFY_POLL_SECONDS,
) -> dict[str, str]:
    """Upload missing artifacts and reconcile every response with PyPI state.

    A timeout or nonzero Twine exit is ambiguous: the server may have accepted
    the immutable file before the client lost its response. This operation
    therefore verifies server state after every attempt. A file is treated as
    published only when both its filename and SHA-256 digest match the local
    artifact, so a retry cannot combine distributions from different builds.
    """
    paths_by_name = {
        Path(distribution_path).name: Path(distribution_path)
        for distribution_path in distribution_paths
    }
    expected = frozenset(expected_filenames)
    if frozenset(paths_by_name) != expected:
        raise ReleaseUploadError(
            "Distribution files do not match the expected release set: "
            f"expected={sorted(expected)}, observed={sorted(paths_by_name)}"
        )
    local_artifacts = {
        filename: file_sha256(path)
        for filename, path in paths_by_name.items()
    }

    published = dict(fetch_release_artifacts())
    for filename in sorted(expected):
        if _has_matching_artifact(
            filename,
            local_artifacts[filename],
            published,
        ):
            print(f"Already published with matching SHA-256: {filename}")
            continue

        upload_error = None
        try:
            upload_file(paths_by_name[filename])
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as error:
            upload_error = error

        if not wait_for_release_file(
            filename,
            local_artifacts[filename],
            fetch_release_artifacts,
            timeout_seconds=verify_timeout_seconds,
            poll_seconds=verify_poll_seconds,
        ):
            message = f"PyPI did not publish {filename} after the upload attempt"
            if upload_error is not None:
                raise ReleaseUploadError(message) from upload_error
            raise ReleaseUploadError(message)
        if upload_error is not None:
            print(f"Reconciled ambiguous upload from PyPI state: {filename}")
        else:
            print(f"Published: {filename}")
        published = dict(fetch_release_artifacts())

    published = wait_for_complete_release(
        local_artifacts,
        fetch_release_artifacts,
        timeout_seconds=verify_timeout_seconds,
        poll_seconds=verify_poll_seconds,
    )
    missing = expected - set(published)
    if missing:
        raise ReleaseUploadError(
            f"PyPI release is missing expected files: {sorted(missing)}"
        )
    return published


def main(argv: list[str] | None = None) -> int:
    """Upload and verify one complete vaxrank release."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", required=True)
    parser.add_argument("--version", required=True)
    parser.add_argument("--json-base-url", default=PYPI_JSON_BASE_URL)
    parser.add_argument("--repository-url")
    parser.add_argument(
        "--upload-timeout-seconds",
        type=float,
        default=float(os.environ.get(
            "VAXRANK_UPLOAD_TIMEOUT_SECONDS",
            DEFAULT_UPLOAD_TIMEOUT_SECONDS,
        )),
    )
    parser.add_argument(
        "--verify-timeout-seconds",
        type=float,
        default=float(os.environ.get(
            "VAXRANK_VERIFY_TIMEOUT_SECONDS",
            DEFAULT_VERIFY_TIMEOUT_SECONDS,
        )),
    )
    parser.add_argument("distributions", nargs="+")
    args = parser.parse_args(argv)

    expected = expected_release_filenames(args.project, args.version)

    def fetch_release_artifacts():
        return pypi_release_artifacts(
            args.project,
            args.version,
            json_base_url=args.json_base_url,
        )

    def upload_file(path):
        upload_distribution(
            path,
            python_executable=sys.executable,
            timeout_seconds=args.upload_timeout_seconds,
            repository_url=args.repository_url,
        )

    published = publish_release(
        args.distributions,
        expected_filenames=expected,
        fetch_release_artifacts=fetch_release_artifacts,
        upload_file=upload_file,
        verify_timeout_seconds=args.verify_timeout_seconds,
    )
    print(
        "Verified PyPI release files: %s"
        % ", ".join(sorted(expected & set(published)))
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
