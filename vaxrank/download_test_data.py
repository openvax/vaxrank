"""Export of the Sid retrieval test cases.

The default exports the 49 retrieval cases from the Sid test data, which is
built from the shared openvax-v1 reads (downloaded once into the osteosarc
cache; ``--offline`` never downloads them). Explicit custom manifests retain the
generic verified-download API for existing callers.
"""

import argparse
import ctypes
import errno
import json
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile
from urllib.parse import urlsplit

from datacache import Cache, FileValidationError, get_data_dir, validate_file
from osteosarc import OsteosarcError


DEFAULT_MANIFEST = None
CACHE_NAMESPACE = "openvax"
CACHE_ENVIRONMENT = "OPENVAX_DATA_CACHE"


def _sid_test_data(offline):
    """The Sid test files; offline, only if openvax-v1 is already cached."""
    import osteosarc
    from .sid_test_data import READS, sid_test_data
    if offline:
        osteosarc.fetch_bundle(READS, offline=True)
    return sid_test_data()


def load_manifest(path=DEFAULT_MANIFEST, *, offline=False):
    """Validate an explicit manifest, or the Sid test data's retrieval manifest."""
    if path is None:
        path = _sid_test_data(offline) / "osteosarc/shared-v1/manifest.json"
    manifest = json.loads(Path(path).read_text())
    if manifest.get("schema_version") != 1 or not manifest.get("dataset") or not manifest.get("data_version"):
        raise ValueError("Unsupported test-data manifest")
    assets = manifest.get("assets")
    if not isinstance(assets, list) or not assets:
        raise ValueError("Manifest must contain assets")
    names = set()
    for asset in assets:
        name = asset["filename"]
        if not isinstance(name, str) or not name or name in (".", "..", "manifest.json") or "/" in name or "\\" in name or ":" in name:
            raise ValueError("Asset filename must be a plain, safe basename")
        if name in names:
            raise ValueError("Duplicate asset filename: " + name)
        names.add(name)
        if not re.fullmatch("[0-9a-f]{64}", asset["sha256"]):
            raise ValueError("Asset needs a lowercase SHA-256 digest")
        if type(asset["size_bytes"]) is not int or asset["size_bytes"] < 0:
            raise ValueError("Asset needs a nonnegative byte size")
        if "bundle_path" in asset:
            if asset["bundle_path"] != "osteosarc/shared-v1/" + name:
                raise ValueError("Invalid bundled asset path")
        elif urlsplit(asset["url"]).scheme not in ("https", "http", "file"):
            raise ValueError("Unsupported asset URL scheme")
    cases = manifest.get("cases", [])
    if len({c["case_id"] for c in cases}) != len(cases):
        raise ValueError("Duplicate test case identity")
    if any(c["bam"] not in names or c["bam"] + ".bai" not in names for c in cases):
        raise ValueError("Test case BAM/index absent from manifest")
    return manifest


def cache_filename(asset):
    """Content identity shared across consumers, URLs and dataset revisions."""
    # Keep the original suffix so datacache never infers decompression from a
    # suffix-less destination for a compressed upstream object.
    return "%s%s" % (asset["sha256"], "".join(Path(asset["filename"]).suffixes))


def verify_dataset(output, manifest=None):
    """Validate a materialized dataset offline, without repair or cache writes."""
    manifest = load_manifest() if manifest is None else manifest
    output = Path(output)
    if output.is_symlink() or not output.is_dir():
        raise ValueError("Dataset must be a real directory: " + str(output))
    expected_names = {a["filename"] for a in manifest["assets"]} | {"manifest.json"}
    if {p.name for p in output.iterdir()} != expected_names:
        raise ValueError("Dataset has missing or unexpected files: " + str(output))
    if (output / "manifest.json").is_symlink():
        raise ValueError("Dataset manifest must not be a symlink")
    if json.loads((output / "manifest.json").read_text()) != manifest:
        raise ValueError("Dataset manifest does not match requested dataset")
    for asset in manifest["assets"]:
        path = output / asset["filename"]
        if path.is_symlink():
            raise ValueError("Dataset assets must not be symlinks")
        validate_file(path, expected_sha256=asset["sha256"], expected_size=asset["size_bytes"])
    return output


def _publish_no_replace(staged, output):
    """Atomically publish a directory, refusing even an empty destination.

    Python's rename replaces empty directories on POSIX. Use the native
    exclusive operation on Linux/macOS, and fail closed if it is unavailable.
    """
    libc = ctypes.CDLL(None, use_errno=True)
    if sys.platform == "linux":
        rename = getattr(libc, "renameat2", None)
        argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_uint]
        # AT_FDCWD = -100; RENAME_NOREPLACE = 1 (Linux renameat2(2)).
        args = (-100, os.fsencode(staged), -100, os.fsencode(output), 1)
    elif sys.platform == "darwin":
        rename = getattr(libc, "renamex_np", None)
        argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_uint]
        # RENAME_EXCL = 0x4 (Darwin sys/stdio.h).
        args = (os.fsencode(staged), os.fsencode(output), 0x4)
    else:
        rename = None
    if rename is None:
        raise OSError(errno.ENOTSUP, "Atomic no-replace directory publication is unavailable", str(output))
    rename.argtypes = argtypes
    rename.restype = ctypes.c_int
    if rename(*args) != 0:
        error = ctypes.get_errno()
        raise OSError(error, os.strerror(error), str(output))


def download_test_data(output=None, *, manifest_path=DEFAULT_MANIFEST,
                       cache_root=None, offline=False, repair_cache=False):
    """Export the Sid retrieval cases, or fetch an explicitly supplied custom manifest.

    Existing output directories are only validated, never overwritten. Failed
    downloads leave verified cache objects reusable but do not publish output.
    ``repair_cache`` explicitly permits replacing invalid cached objects.
    """
    if offline and repair_cache:
        raise ValueError("Offline mode cannot repair the download cache")
    manifest = load_manifest(manifest_path, offline=offline)
    if output is not None:
        output = Path(output).absolute()
        if output.exists() or output.is_symlink():
            return verify_dataset(output, manifest)
    if all("bundle_path" in asset for asset in manifest["assets"]):
        bundled = _sid_test_data(offline) / "osteosarc/shared-v1"
        verify_dataset(bundled, manifest)
        paths = {asset["filename"]: bundled / asset["filename"] for asset in manifest["assets"]}
    else:
        root = cache_root or os.environ.get(CACHE_ENVIRONMENT) or get_data_dir(CACHE_NAMESPACE)
        cache = Cache(CACHE_NAMESPACE, cache_root=Path(root) / "objects" / "sha256")
        paths = {}
        for asset in manifest["assets"]:
            name = cache_filename(asset)
            expected = dict(expected_sha256=asset["sha256"], expected_size=asset["size_bytes"])
            if offline:
                path = cache.local_path(filename=name)
                validate_file(path, **expected)
            else:
                force = False
                if repair_cache:
                    try:
                        validate_file(cache.local_path(filename=name), **expected)
                    except FileValidationError:
                        force = True
                    except FileNotFoundError:
                        pass
                path = cache.fetch(asset["url"], filename=name, timeout=60,
                                   force=force, **expected)
            paths[asset["filename"]] = Path(path)
    if output is None:
        return paths
    output.parent.mkdir(parents=True, exist_ok=True)
    # A one-time output export, not a general mutable bundle registry. Siblings
    # permit atomic exclusive rename; no destination tree is deleted or replaced.
    with tempfile.TemporaryDirectory(prefix=".osteosarc-", dir=output.parent) as temporary:
        staged = Path(temporary) / "dataset"
        staged.mkdir()
        for name, path in paths.items():
            shutil.copyfile(path, staged / name)
        (staged / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        verify_dataset(staged, manifest)
        try:
            _publish_no_replace(staged, output)
        except FileExistsError:
            # A concurrent equivalent publisher is fine; a foreign/modified
            # output is rejected without deleting or repairing it.
            return verify_dataset(output, manifest)
    return verify_dataset(output, manifest)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output", type=Path, help="New directory to export, or an existing identical dataset to verify")
    parser.add_argument("--cache-root", type=Path, help="Shared OpenVax cache root (or OPENVAX_DATA_CACHE)")
    parser.add_argument("--offline", action="store_true", help="Use verified local objects only; never access the network")
    parser.add_argument("--repair-cache", action="store_true", help="Explicitly refetch invalid cached objects")
    parser.add_argument("--verify-only", action="store_true", help="Validate --output without fetching or touching the cache")
    args = parser.parse_args(argv)
    if args.verify_only and (args.output is None or args.repair_cache):
        parser.error("--verify-only requires --output and cannot repair the cache")
    try:
        manifest = load_manifest(args.manifest, offline=args.offline or args.verify_only)
        if args.verify_only:
            result = verify_dataset(args.output, manifest)
        else:
            result = download_test_data(args.output, manifest_path=args.manifest,
                cache_root=args.cache_root, offline=args.offline, repair_cache=args.repair_cache)
    except (OSError, ValueError, FileValidationError, OsteosarcError) as error:
        parser.exit(1, "Test data unavailable: %s\n" % error)
    print(json.dumps(dict(dataset=manifest["dataset"], data_version=manifest["data_version"],
        assets=len(manifest["assets"]), output=str(result) if isinstance(result, Path) else None,
        cached_paths={k: str(v) for k, v in result.items()} if isinstance(result, dict) else None), indent=2))


if __name__ == "__main__":
    main()
