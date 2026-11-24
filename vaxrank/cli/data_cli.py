# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Read-only cache inspection and maintenance subcommands: ``vaxrank data ...``."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path


def _dir_stats(path: str | Path) -> tuple[bool, int, int]:
    root = Path(path)
    if not root.exists():
        return False, 0, 0
    total_bytes = 0
    file_count = 0
    for dirpath, _dirnames, filenames in os.walk(root, followlinks=False):
        for name in filenames:
            try:
                st = os.lstat(os.path.join(dirpath, name))
            except OSError:
                continue
            total_bytes += st.st_size
            file_count += 1
    return True, total_bytes, file_count


def _build_parser(prog: str = "vaxrank data") -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=prog,
        description="Inspect and maintain Vaxrank's local reference-peptide and OpenVax caches.",
    )
    sub = parser.add_subparsers(dest="data_command", required=True, metavar="COMMAND")

    p_dirs = sub.add_parser(
        "dirs",
        help="Show where Vaxrank's cache lives, why it resolved there, and disk usage",
    )
    p_dirs.add_argument("--json", action="store_true", help="Print machine-readable JSON")

    p_status = sub.add_parser(
        "status",
        help="List cached reference-proteome k-mer index files and sizes",
    )
    p_status.add_argument("--json", action="store_true", help="Print machine-readable JSON")

    p_path = sub.add_parser(
        "path",
        help="Print Vaxrank's resolved cache directory without creating it",
    )
    p_path.add_argument("--json", action="store_true", help="Print machine-readable JSON")

    p_prune = sub.add_parser(
        "prune",
        help="Remove stale datacache bundle generations (dry-run by default)",
    )
    p_prune.add_argument("--keep", type=int, default=1, help="Generations to keep per store (default: 1)")
    p_prune.add_argument(
        "--kmer-indexes",
        action="store_true",
        help="Also remove cached reference-proteome k-mer index files (.pkl / .pkl.gz)",
    )
    p_prune.add_argument("--yes", action="store_true", help="Actually delete selected files (default is dry-run)")
    p_prune.add_argument("--json", action="store_true", help="Print machine-readable JSON")

    return parser


def _kmer_index_files(cache_dir: str | Path) -> list[dict]:
    root = Path(cache_dir)
    if not root.is_dir():
        return []
    files = []
    for entry in sorted(os.scandir(root), key=lambda e: e.name):
        if not entry.is_file(follow_symlinks=False):
            continue
        if not (entry.name.endswith(".pkl.gz") or entry.name.endswith(".pkl")):
            continue
        try:
            st = os.lstat(entry.path)
        except OSError:
            continue
        files.append(dict(
            name=entry.name,
            path=entry.path,
            size_bytes=st.st_size,
            mtime=st.st_mtime,
        ))
    return files


def run(argv: list[str] | None = None, prog: str = "vaxrank data") -> int:
    """Execute ``vaxrank data <dirs|status|path|prune>``."""
    parser = _build_parser(prog=prog)
    args = parser.parse_args(argv)

    import datacache
    from datacache.cli import prune_root
    from ..reference_proteome import resolve_cache_location

    loc = resolve_cache_location()

    if args.data_command == "path":
        if args.json:
            print(json.dumps(dict(path=loc.path, origin=loc.origin, detail=loc.detail), indent=2))
        else:
            print(loc.path)
        return 0

    if args.data_command == "dirs":
        exists, total_bytes, file_count = _dir_stats(loc.path)
        offline = datacache.offline_source()
        shared_root = datacache.get_cache_root("openvax", datacache.SHARED_CACHE_ENV_KEY)
        payload = dict(
            cache_dir=loc.path,
            origin=loc.origin,
            detail=loc.detail,
            exists=exists,
            size_bytes=total_bytes,
            file_count=file_count,
            shared_openvax_root=shared_root,
            offline=offline is not None,
            offline_source=offline,
        )
        if args.json:
            print(json.dumps(payload, indent=2))
            return 0
        print(f"Vaxrank cache:     {loc.path}")
        print(f"Origin:            {loc.detail} ({loc.origin})")
        if exists:
            print(f"Size:              {datacache.format_bytes(total_bytes)} ({file_count} file(s))")
        else:
            print("Size:              not created yet")
        print(f"Shared OpenVax:    {shared_root}")
        if offline:
            print(f"Offline:           enabled ({offline})")
        return 0

    if args.data_command == "status":
        indexes = _kmer_index_files(loc.path)
        total_bytes = sum(f["size_bytes"] for f in indexes)
        payload = dict(
            cache_dir=loc.path,
            origin=loc.origin,
            detail=loc.detail,
            kmer_index_count=len(indexes),
            total_bytes=total_bytes,
            kmer_indexes=indexes,
        )
        if args.json:
            print(json.dumps(payload, indent=2))
            return 0
        print(f"Vaxrank cache: {loc.path} ({loc.detail})")
        if not indexes:
            print("No cached reference-proteome k-mer indexes found.")
            return 0
        for item in indexes:
            print(f"  {datacache.format_bytes(item['size_bytes']):>10}  {item['name']}")
        print(f"\nTotal: {datacache.format_bytes(total_bytes)} across {len(indexes)} index file(s).")
        return 0

    if args.data_command == "prune":
        store_results = prune_root(loc.path, keep=args.keep, apply=args.yes)
        removed_bundles = sum(r["removed_count"] for r in store_results)
        reclaimed_bundles = sum(r["reclaimable_bytes"] or 0 for r in store_results)

        kmer_removed: list[str] = []
        kmer_bytes = 0
        if args.kmer_indexes:
            for item in _kmer_index_files(loc.path):
                kmer_bytes += item["size_bytes"]
                kmer_removed.append(item["path"])
                if args.yes:
                    try:
                        os.unlink(item["path"])
                    except FileNotFoundError:
                        pass

        total_reclaimable = reclaimed_bundles + kmer_bytes
        total_items = removed_bundles + len(kmer_removed)
        payload = dict(
            cache_dir=loc.path,
            dry_run=not args.yes,
            keep=args.keep,
            removed_bundle_count=removed_bundles,
            removed_kmer_indexes=kmer_removed,
            reclaimable_bytes=total_reclaimable,
            stores=store_results,
        )
        if args.json:
            print(json.dumps(payload, indent=2))
            return 0
        if total_items == 0:
            print(f"Nothing to prune under {loc.path}.")
            return 0
        verb = "Removed" if args.yes else "Would remove"
        if removed_bundles:
            print(f"{verb} {removed_bundles} stale bundle generation(s) ({datacache.format_bytes(reclaimed_bundles)}).")
        if kmer_removed:
            print(f"{verb} {len(kmer_removed)} k-mer index file(s) ({datacache.format_bytes(kmer_bytes)}).")
        if not args.yes:
            print("\nDry run only. Pass --yes to delete these files.")
        return 0

    return 0
