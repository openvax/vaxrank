"""Pin library-first acquisition choices; retain every source in the inventory."""

import argparse
import json
from pathlib import Path

from build import digest, now, write_json


def select_sources(sources, groups):
    decisions = {}
    known = {s["source_id"] for s in sources}
    for group in groups:
        roles = {"primary": [group["primary"]],
                 "tagged_evidence": group.get("evidence_sources", []),
                 "full_input_cached_separately": group.get("full_input_sources", []),
                 "alternative": group.get("alternatives", []),
                 "unassigned_partition": group.get("unassigned", [])}
        for role, ids in roles.items():
            for sid in ids:
                if sid not in known or sid in decisions:
                    raise ValueError("Unknown or multiply grouped source: " + sid)
                decisions[sid] = dict(source_id=sid, group=group["group"], role=role,
                    acquire=role in ("primary", "tagged_evidence"), reason=group["evidence"])
    for source in sources:
        sid = source["source_id"]
        if sid not in decisions:
            decisions[sid] = dict(source_id=sid, group="ungrouped-" + sid, role="primary",
                acquire=True, reason="No reviewed processing-equivalence group; retain separately to avoid losing a library, provider or DNA control.")
    return [dict(decisions[s["source_id"]], url=s["url"], metadata=s["metadata"]) for s in sources]


def pin_selection(root):
    groups_path = Path(__file__).with_name("library_groups.json")
    sources = json.loads((root / "source_inventory.json").read_text())
    decisions = select_sources(sources, json.loads(groups_path.read_text()))
    policy = dict(created_utc=now(), policy="library_timepoint_coverage_first",
        requested_by="User: prioritize complete library/timepoint coverage; keep alternatives inventoried",
        groups_sha256=digest(groups_path), source_inventory_sha256=digest(root / "source_inventory.json"),
        interpretation="Groups are reviewed acquisition families, not independent-library counts or assertions of identical reads. Preserve conflicting metadata; never pool alternative products.",
        sources=decisions)
    for part in (root, root / "phase-and-sv-context"):
        inventory = json.loads((part / "source_inventory.json").read_text())
        if {s["source_id"]: s["url"] for s in inventory} != {s["source_id"]: s["url"] for s in sources}:
            raise ValueError("Acquisition parts have different source inventories")
    for part in (root, root / "phase-and-sv-context"):
        write_json(part / "source-selection.json", policy)
    print(f"Selected {sum(s['acquire'] for s in decisions)} regional products; retained all {len(sources)} inventory entries")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    pin_selection(parser.parse_args().run)
