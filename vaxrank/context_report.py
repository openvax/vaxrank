# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

"""Whole-context reports with exact-bond, target and all-ligand overlays."""

from pathlib import Path

import jinja2
from mhctools.pred import value_unit

from .context_audit import SequenceContextAudit
from .native_serialization import to_native_json


def _profile_view(audit, profile):
    roles = {b: {"target_internal": [], "target_boundary": [], "ligand_internal": [],
                 "ligand_n": [], "ligand_c": []} for b in range(1, len(audit.context.sequence))}
    for i, target in enumerate(audit.context.targets, 1):
        for b in range(target.start + 1, target.end):
            roles[b]["target_internal"].append(i)
        for b in {target.start, target.end}:
            if b in roles:
                roles[b]["target_boundary"].append(i)
    for i, ligand in enumerate(audit.ligands, 1):
        for b in range(ligand.start + 1, ligand.end):
            roles[b]["ligand_internal"].append(i)
        if ligand.start in roles:
            roles[ligand.start]["ligand_n"].append(i)
        if ligand.end in roles:
            roles[ligand.end]["ligand_c"].append(i)
    return {"profile": profile, "sites": {s.bond: s for s in profile.sites}, "roles": roles}


def write_sequence_context_audits(audits, *, json_path, html_path=None):
    """Export full native graphs plus explicitly derived coverage and overlays.

    The JSON object has ``audits``, ``coverage`` and ``overlays`` entries, each
    round-trippable through the allowlisted native serializer. Derived entries
    are conveniences, never a replacement for original source/model records.
    """
    audits = tuple(audits)
    for audit in audits:
        if not isinstance(audit, SequenceContextAudit):
            raise ValueError("Context report requires typed sequence audits")
        audit.validate()
    payload = to_native_json({
        "audits": audits, "coverage": tuple(a.coverage for a in audits),
        "overlays": tuple(tuple(a.overlays(p) for p in a.profiles) for a in audits)})
    rendered = None
    if html_path is not None:
        environment = jinja2.Environment(
            loader=jinja2.PackageLoader("vaxrank", "templates"),
            autoescape=jinja2.select_autoescape(("html", "xml")))
        records = tuple({"audit": a, "profiles": tuple(_profile_view(a, p) for p in a.profiles)} for a in audits)
        rendered = environment.get_template("sequence_context_audit.html").render(records=records, value_unit=value_unit)
    Path(json_path).write_text(payload, encoding="utf8")
    if rendered is not None:
        Path(html_path).write_text(rendered, encoding="utf8")
