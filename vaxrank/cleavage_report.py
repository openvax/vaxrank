# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Complete bond inventories through small report glue and a Jinja template."""

from pathlib import Path

import jinja2

from .cleavage_profile import CleavageProfile
from .native_serialization import to_native_json


def write_cleavage_profiles(profiles, *, json_path, html_path=None):
    """Export all model observations, missing bonds, settings and source offsets.

    This low-level report does not infer target identities or HLA ligands from
    cleavage data. Those are separate sequence-context overlays, not cleavage
    predictions. JSON retains the complete immutable native object graph.
    """
    profiles = tuple(profiles)
    if any(not isinstance(profile, CleavageProfile) for profile in profiles):
        raise ValueError("Cleavage report requires typed profiles")
    payload = to_native_json(profiles)
    rendered = None
    if html_path is not None:
        environment = jinja2.Environment(
            loader=jinja2.PackageLoader("vaxrank", "templates"),
            autoescape=jinja2.select_autoescape(("html", "xml")))
        records = tuple({"profile": profile,
                         "by_bond": {site.bond: site for site in profile.sites},
                         "padded": frozenset(profile.padded_bonds)} for profile in profiles)
        rendered = environment.get_template("cleavage_profiles.html").render(records=records)
    Path(json_path).write_text(payload, encoding="utf8")
    if rendered is not None:
        Path(html_path).write_text(rendered, encoding="utf8")
