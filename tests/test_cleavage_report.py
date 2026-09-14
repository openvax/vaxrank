"""Packaged report retains all bonds, native scores and unknown coverage."""

from dataclasses import replace

from mhctools.cleavage import (
    CleavageInput, CleavageModel, CleavageResult, CleavageSite,
)
import pytest

from vaxrank.cleavage_profile import profile_from_cleavage_result
from vaxrank.cleavage_report import write_cleavage_profiles
from vaxrank.native_serialization import from_native_json

from .test_cleavage_profile import MODEL, MOTIF_MODEL, profile


def test_complete_profiles_roundtrip_and_render_all_sites(tmp_path):
    p = profile(backend_version="3.41.0", upstream_context=8, downstream_context=8)
    p = replace(p, peptide=replace(p.peptide, source_id="<script>bad()</script>", source_start=40))
    partial = profile_from_cleavage_result(CleavageResult(
        CleavageInput("ACDEFG"), MOTIF_MODEL, (CleavageSite(5, "not_matched", "partial rule"),)))
    json_path, html_path = tmp_path / "sites.json", tmp_path / "sites.html"
    write_cleavage_profiles([p, partial], json_path=json_path, html_path=html_path)
    assert from_native_json(json_path.read_text(), tuple) == (p, partial)
    html = html_path.read_text()
    assert "<script>" not in html and "&lt;script&gt;" in html
    assert html.count("<td><code>") == 10  # five actual bonds for each context
    assert "41</td>" in html and "45</td>" in html
    assert "No observation returned for this bond" in html
    assert "not_matched" in html and "partial rule" in html
    assert "terminal context padded" in html
    assert "endpoint sentinel" in html
    assert "Target masks and patient-HLA ligand overlays are not included" in html
    assert "0.1, 0.2, 0.3, 0.4, 0.5, 0.0" in html


def test_empty_inventory_does_not_look_assessed(tmp_path):
    json_path, html_path = tmp_path / "none.json", tmp_path / "none.html"
    write_cleavage_profiles([], json_path=json_path, html_path=html_path)
    assert from_native_json(json_path.read_text(), tuple) == ()
    assert "processing is unassessed" in html_path.read_text()


def test_report_renders_the_provenance_fields_mhctools_now_requires(tmp_path):
    """scored_endpoint and the motif strictness grade must reach the reader.

    These are safety-relevant provenance: the endpoint says which
    measurement a native score answers, and the grade says how much a motif
    non-match is worth. The JSON export kept them losslessly while the
    human-readable report silently dropped both.
    """
    quantitative = profile()
    motif = profile_from_cleavage_result(CleavageResult(
        CleavageInput("ACDEFG"), MOTIF_MODEL,
        (CleavageSite(5, "not_matched", "partial rule"),)))
    json_path, html_path = tmp_path / "sites.json", tmp_path / "sites.html"
    write_cleavage_profiles([quantitative, motif], json_path=json_path, html_path=html_path)
    html = html_path.read_text()
    assert MODEL.scored_endpoint in html
    assert MOTIF_MODEL.motif_strictness in html
    assert MOTIF_MODEL.strictness_basis in html
    # The grade is rendered with what it means, not as a bare word.
    assert "a match constrains little" in html


def test_archived_quantitative_model_without_endpoint_is_migrated_or_refused():
    """mhctools 3.44.0 made scored_endpoint required for quantitative models.

    Archives written before that have no such field. vaxrank's own Pepsickle
    model has a known endpoint and is migrated; anything else must fail with
    a message naming the problem instead of a bare validation error.
    """
    from vaxrank.native_serialization import _migrate_state

    pepsickle = {"name": "pepsickle-epitope-mammal", "evidence": "quantitative_model"}
    assert _migrate_state(CleavageModel, pepsickle)["scored_endpoint"] == "site_cleavage"

    with pytest.raises(ValueError, match="cannot be inferred"):
        _migrate_state(CleavageModel, {"name": "third-party", "evidence": "quantitative_model"})

    # Untouched where no migration applies.
    motif = {"name": "whatever", "evidence": "motif_rule"}
    assert _migrate_state(CleavageModel, motif) == motif
    already = {"name": "x", "evidence": "quantitative_model", "scored_endpoint": "substrate_depletion"}
    assert _migrate_state(CleavageModel, already) == already
