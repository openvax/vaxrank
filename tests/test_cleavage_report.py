"""Packaged report retains all bonds, native scores and unknown coverage."""

from dataclasses import replace

from mhctools.cleavage import CleavageInput, CleavageResult, CleavageSite

from vaxrank.cleavage_profile import profile_from_cleavage_result
from vaxrank.cleavage_report import write_cleavage_profiles
from vaxrank.native_serialization import from_native_json

from .test_cleavage_profile import MODEL, profile


def test_complete_profiles_roundtrip_and_render_all_sites(tmp_path):
    p = profile(backend_version="3.41.0", upstream_context=8, downstream_context=8)
    p = replace(p, peptide=replace(p.peptide, source_id="<script>bad()</script>", source_start=40))
    motif = replace(MODEL, evidence="motif_rule", score_name=None, score_units=None)
    partial = profile_from_cleavage_result(CleavageResult(
        CleavageInput("ACDEFG"), motif, (CleavageSite(5, "not_matched", "partial rule"),)))
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
