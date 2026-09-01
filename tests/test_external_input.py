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

"""End-to-end tests covering the LENS / pVACseq → vaccine-design path.

Pre-PR-#253: external-input mode hard-short-circuited and the
construct writers were never reached. These tests pin the new
architecture: a single shared ``ranked_variants_with_vaccine_peptides``
intermediate that *both* the VCF/BAM pipeline and the LENS / pVACseq
loaders produce, then a shared dispatch (``emit_outputs``) into
modality-specific construct writers.
"""

import glob
import json
import os

import pytest

from vaxrank.epitope_io import read_lens_report
from vaxrank.external_input import (
    ExternalConstructOptions,
    ExternalVariantEntry,
    peptide_offsets_in_context,
    parse_lens_variant_coordinates,
    lens_ranking_result,
)

DATA_DIR = os.path.join(
    os.path.dirname(__file__), "data", "epitope_fixtures")


def test_external_variant_entry_defaults_describe_an_empty_result():
    entry = ExternalVariantEntry()

    assert entry.variant is None
    assert entry.vaccine_peptide is None
    assert entry.had_transcript_ids is False
    assert entry.resolved_transcript is False
    assert entry.annotation is None
    assert entry.dna_vaf is None
    assert entry.has_rna_support is False
    assert entry.unparseable is False


def test_parse_lens_variant_coordinates_extracts_chr_pos():
    """The parser extracts only (contig, pos), stripping
    the ``chr`` prefix LENS emits — pyensembl uses bare contigs, and
    keeping the prefix breaks ``variant.effect_on_transcript`` and
    renders the variant short_description as ``chrchr3 …``. Ref/alt
    come from the dedicated per-antigen-source columns
    (snv_ref_allele / snv_alt_allele etc.), looked up by
    ``variant_from_lens_row``."""
    assert parse_lens_variant_coordinates(
        "chr17:7675088:C:T") == ("17", 7675088)
    assert parse_lens_variant_coordinates("chr1:26780312") == ("1", 26780312)
    assert parse_lens_variant_coordinates(
        "chr1:26780312:C") == ("1", 26780312)
    # Already-bare contig pass-through:
    assert parse_lens_variant_coordinates("17:7675088") == ("17", 7675088)


def test_parse_lens_variant_coordinates_returns_none_on_malformed_input():
    for bad in ("garbage", "chr1:notapos:A:T", None, "", "   ",
                "nan", "NaN"):
        assert parse_lens_variant_coordinates(bad) is None


def test_normalize_lens_allele_handles_bracket_form():
    """LENS records alt alleles as bracketed strings: '[T]', '[CA]'.
    The helper strips brackets and returns the inner sequence."""
    from vaxrank.external_input import normalize_lens_allele
    assert normalize_lens_allele("[T]") == "T"
    assert normalize_lens_allele("[CA]") == "CA"
    assert normalize_lens_allele("T") == "T"   # already unbracketed
    assert normalize_lens_allele("") is None
    assert normalize_lens_allele(None) is None
    assert normalize_lens_allele("nan") is None


# ---- resolve_external_transcripts ---------------------------------------

class _StubGenome:
    """Minimal stand-in for pyensembl.EnsemblRelease used by tests.

    ``transcript_by_id(bare_id)`` returns the configured Transcript
    object for that ID, or raises ``ValueError`` (the same shape
    pyensembl 2.x raises for unknown IDs). Lets the resolution helper
    exercise both paths without pulling in a real release.
    """
    def __init__(self, mapping):
        self._mapping = mapping
        self.lookups = []

    def transcript_by_id(self, tid):
        self.lookups.append(tid)
        if tid in self._mapping:
            return self._mapping[tid]
        raise ValueError("Transcript not found: %s" % tid)


def test_resolve_external_transcripts_strips_version_suffix():
    """LENS IDs carry version suffixes (``ENST00000312960.4``); pyensembl
    2.x doesn't strip them. The helper must, or every LENS lookup fails."""
    from vaxrank.external_input import resolve_external_transcripts
    sentinel = object()
    g = _StubGenome({'ENST00000312960': sentinel})
    out = resolve_external_transcripts(['ENST00000312960.4'], g)
    assert out == [sentinel]
    assert g.lookups == ['ENST00000312960']  # bare form was sent


def test_resolve_external_transcripts_drops_unresolvable_ids():
    """A release-mismatch yields some IDs the configured release
    doesn't know. Drop quietly (DEBUG-logged) rather than crash; the
    caller's aggregate WARN summarizes the count."""
    from vaxrank.external_input import resolve_external_transcripts
    known = object()
    g = _StubGenome({'ENST00000000001': known})
    out = resolve_external_transcripts(
        ['ENST00000000001.1', 'ENST99999999999.7', ''],
        g)
    assert out == [known]


def test_infer_genome_build_from_lens_grch38(tmp_path):
    """Hsap38 markers in ``origin_descriptor`` (LENS ERV-row format)
    identify the build as GRCh38."""
    import pandas as pd
    from vaxrank.external_input import infer_genome_build_from_lens
    p = tmp_path / "lens.tsv"
    df = pd.DataFrame({
        'origin_descriptor': [
            'ENSG00000004399.13:PLXND1',
            'Hsap38.chr2.156963765.156964472.-',
        ],
    })
    df.to_csv(p, sep='\t', index=False)
    assert infer_genome_build_from_lens(str(p)) == 'GRCh38'


def test_infer_genome_build_from_lens_grch37(tmp_path):
    import pandas as pd
    from vaxrank.external_input import infer_genome_build_from_lens
    p = tmp_path / "lens.tsv"
    df = pd.DataFrame({
        'origin_descriptor': ['Hsap37.chrX.111111.222222.+'],
    })
    df.to_csv(p, sep='\t', index=False)
    assert infer_genome_build_from_lens(str(p)) == 'GRCh37'


def test_installed_ensembl_releases_walks_pyensembl_cache(
        tmp_path, monkeypatch):
    """``installed_ensembl_releases_for_build`` lists release numbers
    by walking the pyensembl cache directory layout
    (``<cache>/pyensembl/<build>/ensembl<N>``). Important so the
    pre-flight ``--ensembl-release`` hint suggests releases the user
    actually has installed, not arbitrary numbers."""
    import platformdirs
    from vaxrank.external_input import installed_ensembl_releases_for_build

    fake_cache = tmp_path / 'pyensembl'
    (fake_cache / 'GRCh38' / 'ensembl100').mkdir(parents=True)
    (fake_cache / 'GRCh38' / 'ensembl112').mkdir(parents=True)
    (fake_cache / 'GRCh37' / 'ensembl75').mkdir(parents=True)
    (fake_cache / 'GRCh38' / 'not_an_ensembl_dir').mkdir(parents=True)

    monkeypatch.setattr(
        platformdirs, 'user_cache_dir',
        lambda app: str(fake_cache) if app == 'pyensembl' else '')

    assert installed_ensembl_releases_for_build('GRCh38') == [100, 112]
    assert installed_ensembl_releases_for_build('GRCh37') == [75]
    assert installed_ensembl_releases_for_build('Unknown_Build_42') == []


def test_infer_genome_build_from_lens_unknown_returns_none(tmp_path):
    """No Hsap markers (SNV-only file) → None; caller falls back to
    a generic warning rather than guessing wrong."""
    import pandas as pd
    from vaxrank.external_input import infer_genome_build_from_lens
    p = tmp_path / "lens.tsv"
    df = pd.DataFrame({
        'origin_descriptor': [
            'ENSG00000004399.13:PLXND1',
            'ENSG00000005339.15:CREBBP',
        ],
    })
    df.to_csv(p, sep='\t', index=False)
    assert infer_genome_build_from_lens(str(p)) is None


def test_resolve_external_transcripts_returns_empty_without_genome():
    """No --ensembl-release set → genome=None; resolution is a no-op
    and downstream code falls back to the empty-transcript path."""
    from vaxrank.external_input import resolve_external_transcripts
    assert resolve_external_transcripts(['ENST00000000001'], None) == []
    assert resolve_external_transcripts([], None) == []


def test_variant_from_lens_row_uses_real_snv_alleles():
    """SNV rows: ref/alt come from snv_ref_allele / snv_alt_allele.
    No placeholder genotype. The synthesized Variant carries the real
    biology so it can be fed to varcode-effect annotation safely."""
    from vaxrank.external_input import variant_from_lens_row
    row = {
        'variant': 'chr1:1624824',
        'antigen_source': 'SNV',
        'snv_ref_allele': 'C',
        'snv_alt_allele': '[T]',  # LENS's bracketed format
    }
    v = variant_from_lens_row(row)
    assert v is not None
    # ``chr`` prefix stripped during parse so pyensembl-style lookups
    # work; see ``parse_lens_variant_coordinates``.
    assert v.contig == "1"
    assert v.start == 1624824
    assert v.ref == "C"
    assert v.alt == "T"


def test_variant_from_lens_row_uses_real_indel_alleles():
    """INDEL rows: ref/alt come from indel_ref_allele / indel_alt_allele.
    varcode normalizes indels (strips shared prefix), so ``CA → C``
    becomes ref='A' alt='' at start+1 — that's the canonical
    representation a downstream caller would expect."""
    from vaxrank.external_input import variant_from_lens_row
    row = {
        'variant': 'chr3:150742445',
        'antigen_source': 'INDEL',
        'indel_ref_allele': 'CA',
        'indel_alt_allele': '[C]',
    }
    v = variant_from_lens_row(row)
    assert v is not None
    # varcode's canonical indel representation: shared prefix
    # stripped, position advanced to the differing base. Contig
    # also has its ``chr`` prefix stripped for pyensembl compatibility.
    assert v.contig == "3"
    assert v.start == 150742446
    assert v.ref == "A"
    assert v.alt == ""


def test_variant_from_lens_row_skips_non_snv_indel():
    """SPLICE / FUSION / CTA-SELF / ERV rows don't have variant_coords
    populated (NaN); the row-level helper returns None for these,
    and the caller skips them earlier on the empty-coords path."""
    from vaxrank.external_input import variant_from_lens_row
    for src in ('SPLICE', 'FUSION', 'CTA/SELF', 'ERV'):
        row = {
            'variant': None,
            'antigen_source': src,
        }
        assert variant_from_lens_row(row) is None


def test_variant_from_lens_row_skips_when_alleles_missing():
    """If the SNV row's allele columns are missing/NaN, return None
    rather than fabricate a placeholder."""
    from vaxrank.external_input import variant_from_lens_row
    row = {
        'variant': 'chr1:1000',
        'antigen_source': 'SNV',
        'snv_ref_allele': None,
        'snv_alt_allele': None,
    }
    assert variant_from_lens_row(row) is None


def test_real_lens_v19_subset_produces_ranked_entries():
    """Regression test for the real-LENS-v1.9 path. The subset fixture
    came from a Hugo IPRES Pt10 dump and tripped both #259 (variant_coords
    parser) and #260 (duplicate-row handling). Pin that the path now
    produces a non-empty ranked list end-to-end."""
    path = os.path.join(
        DATA_DIR, "real_lens_subsets", "lens_v1.9_real_subset.tsv")
    _loaded = read_lens_report(path)
    predictions = list(_loaded.epitopes)
    _ranking = lens_ranking_result(_loaded, predictions)
    ranked, _dna_vaf = _ranking.ranked, _ranking.dna_vaf_by_variant
    # The fixture has at least one parseable variant_coords row.
    # Pre-fix: 0. Now: > 0.
    assert len(ranked) > 0, (
        "Expected real LENS v1.9 fixture to produce a non-empty ranked "
        "list after the variant_coords parser fix")


# ---- Parametrized end-to-end coverage of every real LENS fixture --------

_REAL_LENS_FIXTURES = sorted(glob.glob(
    os.path.join(DATA_DIR, "real_lens_subsets", "*.tsv")))


@pytest.mark.parametrize("fixture_path", _REAL_LENS_FIXTURES)
def test_real_lens_fixture_runs_end_to_end_via_cli(fixture_path, tmp_path):
    """Drive the full CLI ``main()`` against every LENS fixture in
    ``tests/data/epitope_fixtures/real_lens_subsets/``. Any new LENS
    field-format quirk we encounter in the wild gets added as a
    fixture in that directory and this test surfaces breakage
    automatically — no more "shouldn't have to run these manually."
    """
    from vaxrank.cli.entry_point import main
    csv_path = tmp_path / "out.csv"
    xlsx_path = tmp_path / "out.xlsx"
    main([
        "--input-lens", fixture_path,
        "--output-csv", str(csv_path),
        "--output-neoepitope-report", str(xlsx_path),
    ])
    # Both outputs landed; the path didn't crash on any LENS quirk.
    assert csv_path.exists(), f"CSV not written for {fixture_path}"
    assert xlsx_path.exists(), f"XLSX not written for {fixture_path}"


def test_cli_applies_external_score_config_before_ranking(
        tmp_path, monkeypatch):
    """Configured evidence-only scores drive ranking, not just CSV output."""
    from vaxrank.cli import entry_point

    lens_path = tmp_path / "processing-only.tsv"
    lens_path.write_text(
        "peptide\tallele\tpep_context\tmhcflurry_2.1.1.aff\t"
        "mhcflurry_2.1.1.proc_score\tantigen_source\tvariant_coords\t"
        "snv_ref_allele\tsnv_alt_allele\tgene_name\n"
        "SIINFEKL\tHLA-A02:01\tAASIINFEKLLL\t\t0.8\tSNV\t"
        "chr1:1000\tA\t[T]\tGENE1\n")
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        "epitopes:\n"
        "  score_expr: processing[mhcflurry].score\n")
    csv_path = tmp_path / "neoepitopes.csv"
    captured = {}

    def capture_ranked(_args, ranked, source):
        captured["ranked"] = ranked
        captured["source"] = source

    monkeypatch.setattr(entry_point, "emit_outputs", capture_ranked)
    entry_point.main([
        "--input-lens", str(lens_path),
        "--config", str(config_path),
        "--output-csv", str(csv_path),
        "--no-processing-aware-annotation",
    ])

    assert captured["source"] == "external"
    vaccine_peptide = captured["ranked"][0][1][0]
    assert vaccine_peptide.target_epitope_score == pytest.approx(0.8)
    assert vaccine_peptide.target_epitopes[0].per_allele_scores == {
        "HLA-A*02:01": pytest.approx(0.8)}


def test_external_ranking_applies_default_minimum_epitope_score(tmp_path):
    """Zero-score evidence stays auditable but cannot become a target."""
    from types import SimpleNamespace

    from vaxrank.external_input import (
        LENS_PROVENANCE_MARKER,
        load_external_ranked,
    )

    lens_path = tmp_path / "processing-only.tsv"
    lens_path.write_text(
        "peptide\tallele\tpep_context\tmhcflurry_2.1.1.aff\t"
        "mhcflurry_2.1.1.proc_score\tantigen_source\tvariant_coords\t"
        "snv_ref_allele\tsnv_alt_allele\tgene_name\n"
        "SIINFEKL\tHLA-A02:01\tAASIINFEKLLL\t\t0.8\tSNV\t"
        "chr1:1000\tA\t[T]\tGENE1\n")
    args = SimpleNamespace(
        input_lens=str(lens_path),
        input_pvacseq=None,
        output_patient_id="",
        vaccine_peptide_length=25,
        genome=None,
    )

    ranked, report, loaded, patient_info, _vafs = load_external_ranked(args)

    assert ranked == []
    assert len(report) == 1
    assert loaded[0].per_allele_scores == {"HLA-A*02:01": 0.0}
    assert patient_info.mhc_alleles == [
        "HLA-A*02:01", LENS_PROVENANCE_MARKER]


def test_external_filter_prunes_ranked_groups_but_preserves_patient_genotype(
        tmp_path):
    """Topiary filtering controls targets without redefining the patient."""
    from types import SimpleNamespace

    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.external_input import (
        LENS_PROVENANCE_MARKER,
        load_external_ranked,
    )

    lens_path = tmp_path / "partially-filtered.tsv"
    lens_path.write_text(
        "peptide\tallele\tpep_context\tmhcflurry_2.1.1.aff\t"
        "mhcflurry_2.1.1.proc_score\tantigen_source\tvariant_coords\t"
        "snv_ref_allele\t"
        "snv_alt_allele\tgene_name\n"
        "SIINFEKL\tHLA-A02:01\tAASIINFEKLLL\t50\t0.7\tSNV\t"
        "chr1:1000\tA\t[T]\tGENE1\n"
        "SIINFEKL\tHLA-B07:02\tAASIINFEKLLL\t500\t0.7\tSNV\t"
        "chr1:1000\tA\t[T]\tGENE1\n")
    args = SimpleNamespace(
        input_lens=str(lens_path),
        input_pvacseq=None,
        output_patient_id="",
        vaccine_peptide_length=25,
        genome=None,
    )

    ranked, _report, loaded, patient_info, _vafs = load_external_ranked(
        args,
        epitope_config=EpitopeConfig(filter_expr="affinity <= 100"),
    )

    # Genotype comes from the complete input membership, before filtering.
    assert loaded[0].patient_alleles == (
        "HLA-A*02:01", "HLA-B*07:02")
    assert patient_info.mhc_alleles == [
        "HLA-A*02:01", "HLA-B*07:02", LENS_PROVENANCE_MARKER]

    # The filtered B*07:02 group cannot enter ranking, coverage, or design.
    target = ranked[0][1][0].target_epitopes[0]
    assert target.patient_alleles == ("HLA-A*02:01",)
    assert target.per_allele_scores.keys() == {"HLA-A*02:01"}
    assert {
        p.allele for p in target.predictions_flat() if p.allele
    } == {"HLA-A*02:01"}
    processing = target.predictions_for(
        "antigen_processing", predictor="mhcflurry")
    assert len(processing) == 1
    assert processing[0].allele == ""


def test_external_minimum_score_prunes_all_ranking_state_for_rejected_allele(
        tmp_path):
    """A below-threshold allele cannot survive through the score map."""
    from types import SimpleNamespace

    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.external_input import load_external_ranked

    lens_path = tmp_path / "mixed-presentation.tsv"
    lens_path.write_text(
        "peptide\tallele\tpep_context\tmhcflurry_2.1.1.aff\t"
        "mhcflurry_2.1.1.pres_score\tantigen_source\tvariant_coords\t"
        "snv_ref_allele\tsnv_alt_allele\tgene_name\n"
        "SIINFEKL\tHLA-A02:01\tAASIINFEKLLL\t\t0.8\tSNV\t"
        "chr1:1000\tA\t[T]\tGENE1\n"
        "SIINFEKL\tHLA-B07:02\tAASIINFEKLLL\t\t0.2\tSNV\t"
        "chr1:1000\tA\t[T]\tGENE1\n")
    args = SimpleNamespace(
        input_lens=str(lens_path),
        input_pvacseq=None,
        output_patient_id="",
        vaccine_peptide_length=25,
        genome=None,
    )
    config = EpitopeConfig(
        score_expr="presentation[mhcflurry].score",
        min_epitope_score=0.5,
    )

    ranked, _report, loaded, _patient_info, _vafs = load_external_ranked(
        args, epitope_config=config)

    assert loaded[0].per_allele_scores == {
        "HLA-A*02:01": pytest.approx(0.8),
        "HLA-B*07:02": pytest.approx(0.2),
    }
    target = ranked[0][1][0].target_epitopes[0]
    assert target.patient_alleles == ("HLA-A*02:01",)
    assert target.per_allele_scores == {
        "HLA-A*02:01": pytest.approx(0.8)}
    assert target.epitope_score == pytest.approx(0.8)


def test_external_filter_excluding_every_group_produces_no_ranked_vaccine(
        tmp_path, monkeypatch):
    """A hard Topiary filter cannot leave zero-scored construct targets."""
    from types import SimpleNamespace

    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank import external_input
    from vaxrank.external_input import (
        LENS_PROVENANCE_MARKER,
        load_external_ranked,
    )

    monkeypatch.setattr(
        external_input,
        "resolve_external_transcripts",
        lambda transcript_ids, _genome: [object()] if transcript_ids else [],
    )

    lens_path = tmp_path / "fully-filtered.tsv"
    lens_path.write_text(
        "peptide\tallele\tpep_context\tmhcflurry_2.1.1.aff\t"
        "antigen_source\tvariant_coords\tsnv_ref_allele\t"
        "snv_alt_allele\tgene_name\ttranscript_id\t"
        "rna_reads_covering_genomic_origin\t"
        "rna_reads_covering_genomic_origin_with_peptide_cds\n"
        "SIINFEKL\tHLA-A02:01\tAASIINFEKLLL\t50\tSNV\t"
        "chr1:1000\tA\t[T]\tGENE1\tENST1\t10\t3\n")
    args = SimpleNamespace(
        input_lens=str(lens_path),
        input_pvacseq=None,
        output_patient_id="",
        vaccine_peptide_length=25,
        genome=None,
    )

    ranked, _report, loaded, patient_info, _vafs = load_external_ranked(
        args,
        epitope_config=EpitopeConfig(filter_expr="affinity <= 0"),
    )

    assert ranked == []
    assert loaded[0].per_allele_scores == {}
    assert patient_info.mhc_alleles == [
        "HLA-A*02:01", LENS_PROVENANCE_MARKER]
    assert patient_info.num_somatic_variants == 1
    assert patient_info.num_coding_effect_variants == 1
    assert patient_info.num_variants_with_rna_support == 1
    assert patient_info.num_variants_with_vaccine_peptides == 0


def test_external_filter_does_not_cross_lens_source_contexts(tmp_path):
    """A shared sequence cannot transfer eligibility between variants."""
    from types import SimpleNamespace

    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.external_input import load_external_ranked

    lens_path = tmp_path / "shared-peptide.tsv"
    lens_path.write_text(
        "peptide\tallele\tpep_context\tmhcflurry_2.1.1.aff\t"
        "antigen_source\tvariant_coords\tsnv_ref_allele\t"
        "snv_alt_allele\tgene_name\n"
        "SIINFEKL\tHLA-A02:01\tAASIINFEKLLL\t50\tSNV\t"
        "chr1:1000\tA\t[T]\tGENE1\n"
        "SIINFEKL\tHLA-A02:01\tGGSIINFEKLTT\t500\tSNV\t"
        "chr1:2000\tG\t[C]\tGENE2\n")
    args = SimpleNamespace(
        input_lens=str(lens_path),
        input_pvacseq=None,
        output_patient_id="",
        vaccine_peptide_length=25,
        genome=None,
    )

    ranked, report, loaded, _patient_info, _vafs = load_external_ranked(
        args,
        epitope_config=EpitopeConfig(filter_expr="affinity <= 100"),
    )

    assert len(report) == 2
    assert len(loaded) == 2
    assert len(ranked) == 1
    variant, vaccine_peptides = ranked[0]
    assert variant.start == 1000
    assert vaccine_peptides[0].target_epitopes[0].source_sequence == (
        "AASIINFEKLLL")


def test_external_filter_does_not_cross_identical_lens_contexts(tmp_path):
    """Variant provenance, not sequence context, controls eligibility."""
    from types import SimpleNamespace

    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.external_input import load_external_ranked

    lens_path = tmp_path / "identical-context.tsv"
    lens_path.write_text(
        "peptide\tallele\tpep_context\tmhcflurry_2.1.1.aff\t"
        "antigen_source\tvariant_coords\tsnv_ref_allele\t"
        "snv_alt_allele\tgene_name\tgene_id\ttranscript_id\n"
        "SIINFEKL\tHLA-A02:01\tAASIINFEKLLL\t50\tSNV\t"
        "chr1:1000\tA\t[T]\tGENE1\tENSG1\tENST1\n"
        "SIINFEKL\tHLA-A02:01\tAASIINFEKLLL\t500\tSNV\t"
        "chr1:2000\tG\t[C]\tGENE2\tENSG2\tENST2\n")
    args = SimpleNamespace(
        input_lens=str(lens_path),
        input_pvacseq=None,
        output_patient_id="",
        vaccine_peptide_length=25,
        genome=None,
    )

    ranked, report, loaded, _patient_info, _vafs = load_external_ranked(
        args,
        epitope_config=EpitopeConfig(filter_expr="affinity <= 100"),
    )

    assert len(report) == 2
    assert len(loaded) == 2
    assert loaded[0].prediction_id != loaded[1].prediction_id
    assert len(ranked) == 1
    assert ranked[0][0].start == 1000


def test_lens_representative_uses_configured_dsl_score(tmp_path):
    """The selected construct context is the highest DSL-scored context."""
    from types import SimpleNamespace

    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.external_input import load_external_ranked

    lens_path = tmp_path / "presentation-ranked.tsv"
    lens_path.write_text(
        "peptide\tallele\tpep_context\tmhcflurry_2.1.1.aff\t"
        "mhcflurry_2.1.1.pres_score\tantigen_source\tvariant_coords\t"
        "snv_ref_allele\tsnv_alt_allele\tgene_name\ttranscript_id\n"
        "SIINFEKL\tHLA-A02:01\tAASIINFEKLLL\t50\t0.1\tSNV\t"
        "chr1:1000\tA\t[T]\tGENE1\tENST1\n"
        "LATEKSRWS\tHLA-A02:01\tGGLATEKSRWSTT\t500\t0.9\tSNV\t"
        "chr1:1000\tA\t[T]\tGENE1\tENST1\n")
    args = SimpleNamespace(
        input_lens=str(lens_path),
        input_pvacseq=None,
        output_patient_id="",
        vaccine_peptide_length=25,
        genome=None,
    )

    ranked, _report, _loaded, _patient_info, _vafs = load_external_ranked(
        args,
        epitope_config=EpitopeConfig(
            score_expr="presentation[mhcflurry].score"),
    )

    vaccine_peptide = ranked[0][1][0]
    assert vaccine_peptide.amino_acids == "GGLATEKSRWSTT"
    assert [e.sequence for e in vaccine_peptide.target_epitopes] == [
        "LATEKSRWS"]
    assert vaccine_peptide.target_epitope_score == pytest.approx(0.9)


def test_lens_input_summary_reduces_all_rows_before_construct_checks(
        tmp_path, monkeypatch):
    """Input facts survive an unconstructable affinity-winning row."""
    from types import SimpleNamespace

    from vaxrank import external_input
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.external_input import load_external_ranked

    monkeypatch.setattr(
        external_input,
        "resolve_external_transcripts",
        lambda transcript_ids, _genome: [object()] if transcript_ids else [],
    )
    lens_path = tmp_path / "all-row-summary.tsv"
    lens_path.write_text(
        "peptide\tallele\tpep_context\tmhcflurry_2.1.1.aff\t"
        "antigen_source\tvariant_coords\tsnv_ref_allele\t"
        "snv_alt_allele\tgene_name\ttranscript_id\t"
        "rna_reads_covering_genomic_origin\t"
        "rna_reads_covering_genomic_origin_with_peptide_cds\n"
        "SIINFEKL\tHLA-A02:01\tXXSIINFEKLXX\t50\tSNV\t"
        "chr1:1000\tA\t[T]\tGENE1\t\t0\t0\n"
        "LATEKSRWS\tHLA-A02:01\tGGLATEKSRWSTT\t500\tSNV\t"
        "chr1:1000\tA\t[T]\tGENE1\tENST1\t10\t3\n")
    args = SimpleNamespace(
        input_lens=str(lens_path),
        input_pvacseq=None,
        output_patient_id="",
        vaccine_peptide_length=25,
        genome=None,
    )

    _ranked, _report, _loaded, patient_info, _vafs = load_external_ranked(
        args,
        epitope_config=EpitopeConfig(filter_expr="affinity <= 0"),
    )

    assert patient_info.num_somatic_variants == 1
    assert patient_info.num_coding_effect_variants == 1
    assert patient_info.num_variants_with_rna_support == 1
    assert patient_info.num_variants_with_vaccine_peptides == 0


def test_pvacseq_filter_excluding_every_group_produces_no_ranked_vaccine():
    """The shared external filtering contract also applies to pVACseq."""
    from types import SimpleNamespace

    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.external_input import load_external_ranked

    path = os.path.join(DATA_DIR, "pvacseq_example.tsv")
    args = SimpleNamespace(
        input_lens=None,
        input_pvacseq=path,
        output_patient_id="",
        vaccine_peptide_length=25,
        genome=None,
    )

    ranked, _report, loaded, patient_info, _vafs = load_external_ranked(
        args,
        epitope_config=EpitopeConfig(filter_expr="affinity <= 0"),
    )

    assert ranked == []
    assert loaded
    assert all(epitope.per_allele_scores == {} for epitope in loaded)
    assert set(patient_info.mhc_alleles[:-1]) == {
        allele
        for epitope in loaded
        for allele in epitope.patient_alleles
    }
    assert patient_info.num_somatic_variants == 3
    assert patient_info.num_coding_effect_variants == 0
    assert patient_info.num_variants_with_rna_support == 3
    assert patient_info.num_variants_with_vaccine_peptides == 0


def write_multi_peptide_pvacseq_fixture(path):
    """Write one pVACseq variant with score axes that disagree."""
    header = (
        "ID\tIndex\tA*02:01\tGene\tAA Change\tNum Passing Transcripts\t"
        "Best Peptide\tBest Transcript\tAllele\tIC50 MT\tIC50 WT\t"
        "%ile MT\t%ile WT\tRNA Expr\tRNA VAF\tRNA Depth\tDNA VAF\t"
        "Tier\tRef Match\tEvaluation\n")
    path.write_text(
        header
        + "chr1-1000-1001-A-T\t1.GENE1.ENST1.missense.1A/T\t1\t"
        "GENE1\tA1T\t1\tSIINFEKL\t\tHLA-A*02:01\t50\t5000\t"
        "0.1\t10\t1\t0\t0\t0.4\tPass\tFalse\tPending\n"
        + "chr1-1000-1001-A-T\t1.GENE1.ENST2.missense.1A/T\t1\t"
        "GENE1\tA1T\t1\tLATEKSRWS\tENST2\tHLA-A*02:01\t500\t5000\t"
        "1.0\t10\t1\t0.2\t100\t0.4\tPass\tFalse\tPending\n")
    return path


def test_pvacseq_representative_uses_configured_dsl_score(tmp_path):
    """pVACseq construct selection consumes the same score as ranking."""
    from types import SimpleNamespace

    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.external_input import load_external_ranked

    path = write_multi_peptide_pvacseq_fixture(
        tmp_path / "multi-peptide-pvacseq.tsv")
    args = SimpleNamespace(
        input_lens=None,
        input_pvacseq=str(path),
        output_patient_id="",
        vaccine_peptide_length=25,
        genome=None,
    )

    ranked, _report, _loaded, _patient_info, _vafs = load_external_ranked(
        args,
        epitope_config=EpitopeConfig(score_expr="affinity.value"),
    )

    vaccine_peptide = ranked[0][1][0]
    assert vaccine_peptide.amino_acids == "LATEKSRWS"
    assert [e.sequence for e in vaccine_peptide.target_epitopes] == [
        "LATEKSRWS"]
    assert vaccine_peptide.target_epitope_score == pytest.approx(500.0)


def test_pvacseq_input_summary_reduces_all_rows_before_selection(
        tmp_path, monkeypatch):
    """pVACseq RNA/transcript facts are unions over the raw variant rows."""
    from types import SimpleNamespace

    from vaxrank import external_input
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.external_input import load_external_ranked

    monkeypatch.setattr(
        external_input,
        "resolve_external_transcripts",
        lambda transcript_ids, _genome: [object()] if transcript_ids else [],
    )
    path = write_multi_peptide_pvacseq_fixture(
        tmp_path / "multi-row-summary-pvacseq.tsv")
    args = SimpleNamespace(
        input_lens=None,
        input_pvacseq=str(path),
        output_patient_id="",
        vaccine_peptide_length=25,
        genome=None,
    )

    _ranked, _report, _loaded, patient_info, _vafs = load_external_ranked(
        args,
        epitope_config=EpitopeConfig(filter_expr="affinity <= 0"),
    )

    assert patient_info.num_somatic_variants == 1
    assert patient_info.num_coding_effect_variants == 1
    assert patient_info.num_variants_with_rna_support == 1
    assert patient_info.num_variants_with_vaccine_peptides == 0


def test_real_lens_fixtures_present():
    """Pin that the fixture directory has the expected coverage:
    one fixture per known LENS version + one with stop-codon
    edge cases. Adding a new file in real_lens_subsets/ implicitly
    grows the parametrized end-to-end test above."""
    names = sorted(os.path.basename(p) for p in _REAL_LENS_FIXTURES)
    assert "lens_v1.4_real_subset.tsv" in names
    assert "lens_v1.5_real_subset.tsv" in names
    assert "lens_v1.9_real_subset.tsv" in names
    assert "lens_v1.9_with_stop_codons.tsv" in names, (
        "Expected the stop-codon edge-case fixture; this fixture "
        "exercises rows where pep_context contains '*' (stop codon "
        "mid-context, e.g. stop-loss / readthrough variants).")


def test_lens_pep_context_with_stop_codon_truncates():
    """Real LENS files emit ``*`` in pep_context for stop-loss /
    readthrough variants. Pre-fix this crashed manufacturability
    scoring with KeyError: '*'. Now: truncate at first stop, drop
    the row if the neoepitope itself was past the stop."""
    path = os.path.join(
        DATA_DIR, "real_lens_subsets",
        "lens_v1.9_with_stop_codons.tsv")
    _loaded = read_lens_report(path)
    predictions = list(_loaded.epitopes)
    _ranking = lens_ranking_result(_loaded, predictions)
    ranked, _dna_vaf = _ranking.ranked, _ranking.dna_vaf_by_variant
    # No fragment carries a '*' or a non-standard residue.
    for _, peptides in ranked:
        for vp in peptides:
            assert '*' not in vp.mutant_protein_fragment.amino_acids


def test_truncate_at_stop_codon_helper():
    """Translation stops at the first ``*`` — any AA after is
    non-existent in the cell."""
    from vaxrank.vaccine_library import truncate_at_stop_codon
    assert truncate_at_stop_codon("AAVK*GTRPL") == "AAVK"
    assert truncate_at_stop_codon("KLQGHSAPVL") == "KLQGHSAPVL"  # no stop
    assert truncate_at_stop_codon("") == ""
    assert truncate_at_stop_codon("*") == ""  # all stop
    assert truncate_at_stop_codon("AB*CD*EF") == "AB"  # only first stop


def test_has_only_standard_amino_acids_helper():
    """The 20 canonical AAs pass; non-standard residues (selenocysteine
    U, pyrrolysine O, ambiguous X / B / Z / J, stop *) fail."""
    from vaxrank import has_only_standard_amino_acids
    assert has_only_standard_amino_acids("KLQGHSAPVL")
    assert has_only_standard_amino_acids("")
    for non_standard in ("U", "O", "X", "B", "Z", "J", "*"):
        assert not has_only_standard_amino_acids(
            "KLQ" + non_standard + "VL"), \
            f"Should reject {non_standard!r}"


def test_peptide_offsets_in_context_finds_peptide():
    # AASVVGSSSSSGTR contains SVVGSSSSS at offset 2, length 9
    start, end = peptide_offsets_in_context(
        "SVVGSSSSS", "AASVVGSSSSSGTR")
    assert start == 2
    assert end == 11


def test_peptide_offsets_in_context_returns_none_when_not_found():
    """When the peptide isn't a substring of pep_context, return
    ``(None, None)`` so the caller drops the row instead of falsely
    claiming the entire context is the mutation span."""
    start, end = peptide_offsets_in_context("XYZXYZ", "AASVVGSSSSSGTR")
    assert start is None
    assert end is None
    # Empty inputs also return None
    assert peptide_offsets_in_context("", "AAVK") == (None, None)
    assert peptide_offsets_in_context("AAVK", "") == (None, None)


def test_lens_picks_strongest_binder_when_multiple_rows_per_variant():
    """Multi-row LENS fixture: one variant has three peptides at very
    different IC50s. The representative pick must be the strongest
    binder (lowest IC50), not the file-order first row.

    Pre-fix: ``_pick_representative`` looked for non-existent columns
    (``ic50``, ``percentile_rank``) and tied every row at (2, 0.0),
    so it returned the alphabetically-first allele's row regardless
    of binding strength. This test pins the fix.
    """
    path = os.path.join(DATA_DIR, "lens_multi_row_per_variant.tsv")
    _loaded = read_lens_report(path)
    predictions = list(_loaded.epitopes)
    _ranking = lens_ranking_result(_loaded, predictions)
    ranked, _dna_vaf = _ranking.ranked, _ranking.dna_vaf_by_variant
    assert len(ranked) == 1, (
        "Fixture has one variant with 3 peptides; got %d entries" % len(ranked))
    _, peptides = ranked[0]
    fragment = peptides[0].mutant_protein_fragment
    # Strongest binder is STRNGLVLLL (mhcflurry IC50 = 18.0). Its
    # pep_context is "AASTRNGLVLLLGTR".
    assert "STRNGLVLLL" in fragment.amino_acids, (
        "Expected pep_context to come from the strongest binder row "
        "(STRNGLVLLL, IC50=18); got %r" % fragment.amino_acids)


# ---- patient_info_from_external -----------------------------------------

def test_patient_info_from_external_proxy_counts():
    """The synthesized PatientInfo's variant counts come from the
    ranked output, not from VCF/BAM. Pin each count's definition so a
    later loader change can't silently shift the rendered template
    report headers."""
    from vaxrank.external_input import patient_info_from_external

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    _loaded = read_lens_report(path)
    predictions = list(_loaded.epitopes)
    _ranking = lens_ranking_result(_loaded, predictions)
    ranked, _dna_vaf = _ranking.ranked, _ranking.dna_vaf_by_variant
    info = patient_info_from_external(
        ranked, path, 'Pt-X', _ranking.input_summary,
        input_label='LENS report')
    assert info.patient_id == 'Pt-X'
    # PatientInfo on the external path no longer overloads
    # ``vcf_paths`` (LENS files aren't VCFs); the explicit
    # ``inputs`` list carries the path with an accurate label.
    assert info.vcf_paths == []
    assert info.inputs == [('LENS report', path)]
    assert info.bam_path is None
    # Input composition comes from the parse, not from the ranked output:
    # every variant the file produced antigens for, before filtering.
    assert info.num_somatic_variants == (
        _ranking.input_summary.num_somatic_variants)
    # Without --ensembl-release the test fixture won't resolve
    # transcripts; coding-effect count drops to 0.
    assert info.num_coding_effect_variants == 0
    # ``num_variants_with_vaccine_peptides`` matches the populated
    # ranked entries (any with vps).
    assert info.num_variants_with_vaccine_peptides == sum(
        1 for _, vps in ranked if vps)


def test_patient_info_from_external_empty_ranked():
    """No ranked variants → all counts zero, no crash."""
    from vaxrank.external_input import patient_info_from_external
    from vaxrank.external_input import ExternalInputSummary
    info = patient_info_from_external(
        [], '/tmp/empty.tsv', '', ExternalInputSummary())
    assert info.num_somatic_variants == 0
    assert info.num_coding_effect_variants == 0
    assert info.num_variants_with_rna_support == 0
    assert info.num_variants_with_vaccine_peptides == 0


# ---- predicted_effect None tolerance ------------------------------------

def test_predicted_effect_returns_none_with_no_transcripts():
    """``MutantProteinFragment.predicted_effect`` returns None — not
    raises — when no transcripts resolved. Template renderer relies
    on this for ERV / non-genic antigens and release-mismatch rows."""
    from varcode import Variant
    from vaxrank.mutant_protein_fragment import MutantProteinFragment
    frag = MutantProteinFragment(
        variant=Variant(contig='1', start=100, ref='A', alt='G'),
        gene_name='TEST',
        amino_acids='MASSEQUENCE',
        mutant_amino_acid_start_offset=0,
        mutant_amino_acid_end_offset=11,
        supporting_reference_transcripts=[],
        n_overlapping_reads=0, n_alt_reads=0, n_ref_reads=0,
        n_alt_reads_supporting_protein_sequence=0,
    )
    assert frag.predicted_effect() is None


def test_predicted_effect_returns_none_when_varcode_rejects_all_transcripts():
    """``ReferenceMismatchError`` from varcode means LENS / pVACseq
    was called against a different reference than the configured
    pyensembl release. predicted_effect() must return None for this
    case so report rendering keeps going (with "—" placeholders)."""
    from varcode import Variant
    from varcode.errors import ReferenceMismatchError
    from vaxrank.mutant_protein_fragment import MutantProteinFragment

    real_variant = Variant(contig='1', start=100, ref='A', alt='G')

    class _StubVariant:
        contig = '1'
        start = 100
        ref = 'A'
        alt = 'G'
        def __str__(self):
            return 'StubVariant'
        def effect_on_transcript(self, t):
            raise ReferenceMismatchError(
                variant=real_variant,
                transcript=t,
                expected_ref='A',
                observed_ref='C')

    class _StubTranscript:
        transcript_id = 'ENST_STUB'

    frag = MutantProteinFragment(
        variant=_StubVariant(),
        gene_name='TEST',
        amino_acids='MASSEQUENCE',
        mutant_amino_acid_start_offset=0,
        mutant_amino_acid_end_offset=11,
        supporting_reference_transcripts=[_StubTranscript()],
        n_overlapping_reads=0, n_alt_reads=0, n_ref_reads=0,
        n_alt_reads_supporting_protein_sequence=0,
    )
    assert frag.predicted_effect() is None


def test_lens_to_ranked_vaccine_peptides_round_trip():
    """Loading the LENS fixture and converting to ranked vaccine
    peptides should produce one entry per unique variant, each
    carrying the pep_context as the antigen fragment."""
    path = os.path.join(DATA_DIR, "lens_example.tsv")
    _loaded = read_lens_report(path)
    epitopes = list(_loaded.epitopes)
    _ranking = lens_ranking_result(_loaded, epitopes)
    ranked, _dna_vaf = _ranking.ranked, _ranking.dna_vaf_by_variant
    assert len(ranked) == 3, (
        "LENS fixture has 3 unique variants; got %d" % len(ranked))
    # Each ranked entry: (Variant, [VaccinePeptide])
    for variant, peptides in ranked:
        assert len(peptides) == 1
        vp = peptides[0]
        # Fragment AA = pep_context (the SLP-style window from LENS)
        assert vp.mutant_protein_fragment.amino_acids
        # All three fixture rows have pep_context longer than peptide
        assert (len(vp.mutant_protein_fragment.amino_acids)
                > len(vp.epitopes[0].sequence))
        # Mutation span is non-empty inside the fragment
        assert (vp.mutant_protein_fragment.mutant_amino_acid_end_offset
                > vp.mutant_protein_fragment.mutant_amino_acid_start_offset)


def test_lens_drives_mrna_construct_assembly_end_to_end():
    """LENS file → ranked vaccine peptides → assemble_mrna_constructs
    must produce a non-empty construct. Pre-#253 the construct
    writers were unreachable from external input — this test pins
    that the path is now wired."""
    from vaxrank.mrna import RNAConstructConfig, assemble_mrna_constructs

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    _loaded = read_lens_report(path)
    predictions = list(_loaded.epitopes)
    _ranking = lens_ranking_result(_loaded, predictions)
    ranked, _dna_vaf = _ranking.ranked, _ranking.dna_vaf_by_variant
    assert ranked, "Expected non-empty ranked list from LENS fixture"

    options = RNAConstructConfig(
        signal_peptide=None, include_mitd=False,
        utr_3p='HBB', poly_a_length=10,
        antigens_per_construct=3, max_constructs=1,
        max_antigen_length_aa=14,
        optimize_linkers=False,
    )
    constructs = assemble_mrna_constructs(ranked, options=options)
    assert len(constructs) == 1
    c = constructs[0]
    # Three antigens packed into one construct
    assert len(c.antigen_names) == 3
    # cds_aa contains all three antigen AA strings (or windowed
    # subwindows centered on the mutation)
    for vp_pair, name in zip(ranked, c.antigen_names):
        # The antigen name carries the gene_name; LENS fixture genes
        # are TP53, BRAF, KRAS.
        assert any(g in name for g in ("TP53", "BRAF", "KRAS")), name
    assert c.cds_nt.endswith("TAA")
    assert c.full_nt.endswith("A" * 10)


def test_lens_drives_peptide_construct_assembly_end_to_end():
    """LENS file → peptide-vaccine constructs (SLP mode)."""
    from vaxrank.peptide import (
        PeptideConstructConfig, assemble_peptide_constructs)

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    _loaded = read_lens_report(path)
    predictions = list(_loaded.epitopes)
    _ranking = lens_ranking_result(_loaded, predictions)
    ranked, _dna_vaf = _ranking.ranked, _ranking.dna_vaf_by_variant

    options = PeptideConstructConfig(mode='slp')
    constructs = assemble_peptide_constructs(ranked, options=options)
    assert len(constructs) == 3, (
        "SLP mode: one construct per ranked vaccine peptide; "
        "got %d" % len(constructs))


# ---- pVACseq path coverage ----------------------------------------------

def test_pvacseq_to_ranked_vaccine_peptides_round_trip():
    """Smoke test for the pVACseq path: fixture has 3 unique variants
    (TP53/BRAF/KRAS) with dashed-form IDs (chr-start-end-ref-alt).
    The path was previously degenerate — _PVACSEQ_IC50_COLS was
    untested and the ID parser only handled dotted form, producing
    an empty ranked list end-to-end.
    """
    from vaxrank.epitope_io import read_pvacseq_report
    from vaxrank.external_input import pvacseq_ranking_result

    path = os.path.join(DATA_DIR, "pvacseq_example.tsv")
    _loaded = read_pvacseq_report(path)
    predictions = list(_loaded.epitopes)
    _ranking = pvacseq_ranking_result(_loaded, predictions)
    ranked, _dna_vaf = _ranking.ranked, _ranking.dna_vaf_by_variant
    assert len(ranked) == 3, (
        "pVACseq fixture has 3 unique variants; got %d" % len(ranked))
    genes = sorted(
        peps[0].mutant_protein_fragment.gene_name for _, peps in ranked)
    assert genes == ['BRAF', 'KRAS', 'TP53']
    # IDs are dashed-form (chr1-100000-100001-A-T); the parser strips
    # the ``chr`` prefix so pyensembl-style downstream lookups work.
    contigs = sorted(v.contig for v, _ in ranked)
    assert contigs == ['1', '2', '3']


def test_pvacseq_id_parser_handles_dashed_and_dotted():
    """The parser accepts both modern (chr1-100000-100001-A-T) and
    legacy (1.123.A.T) ID forms, plus the 4-part dashed variant
    (chr1-100000-A-T) without an explicit end position. Returns
    The returned Variant carries the real ref/alt nucleotides from the ID.

    The parser strips ``chr`` prefixes so pyensembl-style downstream
    lookups (``variant.effect_on_transcript``) work — same rationale
    as ``parse_lens_variant_coordinates`` on the LENS path."""
    from vaxrank.external_input import parse_pvacseq_variant

    for vid, expected in [
        ("chr1-100000-100001-A-T", ("1", 100000, "A", "T")),
        ("chr1-100000-A-T", ("1", 100000, "A", "T")),
        ("1.123.A.T", ("1", 123, "A", "T")),
    ]:
        v = parse_pvacseq_variant(vid)
        assert v is not None
        assert (v.contig, v.start, v.ref, v.alt) == expected
    # Garbage / wrong shape returns None instead of raising
    for bad in ("garbage", "chr1-notapos-A-T", None, ""):
        assert parse_pvacseq_variant(bad) is None


def test_pvacseq_drives_mrna_construct_assembly_end_to_end():
    """End-to-end: pVACseq report → ranked → mRNA constructs (mirror
    of the LENS end-to-end test). Pre-fix this produced zero
    constructs because the ranked list was empty."""
    from vaxrank.epitope_io import read_pvacseq_report
    from vaxrank.external_input import pvacseq_ranking_result
    from vaxrank.mrna import RNAConstructConfig, assemble_mrna_constructs

    path = os.path.join(DATA_DIR, "pvacseq_example.tsv")
    _loaded = read_pvacseq_report(path)
    predictions = list(_loaded.epitopes)
    _ranking = pvacseq_ranking_result(_loaded, predictions)
    ranked, _dna_vaf = _ranking.ranked, _ranking.dna_vaf_by_variant
    assert ranked

    options = RNAConstructConfig(
        signal_peptide=None, include_mitd=False,
        utr_3p='HBB', poly_a_length=10,
        antigens_per_construct=3, max_constructs=1,
        max_antigen_length_aa=12, min_antigen_length_aa=5,
        optimize_linkers=False)
    constructs = assemble_mrna_constructs(ranked, options=options)
    assert len(constructs) == 1
    c = constructs[0]
    assert len(c.antigen_names) == 3
    assert c.cds_nt.endswith("TAA")


def test_pvacseq_all_epitopes_to_ranked_vaccine_peptides(tmp_path):
    """all_epitopes flavor uses MT Epitope Seq + chr/ref/alt columns."""
    from tests.test_epitope_io import _write_pvacseq_all_epitopes_fixture
    from vaxrank.epitope_io import read_pvacseq_report
    from vaxrank.external_input import pvacseq_ranking_result

    path = _write_pvacseq_all_epitopes_fixture(tmp_path)
    _loaded = read_pvacseq_report(path)
    predictions = list(_loaded.epitopes)
    _ranking = pvacseq_ranking_result(_loaded, predictions)
    ranked, dna_vaf_by_variant = _ranking.ranked, _ranking.dna_vaf_by_variant

    assert len(ranked) == 1
    variant, vaccine_peptides = ranked[0]
    assert (variant.contig, variant.start, variant.ref, variant.alt) == (
        "1", 154590262, "T", "A")
    assert vaccine_peptides[0].mutant_protein_fragment.gene_name == "ADAR"
    # pVACseq has no longer source window that can honestly contain both
    # minimal epitopes. The default DSL chooses AERMGFTVV; evidence from the
    # other peptide cannot inflate this construct's target score.
    assert [e.sequence for e in vaccine_peptides[0].epitopes] == [
        "AERMGFTVV"]
    assert dna_vaf_by_variant[variant] == pytest.approx(0.302)


def test_pvacseq_duplicate_peptides_stay_variant_scoped(tmp_path):
    """Identical peptides from different pVACseq variants must not mix."""
    from tests.test_epitope_io import _write_pvacseq_duplicate_peptide_fixture
    from vaxrank.epitope_io import read_pvacseq_report
    from vaxrank.external_input import pvacseq_ranking_result

    path = _write_pvacseq_duplicate_peptide_fixture(tmp_path)
    _loaded = read_pvacseq_report(path)
    predictions = list(_loaded.epitopes)
    _ranking = pvacseq_ranking_result(_loaded, predictions)
    ranked, _dna_vaf_by_variant = _ranking.ranked, _ranking.dna_vaf_by_variant

    assert len(ranked) == 2
    expected_source_by_variant = {
        ("1", 154590262, "T", "A"): "chr1-154590262-T-A",
        ("2", 200000, "G", "C"): "chr2-200000-G-C",
    }
    for variant, vaccine_peptides in ranked:
        key = (variant.contig, variant.start, variant.ref, variant.alt)
        expected_source = expected_source_by_variant[key]
        epitopes = vaccine_peptides[0].epitopes
        assert len(epitopes) == 1
        assert epitopes[0].sequence == "AERMGFTVV"
        # Variant provenance lives in ``source_name``. ``source_sequence`` is
        # an amino-acid window and must stay one: pVACseq ships no source
        # context, so the peptide is its own window. Storing the variant
        # string there made construct assembly slice a variant *name*.
        assert epitopes[0].source_name == expected_source
        assert epitopes[0].source_sequence == "AERMGFTVV"
        assert (
            vaccine_peptides[0].mutant_protein_fragment.amino_acids
            == epitopes[0].source_sequence)


# ---- LENS scoring-column edge cases --------------------------------------

def test_lens_without_affinity_uses_dsl_score_without_fallback_warning(
        tmp_path, caplog):
    """Non-affinity evidence no longer triggers a raw-affinity fallback."""
    import logging
    from vaxrank.epitope_io import read_lens_report
    from vaxrank.external_input import lens_ranking_result

    no_aff = tmp_path / "no_affinity.tsv"
    no_aff.write_text(
        "allele\tpeptide\tnetmhcstabpan_1.0.halflife_hours\t"
        "netmhcstabpan_1.0.perc_rank_stab\tantigen_source\tmut_aa_pos\t"
        "variant_coords\tsnv_ref_allele\tsnv_alt_allele\t"
        "gene_name\ttpm\tpep_context\n"
        "HLA-A02:01\tSVVGSSSSS\t12.5\t0.5\tSNV\t245\t"
        "chr17:7675088\tC\t[T]\tTP53\t42.5\tAASVVGSSSSSGTR\n")

    _loaded = read_lens_report(str(no_aff))
    predictions = list(_loaded.epitopes)
    with caplog.at_level(logging.WARNING):
        _ranking = lens_ranking_result(_loaded, predictions)
        ranked, _dna_vaf = _ranking.ranked, _ranking.dna_vaf_by_variant
    # The configured/default DSL is now the sole representative selector.
    assert len(ranked) == 1
    assert not any(
        "no pMHC_affinity scoring columns" in record.message
        for record in caplog.records)


# ---- emit_outputs observability -------------------------------------------

def _vaccine_args(output_dir, vaccine_type='mrna'):
    """Build a SimpleNamespace covering both peptide and mrna design
    knobs so tests can pass any ``vaccine_type`` (single or list) and
    have all required attrs present. Centralizes defaults — when a
    knob's CLI default changes, only this helper needs to follow."""
    from types import SimpleNamespace
    types = vaccine_type if isinstance(vaccine_type, list) else [vaccine_type]
    return SimpleNamespace(
        output_csv='', output_xlsx_report='',
        output_neoepitope_report='',
        vaccine_type=types,
        output_dir=output_dir,
        # Shared antigen-design axes (resolved per modality at writer time)
        antigen_content=None, epitopes_per_antigen=None,
        # mRNA knobs
        mrna_csv_full_rows=True,
        mrna_signal_peptide=None, mrna_linker='(G4S)2',
        mrna_include_mitd=False, mrna_mitd='HLA_A',
        mrna_5p_utr='HBB', mrna_3p_utr='HBB',
        mrna_codon_species='h_sapiens',
        mrna_codon_method='use_best_codon',
        mrna_min_antigen_length_aa=5, mrna_max_antigen_length_aa=14,
        mrna_antigens_per_construct=3, mrna_max_constructs=1,
        mrna_candidates_per_slot=1, mrna_max_length_nt=4000,
        mrna_optimize_linkers=False, mrna_junction_candidates='',
        mrna_junction_rank_strong=0.5, mrna_junction_rank_mild=2.0,
        mrna_poly_a_length=20, mrna_poly_a_segmented=False,
        mrna_poly_a_first_segment=30,
        mrna_poly_a_segment_linker='GCATATGACT',
        mrna_antigen_content=None, mrna_epitopes_per_antigen=None,
        # Peptide knobs
        peptide_mode='slp',
        peptide_antigen_content=None, peptide_epitopes_per_antigen=None,
        peptide_linker='G4S3',
        peptide_min_antigen_length_aa=5,
        peptide_max_antigen_length_aa=14,
        peptide_antigens_per_construct=1,
        peptide_max_constructs=5,
        peptide_candidates_per_slot=1,
        peptide_n_terminal_acetyl=False,
        peptide_c_terminal_amide=False,
        peptide_scale_mg=5.0,
        peptide_purity_percent=95.0,
        peptide_counterion='TFA',
    )


# Backward-compat alias for tests that already say ``_mrna_args``.
_mrna_args = _vaccine_args


def test_emit_outputs_logs_active_and_fired_dispatch(caplog, tmp_path):
    """The dispatch summary log line is part of the contract — pins
    the active vaccine type(s) and which writers fired."""
    import logging
    from vaxrank.cli.entry_point import emit_outputs

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    _loaded = read_lens_report(path)
    predictions = list(_loaded.epitopes)
    _ranking = lens_ranking_result(_loaded, predictions)
    ranked, _dna_vaf = _ranking.ranked, _ranking.dna_vaf_by_variant

    out_dir = str(tmp_path / "out")
    args = _mrna_args(out_dir)
    args.input_lens = path  # source-specific label derives from this
    with caplog.at_level(logging.INFO):
        emit_outputs(args, ranked, source='external')
    dispatch = [r.message for r in caplog.records
                if 'Vaccine-type dispatch' in r.message]
    assert len(dispatch) == 1
    line = dispatch[0]
    assert "[LENS import]" in line
    assert "['mrna']" in line  # active types and fired list both ['mrna']


def test_resolve_vaccine_types_default_and_unknown():
    """``--vaccine-type`` is multi-valued (default
    ``['peptide', 'mrna']`` — both modalities). Unknown values
    raise so future types stub-declared in user scripts fail loudly
    instead of silently no-opping."""
    import pytest as _pytest
    from types import SimpleNamespace
    from vaxrank.cli.entry_point import resolve_vaccine_types

    # Defaults: both modalities active when --vaccine-type isn't set
    assert resolve_vaccine_types(SimpleNamespace()) == ['peptide', 'mrna']
    assert resolve_vaccine_types(
        SimpleNamespace(vaccine_type=None)) == ['peptide', 'mrna']
    # Single-string form (caller bypassing argparse)
    assert resolve_vaccine_types(
        SimpleNamespace(vaccine_type='mrna')) == ['mrna']
    # Single-element list narrows to that one type
    assert resolve_vaccine_types(
        SimpleNamespace(vaccine_type=['peptide'])) == ['peptide']
    # Multi-mode + dedup
    assert resolve_vaccine_types(
        SimpleNamespace(vaccine_type=['peptide', 'mrna'])
    ) == ['peptide', 'mrna']
    assert resolve_vaccine_types(
        SimpleNamespace(vaccine_type=['mrna', 'peptide', 'mrna'])
    ) == ['mrna', 'peptide']  # first-occurrence order preserved
    # Unknown
    with _pytest.raises(ValueError, match="Unknown vaccine type"):
        resolve_vaccine_types(SimpleNamespace(vaccine_type=['dna']))


def test_vaccine_target_dir_layout(tmp_path):
    """Single-mode → files land directly in --output-dir; multi-mode
    → per-modality subdirs (so peptide's manifest.json doesn't
    overwrite mrna's)."""
    from vaxrank.cli.entry_point import vaccine_target_dir
    base = str(tmp_path / "out")
    # Single-mode flat
    assert vaccine_target_dir(base, 'peptide', ['peptide']) == base
    # Multi-mode subdirs
    assert vaccine_target_dir(base, 'peptide', ['peptide', 'mrna']) == \
        os.path.join(base, 'peptide')
    assert vaccine_target_dir(base, 'mrna', ['peptide', 'mrna']) == \
        os.path.join(base, 'mrna')
    # No --output-dir → no writer fires
    assert vaccine_target_dir('', 'peptide', ['peptide']) is None


def test_lens_with_vaccine_type_mrna_writes_three_fastas(tmp_path):
    """End-to-end: --input-lens + --vaccine-type=mrna +
    --output-dir=DIR writes three FASTAs + manifest + mrna-sequence-parts.csv
    directly into DIR (single-mode flat layout)."""
    from vaxrank.cli.entry_point import emit_outputs

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    _loaded = read_lens_report(path)
    predictions = list(_loaded.epitopes)
    _ranking = lens_ranking_result(_loaded, predictions)
    ranked, _dna_vaf = _ranking.ranked, _ranking.dna_vaf_by_variant

    out_dir = str(tmp_path / "out")
    args = _mrna_args(out_dir)
    emit_outputs(args, ranked, source='external')
    assert os.path.isfile(os.path.join(out_dir, "cds.fasta"))
    assert os.path.isfile(os.path.join(out_dir, "no_polyA.fasta"))
    assert os.path.isfile(os.path.join(out_dir, "full.fasta"))
    assert os.path.isfile(os.path.join(out_dir, "manifest.json"))
    with open(os.path.join(out_dir, "manifest.json")) as f:
        manifest = json.load(f)
    assert len(manifest) == 1
    assert manifest[0]['modality'] == 'mrna'
    # full ends with polyA
    full_body = ''.join(
        line for line in open(
            os.path.join(out_dir, "full.fasta")).read().split("\n")
        if not line.startswith(">"))
    assert full_body.endswith("A" * 20)


def test_lens_with_multi_vaccine_type_uses_subdirs(tmp_path):
    """Multi-mode (peptide + mrna) writes into per-modality subdirs
    so canonical filenames don't collide on manifest.json."""
    from vaxrank.cli.entry_point import emit_outputs

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    _loaded = read_lens_report(path)
    predictions = list(_loaded.epitopes)
    _ranking = lens_ranking_result(_loaded, predictions)
    ranked, _dna_vaf = _ranking.ranked, _ranking.dna_vaf_by_variant

    out_dir = str(tmp_path / "out")
    args = _vaccine_args(out_dir, vaccine_type=['peptide', 'mrna'])
    emit_outputs(args, ranked, source='external')
    # mRNA outputs in DIR/mrna/
    assert os.path.isfile(os.path.join(out_dir, "mrna", "cds.fasta"))
    assert os.path.isfile(os.path.join(out_dir, "mrna", "manifest.json"))
    # Peptide outputs in DIR/peptide/ — canonical filenames don't
    # collide with mrna's manifest.json
    assert os.path.isfile(os.path.join(out_dir, "peptide", "vaccine.fasta"))
    assert os.path.isfile(os.path.join(out_dir, "peptide", "manifest.json"))


def test_emit_outputs_skips_writer_when_no_output_dir(tmp_path):
    """No --output-dir → writer no-ops (ranking-only run)."""
    from vaxrank.cli.entry_point import emit_outputs

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    _loaded = read_lens_report(path)
    predictions = list(_loaded.epitopes)
    _ranking = lens_ranking_result(_loaded, predictions)
    ranked, _dna_vaf = _ranking.ranked, _ranking.dna_vaf_by_variant

    out_dir = str(tmp_path / "out")
    args = _mrna_args(out_dir)
    args.output_dir = ''  # no destination → no writer fires
    emit_outputs(args, ranked, source='external')
    assert not os.path.exists(out_dir), \
        "Without --output-dir, no writer should have fired"


def test_variant_is_frameshift():
    # varcode trims the shared prefix, so the helper compares the
    # *normalized* allele lengths: e.g. CAT/C -> AT/'' (delta 2).
    from varcode import Variant
    from vaxrank.external_input import variant_is_frameshift
    assert variant_is_frameshift(Variant('1', 100, 'C', 'T')) is False    # SNV, delta 0
    assert variant_is_frameshift(Variant('1', 100, 'CA', 'C')) is True     # 1bp del, delta 1
    assert variant_is_frameshift(Variant('1', 100, 'C', 'CG')) is True     # 1bp ins, delta 1
    assert variant_is_frameshift(Variant('1', 100, 'CAT', 'C')) is True    # delta 2
    assert variant_is_frameshift(Variant('1', 100, 'CATG', 'C')) is False  # delta 3, inframe
    assert variant_is_frameshift(Variant('1', 100, 'CATGA', 'C')) is True  # delta 4


def test_maximal_mutant_span_combines_rows_for_frameshift():
    from vaxrank.external_input import maximal_mutant_span
    ctx = "WTWTWTNOVELAAANOVELBBBNOVELCCC"  # 30 chars
    # non-frameshift: representative span unchanged
    assert maximal_mutant_span(6, 11, ["NOVELAAA"], ctx, False) == (6, 11)
    # frameshift: start = earliest located neoepitope, end = len(ctx)
    peptides = ["NOVELBBB", "NOVELAAA", "NOVELCCC"]
    start, end = maximal_mutant_span(14, 22, peptides, ctx, True)
    assert start == ctx.find("NOVELAAA")   # earliest
    assert end == len(ctx)                 # extends to context end


def test_check_varcode_annotation_with_stubs():
    """Pure comparison logic (no varcode/genome): gene + effect-class
    agreement on the provider's transcript, None when uncomputable."""
    from types import SimpleNamespace
    from vaxrank.external_input import check_varcode_annotation

    class FrameShift:
        gene_name = 'FYN'

    class Substitution:
        gene_name = 'TPMT'

    # No resolved transcript -> skip (None, None).
    assert check_varcode_annotation(
        SimpleNamespace(), None, 'FYN', True) == (None, None)
    # Frameshift + gene agree.
    v_fs = SimpleNamespace(effect_on_transcript=lambda t: FrameShift())
    assert check_varcode_annotation(v_fs, object(), 'FYN', True) == (True, True)
    # Gene mismatch.
    assert check_varcode_annotation(
        v_fs, object(), 'OTHER', True) == (False, True)
    # Effect-class mismatch (provider says frameshift, varcode substitution).
    v_sub = SimpleNamespace(effect_on_transcript=lambda t: Substitution())
    assert check_varcode_annotation(
        v_sub, object(), 'TPMT', True) == (True, False)

    # Exception during effect computation -> skipped, not raised.
    def _boom(_t):
        raise RuntimeError("effect blew up")
    assert check_varcode_annotation(
        SimpleNamespace(effect_on_transcript=_boom),
        object(), 'X', False) == (None, None)


def test_log_varcode_agreement(caplog):
    import logging
    from vaxrank.external_input import log_varcode_agreement

    # All agree -> single INFO; the None entry (uncomputable) is ignored.
    with caplog.at_level(logging.INFO, logger='vaxrank.external_input'):
        log_varcode_agreement(
            [(True, True, 'v1'), (True, True, 'v2'), (None, None, 'v3')],
            'LENS')
    assert any('agrees with LENS on all 2 checked' in r.getMessage()
               for r in caplog.records)

    # Mismatches -> WARNING; headline counts DISTINCT bad variants (v1
    # has both, v2 effect-only -> 2 distinct, not max(1,2)=2 by luck;
    # use a case where it differs: v1 gene-only, v2 effect-only -> 2).
    caplog.clear()
    with caplog.at_level(logging.WARNING, logger='vaxrank.external_input'):
        log_varcode_agreement(
            [(False, True, 'v1'), (True, False, 'v2'), (True, True, 'v3')],
            'pVACseq')
    msg = [r.getMessage() for r in caplog.records
           if 'disagrees' in r.getMessage()][0]
    assert 'disagrees with pVACseq on 2 / 3' in msg
    assert '1 gene' in msg and '1 effect-class' in msg

    # Nothing checked -> silent.
    caplog.clear()
    with caplog.at_level(logging.INFO, logger='vaxrank.external_input'):
        log_varcode_agreement([(None, None, 'x')], 'LENS')
    assert not any('varcode' in r.getMessage() for r in caplog.records)


# ---- identity / evidence separation regressions --------------------------

def _write_pvacseq_all_epitopes_with_index(tmp_path):
    """all_epitopes flavor carrying BOTH an Index and genomic coordinates."""
    path = tmp_path / "pvacseq_all_epitopes_indexed.tsv"
    columns = [
        "Index", "Chromosome", "Start", "Stop", "Reference", "Variant",
        "Transcript", "Variant Type", "Mutation", "Gene Name", "HLA Allele",
        "Mutation Position", "MT Epitope Seq", "WT Epitope Seq",
        "Median MT IC50 Score", "Median WT IC50 Score",
        "Median MT Percentile", "Median WT Percentile",
        "Tumor DNA Depth", "Tumor DNA VAF", "Tumor RNA Depth",
        "Tumor RNA VAF", "Gene Expression", "Transcript Expression",
    ]
    row = [
        "ADAR.ENST00000368474.9.missense.806E/V",
        "chr1", "154590262", "154590263", "T", "A",
        "ENST00000368474.9", "missense", "E806V", "ADAR",
        "HLA-B*45:01", "5", "AERMGFTVV", "AERMGFTEV",
        "76.11", "61.80", "0.10", "0.12",
        "1233", "0.302", "1233", "0.348", "131.832", "84.961",
    ]
    path.write_text("\t".join(columns) + "\n" + "\t".join(row) + "\n")
    return path


def test_pvacseq_index_and_coordinates_both_yield_a_ranked_construct(tmp_path):
    """An ``Index`` column must not cost the run its vaccine peptides.

    topiary normalizes ``variant`` from ``Index`` when that column exists, so
    an identity derived from genomic coordinates never joins. Aligning the
    identity is only half the fix: the ``Index`` string carries no
    coordinates, so the genomic variant has to come from the coordinate
    columns independently.
    """
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import epitopes_for_ranking
    from vaxrank.epitope_io import read_pvacseq_report
    from vaxrank.external_input import pvacseq_ranking_result

    path = _write_pvacseq_all_epitopes_with_index(tmp_path)
    config = EpitopeConfig()
    report = read_pvacseq_report(path, epitope_config=config)
    assert len(report.epitopes) == 1

    result = pvacseq_ranking_result(
        report, epitopes_for_ranking(list(report.epitopes), config))
    assert len(result.ranked) == 1
    variant, vaccine_peptides = result.ranked[0]
    # Genomic placement comes from the coordinate columns, not the Index.
    assert (variant.contig, variant.start, variant.ref, variant.alt) == (
        "1", 154590262, "T", "A")
    assert vaccine_peptides[0].amino_acids == "AERMGFTVV"


def test_lens_construct_excludes_epitopes_outside_its_window(tmp_path):
    """Only epitopes inside the synthesized peptide may score it.

    A LENS ``pep_context`` can be far longer than the vaccine peptide. The
    construct window is centered on the mutation, so a neoepitope elsewhere in
    the context is not in the peptide that would be synthesized and must not
    contribute to ``target_epitope_score`` or appear in its report table.
    """
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import epitopes_for_ranking
    from vaxrank.epitope_io import read_lens_report
    from vaxrank.external_input import lens_ranking_result

    context = (
        "MAAAAAAAAA" "SIINFEKLQQ" "CCCCCCCCCC" "DDDDDDDDDD" "EEEEEEEEEE"
        "FFFFFFFFFF" "GGGGGGGGGG" "YLLPAIVHIW" "HHHHHHHHHH")
    path = tmp_path / "lens_long_context.tsv"
    path.write_text(
        "peptide\tallele\tpep_context\tantigen_source\tvariant_coords\t"
        "snv_ref_allele\tsnv_alt_allele\tgene_name\ttranscript_id\t"
        "mhcflurry_2.1.1.aff\n"
        f"SIINFEKL\tHLA-A02:01\t{context}\tSNV\tchr1:100\tC\t[T]\tG\tENST1\t20\n"
        f"YLLPAIVHI\tHLA-A02:01\t{context}\tSNV\tchr1:100\tC\t[T]\tG\tENST1\t15\n")

    config = EpitopeConfig()
    report = read_lens_report(path, epitope_config=config)
    result = lens_ranking_result(
        report,
        epitopes_for_ranking(list(report.epitopes), config),
        options=ExternalConstructOptions(vaccine_peptide_length=25))

    assert len(result.ranked) == 1
    vaccine_peptide = result.ranked[0][1][0]
    assert len(vaccine_peptide.amino_acids) == 25
    for epitope in vaccine_peptide.epitopes:
        assert epitope.sequence in vaccine_peptide.amino_acids
        # Offsets must be rebased into the window they now describe.
        assert epitope.source_sequence == vaccine_peptide.amino_acids
        assert (
            vaccine_peptide.amino_acids[
                epitope.offset:epitope.offset + len(epitope.sequence)]
            == epitope.sequence)
    assert vaccine_peptide.target_epitope_score == pytest.approx(
        sum(e.epitope_score for e in vaccine_peptide.target_epitopes))


def test_lens_whitespace_in_pep_context_does_not_split_an_identity(tmp_path):
    """One normalization vocabulary, so a stray space cannot fork a group."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import read_lens_report

    path = tmp_path / "lens_whitespace.tsv"
    path.write_text(
        "peptide\tallele\tpep_context\tantigen_source\tvariant_coords\t"
        "snv_ref_allele\tsnv_alt_allele\tgene_name\tmhcflurry_2.1.1.aff\n"
        "SIINFEKL\tHLA-A02:01\tXXSIINFEKLXX\tSNV\tchr1:100\tC\t[T]\tG\t20\n"
        "SIINFEKL\tHLA-B07:02\tXXSIINFEKLXX \tSNV\tchr1:100\tC\t[T]\tG\t25\n")

    report = read_lens_report(path, epitope_config=EpitopeConfig())
    assert len(report.epitopes) == 1
    assert report.epitopes[0].patient_alleles == (
        "HLA-A*02:01", "HLA-B*07:02")
    # The identity index is what construct ranking joins against; it must not
    # raise on rows that only differ by whitespace.
    assert len(report.records_with_epitopes(list(report.epitopes))) == 2


def test_external_prediction_key_rejects_an_impossible_position():
    """The key is a join target, so its position fields must be real."""
    from vaxrank.external_prediction import ExternalPredictionKey

    with pytest.raises(ValueError, match="not at offset"):
        ExternalPredictionKey(
            source_format="lens", peptide="SIINFEKL",
            source_sequence="XXSIINFEKLXX", offset=0)


def test_vaccine_config_reaches_external_construct_assembly(tmp_path):
    """``vaccine_peptides:`` config must apply to LENS runs, not just VCF ones.

    It was resolved only inside ``run_vaxrank_from_parsed_args`` — the
    pipeline-only entry point — so external runs silently used a hard-coded
    25-aa window and kept every epitope regardless of configuration.
    """
    from types import SimpleNamespace

    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.external_input import load_external_ranked
    from vaxrank.vaccine_config import VaccineConfig

    context = (
        "MAAAAAAAAA" "SIINFEKLQQ" "CCCCCCCCCC" "DDDDDDDDDD" "EEEEEEEEEE"
        "FFFFFFFFFF" "GGGGGGGGGG" "YLLPAIVHIW" "HHHHHHHHHH")
    path = tmp_path / "config_lens.tsv"
    path.write_text(
        "peptide\tallele\tpep_context\tantigen_source\tvariant_coords\t"
        "snv_ref_allele\tsnv_alt_allele\tgene_name\ttranscript_id\t"
        "mhcflurry_2.1.1.aff\n"
        f"YLLPAIVHI\tHLA-A02:01\t{context}\tSNV\tchr1:100\tC\t[T]\tG\tENST1\t15\n"
        f"YLLPAIVHIW\tHLA-B07:02\t{context}\tSNV\tchr1:100\tC\t[T]\tG\tENST1\t30\n")

    def run(length, keep):
        args = SimpleNamespace(
            input_lens=str(path), input_pvacseq=None,
            output_patient_id="", genome=None)
        ranked, *_ = load_external_ranked(
            args,
            epitope_config=EpitopeConfig(),
            vaccine_config=VaccineConfig(
                preferred_peptide_length=length,
                min_peptide_length=length,
                max_peptide_length=length,
                num_target_epitopes_to_keep=keep))
        return ranked[0][1][0]

    assert len(run(25, 0).amino_acids) == 25
    assert len(run(41, 0).amino_acids) == 41

    # num_target_epitopes_to_keep was never forwarded at all.
    assert len(run(41, 0).target_epitopes) == 2
    assert [e.sequence for e in run(41, 1).target_epitopes] == ["YLLPAIVHI"]


def test_external_input_parser_accepts_vaccine_peptide_flags():
    """The flags were rejected as unknown, not merely ignored."""
    from vaxrank.cli.arg_parser import external_input_arg_parser

    args = external_input_arg_parser().parse_args([
        "--input-lens", "x.tsv",
        "--vaccine-peptide-length", "31",
        "--num-epitopes-per-vaccine-peptide", "2",
    ])
    assert args.vaccine_peptide_length == 31
    assert args.num_epitopes_per_vaccine_peptide == 2


# ---- code-review regressions ---------------------------------------------

def test_construct_is_labelled_with_the_row_s_own_gene_not_the_alphabetical_one(
        tmp_path):
    """A vaccine peptide must carry the gene the report named for its row.

    ``gene_names`` is sorted so two rows listing the same genes in different
    orders resolve to one identity — which makes ``gene_names[0]``
    alphabetical, not primary. Labelling constructs from it put the wrong
    gene on every template report and produced spurious varcode
    gene-disagreement warnings.
    """
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import epitopes_for_ranking
    from vaxrank.epitope_io import read_lens_report
    from vaxrank.external_input import lens_ranking_result

    path = tmp_path / "gene_precedence.tsv"
    path.write_text(
        "peptide\tallele\tpep_context\tantigen_source\tvariant_coords\t"
        "snv_ref_allele\tsnv_alt_allele\tgene_name\t"
        "all_gene_names_encoding_peptide\ttranscript_id\t"
        "all_transcript_ids_encoding_peptide\tmhcflurry_2.1.1.aff\n"
        "SIINFEKL\tHLA-A02:01\tAASIINFEKLAA\tSNV\tchr17:7577120\tC\t[T]\t"
        "TP53\tATXN7,TP53\tENST00000269305\t"
        "ENST00000000233,ENST00000269305\t20\n")

    config = EpitopeConfig()
    report = read_lens_report(path, epitope_config=config)
    key = report.records[0].key
    # The identity keeps both, sorted, so row order cannot fork it ...
    assert key.gene_names == ("ATXN7", "TP53")
    # ... while the row's own naming is preserved separately.
    assert key.primary_gene_name == "TP53"
    assert key.ordered_transcript_ids[0] == "ENST00000269305"
    # Annotation must stay out of the identity, or two rows naming the same
    # genes in different orders would become different candidates.
    assert "primary_gene_name" not in key.identifier

    result = lens_ranking_result(
        report, epitopes_for_ranking(list(report.epitopes), config))
    fragment = result.ranked[0][1][0].mutant_protein_fragment
    assert fragment.gene_name == "TP53"


def test_lens_read_counts_come_from_one_row(tmp_path):
    """Counts must describe one observation, not per-field maxima.

    Maximizing each field independently can assemble a tuple no row ever
    reported — total from a deep row, alt from a shallow one — and
    n_alt_reads feeds the combined-score DSL, so an incoherent count
    reorders the construct ranking.

    The alt count is depth x VAF, which is why both rows carry a VAF. It
    used to be the CDS-overlap column, which counts reads overlapping the
    peptide's coding sequence rather than reads supporting the variant.
    """
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import epitopes_for_ranking
    from vaxrank.epitope_io import read_lens_report
    from vaxrank.external_input import lens_ranking_result

    path = tmp_path / "counts.tsv"
    path.write_text(
        "peptide\tallele\tpep_context\tantigen_source\tvariant_coords\t"
        "snv_ref_allele\tsnv_alt_allele\tgene_name\t"
        "rna_reads_covering_genomic_origin\t"
        "rna_reads_covering_genomic_origin_with_peptide_cds\tvaf\t"
        "mhcflurry_2.1.1.aff\n"
        # deep coverage, few supporting reads: 100 x 0.05 = 5 alt
        "SIINFEKL\tHLA-A02:01\tAASIINFEKLAA\tSNV\tchr1:100\tC\t[T]\tG\t"
        "100\t5\t0.05\t20\n"
        # shallow coverage, many supporting reads: 10 x 0.9 = 9 alt
        "SIINFEKL\tHLA-B07:02\tAASIINFEKLAA\tSNV\tchr1:100\tC\t[T]\tG\t"
        "10\t9\t0.9\t25\n")

    config = EpitopeConfig()
    report = read_lens_report(path, epitope_config=config)
    result = lens_ranking_result(
        report, epitopes_for_ranking(list(report.epitopes), config))
    fragment = result.ranked[0][1][0].mutant_protein_fragment

    # The best-supported row wins outright: (10, 9), never the cross-product
    # (100, 9) that no row reported.
    assert (fragment.n_overlapping_reads, fragment.n_alt_reads) == (10, 9)
    assert fragment.n_ref_reads == 1


def _lens_counts_file(tmp_path, name, *, depth, cds_overlap, vaf):
    """One LENS row with the two read columns and an optional VAF."""
    path = tmp_path / name
    vaf_header = "\tvaf" if vaf is not None else ""
    vaf_cell = "\t%s" % vaf if vaf is not None else ""
    path.write_text(
        "peptide\tallele\tpep_context\tantigen_source\tvariant_coords\t"
        "snv_ref_allele\tsnv_alt_allele\tgene_name\t"
        "rna_reads_covering_genomic_origin\t"
        "rna_reads_covering_genomic_origin_with_peptide_cds"
        + vaf_header + "\tmhcflurry_2.1.1.aff\n"
        "SIINFEKL\tHLA-A02:01\tAASIINFEKLAA\tSNV\tchr1:100\tC\t[T]\tG\t"
        "%s\t%s" % (depth, cds_overlap) + vaf_cell + "\t20\n")
    return path


def _fragment_for(path):
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import epitopes_for_ranking
    from vaxrank.epitope_io import read_lens_report
    from vaxrank.external_input import lens_ranking_result

    config = EpitopeConfig()
    report = read_lens_report(path, epitope_config=config)
    result = lens_ranking_result(
        report, epitopes_for_ranking(list(report.epitopes), config))
    return result.ranked[0][1][0].mutant_protein_fragment


def test_alt_reads_are_depth_times_vaf_not_the_cds_overlap_count(tmp_path):
    """The CDS-overlap column is not variant support.

    ``rna_reads_covering_genomic_origin_with_peptide_cds`` counts reads
    overlapping the peptide's coding sequence, not reads carrying the
    variant allele. vaxrank assigned it straight to ``n_alt_reads``, which
    feeds the combined-score DSL — ``sqrt(n_alt_reads) *
    target_epitope_score`` is the documented canonical form.

    The two disagree in both directions on real files: 2337x0.022 gives 51
    alt reads where the CDS count says 2, and 291x0.258 gives 75 where it
    says 166.
    """
    path = _lens_counts_file(
        tmp_path, "split.tsv", depth=2337, cds_overlap=2, vaf=0.022)

    fragment = _fragment_for(path)

    assert fragment.n_overlapping_reads == 2337
    assert fragment.n_alt_reads == 51        # 2337 x 0.022, not 2
    assert fragment.n_ref_reads == 2286      # and the remainder, not 2335
    # The CDS count is kept, under the name that describes it.
    assert fragment.n_alt_reads_supporting_protein_sequence == 2


def test_no_vaf_means_no_alt_estimate_rather_than_the_cds_count(tmp_path):
    """A source that cannot answer must not be answered for.

    An estimate needs both depth and VAF. Substituting the CDS-overlap
    count because it is the only other number present is how a missing
    value becomes a number nobody measured.
    """
    path = _lens_counts_file(
        tmp_path, "no-vaf.tsv", depth=100, cds_overlap=40, vaf=None)

    fragment = _fragment_for(path)

    assert fragment.n_overlapping_reads == 100
    assert fragment.n_alt_reads == 0
    assert fragment.n_alt_reads_supporting_protein_sequence == 40


def test_read_counts_arrive_with_the_derivation_that_produced_them(tmp_path):
    """An estimate and a measurement are the same integer, not the same claim.

    After #374 a LENS ``n_alt_reads`` is depth x VAF, rounded. On another
    source the same field may be counted from an alignment. The report
    prints the same number either way and ``sqrt(n_alt_reads)`` weights
    them identically, so the label is what lets a reader tell them apart.
    """
    path = _lens_counts_file(
        tmp_path, "labelled.tsv", depth=2337, cds_overlap=2, vaf=0.022)

    fragment = _fragment_for(path)

    assert fragment.n_alt_reads == 51
    # rna_depth_x_source_vaf, not rna_depth_x_vaf: LENS names its read
    # columns rna_* explicitly and leaves `vaf` bare, so the fraction's
    # assay is unstated. The method says the split used the source's own
    # VAF rather than asserting it was an RNA one.
    assert fragment.rna_evidence_method == "rna_depth_x_source_vaf"
    assert fragment.sequence_source == "lens_pep_context"


def test_no_estimate_means_no_method_rather_than_a_confident_zero(tmp_path):
    """A row that could not be split must not claim a derivation.

    Labelling a zero ``rna_depth_x_vaf`` would assert that the estimate was
    made and came out zero, which is a different statement from "the source
    did not carry a VAF".
    """
    path = _lens_counts_file(
        tmp_path, "unlabelled.tsv", depth=100, cds_overlap=40, vaf=None)

    fragment = _fragment_for(path)

    assert fragment.n_alt_reads == 0
    assert fragment.rna_evidence_method == ""
    # The sequence still has a known source; only the split is missing.
    assert fragment.sequence_source == "lens_pep_context"


def test_the_method_describes_the_row_whose_counts_were_used(tmp_path):
    """Provenance travels with the counts, from one observation.

    The counts already come from a single best-supported row rather than
    per-field maxima. Reading the label separately would reintroduce that
    problem one field over: a method describing a row whose numbers were
    not the ones used.
    """
    path = tmp_path / "two-rows.tsv"
    path.write_text(
        "peptide\tallele\tpep_context\tantigen_source\tvariant_coords\t"
        "snv_ref_allele\tsnv_alt_allele\tgene_name\t"
        "rna_reads_covering_genomic_origin\t"
        "rna_reads_covering_genomic_origin_with_peptide_cds\tvaf\t"
        "mhcflurry_2.1.1.aff\n"
        # deep row, no VAF: no estimate, so no method
        "SIINFEKL\tHLA-A02:01\tAASIINFEKLAA\tSNV\tchr1:100\tC\t[T]\tG\t"
        "1000\t5\t\t20\n"
        # shallow row with a VAF: 10 x 0.9 = 9 alt, and it wins
        "SIINFEKL\tHLA-B07:02\tAASIINFEKLAA\tSNV\tchr1:100\tC\t[T]\tG\t"
        "10\t9\t0.9\t25\n")

    fragment = _fragment_for(path)

    # The winning row is the shallow one, so its depth and its label.
    assert (fragment.n_overlapping_reads, fragment.n_alt_reads) == (10, 9)
    assert fragment.rna_evidence_method == "rna_depth_x_source_vaf"


def test_pvacseq_estimates_are_labelled_too(tmp_path):
    """pVACseq reports a depth and a fraction, never an alt-read count.

    vaxrank multiplies them, which is the same estimate the LENS path
    makes — so it carries the same label. It was computing it inline and
    unlabelled, which is the bug #375 describes, one path over.
    """
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import epitopes_for_ranking
    from vaxrank.epitope_io import read_pvacseq_report
    from vaxrank.external_input import pvacseq_ranking_result

    config = EpitopeConfig()
    report = read_pvacseq_report(
        os.path.join(DATA_DIR, "pvacseq_example.tsv"), epitope_config=config)
    result = pvacseq_ranking_result(
        report, epitopes_for_ranking(list(report.epitopes), config))

    fragment = result.ranked[0][1][0].mutant_protein_fragment

    assert fragment.n_alt_reads > 0
    assert fragment.rna_evidence_method == "rna_depth_x_vaf"


def test_a_pvacseq_row_with_no_rna_vaf_claims_no_derivation():
    """Absent is not zero.

    ``pvacseq_rna_vaf`` used to default to 0.0, which made "the variant was
    looked for and not seen" indistinguishable from "no RNA data" — and
    would have labelled the second one ``rna_depth_x_vaf``, asserting an
    estimate that was never made.
    """
    from vaxrank.external_input import pvacseq_rna_vaf

    assert pvacseq_rna_vaf({"rna_vaf": 0.0}) == 0.0
    assert pvacseq_rna_vaf({}) is None
    assert pvacseq_rna_vaf({"rna_vaf": ""}) is None
