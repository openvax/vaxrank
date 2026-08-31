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

"""
Tests for epitope_io: serialization, deserialization, and external imports.
"""

import dataclasses
import os

import pytest
from mhctools.pred import Prediction

from vaxrank.candidate_epitope import COMPARATOR_WT, CandidateEpitope, Peptide
from vaxrank.epitope_io import (
    epitope_prediction_rows,
    predictions_to_dataframe,
    save_predictions,
    lens_epitope_position,
    load_predictions,
    read_pvacseq_report,
    read_lens_report,
    write_neoepitope_report,
)
from vaxrank.external_prediction import ExternalPredictionKey

DATA_DIR = os.path.join(os.path.dirname(__file__), "data", "epitope_fixtures")


def test_external_prediction_key_is_legible_stable_and_round_trippable():
    key = ExternalPredictionKey(
        source_format="lens",
        variant_id="chr1:1000",
        antigen_source="SNV",
        gene_ids=("ENSG2", "ENSG1"),
        gene_names=("GENE2", "GENE1"),
        transcript_ids=("ENST2", "ENST1"),
        peptide="SIINFEKL",
        source_sequence="AASIINFEKLLL",
        offset=2,
    )

    identifier = key.identifier

    assert '"variant_id":"chr1:1000"' in identifier
    assert '"gene_ids":["ENSG1","ENSG2"]' in identifier
    assert ExternalPredictionKey.from_identifier(identifier) == key
    assert key.construct_identifier != key.identifier

    # Annotation is carried but deliberately excluded from the identity:
    # which gene a report named first is presentation, and folding it in
    # would make two rows describing one candidate into two candidates.
    annotated = dataclasses.replace(
        key, primary_gene_name="GENE2", primary_transcript_id="ENST2")
    assert annotated.identifier == identifier
    # So a reconstructed key carries the identity but not that annotation —
    # equality holds against the un-annotated key, not the annotated one.
    assert ExternalPredictionKey.from_identifier(identifier) != annotated


def test_native_round_trip_preserves_external_prediction_identity(tmp_path):
    key = ExternalPredictionKey(
        source_format="lens",
        variant_id="chr1:1000",
        gene_names=("GENE1",),
        transcript_ids=("ENST1",),
        peptide="SIINFEKL",
        source_sequence="AASIINFEKLLL",
        offset=2,
    )
    original = _make_prediction()
    original = CandidateEpitope.from_peptide(
        original,
        comparators=original.comparators,
        source_class=original.source_class,
        overlaps_mutation=original.overlaps_mutation,
        prediction_id=key.identifier,
    )

    path = tmp_path / "external-provenance.csv"
    save_predictions([original], path)
    reloaded = load_predictions(path)[0]

    assert reloaded.prediction_id == key.identifier


def test_lens_epitope_position_uses_source_context_and_offset():
    assert lens_epitope_position(
        "SIINFEKL", "AASIINFEKLLL") == (
            "SIINFEKL", "AASIINFEKLLL", 2)
    assert lens_epitope_position("SIINFEKL", "") == (
        "SIINFEKL", "", 0)
    assert lens_epitope_position("SIINFEKL", "UNRELATED") is None


def _make_prediction(allele="HLA-A*02:01", peptide_sequence="SIINFEKL",
                     wt_peptide_sequence="SIINFEKL", ic50=50.0,
                     wt_ic50=5000.0, percentile_rank=0.5,
                     prediction_method_name="test",
                     predictor_version="",
                     overlaps_mutation=True,
                     source_sequence="XXSIINFEKLXX", offset=2,
                     occurs_in_reference=False):
    """Build a minimal single-allele single-predictor ``CandidateEpitope`` for
    epitope_io roundtrip tests."""
    mutant_pred = Prediction(
        kind='pMHC_affinity', predictor_name=prediction_method_name,
        predictor_version=predictor_version, allele=allele,
        peptide=peptide_sequence, value=ic50, score=0.0,
        percentile_rank=percentile_rank)
    comparators = {}
    if wt_peptide_sequence and wt_ic50 is not None:
        wt_pred = Prediction(
            kind='pMHC_affinity', predictor_name=prediction_method_name,
            predictor_version=predictor_version, allele=allele,
            peptide=wt_peptide_sequence, value=wt_ic50, score=0.0,
            percentile_rank=None)
        comparators[COMPARATOR_WT] = Peptide(
            sequence=wt_peptide_sequence,
            source_sequence=source_sequence, offset=offset,
            predictions=(wt_pred,))
    return CandidateEpitope.from_peptide(
        Peptide(
            sequence=peptide_sequence,
            source_sequence=source_sequence, offset=offset,
            predictions=(mutant_pred,)),
        comparators=comparators,
        overlaps_mutation=overlaps_mutation,
        occurs_in_reference=occurs_in_reference)


def _first_leaf(epitope, kind='pMHC_affinity'):
    """Helper for tests asserting on per-allele record fields. Returns
    the first leaf ``Prediction`` of the requested kind in the
    epitope's mutant context."""
    return epitope.predictions_for(kind)[0]


def _first_wt_leaf(epitope, kind='pMHC_affinity'):
    if epitope.wt is None:
        return None
    leaves = epitope.wt.predictions_for(kind)
    return leaves[0] if leaves else None


# ── Roundtrip ser/de ─────────────────────────────────────────────────────────

def test_predictions_to_dataframe():
    preds = [_make_prediction(), _make_prediction(allele="HLA-B*07:02", ic50=200.0)]
    df = predictions_to_dataframe(preds)
    assert len(df) == 2
    assert "allele" in df.columns
    assert "peptide_sequence" in df.columns


def test_save_and_load_csv_roundtrip(tmp_path):
    # Two distinct Epitopes via differing allele AND source so the
    # adapter doesn't collapse them; the second has no WT.
    preds = [
        _make_prediction(),
        _make_prediction(allele="HLA-B*07:02", peptide_sequence="DIFFRENT",
                         source_sequence="XXDIFFRENTXX", ic50=200.0,
                         wt_ic50=None),
    ]
    path = tmp_path / "test.csv"
    save_predictions(preds, path)
    loaded = load_predictions(path)
    assert len(loaded) == 2
    e0 = next(e for e in loaded if e.sequence == "SIINFEKL")
    e1 = next(e for e in loaded if e.sequence == "DIFFRENT")
    assert _first_leaf(e0).allele == "HLA-A*02:01"
    assert _first_leaf(e0).value == 50.0
    assert e1.wt is None


def test_save_and_load_tsv_roundtrip(tmp_path):
    preds = [_make_prediction()]
    path = tmp_path / "test.tsv"
    save_predictions(preds, path)
    loaded = load_predictions(path)
    assert len(loaded) == 1
    assert loaded[0].sequence == "SIINFEKL"


def test_roundtrip_preserves_values(tmp_path):
    original = _make_prediction(
        percentile_rank=1.23, occurs_in_reference=True, offset=5)
    path = tmp_path / "test.csv"
    save_predictions([original], path)
    loaded = load_predictions(path)[0]
    assert _first_leaf(loaded).percentile_rank == pytest.approx(1.23)
    assert loaded.occurs_in_reference is True
    assert loaded.offset == 5


def test_native_roundtrip_preserves_evidence_only_predictions(tmp_path):
    """Native exports retain every evidence kind and its scoring provenance."""
    predictions = (
        Prediction(
            kind="pMHC_presentation",
            predictor_name="mhcflurry",
            predictor_version="2.1.1",
            allele="HLA-A*02:01",
            peptide="SIINFEKL",
            value=None,
            score=0.85,
            percentile_rank=0.28,
            tcr="CASSIRSSYEQYF",
            n_flank="AA",
            c_flank="LL",
            source_sequence_name="GENE1",
            offset=2,
        ),
        Prediction(
            kind="antigen_processing",
            predictor_name="mhcflurry",
            predictor_version="2.1.1",
            allele="",
            peptide="SIINFEKL",
            value=None,
            score=0.73,
            percentile_rank=None,
        ),
    )
    original = CandidateEpitope(
        sequence="SIINFEKL",
        source_sequence="XXSIINFEKLXX",
        offset=2,
        predictions=predictions,
        patient_alleles=("HLA-A*02:01", "HLA-B*07:02"),
        per_allele_scores={
            "HLA-A*02:01": 1.58,
            "HLA-B*07:02": 0.73,
        },
    )

    rows = epitope_prediction_rows(original)
    assert [row["prediction_kind"] for row in rows] == [
        "antigen_processing", "pMHC_presentation"]

    path = tmp_path / "evidence-only.csv"
    save_predictions([original], path)
    exported = predictions_to_dataframe([original])
    assert len(exported) == 2
    assert set(exported["prediction_kind"]) == {
        "antigen_processing", "pMHC_presentation"}

    loaded = load_predictions(path)
    assert len(loaded) == 1
    reloaded = loaded[0]
    assert reloaded.patient_alleles == original.patient_alleles
    assert reloaded.per_allele_scores == original.per_allele_scores
    assert reloaded.predictions_flat() == original.predictions_flat()


def test_native_reload_rejects_missing_explicit_prediction_score(tmp_path):
    """New native rows never turn missing scientific evidence into zero."""
    path = tmp_path / "missing-score.csv"
    dataframe = predictions_to_dataframe([_make_prediction()])
    dataframe.loc[0, "prediction_score"] = None
    dataframe.to_csv(path, index=False)

    with pytest.raises(ValueError, match="missing prediction_score"):
        load_predictions(path)


def test_native_roundtrip_preserves_complete_wt_prediction(tmp_path):
    """WT leaves retain the same complete provenance as mutant leaves."""
    mutant = Prediction(
        kind="pMHC_affinity",
        predictor_name="mhcflurry",
        predictor_version="2.1.1",
        allele="HLA-A*02:01",
        peptide="SIINFEKL",
        value=50.0,
        score=0.91,
        percentile_rank=0.5,
    )
    wt_prediction = Prediction(
        kind="pMHC_affinity",
        predictor_name="mhcflurry",
        predictor_version="2.1.1",
        allele="HLA-A*02:01",
        peptide="SIINFEKM",
        value=500.0,
        score=0.42,
        percentile_rank=4.8,
        tcr="CASSWTYEQYF",
        n_flank="NN",
        c_flank="CC",
        source_sequence_name="WT_TRANSCRIPT_1",
        offset=17,
    )
    original = CandidateEpitope(
        sequence="SIINFEKL",
        source_sequence="AASIINFEKLLL",
        offset=2,
        predictions=(mutant,),
        comparators={
            COMPARATOR_WT: Peptide(
                sequence="SIINFEKM",
                predictions=(wt_prediction,),
            ),
        },
        patient_alleles=("HLA-A*02:01",),
        per_allele_scores={"HLA-A*02:01": 0.91},
    )

    path = tmp_path / "complete-wt.csv"
    save_predictions([original], path)
    reloaded = load_predictions(path)[0]

    assert reloaded.wt is not None
    assert reloaded.wt.predictions_flat() == (wt_prediction,)


# ── pVACseq import ───────────────────────────────────────────────────────────

def test_load_pvacseq():
    path = os.path.join(DATA_DIR, "pvacseq_example.tsv")
    _loaded = read_pvacseq_report(path)
    report_df, epitopes = _loaded.report_df, list(_loaded.epitopes)
    # One CandidateEpitope per pVACseq row (each row has its own peptide).
    assert len(epitopes) == 3
    assert len(report_df) == 3

    e = next(
        ep for ep in epitopes
        if ep.sequence == "SVVGSSSSS")
    leaf = _first_leaf(e)
    assert leaf.allele == "HLA-A*02:01"
    assert e.sequence == "SVVGSSSSS"
    assert leaf.value == pytest.approx(76.11)
    wt_leaf = _first_wt_leaf(e)
    assert wt_leaf is not None
    assert wt_leaf.value == pytest.approx(4520.30)
    assert leaf.percentile_rank == pytest.approx(0.5)
    assert e.overlaps_mutation is True

    # Check report DataFrame has expected columns
    assert "Gene name" in report_df.columns
    assert "Genomic variant" in report_df.columns
    assert "Tier" in report_df.columns

    # Check the Ref-Match row landed with occurs_in_reference=True.
    assert any(ep.occurs_in_reference for ep in epitopes)


def test_load_pvacseq_populates_per_allele_scores():
    """Regression: the 3.1.0 single-mechanism refactor moved
    per-(peptide, allele) scoring to the DSL and stored the result on
    ``CandidateEpitope.per_allele_scores``. The pVACseq loader path
    must populate it too — otherwise every loaded epitope has
    ``epitope_score=0`` and every variant gets dropped during
    ranking as 'no epitopes for peptide'. Pin both ends: scores are
    non-empty and (for binding-affinity rows) non-zero."""
    path = os.path.join(DATA_DIR, "pvacseq_example.tsv")
    _loaded = read_pvacseq_report(path)
    epitopes = list(_loaded.epitopes)
    assert epitopes
    # Every epitope has scores keyed by its actual allele(s).
    for e in epitopes:
        assert e.per_allele_scores, (
            f"Epitope {e.sequence!r} loaded with empty per_allele_scores"
        )
        for p in e.predictions_flat():
            assert p.allele in e.per_allele_scores
    # The strong binder in the fixture (SVVGSSSSS, IC50≈76 nM) must
    # produce a non-zero score under the default logistic.
    e_strong = next(
        ep for ep in epitopes if ep.sequence == "SVVGSSSSS")
    assert e_strong.epitope_score > 0.0


def test_load_pvacseq_missing_columns(tmp_path):
    bad_file = tmp_path / "bad.tsv"
    bad_file.write_text("Gene\tAA Change\nTP53\tG245S\n")
    with pytest.raises(ValueError, match="Could not detect pVACseq format"):
        read_pvacseq_report(bad_file)


def _write_pvacseq_all_epitopes_fixture(tmp_path):
    path = tmp_path / "pvacseq_all_epitopes.tsv"
    columns = [
        "Chromosome", "Start", "Stop", "Reference", "Variant",
        "Transcript", "Variant Type", "Mutation", "Gene Name", "HLA Allele",
        "Mutation Position", "MT Epitope Seq", "WT Epitope Seq",
        "Median MT IC50 Score", "Median WT IC50 Score",
        "Median MT Percentile", "Median WT Percentile",
        "NetMHCpan MT IC50 Score", "NetMHCpan WT IC50 Score",
        "NetMHCpan MT Percentile", "NetMHCpan WT Percentile",
        "MHCflurry MT IC50 Score", "MHCflurry WT IC50 Score",
        "MHCflurry MT Percentile", "MHCflurry WT Percentile",
        "Tumor DNA Depth", "Tumor DNA VAF", "Tumor RNA Depth",
        "Tumor RNA VAF", "Gene Expression", "Transcript Expression",
    ]
    rows = [
        [
            "chr1", "154590262", "154590263", "T", "A",
            "ENST00000368474.9", "missense", "E806V", "ADAR",
            "HLA-B*45:01", "5", "AERMGFTVV", "AERMGFTEV",
            "76.11", "61.80", "0.10", "0.12",
            "20.16", "28.54", "0.02", "0.04",
            "61.51", "58.48", "0.09", "0.07",
            "1233", "0.302", "1233", "0.348", "131.832", "84.961",
        ],
        [
            "chr1", "154590262", "154590263", "T", "A",
            "ENST00000368474.9", "missense", "E806V", "ADAR",
            "HLA-B*45:01", "6", "AERMGFTVVT", "AERMGFTEVT",
            "81.43", "40.97", "0.37", "0.20",
            "NA", "NA", "NA", "NA",
            "81.43", "40.97", "0.37", "0.20",
            "1233", "0.302", "1233", "0.348", "131.832", "84.961",
        ],
        [
            "chr1", "154590262", "154590263", "T", "A",
            "ENST00000368474.9", "missense", "E806V", "ADAR",
            "HLA-A*29:02", "5", "AERMGFTVV", "AERMGFTEV",
            "850.30", "900.10", "2.50", "2.80",
            "840.0", "890.0", "2.4", "2.7",
            "860.6", "910.2", "2.6", "2.9",
            "1233", "0.302", "1233", "0.348", "131.832", "84.961",
        ],
    ]
    lines = ["\t".join(columns)]
    lines.extend("\t".join(row) for row in rows)
    path.write_text("\n".join(lines) + "\n")
    return path


def _write_pvacseq_duplicate_peptide_fixture(tmp_path):
    path = tmp_path / "pvacseq_duplicate_peptide.tsv"
    columns = [
        "Chromosome", "Start", "Stop", "Reference", "Variant",
        "Transcript", "Variant Type", "Mutation", "Gene Name", "HLA Allele",
        "Mutation Position", "MT Epitope Seq", "WT Epitope Seq",
        "Median MT IC50 Score", "Median WT IC50 Score",
        "Median MT Percentile", "Median WT Percentile",
        "Tumor DNA Depth", "Tumor DNA VAF", "Tumor RNA Depth",
        "Tumor RNA VAF", "Gene Expression", "Transcript Expression",
    ]
    rows = [
        [
            "chr1", "154590262", "154590263", "T", "A",
            "ENST00000368474.9", "missense", "E806V", "ADAR",
            "HLA-B*45:01", "5", "AERMGFTVV", "AERMGFTEV",
            "76.11", "61.80", "0.10", "0.12",
            "1233", "0.80", "1233", "0.348", "131.832", "84.961",
        ],
        [
            "chr2", "200000", "200001", "G", "C",
            "ENST00000288602.6", "missense", "V600E", "BRAF",
            "HLA-B*45:01", "5", "AERMGFTVV", "AERMGFTEV",
            "410.00", "900.10", "2.50", "2.80",
            "900", "0.20", "900", "0.100", "70.0", "42.0",
        ],
    ]
    lines = ["\t".join(columns)]
    lines.extend("\t".join(row) for row in rows)
    path.write_text("\n".join(lines) + "\n")
    return path


def test_load_pvacseq_all_epitopes_flavor(tmp_path):
    path = _write_pvacseq_all_epitopes_fixture(tmp_path)
    _loaded = read_pvacseq_report(path)
    report_df, epitopes = _loaded.report_df, list(_loaded.epitopes)

    assert len(report_df) == 3
    # Same source + peptide groups the two AERMGFTVV alleles together.
    assert len(epitopes) == 2
    multi_allele = next(e for e in epitopes if e.sequence == "AERMGFTVV")
    assert {p.allele for p in multi_allele.predictions_flat()} == {
        "HLA-A*29:02", "HLA-B*45:01"}
    assert multi_allele.wt.sequence == "AERMGFTEV"
    assert multi_allele.per_allele_scores

    topiary_df = report_df.attrs["topiary_df"]
    assert "pvacseq_mhcflurry_ic50_mt" in topiary_df.columns
    assert "pvacseq_mhcflurry_ic50_mt" in report_df.columns
    assert set(report_df["MHC class"]) == {"I"}


def test_pvacseq_all_epitopes_per_algorithm_column_scores(tmp_path):
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report
    import pandas as pd

    path = _write_pvacseq_all_epitopes_fixture(tmp_path)
    _loaded = read_pvacseq_report(path)
    report_df, epitopes = _loaded.report_df, list(_loaded.epitopes)
    csv_path = tmp_path / "scored.csv"
    cfg = EpitopeConfig(score_expr="pvacseq_mhcflurry_ic50_mt")

    write_neoepitope_report(
        report_df, epitopes, csv_report_path=str(csv_path),
        epitope_config=cfg)

    result = pd.read_csv(csv_path)
    scored = result[
        (result["Mutant peptide sequence"] == "AERMGFTVV")
        & (result["Allele"] == "HLA-B*45:01")
    ].iloc[0]
    assert scored["vaxrank_score"] == pytest.approx(61.51)


def test_pvacseq_source_keyed_filtered_scores_do_not_broadcast(tmp_path):
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report
    import pandas as pd

    path = _write_pvacseq_duplicate_peptide_fixture(tmp_path)
    _loaded = read_pvacseq_report(path)
    report_df, epitopes = _loaded.report_df, list(_loaded.epitopes)
    csv_path = tmp_path / "filtered.csv"
    cfg = EpitopeConfig(
        filter_expr="tumor_dna_vaf > 0.5",
        score_expr="affinity.value")

    write_neoepitope_report(
        report_df, epitopes, csv_report_path=str(csv_path),
        epitope_config=cfg)

    result = pd.read_csv(csv_path)
    by_source = result.set_index("Source sequence name")
    assert by_source.loc[
        "chr1-154590262-T-A", "vaxrank_score"] == pytest.approx(76.11)
    assert by_source.loc[
        "chr1-154590262-T-A", "vaxrank_filter_passed"]
    assert pd.isna(by_source.loc[
        "chr2-200000-G-C", "vaxrank_score"])
    assert not by_source.loc[
        "chr2-200000-G-C", "vaxrank_filter_passed"]
    assert by_source.loc[
        "chr2-200000-G-C", "vaxrank_exclusion_reason"] == "dsl_filter"


# ── LENS import ──────────────────────────────────────────────────────────────

def test_load_lens_emits_one_epitope_per_position_with_multi_predictor_leaves():
    """load_lens groups all per-(allele, predictor, kind) predictions for one
    ``(peptide, source, offset)`` position into a single CandidateEpitope.
    The fixture has three report rows and each CandidateEpitope carries
    separate MHCflurry affinity/presentation plus netMHCpan affinity leaves."""
    path = os.path.join(DATA_DIR, "lens_example.tsv")
    _loaded = read_lens_report(path)
    report_df, epitopes = _loaded.report_df, list(_loaded.epitopes)
    assert len(epitopes) == 3
    assert len(report_df) == 3

    leaves = [
        p for e in epitopes for p in e.predictions_flat()]
    methods = {p.predictor_name for p in leaves}
    assert methods == {"mhcflurry", "netmhcpan"}

    # Each CandidateEpitope carries both predictors' rows for one (peptide, allele).
    for e in epitopes:
        ms = {p.predictor_name for p in e.predictions_flat()}
        assert ms == {"mhcflurry", "netmhcpan"}
    assert report_df["Predictors used"].iloc[0] == "mhcflurry,netmhcpan"
    assert report_df["%ile rank"].iloc[0] == pytest.approx(0.45)
    assert report_df["mhcflurry affinity percentile rank"].iloc[0] == \
        pytest.approx(0.45)
    assert report_df["mhcflurry presentation score"].iloc[0] == \
        pytest.approx(0.85)
    assert report_df[
        "mhcflurry presentation percentile rank"].iloc[0] == \
        pytest.approx(0.28)

    # Find the SVVGSSSSS CandidateEpitope and inspect its mhcflurry leaf.
    e = next(
        ep for ep in epitopes
        if ep.sequence == "SVVGSSSSS")
    mhcf = e.predictions_for('pMHC_affinity', predictor='mhcflurry')[0]
    assert mhcf.allele == "HLA-A*02:01"  # 'HLA-A02:01' normalized
    assert mhcf.peptide == "SVVGSSSSS"
    assert mhcf.value == pytest.approx(95.4)
    assert mhcf.percentile_rank == pytest.approx(0.45)  # mhcflurry aff_perc
    assert mhcf.predictor_version == "2.1.1"
    presentation = e.predictions_for(
        'pMHC_presentation', predictor='mhcflurry')[0]
    assert presentation.allele == "HLA-A*02:01"
    assert presentation.value is None
    assert presentation.score == pytest.approx(0.85)
    assert presentation.percentile_rank == pytest.approx(0.28)
    # wt_ic50 is now derived from mhcflurry_agretopicity (= MT/WT ratio).
    # Fixture row: mhcflurry IC50 = 95.4, agretopicity = 0.020 → WT ≈ 4770.
    # Only mhcflurry leaf gets a WT pair (matches that tool's IC50 scale).
    wt_mhcf = next(
        p for p in e.wt.predictions_flat()
        if p.predictor_name == "mhcflurry")
    assert wt_mhcf.value == pytest.approx(95.4 / 0.020)
    assert e.overlaps_mutation is True
    assert e.source_sequence == "AASVVGSSSSSGTR"


def test_load_lens_populates_per_allele_scores():
    """Parallel of ``test_load_pvacseq_populates_per_allele_scores`` —
    the LENS loader must populate per_allele_scores so downstream
    ranking doesn't drop every variant as zero-scoring."""
    path = os.path.join(DATA_DIR, "lens_example.tsv")
    _loaded = read_lens_report(path)
    epitopes = list(_loaded.epitopes)
    assert epitopes
    for e in epitopes:
        assert e.per_allele_scores, (
            f"Epitope {e.sequence!r} loaded with empty per_allele_scores"
        )


def test_load_lens_report_display_uses_canonical_predictor():
    """The report DataFrame's single 'Predicted mutant pMHC affinity'
    column shows the canonical pMHC_affinity predictor (mhcflurry >
    netmhcpan > alphabetical). Per-predictor raw values land in
    individual columns alongside."""
    path = os.path.join(DATA_DIR, "lens_example.tsv")
    report_df = read_lens_report(path).report_df
    # Display value = mhcflurry 95.4 nM (canonical)
    assert "95.40 nM" in report_df["Predicted mutant pMHC affinity"].iloc[0]
    # Both raw-value columns present
    assert report_df["mhcflurry value (nM)"].iloc[0] == pytest.approx(95.4)
    assert report_df["netmhcpan value (nM)"].iloc[0] == pytest.approx(76.11)


def test_load_lens_mhcflurry_only_v19():
    """A v1.9-dev style file (no netmhcpan cols) emits only mhcflurry rows."""
    path = os.path.join(DATA_DIR, "lens_v1.9_mhcflurry_only.tsv")
    _loaded = read_lens_report(path)
    report_df, epitopes = _loaded.report_df, list(_loaded.epitopes)
    assert len(epitopes) == 3
    leaves = [p for e in epitopes for p in e.predictions_flat()]
    methods = {p.predictor_name for p in leaves}
    assert methods == {"mhcflurry"}

    e = next(
        ep for ep in epitopes
        if ep.sequence == "TAEFYQRY")
    leaf = _first_leaf(e)
    assert leaf.allele == "HLA-A*01:01"
    assert leaf.value == pytest.approx(254.15)
    assert leaf.predictor_version == "2.1.1"
    # Antigen source captured; variant_coords NA fallback to origin_descriptor
    assert report_df["Antigen source"].iloc[0] == "ERV"
    assert "Hsap38" in report_df["Genomic variant"].iloc[0]
    assert "mhcflurry value (nM)" in report_df.columns
    assert "netmhcpan value (nM)" not in report_df.columns


def test_load_lens_missing_required_columns(tmp_path):
    bad_file = tmp_path / "bad.tsv"
    bad_file.write_text("gene_name\tpos\nTP53\t5\n")
    with pytest.raises(ValueError, match="missing required columns"):
        read_lens_report(bad_file)


def test_load_lens_no_predictor_columns(tmp_path):
    """File has peptide/allele but no recognized affinity column."""
    bad_file = tmp_path / "bad.tsv"
    bad_file.write_text(
        "peptide\tallele\tgene_name\nSIINFEKL\tHLA-A*02:01\tOVA\n")
    with pytest.raises(ValueError, match="no recognized predictor"):
        read_lens_report(bad_file)


def test_load_lens_detects_netmhcstabpan():
    """netmhcstabpan columns are detected as a pMHC_stability predictor."""
    path = os.path.join(DATA_DIR, "lens_v1.4_with_stability.tsv")
    _loaded = read_lens_report(path)
    report_df, epitopes = _loaded.report_df, list(_loaded.epitopes)
    leaves = [p for e in epitopes for p in e.predictions_flat()]
    methods = {p.predictor_name for p in leaves}
    assert methods == {"mhcflurry", "netmhcpan", "netmhcstabpan"}
    # Stability value column is labeled with its unit
    assert "netmhcstabpan value (hours)" in report_df.columns
    # Older LENS output can carry pres_perc without pres_score. Preserve
    # the percentile leaf with a neutral score instead of dropping it.
    presentation = epitopes[0].predictions_for(
        "pMHC_presentation", predictor="mhcflurry")[0]
    assert presentation.score == 0.0
    assert presentation.percentile_rank is not None


def test_load_lens_preserves_mhcflurry_axes_for_coverage_and_dsl():
    """Affinity and presentation remain distinct clinical evidence.

    This pins the #265 failure mode: ``pres_perc`` must not overwrite
    ``aff_perc``, and the presentation axis must independently drive
    coverage and Topiary DSL scoring.
    """
    from vaxrank.coverage import compute_coverage
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import score_predictions

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    _loaded = read_lens_report(path)
    epitopes = list(_loaded.epitopes)
    epitope = next(e for e in epitopes if e.sequence == "SVVGSSSSS")

    coverage = compute_coverage([epitope], ["HLA-A*02:01"])["HLA-A*02:01"]
    # Coverage aggregates the best affinity evidence across models;
    # netMHCpan's 0.30 beats MHCflurry's independently preserved 0.45.
    assert coverage.best_affinity_pct == pytest.approx(0.30)
    assert coverage.best_presentation_pct == pytest.approx(0.28)

    scores = score_predictions(
        [epitope],
        EpitopeConfig(score_expr="presentation[mhcflurry].score"))
    assert list(scores) == [pytest.approx(0.85)]


def test_load_lens_preserves_allele_independent_processing_once():
    """Repeated LENS allele rows yield one canonical processing leaf."""
    path = os.path.join(
        DATA_DIR, "real_lens_subsets", "lens_v1.9_real_subset.tsv")
    _loaded = read_lens_report(path)
    report_df, epitopes = _loaded.report_df, list(_loaded.epitopes)
    epitope = epitopes[0]
    processing = epitope.predictions_for(
        "antigen_processing", predictor="mhcflurry")

    assert len(processing) == 1
    assert processing[0].allele == ""
    assert processing[0].score == pytest.approx(
        report_df["mhcflurry processing score"].iloc[0])
    assert all(epitope.per_allele_scores)


def test_lens_processing_dsl_scores_every_patient_allele(tmp_path):
    """Canonical processing evidence participates in allele-scoped DSL groups."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import (
        attach_per_allele_scores,
        epitopes_to_topiary_df,
    )

    path = tmp_path / "multi-allele-processing.tsv"
    path.write_text(
        "peptide\tallele\tpep_context\tmhcflurry_2.1.1.aff\t"
        "mhcflurry_2.1.1.proc_score\n"
        "SIINFEKL\tHLA-A02:01\tXXSIINFEKLXX\t50\t0.8\n"
        "SIINFEKL\tHLA-B07:02\tXXSIINFEKLXX\t60\t0.8\n")

    _loaded = read_lens_report(path)
    epitopes = list(_loaded.epitopes)
    assert len(epitopes) == 1
    epitope = epitopes[0]
    processing = epitope.predictions_for(
        "antigen_processing", predictor="mhcflurry")
    assert len(processing) == 1
    assert processing[0].allele == ""

    topiary_df = epitopes_to_topiary_df(epitopes)
    processing_rows = topiary_df[
        topiary_df["kind"] == "antigen_processing"]
    assert set(processing_rows["allele"]) == {
        "HLA-A*02:01", "HLA-B*07:02"}

    scored = attach_per_allele_scores(
        epitopes,
        EpitopeConfig(score_expr="processing[mhcflurry].score"))
    assert scored[0].per_allele_scores == {
        "HLA-A*02:01": pytest.approx(0.8),
        "HLA-B*07:02": pytest.approx(0.8),
    }


def test_lens_processing_only_rows_preserve_alleles_scores_and_report_rows(
        tmp_path):
    """Processing-only LENS evidence stays allele-scoped end to end."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import (
        attach_per_allele_scores,
        epitopes_to_topiary_df,
    )
    from vaxrank.report import TemplateDataCreator, epitope_report_row_inputs

    path = tmp_path / "processing-only.tsv"
    path.write_text(
        "peptide\tallele\tpep_context\tmhcflurry_2.1.1.aff\t"
        "mhcflurry_2.1.1.proc_score\n"
        "SIINFEKL\tHLA-A02:01\tXXSIINFEKLXX\t\t0.8\n"
        "SIINFEKL\tHLA-B07:02\tXXSIINFEKLXX\t\t0.8\n")

    _loaded = read_lens_report(path)
    report_df, epitopes = _loaded.report_df, list(_loaded.epitopes)
    assert len(epitopes) == 1
    epitope = epitopes[0]
    assert epitope.patient_alleles == (
        "HLA-A*02:01", "HLA-B*07:02")
    processing = epitope.predictions_for(
        "antigen_processing", predictor="mhcflurry")
    assert len(processing) == 1
    assert processing[0].allele == ""

    topiary_df = epitopes_to_topiary_df(epitopes)
    assert set(topiary_df["allele"]) == {
        "HLA-A*02:01", "HLA-B*07:02"}
    config = EpitopeConfig(
        score_expr="processing[mhcflurry].score")
    scored = attach_per_allele_scores(epitopes, config)
    assert scored[0].per_allele_scores == {
        "HLA-A*02:01": pytest.approx(0.8),
        "HLA-B*07:02": pytest.approx(0.8),
    }

    csv_path = tmp_path / "scored.csv"
    write_neoepitope_report(
        report_df,
        epitopes,
        csv_report_path=str(csv_path),
        epitope_config=config,
    )
    import pandas as pd
    assert list(pd.read_csv(csv_path)["vaxrank_score"]) == [0.8, 0.8]

    creator = TemplateDataCreator.__new__(TemplateDataCreator)
    creator.processing_predictions_by_key = {}
    row_inputs = epitope_report_row_inputs(scored[0])
    rows = [
        creator.epitope_data(
            scored[0], row_input.prediction,
            allele=row_input.allele,
            include_additional_prediction_axes=True)
        for row_input in row_inputs
    ]
    assert [row["Allele"] for row in rows] == ["A*02:01", "B*07:02"]
    assert [row["Score"] for row in rows] == ["0.8", "0.8"]
    assert all(row["Integrated processing score"] == "0.800" for row in rows)


def test_lens_presentation_only_candidate_has_a_template_report_row(tmp_path):
    """Presentation evidence anchors a table row when affinity is absent."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import attach_per_allele_scores
    from vaxrank.report import TemplateDataCreator, epitope_report_row_inputs

    path = tmp_path / "presentation-only.tsv"
    path.write_text(
        "peptide\tallele\tpep_context\tmhcflurry_2.1.1.aff\t"
        "mhcflurry_2.1.1.pres_score\tmhcflurry_2.1.1.pres_perc\n"
        "SIINFEKL\tHLA-A02:01\tXXSIINFEKLXX\t\t0.85\t0.28\n")

    _loaded = read_lens_report(path)
    epitopes = list(_loaded.epitopes)
    config = EpitopeConfig(
        score_expr="presentation[mhcflurry].score")
    epitope = attach_per_allele_scores(epitopes, config)[0]
    row_inputs = epitope_report_row_inputs(epitope)
    assert len(row_inputs) == 1
    assert row_inputs[0].prediction.kind == "pMHC_presentation"
    assert row_inputs[0].allele == "HLA-A*02:01"

    creator = TemplateDataCreator.__new__(TemplateDataCreator)
    creator.processing_predictions_by_key = {}
    row = creator.epitope_data(
        epitope,
        row_inputs[0].prediction,
        allele=row_inputs[0].allele,
        include_additional_prediction_axes=True,
    )
    assert row["Allele"] == "A*02:01"
    assert row["Score"] == "0.85"
    assert row["IC50"] == "No prediction"
    assert row["Presentation score"] == "0.850"
    assert row["Presentation %ile"] == "0.280"


def test_lens_report_anchors_mixed_evidence_per_allele_and_predictor(tmp_path):
    """Each allele/predictor gets its strongest available report anchor."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import attach_per_allele_scores
    from vaxrank.report import TemplateDataCreator, epitope_report_row_inputs

    path = tmp_path / "mixed-evidence.tsv"
    path.write_text(
        "peptide\tallele\tpep_context\tmhcflurry_2.1.1.aff\t"
        "mhcflurry_2.1.1.pres_score\tmhcflurry_2.1.1.proc_score\t"
        "netmhcpan_4.1.aff_nm\n"
        "SIINFEKL\tHLA-A02:01\tXXSIINFEKLXX\t50\t0.8\t0.7\t\n"
        "SIINFEKL\tHLA-B07:02\tXXSIINFEKLXX\t\t0.6\t0.7\t75\n"
        "SIINFEKL\tHLA-C07:02\tXXSIINFEKLXX\t\t\t0.7\t\n")

    _loaded = read_lens_report(path)
    epitopes = list(_loaded.epitopes)
    epitope = attach_per_allele_scores(
        epitopes,
        EpitopeConfig(score_expr="processing[mhcflurry].score"),
    )[0]
    row_inputs = epitope_report_row_inputs(epitope)
    assert [
        (row.allele, row.prediction.predictor_name, row.prediction.kind)
        for row in row_inputs
    ] == [
        ("HLA-A*02:01", "mhcflurry", "pMHC_affinity"),
        ("HLA-B*07:02", "netmhcpan", "pMHC_affinity"),
        ("HLA-B*07:02", "mhcflurry", "pMHC_presentation"),
        ("HLA-C*07:02", "mhcflurry", "antigen_processing"),
    ]

    creator = TemplateDataCreator.__new__(TemplateDataCreator)
    creator.processing_predictions_by_key = {}
    report_rows = [
        creator.epitope_data(
            epitope,
            row_input.prediction,
            allele=row_input.allele,
            include_additional_prediction_axes=True,
        )
        for row_input in row_inputs
    ]
    by_key = {
        (row["Allele"], row["Predictor"]): row for row in report_rows}
    assert by_key[("A*02:01", "mhcflurry")]["IC50"] == "50.00 nM"
    assert by_key[("A*02:01", "mhcflurry")][
        "Presentation score"] == "0.800"
    assert by_key[("B*07:02", "netmhcpan")]["IC50"] == "75.00 nM"
    assert by_key[("B*07:02", "mhcflurry")][
        "Presentation score"] == "0.600"
    assert by_key[("C*07:02", "mhcflurry")][
        "Integrated processing score"] == "0.700"
    assert all(row["Score"] == "0.7" for row in report_rows)


def test_lens_context_scores_merge_by_source_position(tmp_path):
    """Equal peptide/allele rows retain context-specific integrated scores."""
    from vaxrank.epitope_config import EpitopeConfig

    path = tmp_path / "context-specific-scores.tsv"
    path.write_text(
        "peptide\tallele\tpep_context\tmhcflurry_2.1.1.aff\t"
        "mhcflurry_2.1.1.pres_score\tmhcflurry_2.1.1.proc_score\n"
        "SIINFEKL\tHLA-A02:01\tXXSIINFEKLAA\t50\t0.2\t0.3\n"
        "SIINFEKL\tHLA-A02:01\tYYYYSIINFEKLZ\t60\t0.9\t0.8\n")

    _loaded = read_lens_report(path)
    report_df, epitopes = _loaded.report_df, list(_loaded.epitopes)
    assert list(report_df["Source sequence name"]) == [
        "XXSIINFEKLAA", "YYYYSIINFEKLZ"]
    assert list(report_df["Peptide offset"]) == [2, 4]

    csv_path = tmp_path / "scored.csv"
    write_neoepitope_report(
        report_df,
        epitopes,
        csv_report_path=str(csv_path),
        epitope_config=EpitopeConfig(score_expr=(
            "presentation[mhcflurry].score + "
            "processing[mhcflurry].score")),
    )

    import pandas as pd
    result = pd.read_csv(csv_path).set_index("Source sequence name")
    assert result.loc["XXSIINFEKLAA", "vaxrank_score"] == pytest.approx(0.5)
    assert result.loc["YYYYSIINFEKLZ", "vaxrank_score"] == pytest.approx(1.7)


def test_load_lens_reports_conflicting_processing_scores_without_aborting(
        tmp_path, caplog):
    """A disagreement is surfaced, not fatal.

    An allele-independent score repeated across allele rows should agree —
    it depends on peptide and flanks, not on MHC. When it doesn't, aborting
    the load gives the operator nothing to act on and costs them the whole
    run; the usual cause is a file printing one value at two precisions.
    """
    import logging

    path = tmp_path / "conflicting-processing.tsv"
    path.write_text(
        "peptide\tallele\tpep_context\tmhcflurry_2.1.1.aff\t"
        "mhcflurry_2.1.1.proc_score\n"
        "SIINFEKL\tHLA-A02:01\tXXSIINFEKLXX\t50\t0.8\n"
        "SIINFEKL\tHLA-B07:02\tXXSIINFEKLXX\t60\t0.4\n")

    with caplog.at_level(logging.WARNING):
        report = read_lens_report(path)

    [epitope] = report.epitopes
    [processing] = epitope.predictions_for(
        "antigen_processing", predictor="mhcflurry")
    # First-seen wins, so the result is deterministic in file order.
    assert processing.score == pytest.approx(0.8)
    assert epitope.patient_alleles == ("HLA-A*02:01", "HLA-B*07:02")

    warning = "\n".join(
        record.getMessage()
        for record in caplog.records if record.levelno >= logging.WARNING)
    assert "processing" in warning
    # The operator needs the peptide and both values to chase it upstream.
    assert "SIINFEKL" in warning
    assert "0.8" in warning and "0.4" in warning


def test_load_lens_tolerates_reprinted_processing_scores(tmp_path, caplog):
    """The same value at two precisions is not a disagreement."""
    import logging

    path = tmp_path / "reprinted-processing.tsv"
    path.write_text(
        "peptide\tallele\tpep_context\tmhcflurry_2.1.1.aff\t"
        "mhcflurry_2.1.1.proc_score\n"
        "SIINFEKL\tHLA-A02:01\tXXSIINFEKLXX\t50\t0.8\n"
        "SIINFEKL\tHLA-B07:02\tXXSIINFEKLXX\t60\t0.8000\n")

    with caplog.at_level(logging.WARNING):
        report = read_lens_report(path)

    [epitope] = report.epitopes
    [processing] = epitope.predictions_for(
        "antigen_processing", predictor="mhcflurry")
    assert processing.score == pytest.approx(0.8)
    assert not [
        record for record in caplog.records
        if record.levelno >= logging.WARNING
        and "processing score" in record.getMessage()]


def test_epitope_report_surfaces_integrated_mhcflurry_axes():
    """Clinical rows distinguish integrated MHCflurry and pepsickle scores.

    Mixed-predictor tables retain identical columns; the netMHCpan affinity
    row gets explicit placeholders because those integrated axes came from
    MHCflurry rather than netMHCpan.
    """
    from vaxrank.report import TemplateDataCreator

    path = os.path.join(
        DATA_DIR, "real_lens_subsets", "lens_v1.4_real_subset.tsv")
    _loaded = read_lens_report(path)
    epitopes = list(_loaded.epitopes)
    epitope = epitopes[0]
    mhcflurry_affinity = epitope.predictions_for(
        "pMHC_affinity", predictor="mhcflurry")[0]
    netmhcpan_affinity = epitope.predictions_for(
        "pMHC_affinity", predictor="netmhcpan")[0]
    creator = TemplateDataCreator.__new__(TemplateDataCreator)
    creator.processing_predictions_by_key = {}

    mhcflurry_row = creator.epitope_data(
        epitope, mhcflurry_affinity,
        include_additional_prediction_axes=True)
    netmhcpan_row = creator.epitope_data(
        epitope, netmhcpan_affinity,
        include_additional_prediction_axes=True)

    assert mhcflurry_row["Presentation score"] != "—"
    assert mhcflurry_row["Presentation %ile"] != "—"
    assert mhcflurry_row["Integrated processing score"] != "—"
    assert netmhcpan_row["Presentation score"] == "—"
    assert netmhcpan_row["Presentation %ile"] == "—"
    assert netmhcpan_row["Integrated processing score"] == "—"
    assert tuple(mhcflurry_row) == tuple(netmhcpan_row)


# ── DSL integration for external inputs ──────────────────────────────────────

def test_default_methods_resolves_multi_model_affinity(tmp_path):
    """When a file exposes multiple affinity models and no score_expr is
    set, cfg.default_methods picks which one the default DSL node scores
    on. Without default_methods, vaxrank auto-picks canonical."""
    from vaxrank.epitope_config import EpitopeConfig

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    _loaded = read_lens_report(path)
    report_df, preds = _loaded.report_df, list(_loaded.epitopes)

    # Explicit default = netmhcpan → score reflects netmhcpan affinity.
    cfg = EpitopeConfig(default_methods={"pMHC_affinity": "netmhcpan"})
    csv_netmhcpan = tmp_path / "netmhcpan.csv"
    write_neoepitope_report(
        report_df.copy(), preds,
        csv_report_path=str(csv_netmhcpan), epitope_config=cfg)

    # Explicit default = mhcflurry → score reflects mhcflurry affinity.
    cfg2 = EpitopeConfig(default_methods={"pMHC_affinity": "mhcflurry"})
    csv_mhcflurry = tmp_path / "mhcflurry.csv"
    write_neoepitope_report(
        report_df.copy(), preds,
        csv_report_path=str(csv_mhcflurry), epitope_config=cfg2)

    import pandas as pd
    r1 = pd.read_csv(csv_netmhcpan).set_index("Mutant peptide sequence")
    r2 = pd.read_csv(csv_mhcflurry).set_index("Mutant peptide sequence")
    # At least one row should differ (the two predictors disagree enough
    # that logistic-normalized output diverges).
    diffs = (r1["vaxrank_score"] - r2["vaxrank_score"]).abs()
    assert (diffs > 1e-6).any(), (
        "Expected different scores when default_methods picks different "
        "affinity predictor")


def test_default_methods_rejects_unknown_method(tmp_path):
    """cfg.default_methods[kind] must match a method actually in the data."""
    from vaxrank.epitope_config import EpitopeConfig

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    _loaded = read_lens_report(path)
    report_df, preds = _loaded.report_df, list(_loaded.epitopes)
    cfg = EpitopeConfig(default_methods={"pMHC_affinity": "mhcnuggets"})
    with pytest.raises(ValueError, match="default_methods.*mhcnuggets"):
        write_neoepitope_report(
            report_df, preds, csv_report_path=str(tmp_path / "x.csv"),
            epitope_config=cfg)


def test_default_methods_partial_qualification_mixed_kinds(tmp_path):
    """score_expr brackets a method for one Kind but leaves another Kind
    unqualified. The bracketed Kind keeps its explicit method; the other
    Kind falls back to its default (auto-picked or user-specified)."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(DATA_DIR, "lens_v1.4_with_stability.tsv")
    _loaded = read_lens_report(path)
    report_df, preds = _loaded.report_df, list(_loaded.epitopes)
    # Formula: affinity qualified (netmhcpan), stability unqualified.
    # default_methods settles the stability kind.
    cfg = EpitopeConfig(
        score_expr=(
            "affinity['netmhcpan'].logistic(350, 150) * 0.7 + "
            "(stability.value / 24).clip(0, 1) * 0.3"
        ),
        default_methods={"pMHC_stability": "netmhcstabpan"},
    )
    csv_path = tmp_path / "out.csv"
    write_neoepitope_report(
        report_df, preds, csv_report_path=str(csv_path), epitope_config=cfg)
    import pandas as pd
    result = pd.read_csv(csv_path)
    assert len(result) == 2
    assert (result["vaxrank_score"] > 0).all()


def test_default_methods_partial_qualification_auto_picks_stability(tmp_path):
    """Same as above but without default_methods — auto-pick fills in
    the stability default."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(DATA_DIR, "lens_v1.4_with_stability.tsv")
    _loaded = read_lens_report(path)
    report_df, preds = _loaded.report_df, list(_loaded.epitopes)
    # Only one stability method (netmhcstabpan) in this fixture, so
    # auto-pick is trivial but tests the path.
    cfg = EpitopeConfig(
        score_expr=(
            "affinity['mhcflurry'].logistic(350, 150) * 0.7 + "
            "(stability.value / 24).clip(0, 1) * 0.3"
        ),
    )
    csv_path = tmp_path / "out.csv"
    write_neoepitope_report(
        report_df, preds, csv_report_path=str(csv_path), epitope_config=cfg)
    import pandas as pd
    result = pd.read_csv(csv_path)
    assert (result["vaxrank_score"] > 0).all()


def test_filter_expr_unqualified_on_multi_method_data(tmp_path):
    """Unqualified filter_expr ('affinity <= X') on a file with multiple
    affinity methods should resolve via the default method, not raise
    "Ambiguous Kind"."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    _loaded = read_lens_report(path)
    report_df, preds = _loaded.report_df, list(_loaded.epitopes)
    # fixture's mhcflurry affinities: 95.4, 180.2, 310.5. Cut at 200 keeps
    # 2 rows; with netmhcpan as the default, cut at 200 would keep rows
    # with IC50 < 200 from netmhcpan: 76.11, 120.50 (2 rows).
    cfg = EpitopeConfig(
        filter_expr="affinity <= 200",
        # auto-pick selects mhcflurry; test it evaluates without error
    )
    csv_path = tmp_path / "out.csv"
    write_neoepitope_report(
        report_df, preds, csv_report_path=str(csv_path), epitope_config=cfg)
    import pandas as pd
    result = pd.read_csv(csv_path)
    # Filter drops rows whose mhcflurry affinity > 200 (only 95.4 passes),
    # so the nonzero-score rows number ≤ 1. Other rows still land in the
    # report with score 0.
    nonzero = result[result["vaxrank_score"] > 0]
    assert len(nonzero) <= len(result)
    assert len(nonzero) >= 1


def test_score_predictions_on_empty_predictions_list(tmp_path):
    """An empty predictions list should produce an empty score Series
    without blowing up."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import score_predictions
    scores = score_predictions([], EpitopeConfig())
    assert len(scores) == 0


def test_resolve_default_methods_noop_on_single_method_data():
    """Single-model Kinds need no default (topiary has only one choice);
    resolve_default_methods omits them from the result."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import (
        epitopes_to_topiary_df, resolve_default_methods)

    path = os.path.join(DATA_DIR, "lens_v1.9_mhcflurry_only.tsv")
    _loaded = read_lens_report(path)
    preds = list(_loaded.epitopes)
    df = epitopes_to_topiary_df(preds)
    assert resolve_default_methods(EpitopeConfig(), df) == {}


def test_resolve_default_methods_auto_picks_canonical():
    """Multi-model Kind with no user override → auto-pick mhcflurry for
    pMHC_affinity."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import (
        epitopes_to_topiary_df, resolve_default_methods)

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    _loaded = read_lens_report(path)
    preds = list(_loaded.epitopes)
    df = epitopes_to_topiary_df(preds)
    assert resolve_default_methods(EpitopeConfig(), df) == {
        "pMHC_affinity": "mhcflurry"}


def test_resolve_default_methods_user_override_wins():
    """Explicit cfg.default_methods beats auto-pick."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import (
        epitopes_to_topiary_df, resolve_default_methods)

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    _loaded = read_lens_report(path)
    preds = list(_loaded.epitopes)
    df = epitopes_to_topiary_df(preds)
    cfg = EpitopeConfig(default_methods={"pMHC_affinity": "netmhcpan"})
    assert resolve_default_methods(cfg, df) == {"pMHC_affinity": "netmhcpan"}


def test_default_methods_works_on_synthetic_main_pipeline_frame(tmp_path):
    """Multi-model resolution doesn't depend on LENS specifically — it's a
    DataFrame-shape concern. A synthetic topiary frame mimicking what the
    main VCF/BAM pipeline would produce if someone passed multiple models
    (e.g. ``TopiaryPredictor(models=[mhcflurry, netmhcpan])``) must score
    cleanly via ``default_methods``."""
    import pandas as pd
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import score_predictions

    # Build a minimal topiary-shaped frame with two models per group —
    # exactly what we'd get from a multi-model TopiaryPredictor in the
    # main pipeline, no LENS loader involved. ``prediction_id`` is present
    # because ``epitopes_to_topiary_df`` always emits it: the score group is
    # an explicit identity, never inferred from the sequence.
    df = pd.DataFrame([
        # group 1: SIINFEKL × HLA-A*02:01, two affinity predictors
        {"prediction_id": "ovalbumin", "source_sequence_name": "ovalbumin",
         "peptide": "SIINFEKL", "peptide_offset": 5, "peptide_length": 8,
         "allele": "HLA-A*02:01", "n_flank": "", "c_flank": "",
         "prediction_method_name": "mhcflurry", "predictor_version": "2.1.1",
         "kind": "pMHC_affinity", "value": 80.0, "affinity": 80.0,
         "percentile_rank": 0.4, "score": 0.0},
        {"prediction_id": "ovalbumin", "source_sequence_name": "ovalbumin",
         "peptide": "SIINFEKL", "peptide_offset": 5, "peptide_length": 8,
         "allele": "HLA-A*02:01", "n_flank": "", "c_flank": "",
         "prediction_method_name": "netmhcpan", "predictor_version": "4.1b",
         "kind": "pMHC_affinity", "value": 120.0, "affinity": 120.0,
         "percentile_rank": 0.7, "score": 0.0},
        # group 2: another peptide+allele
        {"prediction_id": "p53", "source_sequence_name": "p53",
         "peptide": "STPPPGTRV", "peptide_offset": 10, "peptide_length": 9,
         "allele": "HLA-B*07:02", "n_flank": "", "c_flank": "",
         "prediction_method_name": "mhcflurry", "predictor_version": "2.1.1",
         "kind": "pMHC_affinity", "value": 250.0, "affinity": 250.0,
         "percentile_rank": 1.5, "score": 0.0},
        {"prediction_id": "p53", "source_sequence_name": "p53",
         "peptide": "STPPPGTRV", "peptide_offset": 10, "peptide_length": 9,
         "allele": "HLA-B*07:02", "n_flank": "", "c_flank": "",
         "prediction_method_name": "netmhcpan", "predictor_version": "4.1b",
         "kind": "pMHC_affinity", "value": 400.0, "affinity": 400.0,
         "percentile_rank": 2.5, "score": 0.0},
    ])

    # Default config (no formulas, no default_methods) — auto-pick should
    # still drive scoring without raising "Ambiguous Kind".
    cfg = EpitopeConfig()
    scores = score_predictions([], cfg, topiary_df=df)
    assert len(scores) == 2  # two unique (peptide, allele) groups
    # Every score should be a finite float; auto-picked mhcflurry scores
    # both groups against logistic_normalized(350, 150).
    assert all(0 <= s <= 1 for s in scores.values)

    # Explicit default = netmhcpan should give different (lower) scores
    # since netmhcpan's IC50s are higher in this fake data.
    cfg2 = EpitopeConfig(default_methods={"pMHC_affinity": "netmhcpan"})
    scores2 = score_predictions([], cfg2, topiary_df=df)
    for key in scores.index:
        assert scores2[key] < scores[key], (
            f"Expected netmhcpan score < mhcflurry score for {key}; "
            f"got {scores2[key]} vs {scores[key]}")


def test_resolve_default_methods_stability_plus_affinity():
    """Multi-method affinity + single-method stability: only affinity
    needs a default."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import (
        epitopes_to_topiary_df, resolve_default_methods)

    path = os.path.join(DATA_DIR, "lens_v1.4_with_stability.tsv")
    _loaded = read_lens_report(path)
    preds = list(_loaded.epitopes)
    df = epitopes_to_topiary_df(preds)
    resolved = resolve_default_methods(EpitopeConfig(), df)
    # pMHC_affinity has mhcflurry + netmhcpan → needs default
    assert resolved.get("pMHC_affinity") == "mhcflurry"
    # pMHC_stability has only netmhcstabpan → omitted
    assert "pMHC_stability" not in resolved


def test_parse_is_cached():
    """Repeated expression parses should hit the cache
    rather than re-invoking topiary.ranking.parse."""
    from vaxrank.epitope_dsl import parse_epitope_expression
    cache_info_before = parse_epitope_expression.cache_info()
    expr = "affinity['mhcflurry'].logistic(350, 150)"
    parse_epitope_expression(expr)
    parse_epitope_expression(expr)
    parse_epitope_expression(expr)
    info = parse_epitope_expression.cache_info()
    # At least two cache hits from the three calls above
    assert info.hits - cache_info_before.hits >= 2


def test_default_methods_mixed_valid_and_invalid_fails_fast(tmp_path):
    """If default_methods has one valid and one invalid entry, the invalid
    one must fail fast even though other Kinds resolve cleanly."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(DATA_DIR, "lens_v1.4_with_stability.tsv")
    _loaded = read_lens_report(path)
    report_df, preds = _loaded.report_df, list(_loaded.epitopes)
    cfg = EpitopeConfig(default_methods={
        "pMHC_affinity": "mhcflurry",      # valid
        "pMHC_stability": "doesnotexist",  # invalid
    })
    with pytest.raises(ValueError,
                       match=r"default_methods\['pMHC_stability'\]='doesnotexist'"):
        write_neoepitope_report(
            report_df, preds, csv_report_path=str(tmp_path / "x.csv"),
            epitope_config=cfg)


def test_report_columns_expose_all_detected_predictors(tmp_path):
    """Scoring may target one predictor via ``default_methods`` or a
    bracketed formula, but the report DataFrame must still expose every
    detected predictor's raw value column."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    _loaded = read_lens_report(path)
    report_df, preds = _loaded.report_df, list(_loaded.epitopes)
    cfg = EpitopeConfig(
        score_expr="affinity['mhcflurry'].logistic(350, 150)")
    csv_path = tmp_path / "out.csv"
    write_neoepitope_report(
        report_df.copy(), preds,
        csv_report_path=str(csv_path), epitope_config=cfg)
    import pandas as pd
    result = pd.read_csv(csv_path)
    # Even though only mhcflurry drives scoring, netmhcpan values are
    # still reported alongside.
    assert "mhcflurry value (nM)" in result.columns
    assert "netmhcpan value (nM)" in result.columns
    assert result["netmhcpan value (nM)"].notna().all()


def test_predictor_version_in_epitope_prediction_from_lens():
    """Every LENS-loaded leaf ``Prediction`` should carry the sniffed
    predictor version, so version-qualified DSL refs resolve."""
    path = os.path.join(DATA_DIR, "lens_example.tsv")
    _loaded = read_lens_report(path)
    epitopes = list(_loaded.epitopes)
    leaves = [p for e in epitopes for p in e.predictions_flat()]
    by_method = {p.predictor_name: p.predictor_version for p in leaves}
    assert by_method["mhcflurry"] == "2.1.1"
    assert by_method["netmhcpan"] == "4.1b"


def test_dsl_version_mismatch_validation_error(tmp_path):
    """A version-qualified formula against data without that version
    raises a clear error (not silent NaN)."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    _loaded = read_lens_report(path)
    report_df, preds = _loaded.report_df, list(_loaded.epitopes)
    cfg = EpitopeConfig(
        score_expr="affinity['mhcflurry', '0.0.0'].logistic(350, 150)")
    with pytest.raises(ValueError, match="predictor version '0.0.0'"):
        write_neoepitope_report(
            report_df, preds, csv_report_path=str(tmp_path / "x.csv"),
            epitope_config=cfg)


def test_validate_default_methods_empty_df_is_noop():
    """validate_default_methods on an empty topiary frame returns cleanly
    (we can't check methods we don't have)."""
    import pandas as pd
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import validate_default_methods
    cfg = EpitopeConfig(default_methods={"pMHC_affinity": "mhcflurry"})
    validate_default_methods(cfg, pd.DataFrame())  # no raise


def test_validate_default_methods_no_defaults_is_noop():
    """With no default_methods, the validator is a no-op."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import (
        epitopes_to_topiary_df, validate_default_methods)
    path = os.path.join(DATA_DIR, "lens_example.tsv")
    _loaded = read_lens_report(path)
    preds = list(_loaded.epitopes)
    df = epitopes_to_topiary_df(preds)
    validate_default_methods(EpitopeConfig(), df)  # no raise


def test_validate_default_methods_catches_typo_on_single_method_kind(tmp_path):
    """Even when a Kind has only one model in the data (so subsetting
    wouldn't consult the default), validate_default_methods still
    errors on a typo'd default."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import (
        epitopes_to_topiary_df, validate_default_methods)
    # Single-model fixture (v1.9 mhcflurry-only).
    path = os.path.join(DATA_DIR, "lens_v1.9_mhcflurry_only.tsv")
    _loaded = read_lens_report(path)
    preds = list(_loaded.epitopes)
    df = epitopes_to_topiary_df(preds)
    cfg = EpitopeConfig(default_methods={"pMHC_affinity": "bogus"})
    with pytest.raises(ValueError,
                       match=r"default_methods\['pMHC_affinity'\]='bogus'"):
        validate_default_methods(cfg, df)


def test_filter_plus_score_expr_both_bracketed(tmp_path):
    """filter_expr and score_expr both set, each with bracketed predictor
    refs — the union of referenced methods should survive subsetting."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    _loaded = read_lens_report(path)
    report_df, preds = _loaded.report_df, list(_loaded.epitopes)
    cfg = EpitopeConfig(
        filter_expr="affinity['netmhcpan'] <= 500",
        score_expr="affinity['mhcflurry'].logistic(350, 150)",
    )
    csv_path = tmp_path / "out.csv"
    write_neoepitope_report(
        report_df, preds, csv_report_path=str(csv_path), epitope_config=cfg)
    import pandas as pd
    result = pd.read_csv(csv_path)
    # Filter keeps rows whose netmhcpan IC50 <= 500 (all 3 in fixture:
    # 76.11, 120.50, 250.00); score comes from mhcflurry.
    assert (result["vaxrank_score"] > 0).any()


def test_filter_expr_drops_all_groups_marks_audit_rows_without_fake_scores(
        tmp_path):
    """Filtered audit rows are explicit and never acquire synthetic zeros."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    _loaded = read_lens_report(path)
    report_df, preds = _loaded.report_df, list(_loaded.epitopes)
    # No affinity < 0 anywhere.
    cfg = EpitopeConfig(filter_expr="affinity <= 0")
    csv_path = tmp_path / "out.csv"
    write_neoepitope_report(
        report_df, preds, csv_report_path=str(csv_path), epitope_config=cfg)
    import pandas as pd
    result = pd.read_csv(csv_path)
    assert result["vaxrank_score"].isna().all()
    assert not result["vaxrank_filter_passed"].any()
    assert not result["vaxrank_rank_eligible"].any()
    assert set(result["vaxrank_exclusion_reason"]) == {"dsl_filter"}


def test_report_distinguishes_below_minimum_from_dsl_filter(tmp_path):
    """A real score below the ranking gate remains visible and classified."""
    from vaxrank.epitope_config import EpitopeConfig

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    _loaded = read_lens_report(path)
    report_df, preds = _loaded.report_df, list(_loaded.epitopes)
    csv_path = tmp_path / "below-minimum.csv"
    write_neoepitope_report(
        report_df,
        preds,
        csv_report_path=str(csv_path),
        epitope_config=EpitopeConfig(min_epitope_score=2.0),
    )
    import pandas as pd
    result = pd.read_csv(csv_path)

    assert result["vaxrank_score"].notna().all()
    assert result["vaxrank_filter_passed"].all()
    assert not result["vaxrank_rank_eligible"].any()
    assert set(result["vaxrank_exclusion_reason"]) == {"min_epitope_score"}


def test_mixed_bracketed_and_unqualified_evaluates_cleanly(tmp_path):
    """Topiary 5.10's ``EvalContext(default_methods=...)`` handles
    bracketed + unqualified refs to the same Kind: bracketed refs use
    their specific method, unqualified refs use the default. Previously
    vaxrank had to preempt this with its own error; now it just works."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    _loaded = read_lens_report(path)
    report_df, preds = _loaded.report_df, list(_loaded.epitopes)
    # Brackets netmhcpan AND references affinity unqualified. Default
    # auto-picks mhcflurry for unqualified refs. Result: netmhcpan_value +
    # mhcflurry_value per row, evaluated cleanly.
    cfg = EpitopeConfig(
        score_expr="affinity['netmhcpan'].value + affinity.value")
    csv_path = tmp_path / "out.csv"
    write_neoepitope_report(
        report_df, preds, csv_report_path=str(csv_path), epitope_config=cfg)
    import pandas as pd
    result = pd.read_csv(csv_path)
    assert len(result) == 3
    # Verify the sum: for the SVVGSSSSS row, netmhcpan=76.11 + mhcflurry=95.4
    svv = result[result["Mutant peptide sequence"] == "SVVGSSSSS"].iloc[0]
    assert svv["vaxrank_score"] == pytest.approx(76.11 + 95.4)


def test_pvacseq_single_method_per_group_has_no_ambiguity(tmp_path):
    """pVACseq emits at most one row per (peptide, allele). Even when
    different rows use different methods, each group is
    unambiguous — the default node should score without error."""
    from vaxrank.epitope_io import write_neoepitope_report
    path = os.path.join(DATA_DIR, "pvacseq_example.tsv")
    _loaded = read_pvacseq_report(path)
    report_df, preds = _loaded.report_df, list(_loaded.epitopes)
    csv_path = tmp_path / "out.csv"
    # Default EpitopeConfig, no default_methods set.
    write_neoepitope_report(
        report_df, preds, csv_report_path=str(csv_path))
    import pandas as pd
    result = pd.read_csv(csv_path)
    assert len(result) == 3
    assert "vaxrank_score" in result.columns


def test_unknown_column_in_lens_file_is_ignored(tmp_path):
    """LENS files with extra columns that don't match the predictor
    regex are simply ignored — loader shouldn't error."""
    # Valid mhcflurry columns plus a novel experimental column.
    path = tmp_path / "lens.tsv"
    path.write_text(
        "peptide\tallele\tmhcflurry_2.1.1.aff\t"
        "mhcflurry_2.1.1.pres_perc\tfuture_tool_99.blarg.metric\t"
        "random_non_matching_col\n"
        "SIINFEKL\tHLA-A*02:01\t100.0\t0.5\t42\thello\n"
    )
    _loaded = read_lens_report(path)
    epitopes = list(_loaded.epitopes)
    assert len(epitopes) == 1
    leaf = _first_leaf(epitopes[0])
    assert leaf.predictor_name == "mhcflurry"


def test_default_methods_auto_picked_when_unset(tmp_path, caplog):
    """When multiple models for a Kind are present and default_methods
    doesn't specify one, vaxrank auto-picks canonical and logs which."""
    from vaxrank.epitope_config import EpitopeConfig
    import logging

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    _loaded = read_lens_report(path)
    report_df, preds = _loaded.report_df, list(_loaded.epitopes)
    cfg = EpitopeConfig()  # no default_methods

    with caplog.at_level(logging.INFO, logger="vaxrank.epitope_dsl"):
        write_neoepitope_report(
            report_df, preds,
            csv_report_path=str(tmp_path / "out.csv"), epitope_config=cfg)

    assert any(
        "canonical pick is 'mhcflurry'" in rec.getMessage()
        for rec in caplog.records), (
        f"Expected auto-pick info log; got: "
        f"{[r.getMessage() for r in caplog.records]}")


# ── NaN handling in vaxrank-native roundtrip ─────────────────────────────────

def test_load_predictions_empty_string_fields_become_empty_not_nan(tmp_path):
    """CSV cells that read as NaN for string-typed fields should come
    back as '' (not float NaN) so downstream comparisons work."""
    from vaxrank.epitope_io import load_predictions
    path = tmp_path / "preds.csv"
    # Empty wt_peptide_sequence, source_sequence, prediction_method_name,
    # predictor_version columns all present but blank.
    path.write_text(
        "allele,peptide_sequence,wt_peptide_sequence,ic50,wt_ic50,"
        "percentile_rank,prediction_method_name,overlaps_mutation,"
        "source_sequence,offset,occurs_in_reference,predictor_version\n"
        "HLA-A*02:01,SIINFEKL,,50.0,,0.5,,True,,0,False,\n"
    )
    epitopes = load_predictions(path)
    assert len(epitopes) == 1
    e = epitopes[0]
    # WT comparator is absent because wt_ic50 was blank → no WT signal.
    assert e.wt is None
    assert e.source_sequence == ""
    leaf = _first_leaf(e)
    assert leaf.predictor_name == ""
    assert leaf.predictor_version == ""


def test_write_neoepitope_report_scores_duplicate_peptides_per_source(tmp_path):
    """The same (peptide, allele) from two sources gets two scores, not one.

    Real LENS / pVACseq files emit identical peptide-allele pairs from
    alternative transcripts, homologous regions, and distinct variants. Those
    are different biological candidates with different evidence, so each row
    joins on its own prediction identity. Broadcasting one score across them —
    the pre-identity behavior — silently attributed one source's evidence to
    another.
    """
    import pandas as pd
    from vaxrank.epitope_io import write_neoepitope_report

    report_df = pd.DataFrame([
        {'Allele': 'HLA-A*02:01', 'Mutant peptide sequence': 'SIINFEKL',
         'Genomic variant': 'chr1:100', 'Gene name': 'GENEA',
         'Prediction identity': 'src-a', 'Peptide offset': 2},
        {'Allele': 'HLA-A*02:01', 'Mutant peptide sequence': 'SIINFEKL',
         'Genomic variant': 'chr2:200', 'Gene name': 'GENEB',
         'Prediction identity': 'src-b', 'Peptide offset': 2},
    ])
    # Same peptide and allele, different sources, different affinities.
    preds = [
        dataclasses.replace(
            _make_prediction(peptide_sequence='SIINFEKL',
                             allele='HLA-A*02:01', ic50=50.0),
            prediction_id='src-a'),
        dataclasses.replace(
            _make_prediction(peptide_sequence='SIINFEKL',
                             allele='HLA-A*02:01', ic50=4000.0),
            prediction_id='src-b'),
    ]
    csv_path = tmp_path / "out.csv"
    write_neoepitope_report(report_df, preds, csv_report_path=str(csv_path))
    result = pd.read_csv(csv_path)

    assert len(result) == 2
    assert sorted(result['Genomic variant'].tolist()) == [
        'chr1:100', 'chr2:200']
    by_source = dict(zip(result['Prediction identity'],
                         result['vaxrank_score']))
    # The strong binder scores above the weak one; neither borrows the other's.
    assert by_source['src-a'] > by_source['src-b']


def test_write_neoepitope_report_requires_prediction_identity(tmp_path):
    """A frame without identity columns is refused, not silently broadcast."""
    import pandas as pd
    from vaxrank.epitope_io import write_neoepitope_report

    report_df = pd.DataFrame([
        {'Allele': 'HLA-A*02:01', 'Mutant peptide sequence': 'SIINFEKL',
         'Genomic variant': 'chr1:100', 'Gene name': 'GENEA'},
    ])
    preds = [_make_prediction()]
    with pytest.raises(ValueError, match="Prediction identity"):
        write_neoepitope_report(
            report_df, preds, csv_report_path=str(tmp_path / "out.csv"))


def test_lens_dsl_combines_both_predictors(tmp_path):
    """A score_expr combining two predictors' affinities should average them."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    _loaded = read_lens_report(path)
    report_df, preds = _loaded.report_df, list(_loaded.epitopes)

    # 50/50 blend of the two predictors' logistic-normalized affinity scores
    cfg = EpitopeConfig(
        score_expr=(
            "affinity['mhcflurry'].logistic(350, 150) * 0.5 + "
            "affinity['netmhcpan'].logistic(350, 150) * 0.5"
        )
    )
    csv_path = tmp_path / "report.csv"
    write_neoepitope_report(
        report_df, preds, csv_report_path=str(csv_path), epitope_config=cfg)

    import pandas as pd
    result = pd.read_csv(csv_path)
    assert len(result) == 3
    # Non-zero score means the DSL actually picked up both predictors
    assert (result["vaxrank_score"] > 0).any()


def test_lens_dsl_validation_catches_missing_predictor(tmp_path):
    """Formula referencing a predictor not in the file should error clearly."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(DATA_DIR, "lens_v1.9_mhcflurry_only.tsv")
    _loaded = read_lens_report(path)
    report_df, preds = _loaded.report_df, list(_loaded.epitopes)

    cfg = EpitopeConfig(
        score_expr="affinity['netmhcpan'].logistic(350, 150)"
    )
    with pytest.raises(ValueError, match="predictor 'netmhcpan'"):
        write_neoepitope_report(
            report_df, preds, csv_report_path=str(tmp_path / "unused.csv"),
            epitope_config=cfg)


def test_lens_dsl_validation_catches_missing_kind(tmp_path):
    """Formula using a Kind not present in the data (e.g. stability on a
    file without a stability predictor) should error clearly."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(DATA_DIR, "lens_v1.9_mhcflurry_only.tsv")
    _loaded = read_lens_report(path)
    report_df, preds = _loaded.report_df, list(_loaded.epitopes)
    cfg = EpitopeConfig(score_expr="stability['netmhcstabpan'].value")
    with pytest.raises(ValueError, match="references kind 'pMHC_stability'"):
        write_neoepitope_report(
            report_df, preds, csv_report_path=str(tmp_path / "unused.csv"),
            epitope_config=cfg)


def test_lens_dsl_validation_catches_missing_version(tmp_path):
    """Formula pinned to a specific predictor version not in the data should
    error clearly."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(DATA_DIR, "lens_v1.9_mhcflurry_only.tsv")
    _loaded = read_lens_report(path)
    report_df, preds = _loaded.report_df, list(_loaded.epitopes)
    cfg = EpitopeConfig(
        score_expr="affinity['mhcflurry', '9.9.9'].logistic(350, 150)"
    )
    with pytest.raises(ValueError, match="predictor version '9.9.9'"):
        write_neoepitope_report(
            report_df, preds, csv_report_path=str(tmp_path / "unused.csv"),
            epitope_config=cfg)


def test_lens_dsl_version_qualified_formula(tmp_path):
    """A formula pinned to a specific predictor version should match
    leaf predictions carrying that version from LENS column sniffing."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    _loaded = read_lens_report(path)
    report_df, epitopes = _loaded.report_df, list(_loaded.epitopes)
    leaves = [p for e in epitopes for p in e.predictions_flat()]
    # Every LENS leaf should carry predictor_version from columns.
    assert all(p.predictor_version for p in leaves)
    mhcflurry_version = next(
        p.predictor_version for p in leaves if p.predictor_name == "mhcflurry")
    assert mhcflurry_version == "2.1.1"

    cfg = EpitopeConfig(
        score_expr=(
            "affinity['mhcflurry', '2.1.1'].logistic(350, 150)"
        )
    )
    csv_path = tmp_path / "report.csv"
    write_neoepitope_report(
        report_df, epitopes, csv_report_path=str(csv_path), epitope_config=cfg)
    import pandas as pd
    result = pd.read_csv(csv_path)
    assert (result["vaxrank_score"] > 0).any()


# ── Integration tests against real-LENS subsamples ───────────────────────────
# Fixtures in real_lens_subsets/ are 30-row stratified samples from actual
# LENS reports (v1.4-dev, v1.5.1, v1.9-dev) covering the full column surface
# (94 / 106 / 108 cols). They exercise the real schema including ERV / CTA /
# FUSION / SPLICE / INDEL antigen sources that synthetic fixtures miss.

REAL_LENS_DIR = os.path.join(DATA_DIR, "real_lens_subsets")


def _leaves(epitopes):
    return [p for e in epitopes for p in e.predictions_flat()]


def test_real_lens_v14_detects_all_three_predictors():
    """v1.4-dev files emit mhcflurry + netmhcpan + netmhcstabpan. All
    three should be detected and emit leaf predictions inside the
    loaded Epitopes."""
    path = os.path.join(REAL_LENS_DIR, "lens_v1.4_real_subset.tsv")
    _loaded = read_lens_report(path)
    report_df, epitopes = _loaded.report_df, list(_loaded.epitopes)
    leaves = _leaves(epitopes)
    methods = {p.predictor_name for p in leaves}
    assert methods == {"mhcflurry", "netmhcpan", "netmhcstabpan"}
    # Every real LENS leaf carries a nonempty version.
    assert all(p.predictor_version for p in leaves)
    # Report exposes per-tool columns with correct units.
    assert "mhcflurry value (nM)" in report_df.columns
    assert "netmhcpan value (nM)" in report_df.columns
    assert "netmhcstabpan value (hours)" in report_df.columns


def test_real_lens_v15_emits_both_affinity_predictors():
    """v1.5.1 has both mhcflurry and netmhcpan; load_lens collects all
    per-predictor predictions inside each CandidateEpitope. Canonical
    'mhcflurry' drives the report's display columns."""
    path = os.path.join(REAL_LENS_DIR, "lens_v1.5_real_subset.tsv")
    _loaded = read_lens_report(path)
    report_df, epitopes = _loaded.report_df, list(_loaded.epitopes)
    leaves = _leaves(epitopes)
    methods = {p.predictor_name for p in leaves}
    assert methods == {"mhcflurry", "netmhcpan"}
    # MHCflurry contributes affinity + presentation + processing while
    # netMHCpan contributes affinity. NA cells produce fewer leaves.
    assert len(leaves) <= 4 * len(report_df)
    predictors_used = report_df["Predictors used"].iloc[0]
    assert predictors_used.startswith("mhcflurry")


def test_real_lens_v19_is_mhcflurry_only():
    """v1.9-dev emits only mhcflurry (netMHCpan was dropped)."""
    path = os.path.join(REAL_LENS_DIR, "lens_v1.9_real_subset.tsv")
    _loaded = read_lens_report(path)
    epitopes = list(_loaded.epitopes)
    leaves = _leaves(epitopes)
    methods = {p.predictor_name for p in leaves}
    assert methods == {"mhcflurry"}
    # Real data has HLA alleles in the un-asterisked form; loader normalizes
    # every allele-specific leaf. Processing is intentionally allele-less.
    allele_specific = [p for p in leaves if p.allele]
    assert all(
        p.allele.startswith("HLA-") and "*" in p.allele
        for p in allele_specific)
    assert any(p.kind == "antigen_processing" and not p.allele for p in leaves)


def test_real_lens_dsl_scoring_all_versions(tmp_path):
    """DSL scoring should produce nonzero scores for at least some rows in
    each real LENS fixture (with all predictors emitted). A 50/50 blend of
    mhcflurry + netmhcpan where both present; plain mhcflurry otherwise."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    blend = EpitopeConfig(score_expr=(
        "affinity['mhcflurry'].logistic(350, 150) * 0.5 + "
        "affinity['netmhcpan'].logistic(350, 150) * 0.5"))
    solo = EpitopeConfig(score_expr="affinity['mhcflurry'].logistic(350, 150)")

    cases = [
        ("lens_v1.4_real_subset.tsv", blend),
        ("lens_v1.5_real_subset.tsv", blend),
        ("lens_v1.9_real_subset.tsv", solo),  # no netmhcpan
    ]

    import pandas as pd
    for name, cfg in cases:
        path = os.path.join(REAL_LENS_DIR, name)
        _loaded = read_lens_report(path)
        report_df, preds = _loaded.report_df, list(_loaded.epitopes)
        csv_path = tmp_path / f"{name}.csv"
        write_neoepitope_report(
            report_df, preds, csv_report_path=str(csv_path),
            epitope_config=cfg)
        result = pd.read_csv(csv_path)
        assert len(result) == len(report_df), f"{name}: row count changed"
        assert (result["vaxrank_score"] > 0).any(), (
            f"{name}: no rows scored > 0 under {cfg.score_expr!r}")


def test_real_lens_v14_stability_formula(tmp_path):
    """v1.4 real fixture has netmhcstabpan — verify combined
    affinity + stability formula evaluates successfully."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(REAL_LENS_DIR, "lens_v1.4_real_subset.tsv")
    _loaded = read_lens_report(path)
    report_df, preds = _loaded.report_df, list(_loaded.epitopes)
    cfg = EpitopeConfig(score_expr=(
        "affinity['mhcflurry'].logistic(350, 150) * 0.7 + "
        "(stability['netmhcstabpan'].value / 24).clip(0, 1) * 0.3"))
    csv_path = tmp_path / "report.csv"
    write_neoepitope_report(
        report_df, preds, csv_report_path=str(csv_path), epitope_config=cfg)
    import pandas as pd
    result = pd.read_csv(csv_path)
    assert (result["vaxrank_score"] > 0).any()


def test_real_lens_version_qualified_auto_detects(tmp_path):
    """predictor_version sniffed from real LENS column names matches what
    a version-qualified DSL formula expects."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(REAL_LENS_DIR, "lens_v1.5_real_subset.tsv")
    _loaded = read_lens_report(path)
    report_df, epitopes = _loaded.report_df, list(_loaded.epitopes)
    leaves = _leaves(epitopes)
    versions = {(p.predictor_name, p.predictor_version) for p in leaves}
    # Real v1.5.1 file: mhcflurry 2.1.1 + netmhcpan 4.1b
    assert ("mhcflurry", "2.1.1") in versions
    assert ("netmhcpan", "4.1b") in versions

    cfg = EpitopeConfig(score_expr=(
        "affinity['mhcflurry', '2.1.1'].logistic(350, 150)"))
    csv_path = tmp_path / "report.csv"
    write_neoepitope_report(
        report_df, epitopes, csv_report_path=str(csv_path), epitope_config=cfg)
    import pandas as pd
    result = pd.read_csv(csv_path)
    assert (result["vaxrank_score"] > 0).any()


def test_lens_dsl_stability_and_affinity_combined(tmp_path):
    """Formula combining pMHC_affinity and pMHC_stability Kinds should
    evaluate when both predictors are emitted."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(DATA_DIR, "lens_v1.4_with_stability.tsv")
    _loaded = read_lens_report(path)
    report_df, epitopes = _loaded.report_df, list(_loaded.epitopes)
    leaves = [p for e in epitopes for p in e.predictions_flat()]
    kinds_seen = {p.kind for p in leaves}
    assert kinds_seen == {
        "pMHC_affinity", "pMHC_presentation", "pMHC_stability"}

    cfg = EpitopeConfig(
        score_expr=(
            "affinity['mhcflurry'].logistic(350, 150) * 0.5 + "
            # stability: bound hours / 24h window, clipped to [0, 1]
            "(stability['netmhcstabpan'].value / 24).clip(0, 1) * 0.5"
        )
    )
    csv_path = tmp_path / "report.csv"
    write_neoepitope_report(
        report_df, epitopes, csv_report_path=str(csv_path), epitope_config=cfg)
    import pandas as pd
    result = pd.read_csv(csv_path)
    assert len(result) == 2
    assert (result["vaxrank_score"] > 0).all()


def _assert_default_scores_match_legacy(cfg, path, tmp_path):
    """Shared helper: DSL-routed default scoring should equal the
    legacy per-prediction logistic score on the canonical (mhcflurry)
    leaf records selected by default_methods."""
    from vaxrank.epitope_io import write_neoepitope_report
    from tests._legacy_score_reference import legacy_score_one as _legacy_score_one

    _loaded = read_lens_report(path)
    report_df, epitopes = _loaded.report_df, list(_loaded.epitopes)
    # Auto-picked default for pMHC_affinity is mhcflurry — that's what
    # resolves unqualified Affinity refs in topiary's EvalContext.
    legacy_scores = {}
    for e in epitopes:
        for p in e.predictions_flat():
            if p.kind != 'pMHC_affinity' or p.predictor_name != 'mhcflurry':
                continue
            legacy_scores[(e.sequence, p.allele)] = round(
                _legacy_score_one(
                    ic50=p.value,
                    percentile_rank=p.percentile_rank,
                    midpoint=cfg.logistic_epitope_score_midpoint,
                    width=cfg.logistic_epitope_score_width,
                    ic50_cutoff=cfg.binding_affinity_cutoff,
                    scoring_mode=cfg.scoring_mode,
                    percentile_rank_cutoff=cfg.percentile_rank_cutoff,
                ), 6)
    csv_path = tmp_path / "report.csv"
    write_neoepitope_report(
        report_df, epitopes, csv_report_path=str(csv_path), epitope_config=cfg)
    import pandas as pd
    result = pd.read_csv(csv_path)
    for _, row in result.iterrows():
        key = (row["Mutant peptide sequence"], row["Allele"])
        assert row["vaxrank_score"] == pytest.approx(legacy_scores[key])


def test_lens_percentile_rank_scoring_matches_legacy(tmp_path):
    """scoring_mode='percentile_rank' through the DSL path matches legacy
    logistic_epitope_score computed on the canonical (mhcflurry) rows."""
    from vaxrank.epitope_config import EpitopeConfig
    _assert_default_scores_match_legacy(
        EpitopeConfig(scoring_mode="percentile_rank"),
        os.path.join(DATA_DIR, "lens_example.tsv"),
        tmp_path,
    )


def test_default_scoring_matches_legacy(tmp_path):
    """With no filter_expr/score_expr the DSL default node matches legacy
    logistic_epitope_score byte-for-byte on the canonical predictor rows."""
    from vaxrank.epitope_config import EpitopeConfig
    _assert_default_scores_match_legacy(
        EpitopeConfig(),
        os.path.join(DATA_DIR, "lens_example.tsv"),
        tmp_path,
    )


# ── Report generation from imports ───────────────────────────────────────────

def test_write_neoepitope_report_csv(tmp_path):
    path = os.path.join(DATA_DIR, "pvacseq_example.tsv")
    _loaded = read_pvacseq_report(path)
    report_df, preds = _loaded.report_df, list(_loaded.epitopes)
    csv_path = tmp_path / "report.csv"
    write_neoepitope_report(report_df, preds, csv_report_path=str(csv_path))
    assert csv_path.exists()

    import pandas as pd
    result = pd.read_csv(csv_path)
    assert "rank" in result.columns
    assert "vaxrank_score" in result.columns
    assert len(result) == 3
    # Should be sorted by score descending
    assert result["vaxrank_score"].iloc[0] >= result["vaxrank_score"].iloc[1]


def test_write_neoepitope_report_xlsx(tmp_path):
    path = os.path.join(DATA_DIR, "lens_example.tsv")
    _loaded = read_lens_report(path)
    report_df, preds = _loaded.report_df, list(_loaded.epitopes)
    xlsx_path = tmp_path / "report.xlsx"
    write_neoepitope_report(report_df, preds, excel_report_path=str(xlsx_path))
    assert xlsx_path.exists()


# ── Scoring modes ────────────────────────────────────────────────────────────
# The legacy per-prediction scorer now lives as a free function in
# ``vaxrank.vaccine_peptide`` (``_legacy_score_one``); these tests pin
# its scoring-mode contract directly against IC50 / percentile_rank.

def test_scoring_mode_affinity():
    from tests._legacy_score_reference import legacy_score_one as _legacy_score_one
    assert _legacy_score_one(100.0, 0.5, scoring_mode="affinity") > 0.5


def test_scoring_mode_percentile_rank():
    from tests._legacy_score_reference import legacy_score_one as _legacy_score_one
    score = _legacy_score_one(100.0, 0.5, scoring_mode="percentile_rank")
    assert score == pytest.approx(0.95)


def test_scoring_mode_percentile_rank_weak():
    from tests._legacy_score_reference import legacy_score_one as _legacy_score_one
    score = _legacy_score_one(100.0, 10.0, scoring_mode="percentile_rank")
    assert score == 0.0


def test_scoring_mode_percentile_rank_none():
    from tests._legacy_score_reference import legacy_score_one as _legacy_score_one
    score = _legacy_score_one(100.0, None, scoring_mode="percentile_rank")
    assert score == 0.0


# ── CLI integration test ─────────────────────────────────────────────────────

def test_pvacseq_cli_csv_output(tmp_path):
    """Test the --input-pvacseq flag produces CSV output."""
    from vaxrank.cli.entry_point import main
    pvacseq_path = os.path.join(DATA_DIR, "pvacseq_example.tsv")
    csv_path = str(tmp_path / "output.csv")
    main([
        "--input-pvacseq", pvacseq_path,
        "--output-csv", csv_path,
    ])
    assert os.path.exists(csv_path)
    import pandas as pd
    df = pd.read_csv(csv_path)
    assert len(df) == 3
    assert "vaxrank_score" in df.columns


def test_lens_cli_csv_output(tmp_path):
    """Test the --input-lens flag produces CSV output."""
    from vaxrank.cli.entry_point import main
    lens_path = os.path.join(DATA_DIR, "lens_example.tsv")
    csv_path = str(tmp_path / "output.csv")
    main([
        "--input-lens", lens_path,
        "--output-csv", csv_path,
    ])
    assert os.path.exists(csv_path)
    import pandas as pd
    df = pd.read_csv(csv_path)
    assert len(df) == 3


def test_external_arg_parser_accepts_ensembl_release():
    """``--ensembl-release N`` is reachable from the external-input
    parser (regression: the flag was originally pipeline-only,
    leaving LENS-driven template reports without a way to resolve
    transcript IDs)."""
    from vaxrank.cli.arg_parser import parse_vaxrank_args
    args = parse_vaxrank_args([
        '--input-lens', '/dev/null',
        '--output-csv', '/dev/null',
        '--ensembl-release', '75',
    ])
    assert args.ensembl_release == 75


def test_external_path_resolves_genome_from_ensembl_release():
    """End-to-end plumbing: parse + ``resolve_ensembl_release`` on
    an external-input args namespace must leave ``args.genome`` as a
    real ``pyensembl.EnsemblRelease`` so transcript resolution can
    use it."""
    import pyensembl
    from vaxrank.cli.arg_parser import parse_vaxrank_args
    from vaxrank.cli.entry_point import resolve_ensembl_release
    args = parse_vaxrank_args([
        '--input-lens', '/dev/null',
        '--output-csv', '/dev/null',
        '--ensembl-release', '75',
    ])
    resolve_ensembl_release(args)
    assert isinstance(args.genome, pyensembl.EnsemblRelease)
    assert args.genome.release == 75


def test_lens_cli_writes_ascii_report():
    """ASCII / HTML / PDF template reports work on the external
    (LENS / pVACseq) path now that transcript IDs resolve and a
    ``PatientInfo`` is synthesized from the loaded antigens
    (closes #268). Writes to a tempfile and asserts the patient-info
    section + at least one variant block rendered."""
    import tempfile
    from vaxrank.cli.entry_point import main
    lens_path = os.path.join(DATA_DIR, "lens_example.tsv")
    with tempfile.TemporaryDirectory() as td:
        out = os.path.join(td, "report.txt")
        main([
            "--input-lens", lens_path,
            "--output-ascii-report", out,
            # Required by ``check_args`` pre-flight when external
            # input is paired with template reports — the release
            # value doesn't have to resolve every transcript here,
            # just be present so the run isn't rejected up front.
            "--ensembl-release", "75",
        ])
        assert os.path.exists(out), "ASCII report not written"
        with open(out) as f:
            text = f.read()
    # Patient-info header populated from the synthesized PatientInfo
    assert "Total number of somatic variants:" in text
    # At least one variant block — exact count depends on the fixture
    # but we want the loop to have run, not be silently empty.
    assert "Vaccine antigens:" in text


def test_lens_cli_errors_when_template_report_missing_ensembl_release():
    """Pre-flight check: external input + template-report flag with
    no --ensembl-release must fail fast (was a late warning that
    silently degraded reports to empty effect annotations)."""
    import pytest
    from vaxrank.cli.entry_point import main
    lens_path = os.path.join(DATA_DIR, "lens_example.tsv")
    with pytest.raises(ValueError, match="--ensembl-release"):
        main([
            "--input-lens", lens_path,
            "--output-ascii-report", "/tmp/should-not-be-written.txt",
        ])


def test_lens_cli_errors_when_no_output_flag_set():
    """Running ``vaxrank --input-lens FILE`` with no --output-* flag
    must fail fast — every output is opt-in via its own flag and the
    earlier behavior (run to completion, write nothing, log a quiet
    "wrote=['(none)']") is exactly the surprise the guard prevents."""
    import pytest
    from vaxrank.cli.entry_point import main
    lens_path = os.path.join(DATA_DIR, "lens_example.tsv")
    with pytest.raises(ValueError, match="No output path specified"):
        main(["--input-lens", lens_path])


def test_neoepitope_core_row_shared_contract():
    """The shared core-row builder used by all three input paths
    (VCF pipeline, LENS, pVACseq): consistent columns + blank for
    missing affinity, and Score only when provided."""
    from vaxrank.epitope_io import neoepitope_core_row
    # Missing affinities → blank, present → 'NN.NN nM'.
    row = neoepitope_core_row(
        allele="HLA-A*02:01", mutant_peptide="SIINFEKL",
        mutant_affinity=123.456, wt_peptide="", wt_affinity=None,
        gene_name="TP53", variant="17:7675088 C>T")
    assert row["Allele"] == "HLA-A*02:01"
    assert row["Predicted mutant pMHC affinity"] == "123.46 nM"
    assert row["Predicted wildtype pMHC affinity"] == ""   # blank, not "No prediction"
    assert "Score" not in row                              # omitted when None
    # Score is inserted only when supplied (pipeline path).
    scored = neoepitope_core_row(
        allele="A", mutant_peptide="P", mutant_affinity=None,
        wt_peptide="", wt_affinity=None, gene_name="G", variant="v",
        score=0.5)
    assert scored["Score"] == 0.5
    assert scored["Predicted mutant pMHC affinity"] == ""


def test_native_roundtrip_preserves_targetability_and_self_reference(tmp_path):
    """Safety decisions must survive a save/load cycle.

    ``overlaps_targetable`` and ``self_reference_match`` both fall back to the
    permissive legacy behavior when absent, so dropping them on write did not
    fail loudly — it silently re-admitted a held-out epitope as a vaccine
    target on reload.
    """
    from mhctools.pred import Prediction

    from vaxrank.candidate_epitope import (
        SOURCE_CLASS_MUTATION, CandidateEpitope,
    )
    from vaxrank.epitope_io import load_predictions, save_predictions
    from vaxrank.vaccine_antigen import SelfReferenceMatch

    prediction = Prediction(
        kind="pMHC_affinity", predictor_name="mhcflurry",
        predictor_version="2.1.1", allele="HLA-A*02:01",
        peptide="SIINFEKL", value=50.0, score=0.9, percentile_rank=0.4)
    epitope = CandidateEpitope(
        sequence="SIINFEKL", n_flank="AAA", c_flank="CCC",
        source_sequence="XXSIINFEKLXX", source_name="ENST00000123456",
        offset=2, predictions=(prediction,),
        source_class=SOURCE_CLASS_MUTATION, overlaps_mutation=True,
        overlaps_targetable=False,
        self_reference_match=SelfReferenceMatch(
            peptide="SIINFEKL", occurs=False, antigen_kind="CTA",
            excluded_gene_ids=("ENSG00000141510",),
            genome_release="GRCh38"),
        patient_alleles=("HLA-A*02:01",), prediction_id="pid-1",
        per_allele_scores={"HLA-A*02:01": 0.9})

    path = tmp_path / "native.csv"
    save_predictions([epitope], path)
    [reloaded] = load_predictions(path)

    assert reloaded.overlaps_targetable is False
    assert reloaded.is_targetable is False
    assert reloaded.self_reference_match == epitope.self_reference_match
    assert reloaded.occurs_in_self_reference is False
    # Flanks drive processing predictors; a reload without them cannot
    # reproduce the scores of the file it read.
    assert (reloaded.n_flank, reloaded.c_flank) == ("AAA", "CCC")
    assert reloaded.source_name == "ENST00000123456"


def test_native_roundtrip_keeps_undecided_targetability_undecided(tmp_path):
    """A blank cell reloads as "no decision", never as False."""
    from mhctools.pred import Prediction

    from vaxrank.candidate_epitope import CandidateEpitope
    from vaxrank.epitope_io import load_predictions, save_predictions

    epitope = CandidateEpitope(
        sequence="SIINFEKL", source_sequence="SIINFEKL", offset=0,
        predictions=(Prediction(
            kind="pMHC_affinity", predictor_name="mhcflurry",
            predictor_version="2.1.1", allele="HLA-A*02:01",
            peptide="SIINFEKL", value=50.0, score=0.9),),
        prediction_id="pid-2")
    path = tmp_path / "undecided.csv"
    save_predictions([epitope], path)
    [reloaded] = load_predictions(path)

    assert reloaded.overlaps_targetable is None
    assert reloaded.is_targetable is True
