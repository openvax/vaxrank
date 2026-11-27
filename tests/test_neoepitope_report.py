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
Tests for the minimal neoepitope XLSX report.
"""

from types import SimpleNamespace

import pandas as pd
from mhctools.pred import Prediction

from vaxrank.candidate_epitope import COMPARATOR_WT, CandidateEpitope, Peptide
from vaxrank.report import make_minimal_neoepitope_report


def _dummy_variant():
    return SimpleNamespace(short_description="chr1:1 A>G")


def _dummy_vaccine_peptide(epitopes):
    fragment = SimpleNamespace(gene_name="GENE1", n_alt_reads=10)
    return SimpleNamespace(
        target_epitopes=epitopes,
        target_epitope_score=1.23,
        mutant_protein_fragment=fragment,
    )


def _make_epitope(peptide, ic50, wt_ic50, allele="HLA-A*02:01",
                  wt_peptide=None):
    mutant_pred = Prediction(
        kind='pMHC_affinity', predictor_name='test',
        predictor_version='', allele=allele, peptide=peptide,
        value=ic50, score=0.0, percentile_rank=0.5)
    comparators = {}
    if wt_ic50 is not None:
        wt_pred = Prediction(
            kind='pMHC_affinity', predictor_name='test',
            predictor_version='', allele=allele,
            peptide=wt_peptide or peptide,
            value=wt_ic50, score=0.0, percentile_rank=None)
        comparators[COMPARATOR_WT] = Peptide(
            sequence=wt_peptide or peptide,
            predictions=(wt_pred,))
    return CandidateEpitope.from_peptide(
        Peptide(
            sequence=peptide,
            source_sequence=peptide, offset=0,
            predictions=(mutant_pred,)),
        comparators=comparators,
        overlaps_mutation=True, occurs_in_reference=False)


def test_neoepitope_report_limits_epitopes_per_peptide(tmp_path):
    e1 = _make_epitope("ABCDEFGHI", ic50=100.0, wt_ic50=200.0)
    e2 = _make_epitope("BCDEFGHIJ", ic50=150.0, wt_ic50=250.0)
    vaccine_peptide = _dummy_vaccine_peptide([e1, e2])
    excel_path = tmp_path / "neoepitopes.xlsx"

    make_minimal_neoepitope_report(
        ranked_variants_with_vaccine_peptides=[(_dummy_variant(), [vaccine_peptide])],
        num_epitopes_per_peptide=1,
        excel_report_path=str(excel_path),
    )

    df = pd.read_excel(excel_path, engine="openpyxl")
    assert len(df) == 1
    assert df.iloc[0]["Mutant peptide sequence"] == "ABCDEFGHI"


def test_neoepitope_report_handles_missing_wt_ic50(tmp_path):
    e = _make_epitope("ABCDEFGHI", ic50=100.0, wt_ic50=None)
    vaccine_peptide = _dummy_vaccine_peptide([e])
    excel_path = tmp_path / "neoepitopes_missing_wt.xlsx"

    make_minimal_neoepitope_report(
        ranked_variants_with_vaccine_peptides=[(_dummy_variant(), [vaccine_peptide])],
        num_epitopes_per_peptide=None,
        excel_report_path=str(excel_path),
    )

    df = pd.read_excel(excel_path, engine="openpyxl")
    # Missing WT affinity is blank (read back as NaN) — the same
    # missing-value convention all three report paths now share via
    # neoepitope_core_row / _format_nm.
    assert pd.isna(df.iloc[0]["Predicted wildtype pMHC affinity"])


def test_xlsx_reports_store_numeric_cells_with_explicit_number_formats(tmp_path):
    """#459: Excel workbooks store scores and nM affinities as numeric cells with formats."""
    import openpyxl
    import pytest
    from vaxrank.report import make_csv_report

    e = _make_epitope("SIINFEKL", ic50=123.456, wt_ic50=5000.5)
    vaccine_peptide = _dummy_vaccine_peptide([e])
    variant = _dummy_variant()

    neo_xlsx = tmp_path / "neoepitopes.xlsx"
    make_minimal_neoepitope_report(
        ranked_variants_with_vaccine_peptides=[(variant, [vaccine_peptide])],
        num_epitopes_per_peptide=None,
        excel_report_path=str(neo_xlsx),
    )
    neo_ws = openpyxl.load_workbook(neo_xlsx)["Neoepitopes"]
    neo_headers = {cell.value: idx for idx, cell in enumerate(neo_ws[1], 1)}
    mut_cell = neo_ws.cell(row=2, column=neo_headers["Predicted mutant pMHC affinity"])
    wt_cell = neo_ws.cell(row=2, column=neo_headers["Predicted wildtype pMHC affinity"])
    score_cell = neo_ws.cell(row=2, column=neo_headers["Score"])
    assert mut_cell.data_type == "n"
    assert mut_cell.value == pytest.approx(123.46, abs=0.01)
    assert mut_cell.number_format == '0.00" nM"'
    assert wt_cell.data_type == "n"
    assert wt_cell.value == pytest.approx(5000.5, abs=0.01)
    assert wt_cell.number_format == '0.00" nM"'
    assert score_cell.data_type == "n"
    assert score_cell.number_format == "0.0000"

    from vaxrank.vaccine_peptide import VaccinePeptide
    from .test_assembly_antigens import candidate

    antigen, _ = candidate(0, sequence="SIINFEKL")
    real_vp = VaccinePeptide(antigen=antigen, epitopes=[e])
    ranked_xlsx = tmp_path / "ranked.xlsx"
    make_csv_report(
        ranked_variants_with_vaccine_peptides=[(antigen, [real_vp])],
        excel_report_path=str(ranked_xlsx),
    )
    all_ws = openpyxl.load_workbook(ranked_xlsx)["All"]
    all_headers = [cell.value for cell in all_ws[1]]
    assert all_headers[0] == "source_rank"
    assert None not in all_headers
    col_map = {name: idx for idx, name in enumerate(all_headers, 1)}
    comb_cell = all_ws.cell(row=2, column=col_map["combined_score"])
    gravy_cell = all_ws.cell(row=2, column=col_map["max_7mer_gravy_score"])
    cys_cell = all_ws.cell(row=2, column=col_map["cysteine_count"])
    assert comb_cell.data_type == "n"
    assert comb_cell.number_format == "0.0000"
    assert gravy_cell.data_type == "n"
    assert gravy_cell.number_format == "0.00"
    assert cys_cell.data_type == "n"
    assert cys_cell.number_format == "0"
