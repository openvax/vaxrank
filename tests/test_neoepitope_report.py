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

from vaxrank.epitope_prediction import EpitopePrediction
from vaxrank.report import make_minimal_neoepitope_report


def _dummy_variant():
    return SimpleNamespace(short_description="chr1:1 A>G")


def _dummy_vaccine_peptide(epitope_predictions):
    fragment = SimpleNamespace(gene_name="GENE1", n_alt_reads=10)
    return SimpleNamespace(
        mutant_epitope_predictions=epitope_predictions,
        mutant_epitope_score=1.23,
        mutant_protein_fragment=fragment,
    )


def test_neoepitope_report_limits_epitopes_per_peptide(tmp_path):
    prediction1 = EpitopePrediction(
        allele="HLA-A*02:01",
        peptide_sequence="ABCDEFGHI",
        wt_peptide_sequence="ABCDEFGHI",
        ic50=100.0,
        wt_ic50=200.0,
        percentile_rank=0.5,
        prediction_method_name="test",
        overlaps_mutation=True,
        source_sequence="ABCDEFGHI",
        offset=0,
        occurs_in_reference=False,
    )
    prediction2 = EpitopePrediction(
        allele="HLA-A*02:01",
        peptide_sequence="BCDEFGHIJ",
        wt_peptide_sequence="BCDEFGHIJ",
        ic50=150.0,
        wt_ic50=250.0,
        percentile_rank=0.4,
        prediction_method_name="test",
        overlaps_mutation=True,
        source_sequence="BCDEFGHIJ",
        offset=0,
        occurs_in_reference=False,
    )
    vaccine_peptide = _dummy_vaccine_peptide([prediction1, prediction2])
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
    prediction = EpitopePrediction(
        allele="HLA-A*02:01",
        peptide_sequence="ABCDEFGHI",
        wt_peptide_sequence="ABCDEFGHI",
        ic50=100.0,
        wt_ic50=None,
        percentile_rank=0.5,
        prediction_method_name="test",
        overlaps_mutation=True,
        source_sequence="ABCDEFGHI",
        offset=0,
        occurs_in_reference=False,
    )
    vaccine_peptide = _dummy_vaccine_peptide([prediction])
    excel_path = tmp_path / "neoepitopes_missing_wt.xlsx"

    make_minimal_neoepitope_report(
        ranked_variants_with_vaccine_peptides=[(_dummy_variant(), [vaccine_peptide])],
        num_epitopes_per_peptide=None,
        excel_report_path=str(excel_path),
    )

    df = pd.read_excel(excel_path, engine="openpyxl")
    assert df.iloc[0]["Predicted wildtype pMHC affinity"] == "No prediction"
