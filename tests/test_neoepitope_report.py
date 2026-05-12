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
        mutant_epitopes=epitopes,
        mutant_epitope_score=1.23,
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
    assert df.iloc[0]["Predicted wildtype pMHC affinity"] == "No prediction"
