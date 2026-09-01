# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""What ``n_alt_reads`` holds, and how a consumer can tell.

The field feeds the combined-score DSL, where ``sqrt(n_alt_reads) *
target_epitope_score`` is the documented canonical form. It holds a
different quantity per input path, and these tests pin that each one says
which (openvax/vaxrank#379, #382, #383).
"""

import os
import types

from vaxrank.epitope_config import EpitopeConfig
from vaxrank.epitope_dsl import epitopes_for_ranking
from vaxrank.epitope_io import read_pvacseq_report
from vaxrank.external_input import pvacseq_ranking_result
from topiary import FRAGMENTS, READS

from vaxrank.mutant_protein_fragment import (
    READ_COUNT_METHOD_RNA_READS,
    SEQUENCE_SOURCE_ISOVAR_ASSEMBLY,
    SEQUENCE_SOURCE_VARCODE_TRANSLATION,
    MutantProteinFragment,
)

DATA_DIR = os.path.join(
    os.path.dirname(__file__), "data", "epitope_fixtures")


def _fake_isovar_result():
    """Minimal stand-in exposing only what from_isovar_result reads."""
    protein_sequence = types.SimpleNamespace(
        gene_name="BRAF", amino_acids="SIINFEKLAA",
        mutation_start_idx=2, mutation_end_idx=3,
        num_supporting_fragments=7, transcripts=[])
    return types.SimpleNamespace(
        variant=object(), top_protein_sequence=protein_sequence,
        num_total_fragments=20, num_alt_fragments=8, num_ref_fragments=12)


def test_isovar_counts_are_labelled_as_fragments():
    """The field says reads and the value is fragments, so the label must say.

    For paired-end data a fragment is two reads, so this is roughly half
    the read count for the same evidence — a different quantity from a
    LENS depth x VAF estimate, not a different way of obtaining one.
    Nothing converts, because the conversion needs to know whether the
    library was paired-end and isovar does not say.
    """
    fragment = MutantProteinFragment.from_isovar_result(_fake_isovar_result())

    assert fragment.n_alt_reads == 8            # fragments, per isovar
    # Two axes, not one. rna_reads says the count came from an alignment
    # and deliberately fixes no subject; the subject says what was counted.
    assert fragment.read_count_method == READ_COUNT_METHOD_RNA_READS
    assert fragment.read_count_subject == FRAGMENTS
    assert fragment.sequence_source == SEQUENCE_SOURCE_ISOVAR_ASSEMBLY


def test_the_sequence_source_terms_are_read_from_topiary():
    """No private copy of a controlled vocabulary.

    A rename upstream should break the import rather than silently stamp
    fragments with a term nothing recognizes.
    """
    from topiary import SEQUENCE_SOURCES

    assert SEQUENCE_SOURCE_ISOVAR_ASSEMBLY in SEQUENCE_SOURCES
    assert SEQUENCE_SOURCE_VARCODE_TRANSLATION in SEQUENCE_SOURCES


def test_a_count_cannot_be_read_in_the_wrong_subject():
    """Asking for the other unit returns nothing, never a conversion.

    vaxrank asked topiary for an ``rna_fragments`` *method* and that was
    the wrong shape (openvax/topiary#239): a subject is not a derivation,
    and a term alongside the methods would have left n_alt_reads ambiguous
    on the isovar path while looking resolved. The subject is its own axis
    now, and ``count_in`` is how a caller states which bar it means.
    """
    fragment = MutantProteinFragment.from_isovar_result(_fake_isovar_result())

    assert fragment.count_in("n_alt_reads", FRAGMENTS) == 8
    assert fragment.count_in("n_alt_reads", READS) is None


def test_pvacseq_reference_reads_come_from_the_same_split_as_the_alt():
    """alt and ref are two halves of one estimate, not one and a remainder.

    The reference count was ``max(0, depth - alt)``, which clamps a
    negative to zero and so hides an inconsistent depth/VAF pair instead of
    surfacing it. topiary's split_reads_by_vaf returns both halves, and
    NA when either input is missing.
    """
    config = EpitopeConfig()
    report = read_pvacseq_report(
        os.path.join(DATA_DIR, "pvacseq_example.tsv"), epitope_config=config)
    result = pvacseq_ranking_result(
        report, epitopes_for_ranking(list(report.epitopes), config))

    fragment = result.ranked[0][1][0].mutant_protein_fragment

    assert fragment.n_alt_reads > 0
    assert fragment.n_alt_reads + fragment.n_ref_reads == (
        fragment.n_overlapping_reads)
    assert fragment.read_count_method == "rna_depth_x_vaf"


def test_pvacseq_claims_no_protein_supporting_count():
    """pVACseq does not report reads spanning the assembled sequence.

    That field used to be set to ``n_alt_reads``, asserting that every
    variant-supporting read also spans the peptide — a measurement this
    source never made, and one the LENS path fills with a genuinely
    different count (cds_overlap_reads).
    """
    config = EpitopeConfig()
    report = read_pvacseq_report(
        os.path.join(DATA_DIR, "pvacseq_example.tsv"), epitope_config=config)
    result = pvacseq_ranking_result(
        report, epitopes_for_ranking(list(report.epitopes), config))

    fragment = result.ranked[0][1][0].mutant_protein_fragment

    assert fragment.n_alt_reads_supporting_protein_sequence == 0


def test_delegating_the_split_did_not_move_any_number():
    """Migrating to topiary's arithmetic must be value-preserving.

    The local copy was ``int(round(depth * vaf))``. Python's ``round`` is
    banker's rounding, and nothing had checked that topiary agrees — so
    this pins that it does, at the ``.5`` boundaries where the two could
    differ. The reason to delegate is that a private copy of an upstream
    rule drifts, not that this one had drifted.
    """
    import pandas as pd
    from topiary import split_reads_by_vaf

    cases = [(1, 0.5), (3, 0.5), (5, 0.5), (7, 0.5),
             (2, 0.25), (10, 0.05), (100, 0.005)]
    depth = pd.Series([c[0] for c in cases], dtype="Float64")
    vaf = pd.Series([c[1] for c in cases], dtype="Float64")

    alt, _ref = split_reads_by_vaf(depth, vaf)

    assert [int(a) for a in alt] == [
        int(round(d * v)) for d, v in cases]


def test_an_absent_vaf_yields_no_split_rather_than_zero():
    """Both halves are NA when either input is missing.

    The local copy returned 0 for a missing VAF and derived the reference
    count as ``max(0, depth - 0)`` — the whole depth, as if every read
    supported the reference. topiary returns NA for both, which is what
    lets the caller record "no estimate" instead of a confident zero.
    """
    import pandas as pd
    from topiary import split_reads_by_vaf

    alt, ref = split_reads_by_vaf(
        pd.Series([100.0], dtype="Float64"),
        pd.Series([None], dtype="Float64"))

    assert pd.isna(alt.iloc[0])
    assert pd.isna(ref.iloc[0])


def test_the_external_paths_state_reads():
    """A depth x VAF estimate is in reads, and says so.

    Without this, a threshold written against a native run and copied to a
    LENS run silently changes what it means — the portability harm the
    subject axis exists to prevent. Within one run the unit is consistent
    and the ranking is unaffected either way.
    """
    config = EpitopeConfig()
    report = read_pvacseq_report(
        os.path.join(DATA_DIR, "pvacseq_example.tsv"), epitope_config=config)
    result = pvacseq_ranking_result(
        report, epitopes_for_ranking(list(report.epitopes), config))

    fragment = result.ranked[0][1][0].mutant_protein_fragment

    assert fragment.read_count_subject == READS
    assert fragment.count_in("n_alt_reads", READS) == fragment.n_alt_reads
    assert fragment.count_in("n_alt_reads", FRAGMENTS) is None


def test_a_fragment_that_states_no_subject_answers_either():
    """Silence is not a claim about the unit.

    The varcode path consulted no RNA, so its zeros have no subject. A
    caller asking for one gets the value rather than None, because there is
    no competing claim to contradict — refusing would treat "not stated" as
    "stated to be the other thing".
    """
    fragment = MutantProteinFragment(
        variant=object(), gene_name="BRAF", amino_acids="SIINFEKL",
        mutant_amino_acid_start_offset=0, mutant_amino_acid_end_offset=1,
        supporting_reference_transcripts=[],
        n_overlapping_reads=0, n_alt_reads=0, n_ref_reads=0,
        n_alt_reads_supporting_protein_sequence=0)

    assert fragment.read_count_subject == ""
    assert fragment.count_in("n_alt_reads", READS) == 0
    assert fragment.count_in("n_alt_reads", FRAGMENTS) == 0
