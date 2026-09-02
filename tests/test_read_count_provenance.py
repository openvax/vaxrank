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

"""What RNA evidence fields hold, and how a consumer can tell.

The default combined score uses ``n_rna_alt``, which prefers independent
fragments and falls back to reads. Explicit read-based expressions can still
use ``n_alt_reads``. These tests pin both values and their provenance
(openvax/vaxrank#379, #382, #383, #394).
"""

import os
import types

from vaxrank.epitope_config import EpitopeConfig
from vaxrank.epitope_dsl import epitopes_for_ranking
from vaxrank.epitope_io import read_pvacseq_report
from vaxrank.external_input import pvacseq_ranking_result
from topiary import FRAGMENTS, READS

from vaxrank.mutant_protein_fragment import (
    RNA_EVIDENCE_METHOD_ALIGNMENT,
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
        num_supporting_fragments=7, num_supporting_reads=13, transcripts=[])
    return types.SimpleNamespace(
        variant=object(), top_protein_sequence=protein_sequence,
        num_total_fragments=20, num_alt_fragments=8, num_ref_fragments=12,
        num_total_reads=38, num_alt_reads=15, num_ref_reads=23)


def test_isovar_counts_are_carried_in_both_units():
    """Both units, each in the field named for it.

    vaxrank stored isovar's *fragment* counts in the fields named for
    reads, and then described the discrepancy as a subject a caller had to
    ask about. isovar reports every count in both units — the reads were
    never unavailable, and the accessor that explained their absence was
    describing a lie rather than a limitation (openvax/topiary#239).
    """
    fragment = MutantProteinFragment.from_isovar_result(_fake_isovar_result())

    assert fragment.n_alt_reads == 15           # isovar's num_alt_reads
    assert fragment.n_alt_fragments == 8        # isovar's num_alt_fragments
    assert fragment.rna_evidence_method == RNA_EVIDENCE_METHOD_ALIGNMENT
    assert fragment.rna_evidence_subject == FRAGMENTS
    assert fragment.sequence_source == SEQUENCE_SOURCE_ISOVAR_ASSEMBLY


def test_the_sequence_source_terms_are_read_from_topiary():
    """No private copy of a controlled vocabulary.

    A rename upstream should break the import rather than silently stamp
    fragments with a term nothing recognizes.
    """
    from topiary import SEQUENCE_SOURCES

    assert SEQUENCE_SOURCE_ISOVAR_ASSEMBLY in SEQUENCE_SOURCES
    assert SEQUENCE_SOURCE_VARCODE_TRANSLATION in SEQUENCE_SOURCES


def test_the_rna_accessors_prefer_fragments():
    """``n_rna_alt`` is the number to weight, and it prefers fragments.

    Two mates of one fragment are one piece of evidence, so a read count
    overstates paired-end support by roughly two — the wrong input to a
    diminishing-returns weight like sqrt(). Where a source has only reads,
    reads are what n_rna_alt reports, and the subject says which.
    """
    fragment = MutantProteinFragment.from_isovar_result(_fake_isovar_result())

    assert fragment.n_rna_alt == 8              # fragments win
    assert fragment.n_rna_ref == 12
    assert fragment.n_rna_overlapping == 20
    assert fragment.n_rna_supporting_protein_sequence == 7
    assert fragment.rna_evidence_subject == FRAGMENTS


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
    assert fragment.rna_evidence_method == "rna_depth_x_vaf"


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

    assert fragment.rna_evidence_subject == READS
    # No fragment counts, so the accessors report the reads.
    assert fragment.n_alt_fragments is None
    assert fragment.n_rna_alt == fragment.n_alt_reads


def test_a_fragment_that_states_no_subject_still_reports_its_counts():
    """Silence about the unit does not hide the numbers.

    The varcode path consulted no RNA, so its zeros have no subject. The
    accessors still answer, because there is no competing claim to
    contradict — treating "not stated" as "stated to be the other thing"
    would make an unlabelled zero unreadable.
    """
    fragment = MutantProteinFragment(
        variant=object(), gene_name="BRAF", amino_acids="SIINFEKL",
        mutant_amino_acid_start_offset=0, mutant_amino_acid_end_offset=1,
        supporting_reference_transcripts=[],
        n_overlapping_reads=0, n_alt_reads=0, n_ref_reads=0,
        n_alt_reads_supporting_protein_sequence=0)

    assert fragment.rna_evidence_subject == ""
    assert fragment.n_rna_alt == 0


def test_lens_vaf_is_not_recorded_as_a_dna_fraction(tmp_path):
    """LENS leaves `vaf` unqualified, so nothing may call it DNA.

    LENS names its read columns ``rna_*`` explicitly and leaves ``vaf``
    bare. vaxrank filed it as ``dna_vaf`` and the report printed it under
    "DNA VAF", asserting an assay the file never stated — and topiary was
    using the same column as both an RNA and a DNA fraction, which is how
    the ambiguity surfaced.

    It is kept under its own name now, and the split that uses it says
    ``rna_depth_x_source_vaf`` rather than claiming the fraction was RNA.
    """
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import read_lens_report
    from vaxrank.external_input import lens_ranking_result

    path = tmp_path / "vaf.tsv"
    path.write_text(
        "peptide\tallele\tpep_context\tantigen_source\tvariant_coords\t"
        "snv_ref_allele\tsnv_alt_allele\tgene_name\t"
        "rna_reads_covering_genomic_origin\tvaf\tmhcflurry_2.1.1.aff\n"
        "SIINFEKL\tHLA-A02:01\tAASIINFEKLAA\tSNV\tchr1:100\tC\t[T]\tG\t"
        "100\t0.4\t20\n")

    config = EpitopeConfig()
    report = read_lens_report(path, epitope_config=config)
    result = lens_ranking_result(
        report, epitopes_for_ranking(list(report.epitopes), config))

    assert result.dna_vaf_by_variant == {}
    assert set(result.source_vaf_by_variant.values()) == {0.4}

    fragment = result.ranked[0][1][0].mutant_protein_fragment
    assert fragment.rna_evidence_method == "rna_depth_x_source_vaf"


def test_pvacseq_keeps_its_qualified_dna_vaf():
    """pVACseq names the assay, so the DNA fraction stays a DNA fraction.

    The LENS fix must not collapse a source that *did* say which assay it
    measured into the same unqualified bucket.
    """
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import read_pvacseq_report
    from vaxrank.external_input import pvacseq_ranking_result

    config = EpitopeConfig()
    report = read_pvacseq_report(
        os.path.join(DATA_DIR, "pvacseq_example.tsv"), epitope_config=config)
    result = pvacseq_ranking_result(
        report, epitopes_for_ranking(list(report.epitopes), config))

    assert result.dna_vaf_by_variant
    assert result.source_vaf_by_variant == {}


def test_other_allele_reads_are_absent_when_the_reference_was_derived():
    """A subtraction of a subtraction is not an observation.

    Where a source reports a depth and a fraction, the reference count is
    ``depth - alt``, so ``depth - (ref + alt)`` is structurally zero — and
    the report printed that zero as a fact about the locus. Worse when the
    source reports no fraction: alt and ref are both 0 and the subtraction
    returns the whole depth, rendering a clean locus as one where every
    read supports a third allele.
    """
    config = EpitopeConfig()
    report = read_pvacseq_report(
        os.path.join(DATA_DIR, "pvacseq_example.tsv"), epitope_config=config)
    result = pvacseq_ranking_result(
        report, epitopes_for_ranking(list(report.epitopes), config))

    fragment = result.ranked[0][1][0].mutant_protein_fragment

    assert fragment.n_ref_reads > 0        # derived, not counted
    assert fragment.reference_count_was_measured is False
    assert fragment.n_other_reads is None


def test_other_allele_reads_survive_where_the_reference_was_counted():
    """isovar counts the reference independently, so the number is real."""
    fragment = MutantProteinFragment.from_isovar_result(_fake_isovar_result())

    assert fragment.reference_count_was_measured is True
    # 38 total reads, 15 alt, 23 ref — nothing left for a third allele here,
    # but the zero is an observation rather than an artifact of the method.
    assert fragment.n_other_reads == 0


def test_the_report_omits_the_row_when_there_is_nothing_to_report():
    """An absent measurement is not a zero on the page."""
    from vaxrank.report import TemplateDataCreator

    config = EpitopeConfig()
    report = read_pvacseq_report(
        os.path.join(DATA_DIR, "pvacseq_example.tsv"), epitope_config=config)
    result = pvacseq_ranking_result(
        report, epitopes_for_ranking(list(report.epitopes), config))
    vaccine_peptide = result.ranked[0][1][0]

    creator = TemplateDataCreator.__new__(TemplateDataCreator)
    creator.dna_vaf_by_variant = {}
    creator.source_vaf_by_variant = {}
    creator.gene_pathway_check = None
    variant_data = creator._variant_data(
        vaccine_peptide.mutant_protein_fragment.variant, vaccine_peptide)

    assert 'RNA reads supporting other alleles' not in variant_data
    assert 'RNA reads supporting variant allele' in variant_data


def test_what_counts_as_measured_is_topiarys_rule():
    """No local copy of which methods count.

    Only a direct count qualifies today, so naming ``rna_alignment`` here
    would agree — and would silently exclude any measured method topiary
    adds later, which is the drift this series has spent its length
    removing.
    """
    import dataclasses

    from topiary import READ_COUNT_METHODS, provenance_for_method

    base = MutantProteinFragment.from_isovar_result(_fake_isovar_result())
    for method in sorted(READ_COUNT_METHODS):
        fragment = dataclasses.replace(base, rna_evidence_method=method)
        assert fragment.reference_count_was_measured is (
            provenance_for_method(method) == "measured"), method


def _all_epitopes_fragment():
    import pathlib
    import tempfile

    from tests.test_epitope_io import _write_pvacseq_all_epitopes_fixture

    path = _write_pvacseq_all_epitopes_fixture(
        pathlib.Path(tempfile.mkdtemp()))
    config = EpitopeConfig()
    report = read_pvacseq_report(path, epitope_config=config)
    result = pvacseq_ranking_result(
        report, epitopes_for_ranking(list(report.epitopes), config))
    return result.ranked[0][1][0].mutant_protein_fragment


def test_dna_evidence_mirrors_rna_where_the_source_reports_a_depth():
    """A variant fraction alone cannot say how well a call is backed.

    0.45 from two reads and 0.45 from two hundred are not the same claim,
    and the report carried only the fraction. pVACseq's all_epitopes
    flavour reports a DNA depth, so the same depth x VAF split that
    produces the RNA counts produces these (openvax/topiary#248).
    """
    fragment = _all_epitopes_fragment()

    assert (fragment.n_dna_overlapping, fragment.n_dna_alt,
            fragment.n_dna_ref) == (1233, 372, 861)
    assert fragment.n_dna_alt + fragment.n_dna_ref == fragment.n_dna_overlapping
    assert fragment.dna_evidence_method == "dna_depth_x_vaf"
    assert fragment.dna_evidence_subject == "reads"


def test_a_fraction_without_a_depth_yields_no_dna_counts():
    """The aggregated flavour reports a DNA VAF and no DNA depth.

    Deriving counts from the fraction alone would invent the denominator,
    which is the substitution this whole series exists to prevent.
    """
    config = EpitopeConfig()
    report = read_pvacseq_report(
        os.path.join(DATA_DIR, "pvacseq_example.tsv"), epitope_config=config)
    result = pvacseq_ranking_result(
        report, epitopes_for_ranking(list(report.epitopes), config))

    fragment = result.ranked[0][1][0].mutant_protein_fragment

    assert fragment.dna_vaf is not None      # the fraction is reported
    assert fragment.n_dna_alt is None        # the depth is not
    assert fragment.n_dna_overlapping is None
    assert fragment.dna_evidence_method == ""


def test_lens_gets_no_dna_evidence_at_all():
    """LENS states no assay for its one fraction, so nothing may be DNA.

    Populating DNA counts from an unqualified fraction asserts the assay
    the file declined to name — the same reason its RNA split is labelled
    rna_depth_x_source_vaf rather than rna_depth_x_vaf.
    """
    from vaxrank.epitope_io import read_lens_report
    from vaxrank.external_input import lens_ranking_result

    config = EpitopeConfig()
    report = read_lens_report(
        os.path.join(DATA_DIR, "real_lens_subsets",
                     "lens_v1.4_real_subset.tsv"), epitope_config=config)
    result = lens_ranking_result(
        report, epitopes_for_ranking(list(report.epitopes), config))

    fragment = result.ranked[0][1][0].mutant_protein_fragment

    assert fragment.dna_vaf is None
    assert fragment.n_dna_alt is None
    assert fragment.dna_evidence_method == ""


def test_the_report_shows_dna_counts_only_where_they_exist():
    """Mirrors the RNA rows, and omits rather than zero-fills."""
    from vaxrank.report import TemplateDataCreator

    def variant_data_for(fragment, peptide):
        creator = TemplateDataCreator.__new__(TemplateDataCreator)
        creator.dna_vaf_by_variant = {}
        creator.source_vaf_by_variant = {}
        creator.gene_pathway_check = None
        return creator._variant_data(fragment.variant, peptide)

    import pathlib
    import tempfile

    from tests.test_epitope_io import _write_pvacseq_all_epitopes_fixture

    config = EpitopeConfig()
    path = _write_pvacseq_all_epitopes_fixture(
        pathlib.Path(tempfile.mkdtemp()))
    report = read_pvacseq_report(path, epitope_config=config)
    result = pvacseq_ranking_result(
        report, epitopes_for_ranking(list(report.epitopes), config))
    peptide = result.ranked[0][1][0]

    rows = variant_data_for(peptide.mutant_protein_fragment, peptide)
    assert rows['DNA reads supporting variant allele'] == 372
    assert rows['DNA read count method'] == 'dna_depth_x_vaf'

    aggregated_report = read_pvacseq_report(
        os.path.join(DATA_DIR, "pvacseq_example.tsv"), epitope_config=config)
    aggregated = pvacseq_ranking_result(
        aggregated_report,
        epitopes_for_ranking(list(aggregated_report.epitopes), config))
    aggregated_peptide = aggregated.ranked[0][1][0]

    rows = variant_data_for(
        aggregated_peptide.mutant_protein_fragment, aggregated_peptide)
    assert 'DNA reads supporting variant allele' not in rows


def test_the_reported_dna_vaf_comes_from_the_row_its_counts_came_from():
    """One observation, not a fraction and a depth from different rows.

    The report used to read the DNA VAF from a per-variant map while the
    counts came from a specific row. They agree today, and nothing enforced
    it — a variant group whose rows disagreed would print a fraction beside
    a depth that never produced it, which is the per-field maximum problem
    one field over.
    """
    from vaxrank.report import TemplateDataCreator

    fragment = _all_epitopes_fragment()
    creator = TemplateDataCreator.__new__(TemplateDataCreator)
    # A deliberately contradictory map: the fragment's own value must win.
    creator.dna_vaf_by_variant = {fragment.variant: 0.999}
    creator.source_vaf_by_variant = {}
    creator.gene_pathway_check = None

    class _Peptide:
        mutant_protein_fragment = fragment
        combined_score = 1.0

    rows = creator._variant_data(fragment.variant, _Peptide())

    assert rows['DNA VAF'] == '%.3f' % fragment.dna_vaf
    assert rows['DNA reads covering the position'] == fragment.n_dna_overlapping
