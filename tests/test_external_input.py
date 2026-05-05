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
loaders produce, then a shared dispatch (``_emit_outputs``) into
modality-specific construct writers.
"""

import glob
import json
import os

import pytest

from vaxrank.epitope_io import load_lens
from vaxrank.external_input import (
    _mut_offsets_in_context,
    _parse_variant_coords,
    ranked_from_lens_predictions,
)

DATA_DIR = os.path.join(
    os.path.dirname(__file__), "data", "epitope_fixtures")


def test_parse_variant_coords_extracts_chr_pos():
    """``_parse_variant_coords`` extracts only (contig, pos), stripping
    the ``chr`` prefix LENS emits — pyensembl uses bare contigs, and
    keeping the prefix breaks ``variant.effect_on_transcript`` and
    renders the variant short_description as ``chrchr3 …``. Ref/alt
    come from the dedicated per-antigen-source columns
    (snv_ref_allele / snv_alt_allele etc.), looked up by
    ``_variant_from_lens_row``."""
    assert _parse_variant_coords("chr17:7675088:C:T") == ("17", 7675088)
    assert _parse_variant_coords("chr1:26780312") == ("1", 26780312)
    assert _parse_variant_coords("chr1:26780312:C") == ("1", 26780312)
    # Already-bare contig pass-through:
    assert _parse_variant_coords("17:7675088") == ("17", 7675088)


def test_parse_variant_coords_returns_none_on_missing_or_malformed():
    for bad in ("garbage", "chr1:notapos:A:T", None, "", "   ",
                "nan", "NaN"):
        assert _parse_variant_coords(bad) is None


def test_strip_lens_allele_handles_bracket_form():
    """LENS records alt alleles as bracketed strings: '[T]', '[CA]'.
    The helper strips brackets and returns the inner sequence."""
    from vaxrank.external_input import _strip_lens_allele
    assert _strip_lens_allele("[T]") == "T"
    assert _strip_lens_allele("[CA]") == "CA"
    assert _strip_lens_allele("T") == "T"   # already unbracketed
    assert _strip_lens_allele("") is None
    assert _strip_lens_allele(None) is None
    assert _strip_lens_allele("nan") is None


# ---- _resolve_transcripts ------------------------------------------------

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


def test_resolve_transcripts_strips_version_suffix():
    """LENS IDs carry version suffixes (``ENST00000312960.4``); pyensembl
    2.x doesn't strip them. The helper must, or every LENS lookup fails."""
    from vaxrank.external_input import _resolve_transcripts
    sentinel = object()
    g = _StubGenome({'ENST00000312960': sentinel})
    out = _resolve_transcripts(['ENST00000312960.4'], g)
    assert out == [sentinel]
    assert g.lookups == ['ENST00000312960']  # bare form was sent


def test_resolve_transcripts_drops_unresolvable_ids():
    """A release-mismatch yields some IDs the configured release
    doesn't know. Drop quietly (DEBUG-logged) rather than crash; the
    caller's aggregate WARN summarizes the count."""
    from vaxrank.external_input import _resolve_transcripts
    known = object()
    g = _StubGenome({'ENST00000000001': known})
    out = _resolve_transcripts(
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


def test_resolve_transcripts_returns_empty_when_no_genome():
    """No --ensembl-release set → genome=None; resolution is a no-op
    and downstream code falls back to the empty-transcript path."""
    from vaxrank.external_input import _resolve_transcripts
    assert _resolve_transcripts(['ENST00000000001'], None) == []
    assert _resolve_transcripts([], None) == []


def test_variant_from_lens_row_uses_real_snv_alleles():
    """SNV rows: ref/alt come from snv_ref_allele / snv_alt_allele.
    No placeholder genotype. The synthesized Variant carries the real
    biology so it can be fed to varcode-effect annotation safely."""
    from vaxrank.external_input import _variant_from_lens_row
    row = {
        'variant_coords': 'chr1:1624824',
        'antigen_source': 'SNV',
        'snv_ref_allele': 'C',
        'snv_alt_allele': '[T]',  # LENS's bracketed format
    }
    v = _variant_from_lens_row(row)
    assert v is not None
    # ``chr`` prefix stripped during parse so pyensembl-style lookups
    # work; see ``_parse_variant_coords``.
    assert v.contig == "1"
    assert v.start == 1624824
    assert v.ref == "C"
    assert v.alt == "T"


def test_variant_from_lens_row_uses_real_indel_alleles():
    """INDEL rows: ref/alt come from indel_ref_allele / indel_alt_allele.
    varcode normalizes indels (strips shared prefix), so ``CA → C``
    becomes ref='A' alt='' at start+1 — that's the canonical
    representation a downstream caller would expect."""
    from vaxrank.external_input import _variant_from_lens_row
    row = {
        'variant_coords': 'chr3:150742445',
        'antigen_source': 'INDEL',
        'indel_ref_allele': 'CA',
        'indel_alt_allele': '[C]',
    }
    v = _variant_from_lens_row(row)
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
    from vaxrank.external_input import _variant_from_lens_row
    for src in ('SPLICE', 'FUSION', 'CTA/SELF', 'ERV'):
        row = {
            'variant_coords': None,
            'antigen_source': src,
        }
        assert _variant_from_lens_row(row) is None


def test_variant_from_lens_row_skips_when_alleles_missing():
    """If the SNV row's allele columns are missing/NaN, return None
    rather than fabricate a placeholder."""
    from vaxrank.external_input import _variant_from_lens_row
    row = {
        'variant_coords': 'chr1:1000',
        'antigen_source': 'SNV',
        'snv_ref_allele': None,
        'snv_alt_allele': None,
    }
    assert _variant_from_lens_row(row) is None


def test_real_lens_v19_subset_produces_ranked_entries():
    """Regression test for the real-LENS-v1.9 path. The subset fixture
    came from a Hugo IPRES Pt10 dump and tripped both #259 (variant_coords
    parser) and #260 (duplicate-row handling). Pin that the path now
    produces a non-empty ranked list end-to-end."""
    path = os.path.join(
        DATA_DIR, "real_lens_subsets", "lens_v1.9_real_subset.tsv")
    report_df, predictions = load_lens(path)
    ranked, _dna_vaf = ranked_from_lens_predictions(predictions, path)
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
    report_df, predictions = load_lens(path)
    ranked, _dna_vaf = ranked_from_lens_predictions(predictions, path)
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
    from vaxrank.vaccine_library import has_only_standard_amino_acids
    assert has_only_standard_amino_acids("KLQGHSAPVL")
    assert has_only_standard_amino_acids("")
    for non_standard in ("U", "O", "X", "B", "Z", "J", "*"):
        assert not has_only_standard_amino_acids(
            "KLQ" + non_standard + "VL"), \
            f"Should reject {non_standard!r}"


def test_mut_offsets_in_context_finds_peptide():
    # AASVVGSSSSSGTR contains SVVGSSSSS at offset 2, length 9
    start, end = _mut_offsets_in_context("SVVGSSSSS", "AASVVGSSSSSGTR")
    assert start == 2
    assert end == 11


def test_mut_offsets_in_context_returns_none_when_not_found():
    """When the peptide isn't a substring of pep_context, return
    ``(None, None)`` so the caller drops the row instead of falsely
    claiming the entire context is the mutation span."""
    start, end = _mut_offsets_in_context("XYZXYZ", "AASVVGSSSSSGTR")
    assert start is None
    assert end is None
    # Empty inputs also return None
    assert _mut_offsets_in_context("", "AAVK") == (None, None)
    assert _mut_offsets_in_context("AAVK", "") == (None, None)


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
    report_df, predictions = load_lens(path)
    ranked, _dna_vaf = ranked_from_lens_predictions(predictions, path)
    assert len(ranked) == 1, (
        "Fixture has one variant with 3 peptides; got %d entries" % len(ranked))
    _, peptides = ranked[0]
    fragment = peptides[0].mutant_protein_fragment
    # Strongest binder is STRNGLVLLL (mhcflurry IC50 = 18.0). Its
    # pep_context is "AASTRNGLVLLLGTR".
    assert "STRNGLVLLL" in fragment.amino_acids, (
        "Expected pep_context to come from the strongest binder row "
        "(STRNGLVLLL, IC50=18); got %r" % fragment.amino_acids)


# ---- _patient_info_from_external ----------------------------------------

def test_patient_info_from_external_proxy_counts():
    """The synthesized PatientInfo's variant counts come from the
    ranked output, not from VCF/BAM. Pin each count's definition so a
    later loader change can't silently shift the rendered template
    report headers."""
    from vaxrank.external_input import _patient_info_from_external

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    report_df, predictions = load_lens(path)
    ranked, _dna_vaf = ranked_from_lens_predictions(predictions, path)
    info = _patient_info_from_external(
        ranked, path, patient_id='Pt-X', input_label='LENS report')
    assert info.patient_id == 'Pt-X'
    # PatientInfo on the external path no longer overloads
    # ``vcf_paths`` (LENS files aren't VCFs); the explicit
    # ``inputs`` list carries the path with an accurate label.
    assert info.vcf_paths == []
    assert info.inputs == [('LENS report', path)]
    assert info.bam_path is None
    # All ranked variants are counted as somatic (LENS antigen-only
    # files don't carry pre-pipeline silent variants).
    assert info.num_somatic_variants == len(ranked)
    # Without --ensembl-release the test fixture won't resolve
    # transcripts; coding-effect count drops to 0.
    assert info.num_coding_effect_variants == 0
    # ``num_variants_with_vaccine_peptides`` matches the populated
    # ranked entries (any with vps).
    assert info.num_variants_with_vaccine_peptides == sum(
        1 for _, vps in ranked if vps)


def test_patient_info_from_external_empty_ranked():
    """No ranked variants → all counts zero, no crash."""
    from vaxrank.external_input import _patient_info_from_external
    info = _patient_info_from_external([], '/tmp/empty.tsv', patient_id='')
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
    """Loading the LENS fixture and converting to ranked vaccine peptides
    should produce one entry per unique variant, each carrying the
    pep_context as the antigen fragment."""
    path = os.path.join(DATA_DIR, "lens_example.tsv")
    report_df, predictions = load_lens(path)
    ranked, _dna_vaf = ranked_from_lens_predictions(predictions, path)
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
                > len(vp.epitope_predictions[0].peptide_sequence))
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
    report_df, predictions = load_lens(path)
    ranked, _dna_vaf = ranked_from_lens_predictions(predictions, path)
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
    report_df, predictions = load_lens(path)
    ranked, _dna_vaf = ranked_from_lens_predictions(predictions, path)

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
    from vaxrank.epitope_io import load_pvacseq
    from vaxrank.external_input import ranked_from_pvacseq_predictions

    path = os.path.join(DATA_DIR, "pvacseq_example.tsv")
    report_df, predictions = load_pvacseq(path)
    ranked, _dna_vaf = ranked_from_pvacseq_predictions(predictions, path)
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
    ``(Variant, alleles_real)`` matching the LENS-path shape;
    pVACseq IDs always carry real ref/alt so alleles_real is True.

    The parser strips ``chr`` prefixes so pyensembl-style downstream
    lookups (``variant.effect_on_transcript``) work — same rationale
    as ``_parse_variant_coords`` on the LENS path."""
    from vaxrank.external_input import _parse_pvacseq_id

    for vid, expected in [
        ("chr1-100000-100001-A-T", ("1", 100000, "A", "T")),
        ("chr1-100000-A-T", ("1", 100000, "A", "T")),
        ("1.123.A.T", ("1", 123, "A", "T")),
    ]:
        v, alleles_real = _parse_pvacseq_id(vid)
        assert v is not None
        assert (v.contig, v.start, v.ref, v.alt) == expected
        assert alleles_real is True, (
            "pVACseq IDs carry real ref/alt; alleles_real must be True")
    # Garbage / wrong shape returns (None, False) instead of raising
    for bad in ("garbage", "chr1-notapos-A-T", None, ""):
        v, alleles_real = _parse_pvacseq_id(bad)
        assert v is None
        assert alleles_real is False


def test_pvacseq_drives_mrna_construct_assembly_end_to_end():
    """End-to-end: pVACseq report → ranked → mRNA constructs (mirror
    of the LENS end-to-end test). Pre-fix this produced zero
    constructs because the ranked list was empty."""
    from vaxrank.epitope_io import load_pvacseq
    from vaxrank.external_input import ranked_from_pvacseq_predictions
    from vaxrank.mrna import RNAConstructConfig, assemble_mrna_constructs

    path = os.path.join(DATA_DIR, "pvacseq_example.tsv")
    _, predictions = load_pvacseq(path)
    ranked, _dna_vaf = ranked_from_pvacseq_predictions(predictions, path)
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


# ---- LENS scoring-column edge cases --------------------------------------

def test_lens_warns_when_no_affinity_columns_detected(tmp_path, caplog):
    """When no pMHC_affinity predictor is detected, the
    representative-peptide pick is degenerate (every row ties at
    (2, 0.0); file-order first row "wins"). Pin that we warn
    instead of silently picking arbitrarily.

    Today this can be triggered by a stability-only LENS file
    (netmhcstabpan present without mhcflurry/netmhcpan). It can also
    fire if LENS ever adds tools whose ``kind`` is something other
    than ``pMHC_affinity`` (presentation-only, cleavage,
    immunogenicity, etc.) and the registry classifies them
    accordingly. The warning is kind-agnostic — it reports the
    detected kinds in its message.
    """
    import logging
    from vaxrank.epitope_io import load_lens
    from vaxrank.external_input import ranked_from_lens_predictions

    # Stability-only fixture (the only "no affinity" case the
    # current registry can produce). If new predictor kinds are added
    # to _LENS_PREDICTOR_REGISTRY, add fixtures here covering them.
    no_aff = tmp_path / "no_affinity.tsv"
    no_aff.write_text(
        "allele\tpeptide\tnetmhcstabpan_1.0.halflife_hours\t"
        "netmhcstabpan_1.0.perc_rank_stab\tantigen_source\tmut_aa_pos\t"
        "variant_coords\tsnv_ref_allele\tsnv_alt_allele\t"
        "gene_name\ttpm\tpep_context\n"
        "HLA-A02:01\tSVVGSSSSS\t12.5\t0.5\tSNV\t245\t"
        "chr17:7675088\tC\t[T]\tTP53\t42.5\tAASVVGSSSSSGTR\n")

    _, predictions = load_lens(str(no_aff))
    with caplog.at_level(logging.WARNING):
        ranked, _dna_vaf = ranked_from_lens_predictions(predictions, str(no_aff))
    # Still produces output (the file-order representative)
    assert len(ranked) == 1
    # Warning was emitted, mentions the detected non-affinity kinds
    msg = next(
        (r.message for r in caplog.records
         if "no pMHC_affinity scoring columns" in r.message), None)
    assert msg is not None, \
        "Expected a warn-once for LENS input lacking affinity predictors"
    assert "pMHC_stability" in msg, \
        "Warning should report the detected non-affinity kinds; got %r" % msg


# ---- _emit_outputs observability -----------------------------------------

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
    from vaxrank.cli.entry_point import _emit_outputs

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    report_df, predictions = load_lens(path)
    ranked, _dna_vaf = ranked_from_lens_predictions(predictions, path)

    out_dir = str(tmp_path / "out")
    args = _mrna_args(out_dir)
    with caplog.at_level(logging.INFO):
        _emit_outputs(args, ranked, source='external')
    dispatch = [r.message for r in caplog.records
                if 'Vaccine-type dispatch' in r.message]
    assert len(dispatch) == 1
    line = dispatch[0]
    assert "[external]" in line
    assert "['mrna']" in line  # active types and fired list both ['mrna']


def test_resolve_vaccine_types_default_and_unknown():
    """``--vaccine-type`` is multi-valued (default
    ``['peptide', 'mrna']`` — both modalities). Unknown values
    raise so future types stub-declared in user scripts fail loudly
    instead of silently no-opping."""
    import pytest as _pytest
    from types import SimpleNamespace
    from vaxrank.cli.entry_point import _resolve_vaccine_types

    # Defaults: both modalities active when --vaccine-type isn't set
    assert _resolve_vaccine_types(SimpleNamespace()) == ['peptide', 'mrna']
    assert _resolve_vaccine_types(
        SimpleNamespace(vaccine_type=None)) == ['peptide', 'mrna']
    # Single-string form (caller bypassing argparse)
    assert _resolve_vaccine_types(
        SimpleNamespace(vaccine_type='mrna')) == ['mrna']
    # Single-element list narrows to that one type
    assert _resolve_vaccine_types(
        SimpleNamespace(vaccine_type=['peptide'])) == ['peptide']
    # Multi-mode + dedup
    assert _resolve_vaccine_types(
        SimpleNamespace(vaccine_type=['peptide', 'mrna'])
    ) == ['peptide', 'mrna']
    assert _resolve_vaccine_types(
        SimpleNamespace(vaccine_type=['mrna', 'peptide', 'mrna'])
    ) == ['mrna', 'peptide']  # first-occurrence order preserved
    # Unknown
    with _pytest.raises(ValueError, match="Unknown vaccine type"):
        _resolve_vaccine_types(SimpleNamespace(vaccine_type=['dna']))


def test_vaccine_target_dir_layout(tmp_path):
    """Single-mode → files land directly in --output-dir; multi-mode
    → per-modality subdirs (so peptide's manifest.json doesn't
    overwrite mrna's)."""
    from vaxrank.cli.entry_point import _vaccine_target_dir
    base = str(tmp_path / "out")
    # Single-mode flat
    assert _vaccine_target_dir(base, 'peptide', ['peptide']) == base
    # Multi-mode subdirs
    assert _vaccine_target_dir(base, 'peptide', ['peptide', 'mrna']) == \
        os.path.join(base, 'peptide')
    assert _vaccine_target_dir(base, 'mrna', ['peptide', 'mrna']) == \
        os.path.join(base, 'mrna')
    # No --output-dir → no writer fires
    assert _vaccine_target_dir('', 'peptide', ['peptide']) is None


def test_lens_with_vaccine_type_mrna_writes_three_fastas(tmp_path):
    """End-to-end: --input-lens + --vaccine-type=mrna +
    --output-dir=DIR writes three FASTAs + manifest + layers.csv
    directly into DIR (single-mode flat layout)."""
    from vaxrank.cli.entry_point import _emit_outputs

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    report_df, predictions = load_lens(path)
    ranked, _dna_vaf = ranked_from_lens_predictions(predictions, path)

    out_dir = str(tmp_path / "out")
    args = _mrna_args(out_dir)
    _emit_outputs(args, ranked, source='external')
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
    from vaxrank.cli.entry_point import _emit_outputs

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    report_df, predictions = load_lens(path)
    ranked, _dna_vaf = ranked_from_lens_predictions(predictions, path)

    out_dir = str(tmp_path / "out")
    args = _vaccine_args(out_dir, vaccine_type=['peptide', 'mrna'])
    _emit_outputs(args, ranked, source='external')
    # mRNA outputs in DIR/mrna/
    assert os.path.isfile(os.path.join(out_dir, "mrna", "cds.fasta"))
    assert os.path.isfile(os.path.join(out_dir, "mrna", "manifest.json"))
    # Peptide outputs in DIR/peptide/ — canonical filenames don't
    # collide with mrna's manifest.json
    assert os.path.isfile(os.path.join(out_dir, "peptide", "vaccine.fasta"))
    assert os.path.isfile(os.path.join(out_dir, "peptide", "manifest.json"))


def test_emit_outputs_skips_writer_when_no_output_dir(tmp_path):
    """No --output-dir → writer no-ops (ranking-only run)."""
    from vaxrank.cli.entry_point import _emit_outputs

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    report_df, predictions = load_lens(path)
    ranked, _dna_vaf = ranked_from_lens_predictions(predictions, path)

    out_dir = str(tmp_path / "out")
    args = _mrna_args(out_dir)
    args.output_dir = ''  # no destination → no writer fires
    _emit_outputs(args, ranked, source='external')
    assert not os.path.exists(out_dir), \
        "Without --output-dir, no writer should have fired"
