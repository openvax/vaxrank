"""Default reports show source biology without changing selection evidence."""

from dataclasses import replace
import json

import pytest

from vaxrank import antigens_from_isovar
from vaxrank.antigen_biology import antigen_biology, sequence_difference
from vaxrank.antigen_allocation import AntigenAllocation
from vaxrank.peptide import assemble_peptide_constructs, PeptideConstructConfig
from vaxrank.vaccine_peptide import VaccinePeptide

from .test_isovar_inputs import result
from .test_product_reports import template, render
from .test_assembly_antigens import candidate


def test_ncorf_somatic_declaration_reaches_default_report_with_exact_comparison():
    data = antigens_from_isovar(result('ncorf'), somatic=True)
    antigen, = data.antigens.values()
    assert antigen.tumor_specificity.admits_construct
    peptide = VaccinePeptide(antigen=antigen, epitopes=[])
    view = antigen_biology(antigen, source_evidence=data.source_evidence)
    comparison, = view['comparisons']
    assert comparison['tumor'] == 'MAAAAAVAAAA'
    assert comparison['reference'] == 'MAAAAAAAAAA'
    assert comparison['differences'] == '      ^    '
    assert dict(view['facts'])['Patient translation'] == 'Not measured by RNA reconstruction'
    for extension in ('html', 'txt'):
        report = render(template([(antigen, [peptide])], source_evidence=data.source_evidence), extension)
        assert 'Antigen biology and targeting' in report
        assert 'Focal allele reverted' in report and 'MAAAAAAAAAA' in report


def test_somatic_declaration_does_not_override_sequence_or_support_gates():
    synonymous = antigens_from_isovar(result('synonymous'), somatic=True)
    assert not any(a.tumor_specificity.admits_construct for a in synonymous.antigens.values())
    absent = antigens_from_isovar(result('ncorf'))
    assert all(d['reasons'] for source in absent.source_evidence.values() for d in source['decisions'])
    with pytest.raises(ValueError, match='only to mutations'):
        antigens_from_isovar(result('viral'), somatic=True)


def test_sequence_alignment_crops_by_source_coordinates_with_repeated_sequence_and_indel():
    alignment = sequence_difference('ACDACDVGH', 'ACDACDGH', 3, 8)
    assert alignment['tumor'] == 'ACDVG'
    assert alignment['reference'] == 'ACD-G'
    assert alignment['tumor_start'] == 4 and alignment['reference_start'] == 4


def test_sv_junction_offsets_follow_nested_slices_and_original_comparison():
    data = antigens_from_isovar(result('sv'), somatic=True,
        options={'admission': {'allow_ambiguous': True, 'min_fragments': 1}})
    antigen = next(a for a in data.antigens.values()
        if json.loads(dict(a.source_metadata)['isovar_candidate']).get('reference_comparisons'))
    crossing = antigen.targetable_mask.intervals[0].start
    window = antigen.sliced(crossing - 4, crossing + 6).sliced(1, 9)
    view = antigen_biology(window, source_evidence=data.source_evidence)
    assert window.source_protein_offset == crossing - 3
    assert view['junctions'][0]['boundaries'] == [(3, 4)]
    assert view['comparisons'][0]['tumor'].replace('-', '') == window.amino_acids


def test_ncorf_switch_is_independent_of_ordinary_mutations():
    ordinary, peptides = candidate(0)
    ncorf = replace(ordinary, source_identifier='ncorf', source_metadata=(('isovar_candidate',
        json.dumps({'coding_reference': {'kind': 'extra_orf', 'orf_type': 'uORF'}})),))
    ncorf_peptide = VaccinePeptide(antigen=ncorf, epitopes=peptides[0].epitopes)
    # Different windows avoid sequence deduplication attributing one product
    # to its first source, which is intentionally independent of the switch.
    ordinary2, peptides2 = candidate(1, 'GILGFVFTL')
    ranked = [(ordinary2, peptides2), (ncorf, [ncorf_peptide])]
    config = PeptideConstructConfig(antigen_allocation=AntigenAllocation(max_per_kind={'ncORF': 0}))
    products = assemble_peptide_constructs(ranked, config)
    assert len(products) == 1 and 'GILGFVFTL' in products[0].sequence
    config = PeptideConstructConfig(antigen_allocation=AntigenAllocation(max_per_kind={'mutation': 0}))
    products = assemble_peptide_constructs(ranked, config)
    assert len(products) == 1 and 'SIINFEKLA' in products[0].sequence


def test_cli_can_select_only_vanilla_mutations_or_add_ncorfs():
    from vaxrank.cli.arg_parser import parse_vaxrank_args
    from vaxrank.cli.vaccine_config_args import antigen_allocation_from_args
    args = parse_vaxrank_args(['--input-epitopes', 'saved.tsv', '--antigen-types', 'mutation', 'ncorf'])
    limits = antigen_allocation_from_args(args)['max_per_kind']
    assert limits['CTA'] == limits['SV'] == limits['fusion'] == 0
    assert 'mutation' not in limits and 'ncORF' not in limits
