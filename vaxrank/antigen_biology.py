"""Compact, source-grounded views of antigen biology for every report format."""

from Bio.Align import PairwiseAligner

from .coding_context import antigen_source_label, isovar_candidates
from .report_display import marked_sequence


def sequence_difference(tumor, reference, start=0, end=None):
    """Display a sequence alignment; no inference of ancestry or read phase."""
    if tumor and reference:
        aligner = PairwiseAligner(mode='global', match_score=2, mismatch_score=-1,
                                 open_gap_score=-3, extend_gap_score=-0.5)
        alignment = aligner.align(tumor, reference)[0]
        left, right = alignment[0], alignment[1]
    else:
        width = max(len(tumor), len(reference))
        left, right = tumor.ljust(width, '-'), reference.ljust(width, '-')
    marks = ''.join(' ' if a == b else '^' for a, b in zip(left, right))
    end = len(tumor) if end is None else end
    # Map source coordinates through alignment gaps, never by substring search.
    boundaries = [0]
    for index, residue in enumerate(left):
        if residue != '-':
            boundaries.append(index + 1)
    low, high = boundaries[start], boundaries[end]
    if end == len(tumor):
        high = len(left)
    return dict(tumor=left[low:high], reference=right[low:high], differences=marks[low:high],
                tumor_start=start + 1, tumor_end=end,
                reference_start=len(right[:low].replace('-', '')) + 1,
                reference_end=len(right[:high].replace('-', '')))


def antigen_biology(antigen, *, source_evidence=None, peptides=()):
    """Summarize recorded source facts, keeping RNA and translation separate.

    Comparison sequences come only from producer exports or retained peptide
    comparators. Source offsets, rather than sequence search, locate a window
    even when the original protein contains repeated sequences.
    """
    from .antigen_origin import biology_origin
    candidates = isovar_candidates(antigen)
    facts, comparisons, junctions = [], [], []
    identities = [value for key, value in antigen.source_metadata if key == 'isovar_evidence_id']
    evidence = (source_evidence or {}).get(identities[0], {}) if len(identities) == 1 else {}
    result = evidence.get('isovar_result', {})
    origin = biology_origin(antigen, result=result, candidates=candidates, peptides=peptides)
    if origin is not None:
        for variant in origin.variants:
            facts.append(('Variant type', variant.variant_type + (':' + variant.sv_type if variant.sv_type else '')))
            if variant.contig and variant.start is not None:
                facts.append(('Variant', '%s:%s %s → %s' % (
                    variant.contig, variant.start, variant.ref or '∅', variant.alt or '∅')))
            if variant.end is not None and (variant.sv_type != 'BND' or variant.end != variant.start):
                facts.append(('End coordinate', str(variant.end)))
            if variant.mate_contig:
                mate = '%s:%s' % (variant.mate_contig, variant.mate_start if variant.mate_start is not None else '?')
                facts.append(('Breakend mate', mate + (' ' + variant.mate_orientation if variant.mate_orientation else '')))
        for junction in origin.junctions:
            boundaries = [(a - antigen.source_protein_offset, b - antigen.source_protein_offset)
                for a, b in junction.protein_intervals
                if antigen.source_protein_offset <= a <= b <= antigen.source_protein_offset + len(antigen.amino_acids)]
            if boundaries or junction.left.label() != 'Not recorded' or junction.right.label() != 'Not recorded':
                junctions.append(dict(left=junction.left.label(), right=junction.right.label(), boundaries=boundaries))
            if junction.description:
                facts.append(('RNA junction', junction.description))
        if origin.rna_product:
            facts.append(('RNA product', origin.rna_product))
            if not origin.variants:
                facts.append(('Causal DNA variant', 'Not recorded'))
    for candidate in candidates:
        original = candidate.get('protein_sequence')
        if original is None and len(candidates) == 1:
            original = next((row.get('protein_sequence') for row in result.get('protein_candidates', [])
                             if row.get('candidate_id') == antigen.source_identifier), None)
        low = antigen.source_protein_offset
        high = low + len(antigen.amino_acids)
        if original is not None and original[low:high] != antigen.amino_acids:
            original = None  # Merged origins have no uniquely recorded window mapping.
        support = candidate.get('occurrence_rna_support', candidate.get('rna_support')) or {}
        if support.get('fragments') is not None:
            facts.append(('Patient RNA', '%s fragments%s' % (support['fragments'],
                ' (' + candidate['support_scope'].replace('_', ' ') + ')' if candidate.get('support_scope') else '')))
        recorded = candidate.get('translation_evidence') or {}
        reference_model = candidate.get('coding_reference') or {}
        if reference_model.get('orf_id'):
            facts.append(('Reference ORF', ':'.join(str(reference_model[key])
                for key in ('source', 'orf_id') if reference_model.get(key))))
        ribo = recorded.get('ribosome_studies', [])
        ms = recorded.get('peptides', [])
        if recorded:
            facts.append(('Reference translation evidence',
                '%d ribosome-profiling records; %d peptide records' % (len(ribo), len(ms))))
        if 'translation_observed' in candidate:
            observed = candidate['translation_observed']
            facts.append(('Patient translation', 'Observed in source evidence' if observed is True else
                          'Not measured by RNA reconstruction' if observed is False else 'Not recorded'))
        if candidate.get('ambiguous_coding_candidate'):
            facts.append(('Coding interpretation', 'Competing or unresolved reading frame'))
        if candidate.get('partial_coding_candidate'):
            facts.append(('Sequence extent', 'Local coding window; full protein not reconstructed'))
        reference = candidate.get('reference_allele_counterfactual_protein')
        label = 'Focal allele reverted on the same RNA background'
        if reference is None:
            reference = candidate.get('reference_protein_window')
            label = 'Annotation reference'
        if reference is not None and original is not None:
            comparisons.append(dict(label=label, **sequence_difference(original, reference,
                antigen.source_protein_offset, antigen.source_protein_offset + len(antigen.amino_acids))))
        for comparison in candidate.get('reference_comparisons', []):
            reference = comparison.get('amino_acids')
            if reference and original is not None:
                comparisons.append(dict(label='Reference ORF: ' + comparison.get('transcript_id', ''),
                    **sequence_difference(original, reference, antigen.source_protein_offset,
                                          antigen.source_protein_offset + len(antigen.amino_acids))))
    if not comparisons:
        for peptide in peptides:
            for epitope in peptide.epitopes:
                comparator = epitope.reference_comparator
                if comparator is not None and comparator.sequence:
                    comparison = dict(label='Patient background' if epitope.comparison_evidence else 'Reference peptide',
                                      **sequence_difference(epitope.sequence, comparator.sequence))
                    if comparison not in comparisons:
                        comparisons.append(comparison)
    targeting = []
    from .cta_targeting import targeting_evidence, matching_gene_summary
    for peptide in peptides:
        for epitope in peptide.epitopes:
            target = targeting_evidence(epitope)
            if target is None:
                continue
            sources = sorted({(source.get('gene_name') or source.get('gene_id') or match['record_id']) +
                              ' / ' + (source.get('protein_id') or match['record_id'])
                              for match in target['matches'] if match['role'] in ('exclude', 'review')
                              for source in match['occurrences']})
            row = dict(peptide=epitope.sequence,
                       start=peptide.antigen.source_protein_offset + epitope.offset,
                       end=peptide.antigen.source_protein_offset + epitope.offset + len(epitope.sequence),
                       status={'eligible': 'Eligible under OncoRef sequence policy',
                               'excluded': 'Do not target: contained in Unsafe protein',
                               'review': 'Review: unresolved reference occurrence',
                               'unassessed': 'Unassessed: no qualified reference occurrence'}[target['status']],
                       sources='; '.join(sources) or 'No Unsafe/review match in this bundle')
            row['source_summary'] = matching_gene_summary(target['matches'],
                roles=('exclude',) if target['status'] == 'excluded' else ('review',)) or 'No Unsafe/review match'
            row['status_code'] = target['status']
            if row not in targeting:
                targeting.append(row)
    if antigen.kind == 'CTA' and not targeting:
        facts.append(('Unsafe containment', 'Not assessed with the current OncoRef sequence bundle'))
    if antigen.kind == 'CTA':
        meta = dict(antigen.source_metadata)
        ts_cat = meta.get('oncoref_target_selection_category')
        ts_policy = meta.get('oncoref_target_selection_policy_id')
        if ts_cat:
            label = ts_cat + (f' ({ts_policy})' if ts_policy else '')
            facts.append(('OncoRef target class', label))
        elif ts_policy:
            facts.append(('OncoRef target policy', ts_policy))
    return dict(kind=antigen_source_label(antigen, source_evidence),
                gene=antigen.display_gene_name, identifier=antigen.display_identifier,
                facts=list(dict.fromkeys(facts)), junctions=junctions, comparisons=comparisons,
                sequence=antigen.amino_acids,
                sequence_html=marked_sequence(antigen.amino_acids,
                    [] if antigen.kind == 'CTA' else [(i.start, i.end) for i in antigen.targetable_mask.intervals]),
                sequence_start=antigen.source_protein_offset + 1,
                sequence_end=antigen.source_protein_offset + len(antigen.amino_acids),
                targeting=targeting,
                targeting_preview=sorted(targeting, key=lambda row: row['status_code'] == 'eligible')[:8],
                admission=antigen.tumor_specificity.status)
