"""Compact, source-grounded views of antigen biology for every report format."""

import hashlib

from Bio.Align import PairwiseAligner

from .coding_context import antigen_source_label, isovar_candidates
from .report_display import (
    _clip_intervals,
    _merge_intervals,
    _subtract_intervals,
    effective_target_intervals,
    marked_sequence,
    protein_context_key,
)


def resolve_cta_protein_contexts(ranked_variants_with_vaccine_peptides, source_evidence=None, *,
                                 include_keys=None):
    """Resolve SHA-256-verified full CTA proteins and exact 8-mer+ Unsafe overlaps.

    Runs one batched OncoRef containment pass only when CTAs carry a matching
    ``oncoref_sequence_bundle`` record and ``oncoref`` is installed. Fails open
    to an empty dict when offline or when no matching bundle is present.
    """
    cta_items = []
    bundle_ids = set()
    policies = set()
    for _, peptides in ranked_variants_with_vaccine_peptides:
        if not peptides:
            continue
        for peptide in peptides:
            antigen = peptide.antigen
            if antigen.kind != 'CTA':
                continue
            key = protein_context_key(antigen)
            if include_keys is not None and key not in include_keys:
                continue
            meta = dict(antigen.source_metadata)
            bundle_id = meta.get('oncoref_sequence_bundle')
            if not bundle_id:
                continue
            bundle_ids.add(bundle_id)
            if meta.get('oncoref_target_policy'):
                policies.add(meta['oncoref_target_policy'])
            cta_items.append((key, antigen))
    if not cta_items or len(bundle_ids) != 1:
        return {}
    recorded_bundle_id = next(iter(bundle_ids))
    if not policies and source_evidence:
        for rec in source_evidence.values():
            if (isinstance(rec, dict)
                    and rec.get('schema') == 'vaxrank.oncoref_target_reference.v1'
                    and rec.get('target_policy')):
                policies.add(rec['target_policy'])
    target_policy = next(iter(policies)) if len(policies) == 1 else None

    try:
        from .cta_targeting import load_bundle, load_target_selection, peptide_containment
        bundle = load_bundle('current')
    except Exception:
        return {}
    if getattr(bundle, 'identity', None) != recorded_bundle_id:
        return {}

    wanted_pids = {pid for _, antigen in cta_items for pid in antigen.protein_ids}
    wanted_shas = {antigen.source_protein_sha256 for _, antigen in cta_items if antigen.source_protein_sha256}
    seq_by_pid = {}
    seq_by_sha = {}
    for rec in bundle.records():
        rid = rec.get('record_id')
        seq = rec.get('sequence') or ''
        sha = rec.get('sequence_sha256') or hashlib.sha256(seq.encode()).hexdigest()
        if rid in wanted_pids:
            seq_by_pid[rid] = (seq, sha)
        if sha in wanted_shas and sha not in seq_by_sha:
            seq_by_sha[sha] = seq

    resolved = {}
    kmers = set()
    for key, antigen in cta_items:
        if key in resolved:
            continue
        full_seq = None
        for pid in antigen.protein_ids:
            cand = seq_by_pid.get(pid)
            if cand is None:
                continue
            seq, sha = cand
            if antigen.source_protein_sha256 and sha != antigen.source_protein_sha256:
                continue
            lo = antigen.source_protein_offset
            hi = lo + len(antigen.amino_acids)
            if seq[lo:hi] == antigen.amino_acids:
                full_seq = seq
                break
        if full_seq is None and antigen.source_protein_sha256:
            seq = seq_by_sha.get(antigen.source_protein_sha256)
            lo = antigen.source_protein_offset
            hi = lo + len(antigen.amino_acids)
            if seq is not None and seq[lo:hi] == antigen.amino_acids:
                full_seq = seq
        resolved[key] = (antigen, full_seq)
        scan_seq = full_seq if full_seq is not None else antigen.amino_acids
        for i in range(max(0, len(scan_seq) - 7)):
            kmers.add(scan_seq[i:i + 8])

    if not kmers:
        return {}
    try:
        ts = load_target_selection() if target_policy in ('target_selection', 'target_selection_cancer_antigens') else None
        hits = peptide_containment(kmers, bundle, target_policy=target_policy, target_selection_index=ts)
    except Exception:
        return {}

    contexts = {}
    for key, (antigen, full_seq) in resolved.items():
        scan_seq = full_seq if full_seq is not None else antigen.amino_acids
        win_start = antigen.source_protein_offset if full_seq is not None else 0
        win_end = win_start + len(antigen.amino_acids)
        full_excl = []
        full_genes = set()
        win_excl = []
        win_genes = set()
        for i in range(max(0, len(scan_seq) - 7)):
            kmer = scan_seq[i:i + 8]
            matches = hits.get(kmer, ())
            excl_matches = [m for m in matches if m.get('role') == 'exclude']
            if not excl_matches:
                continue
            iv_s, iv_e = i, i + 8
            full_excl.append((iv_s, iv_e))
            genes = {
                occ.get('gene_name') or occ.get('gene_id')
                for m in excl_matches
                for occ in m.get('occurrences', ())
                if occ.get('shared_sequence', m.get('role')) == 'exclude'
                and (occ.get('gene_name') or occ.get('gene_id'))
            }
            full_genes.update(genes)
            if iv_s >= win_start and iv_e <= win_end:
                win_excl.append((iv_s - win_start, iv_e - win_start))
                win_genes.update(genes)
        contexts[key] = dict(
            full_sequence=full_seq,
            full_excluded_intervals=_merge_intervals(full_excl) if full_seq is not None else (),
            full_excluded_genes=tuple(sorted(full_genes)) if full_seq is not None else (),
            window_excluded_intervals=_merge_intervals(win_excl),
            window_excluded_genes=tuple(sorted(win_genes)),
        )
    return contexts


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


def _format_excluded_title(genes):
    if genes:
        return 'Excluded: overlaps Unsafe/self-reference protein (' + ', '.join(genes[:6]) + (
            ', ...' if len(genes) > 6 else '') + ')'
    return 'Excluded: overlaps Unsafe or self-reference protein'


def _interval_aa_count(intervals):
    return sum(end - start for start, end in intervals if end > start)


def antigen_biology(antigen, *, source_evidence=None, peptides=(), protein_context=None):
    """Summarize recorded source facts, keeping RNA and translation separate.

    Comparison sequences come only from producer exports or retained peptide
    comparators. Source offsets, rather than sequence search, locate a window
    even when the original protein contains repeated sequences.
    """
    from .antigen_origin import biology_origin
    candidates = isovar_candidates(antigen)
    facts, comparisons, junctions = [], [], []
    source_full_sequence = None
    identities = [value for key, value in antigen.source_metadata if key == 'isovar_evidence_id']
    evidence = (source_evidence or {}).get(identities[0], {}) if len(identities) == 1 else {}
    result = evidence.get('isovar_result', {})
    origin = biology_origin(antigen, result=result, candidates=candidates, peptides=peptides)
    if origin is not None:
        from .report_display import compact_allele
        for variant in origin.variants:
            facts.append(('Variant type', variant.variant_type + (':' + variant.sv_type if variant.sv_type else '')))
            if variant.contig and variant.start is not None:
                facts.append(('Variant', '%s:%s %s → %s' % (
                    variant.contig,
                    variant.start,
                    compact_allele(variant.ref) or '∅',
                    compact_allele(variant.alt) or '∅',
                )))
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
        if original is not None and source_full_sequence is None:
            source_full_sequence = original
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
        # Merge overlapping equal-length recorded epitope comparators into one
        # contiguous window comparison per label (zero annotation queries; uses
        # only recorded comparators) instead of emitting a separate grey box
        # for every overlapping 8-11-mer.
        by_label = {}
        for peptide in peptides:
            for epitope in peptide.epitopes:
                comparator = epitope.reference_comparator
                if comparator is not None and comparator.sequence:
                    label = 'Patient background' if epitope.comparison_evidence else 'Reference peptide'
                    by_label.setdefault(label, []).append((epitope, comparator.sequence))
        for label, pairs in by_label.items():
            ref_by_pos = {}
            consistent = True
            for epitope, comp_seq in pairs:
                off = getattr(epitope, 'offset', None)
                if (
                    off is None
                    or len(comp_seq) != len(epitope.sequence)
                    or antigen.amino_acids[off:off + len(epitope.sequence)] != epitope.sequence
                ):
                    consistent = False
                    break
                for i, ch in enumerate(comp_seq):
                    pos = off + i
                    prev = ref_by_pos.get(pos)
                    if prev is not None and prev != ch:
                        consistent = False
                        break
                    ref_by_pos[pos] = ch
                if not consistent:
                    break
            if consistent and ref_by_pos:
                cov_start = min(ref_by_pos)
                cov_end = max(ref_by_pos) + 1
                if len(ref_by_pos) == cov_end - cov_start:
                    ref_slice = ''.join(ref_by_pos[p] for p in range(cov_start, cov_end))
                    tum_slice = antigen.amino_acids[cov_start:cov_end]
                    diff = sequence_difference(tum_slice, ref_slice)
                    diff['tumor_start'] = antigen.source_protein_offset + cov_start + 1
                    diff['tumor_end'] = antigen.source_protein_offset + cov_end
                    diff['reference_start'] = diff['tumor_start']
                    diff['reference_end'] = diff['tumor_end']
                    comparisons.append(dict(label=label, **diff))
                    continue
            # Fallback when comparators differ in length or are disjoint: keep
            # the single highest-scoring representative epitope comparator.
            best_ep, best_comp = max(
                pairs,
                key=lambda item: (getattr(item[0], 'score', 0.0) or 0.0),
            )
            comparisons.append(dict(
                label=label,
                **sequence_difference(best_ep.sequence, best_comp),
            ))
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

    all_epitopes = [e for p in peptides for e in p.epitopes]
    ctx = protein_context or {}
    eff = effective_target_intervals(
        antigen,
        all_epitopes,
        extra_excluded=ctx.get('window_excluded_intervals', ()),
        extra_genes=ctx.get('window_excluded_genes', ()),
    )
    css_class = 'targetable' if antigen.kind == 'CTA' else 'mutant'
    excluded_title = _format_excluded_title(eff['excluded_genes'])
    sequence_html = marked_sequence(
        antigen.amino_acids,
        eff['target'],
        css_class=css_class,
        excluded_intervals=eff['excluded'],
        excluded_title=excluded_title,
    )

    full_seq = ctx.get('full_sequence') or source_full_sequence
    full_protein_view = None
    if full_seq and len(full_seq) > len(antigen.amino_acids):
        win_start = antigen.source_protein_offset
        win_end = win_start + len(antigen.amino_acids)
        if full_seq[win_start:win_end] == antigen.amino_acids:
            shifted_win_excl = tuple((s + win_start, e + win_start) for s, e in eff['excluded'])
            full_excl = _merge_intervals(
                tuple(ctx.get('full_excluded_intervals', ())) + shifted_win_excl
            )
            full_genes = tuple(sorted(set(ctx.get('full_excluded_genes', ())) | set(eff['excluded_genes'])))
            if antigen.kind == 'CTA':
                full_target = _subtract_intervals([(0, len(full_seq))], full_excl)
            else:
                shifted_win_target = tuple((s + win_start, e + win_start) for s, e in eff['target'])
                full_target = _subtract_intervals(shifted_win_target, full_excl)
            cropped = len(full_seq) > 450
            disp_start = max(0, win_start - 100) if cropped else 0
            disp_end = min(len(full_seq), win_end + 100) if cropped else len(full_seq)
            full_protein_view = dict(
                length=len(full_seq),
                display_start=disp_start + 1,
                display_end=disp_end,
                cropped=cropped,
                sequence_html=marked_sequence(
                    full_seq[disp_start:disp_end],
                    _clip_intervals(full_target, disp_start, disp_end),
                    css_class=css_class,
                    excluded_intervals=_clip_intervals(full_excl, disp_start, disp_end),
                    excluded_title=_format_excluded_title(full_genes),
                    window_intervals=[(win_start - disp_start, win_end - disp_start)],
                ),
                full_excluded_aa=_interval_aa_count(full_excl),
                full_excluded_genes=full_genes,
                window_excluded_aa=_interval_aa_count(eff['excluded']),
                window_excluded_genes=eff['excluded_genes'],
            )

    return dict(kind=antigen_source_label(antigen, source_evidence),
                gene=antigen.display_gene_name, identifier=antigen.display_identifier,
                facts=list(dict.fromkeys(facts)), junctions=junctions, comparisons=comparisons,
                sequence=antigen.amino_acids,
                sequence_html=sequence_html,
                target_intervals=eff['target'],
                excluded_intervals=eff['excluded'],
                excluded_genes=eff['excluded_genes'],
                excluded_aa=_interval_aa_count(eff['excluded']),
                full_protein=full_protein_view,
                sequence_start=antigen.source_protein_offset + 1,
                sequence_end=antigen.source_protein_offset + len(antigen.amino_acids),
                targeting=targeting,
                targeting_preview=sorted(targeting, key=lambda row: row['status_code'] == 'eligible')[:8],
                admission=antigen.tumor_specificity.status)
