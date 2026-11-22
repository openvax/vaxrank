"""Small, evidence-preserving views for the human-readable report."""

from markupsafe import Markup, escape


_RUN_PROVENANCE_ARGUMENT_NAMES = frozenset({
    '_declared_mhc_alleles',
    '_inferred_mhc_alleles_from_lens',
    '_input_species_scope',
})


def run_provenance_arguments(arguments):
    """Return non-empty resolved scope that belongs in human-readable output.

    Runtime state is stored on the argument namespace with leading underscores,
    but an underscore alone does not make a value either reportable or internal.
    Keep this allow-list as the one classification rule used by the console and
    report renderers.
    """
    return {
        name: value
        for name, value in arguments.items()
        if name in _RUN_PROVENANCE_ARGUMENT_NAMES and bool(value)
    }


def _merge_intervals(intervals):
    """Merge overlapping/adjacent half-open intervals, keeping zero-width junctions."""
    junctions = sorted({(start, end) for start, end in intervals if start == end})
    spans = sorted((start, end) for start, end in intervals if end > start)
    merged = []
    for start, end in spans:
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return tuple(sorted(merged + junctions))


def _subtract_intervals(intervals, excluded):
    """Strictly subtract half-open excluded spans from targetable intervals."""
    excl = [iv for iv in _merge_intervals(excluded) if iv[1] > iv[0]]
    if not excl:
        return _merge_intervals(intervals)
    result = []
    for start, end in _merge_intervals(intervals):
        if start == end:
            if not any(ex_s < start < ex_e for ex_s, ex_e in excl):
                result.append((start, end))
            continue
        cur = start
        for ex_s, ex_e in excl:
            if ex_e <= cur:
                continue
            if ex_s >= end:
                break
            if ex_s > cur:
                result.append((cur, min(ex_s, end)))
            cur = max(cur, ex_e)
            if cur >= end:
                break
        if cur < end:
            result.append((cur, end))
    return tuple(sorted(result))


def _clip_intervals(intervals, start, end):
    """Clip and rebase intervals onto slice [start, end)."""
    clipped = []
    for iv_s, iv_e in intervals:
        if iv_s == iv_e:
            if start <= iv_s <= end:
                clipped.append((iv_s - start, iv_e - start))
            continue
        lo = max(iv_s, start)
        hi = min(iv_e, end)
        if hi > lo:
            clipped.append((lo - start, hi - start))
    return _merge_intervals(clipped)


def excluded_self_intervals(antigen, epitopes=(), *, extra_excluded=(), extra_genes=()):
    """Intervals on ``antigen.amino_acids`` matching Unsafe or self-reference proteins.

    For CTAs, any self-reference or OncoRef-excluded epitope marks its span.
    For non-CTAs (mutations, SVs, ncORFs), only epitopes that overlap the
    targetable mask and still match self-reference are marked as excluded so
    expected wild-type flanking residues are not painted as hazards.
    """
    from .cta_targeting import targeting_evidence

    source = antigen.amino_acids
    is_cta = antigen.kind == 'CTA'
    spans = list(extra_excluded or ())
    genes = set(extra_genes or ())

    for epitope in epitopes or ():
        offset = getattr(epitope, 'offset', -1)
        seq = getattr(epitope, 'sequence', '')
        if (
            getattr(epitope, 'source_sequence', None) != source
            or offset < 0
            or source[offset:offset + len(seq)] != seq
        ):
            continue
        end = offset + len(seq)
        te = targeting_evidence(epitope)
        cta_excluded = te is not None and te.get('status') == 'excluded'
        is_self = bool(getattr(epitope, 'occurs_in_self_reference', False)) or cta_excluded
        if not is_self:
            continue
        if not is_cta:
            overlaps = getattr(epitope, 'overlaps_targetable', None)
            if overlaps is None:
                overlaps = antigen.interval_is_targetable(offset, end)
            if not overlaps:
                continue
        spans.append((offset, end))
        match = getattr(epitope, 'self_reference_match', None)
        if match is not None:
            for src in getattr(match, 'sources', ()) or ():
                label = getattr(src, 'gene_name', '') or getattr(src, 'gene_id', '')
                if label:
                    genes.add(label)
        if te is not None:
            for m in te.get('matches', ()) or ():
                if m.get('role') != 'exclude':
                    continue
                for occ in m.get('occurrences', ()) or ():
                    if occ.get('shared_sequence', m.get('role')) == 'exclude':
                        label = occ.get('gene_name') or occ.get('gene_id') or ''
                        if label:
                            genes.add(label)

    return _merge_intervals(spans), tuple(sorted(genes))


def effective_target_intervals(antigen, epitopes=(), *, extra_excluded=(), extra_genes=()):
    """Compute targetable and self/Unsafe-excluded intervals with strict non-overlap."""
    from .cta_targeting import targeting_evidence

    excluded, excluded_genes = excluded_self_intervals(
        antigen, epitopes, extra_excluded=extra_excluded, extra_genes=extra_genes,
    )
    raw_target = tuple((i.start, i.end) for i in antigen.targetable_mask.intervals)
    if antigen.kind == 'CTA':
        source = antigen.amino_acids
        eligible_spans = []
        for epitope in epitopes or ():
            offset = getattr(epitope, 'offset', -1)
            seq = getattr(epitope, 'sequence', '')
            if (
                getattr(epitope, 'source_sequence', None) != source
                or offset < 0
                or source[offset:offset + len(seq)] != seq
            ):
                continue
            te = targeting_evidence(epitope)
            if te is not None and te.get('status') != 'eligible':
                continue
            if getattr(epitope, 'is_targetable', False) and not getattr(epitope, 'occurs_in_self_reference', False):
                eligible_spans.append((offset, offset + len(seq)))
        base = _merge_intervals(eligible_spans) if eligible_spans else raw_target
        target = _subtract_intervals(base, excluded)
    else:
        target = _subtract_intervals(raw_target, excluded)
    return dict(
        target=target,
        excluded=excluded,
        excluded_genes=excluded_genes,
    )


def marked_sequence(
    sequence,
    intervals=(),
    css_class="mutant",
    *,
    excluded_intervals=(),
    excluded_title="Excluded: matches Unsafe or self-reference protein",
    window_intervals=(),
):
    """Render exact targetable, excluded, window, and junction intervals cleanly.

    Excluded intervals take precedence over targetable intervals so residues
    overlapping Unsafe/self-reference proteins are never painted as targetable.
    """
    starts = {start for start, end in intervals if start == end}
    excl_spans = [iv for iv in excluded_intervals if iv[1] > iv[0]]
    target_spans = [iv for iv in intervals if iv[1] > iv[0]]
    win_spans = [iv for iv in window_intervals if iv[1] > iv[0]]

    def _state(idx):
        in_win = any(s <= idx < e for s, e in win_spans)
        if any(s <= idx < e for s, e in excl_spans):
            return ('self-excluded', in_win)
        if any(s <= idx < e for s, e in target_spans):
            return (css_class, in_win)
        return ('', in_win)

    def _flush(buf, state):
        if not buf:
            return
        text = escape(''.join(buf))
        cls, in_win = state
        if cls == 'self-excluded':
            node = Markup('<span class="self-excluded" title="{}">{}</span>').format(
                excluded_title, text)
        elif cls:
            node = Markup('<span class="{}">{}</span>').format(cls, text)
        else:
            node = text
        if in_win:
            node = Markup('<span class="vaccine-window" title="Selected vaccine window">{}</span>').format(node)
        pieces.append(node)

    pieces = []
    buf = []
    cur_state = None
    for index, residue in enumerate(sequence):
        if index in starts:
            _flush(buf, cur_state)
            buf = []
            cur_state = None
            pieces.append(Markup('<span class="junction" title="Target junction">|</span>'))
        if index and index % 20 == 0:
            _flush(buf, cur_state)
            buf = []
            cur_state = None
            pieces.append(Markup('<wbr>'))
        state = _state(index)
        if state != cur_state:
            _flush(buf, cur_state)
            buf = [residue]
            cur_state = state
        else:
            buf.append(residue)
    _flush(buf, cur_state)
    if len(sequence) in starts:
        pieces.append(Markup('<span class="junction" title="Target junction">|</span>'))
    return Markup('').join(pieces)


def reference_affinity_values(epitope, allele, predictor, version):
    """Distinct recorded affinities for the exact displayed comparison model.

    No values means unmeasured; multiple values means conflicting evidence.
    Both human-readable and spreadsheet reports preserve that distinction
    rather than choosing a leaf by sort order or another predictor version.
    """
    comparator = epitope.reference_comparator
    return {p.value for p in comparator.predictions_flat()
            if p.kind == 'pMHC_affinity' and p.allele == allele
            and p.predictor_name == predictor and p.predictor_version == version
            and p.value is not None} if comparator is not None else set()


def epitope_display(epitope, antigen, row, prediction):
    """Pair recorded sequences; highlight only coordinates supported by the source."""
    offset = epitope.offset
    source = antigen.amino_acids
    located = (
        epitope.source_sequence == source
        and offset >= 0
        and source[offset:offset + len(epitope.sequence)] == epitope.sequence)
    intervals = ([(i.start, i.end) for i in antigen.targetable_mask.sliced(
        offset, offset + len(epitope.sequence)).intervals] if located else [])
    comparator = epitope.reference_comparator
    wt = comparator.sequence if comparator is not None else ''
    evidence = epitope.comparison_evidence
    comparison_label = ('Patient background (' + evidence['baseline_policy'] + ')'
                        if evidence else 'Annotation / supplied WT')
    comparison_notes = []
    if evidence:
        comparison_notes.append(evidence['comparison'].replace('_', ' ') + '.')
        comparison_notes.extend(sorted({a['reason'].replace('_', ' ') + '.'
                                        for a in evidence.get('alternatives', ()) if a.get('reason')}))
        if evidence.get('reason'):
            comparison_notes.append(evidence['reason'].replace('_', ' ') + '.')
        prediction_state = evidence.get('prediction_state')
        if prediction_state == 'outside_requested_peptide_lengths':
            comparison_notes.append('Comparator unmeasured: length outside the requested prediction lengths.')
        elif prediction_state == 'no_measurements_recorded':
            comparison_notes.append('Comparator sequence available; no predictions recorded.')
    if wt and len(wt) != len(epitope.sequence):
        comparison_notes.append('Different lengths; no residue alignment inferred.')
    differences = [(i, i + 1) for i, (target, reference) in enumerate(zip(epitope.sequence, wt))
                   if target != reference] if (len(wt) == len(epitope.sequence)
                       and (evidence is None or evidence.get('position_aligned') is True)) else []
    return dict(
        values=row, version=prediction.predictor_version,
        target=marked_sequence(epitope.sequence, intervals),
        reference=marked_sequence(wt, differences, "reference-change"),
        location_known=located,
        comparison_label=comparison_label,
        comparison_note=' '.join(comparison_notes),
        presentation_label=('NetMHCpan Score_EL' if prediction.predictor_name.lower().startswith('netmhcpan')
                            else prediction.predictor_name + ' presentation'),
    )


def model_inventory(ranked, processing_predictions, processing_enabled):
    """Describe available evidence, without treating absent predictions as failures."""
    labels = {'pMHC_affinity': 'binding affinity', 'pMHC_presentation': 'presentation',
              'antigen_processing': 'antigen processing'}
    models = {}
    for _, peptides in ranked:
        for peptide in peptides:
            for epitope in peptide.epitopes:
                for prediction in epitope.predictions_flat():
                    key = (prediction.predictor_name, prediction.predictor_version)
                    models.setdefault(key, set()).add(labels.get(prediction.kind, prediction.kind))
    rows = [dict(name=name or 'Unnamed predictor', version=version,
                 status='Recorded: ' + ', '.join(sorted(kinds)))
            for (name, version), kinds in sorted(models.items())]
    for name in ('NetMHCpan', 'MHCflurry'):
        if not any(model.lower().startswith(name.lower()) for model, _ in models):
            rows.append(dict(name=name, version='', status='No results in the reported candidates.'))
    processing_versions = sorted({p.predictor_version for p in processing_predictions.values()
                                  if p.predictor_version})
    rows.append(dict(name='Pepsickle', version=', '.join(processing_versions), status=(
        'Recorded cleavage annotations; coverage can vary by epitope.' if processing_predictions else
        'Disabled for this run; cleavage unassessed.' if processing_enabled is False else
        'No cleavage results recorded; cleavage unassessed.')))
    return rows


def protein_context_key(antigen):
    """Stable lookup key for optional full-protein Unsafe/self context."""
    return (antigen.source_identifier, antigen.source_protein_offset, antigen.amino_acids)


def product_sequence_highlights(product, ranked, *, protein_contexts=None):
    """Map effective targetable and Unsafe/self-excluded intervals to product antigens."""
    from .assembly_antigens import source_matches_candidate

    protein_contexts = protein_contexts or {}
    masks = {antigen['name']: [[] for _ in antigen.get('sources', ())] for antigen in product.antigens}
    excluded_by_name = {antigen['name']: [] for antigen in product.antigens}
    genes_by_name = {antigen['name']: set() for antigen in product.antigens}
    kinds_by_name = {antigen['name']: set() for antigen in product.antigens}

    for _, peptides in ranked:
        for peptide in peptides:
            prediction_ids = {e.prediction_group_source for e in peptide.epitopes}
            input_ids = {e.input_provenance.input_id for e in peptide.epitopes if e.input_provenance}
            ctx = protein_contexts.get(protein_context_key(peptide.antigen)) or {}
            eff = effective_target_intervals(
                peptide.antigen,
                peptide.epitopes,
                extra_excluded=ctx.get('window_excluded_intervals', ()),
                extra_genes=ctx.get('window_excluded_genes', ()),
            )
            for antigen in product.antigens:
                sequence = antigen['aa']
                offset = peptide.amino_acids.find(sequence)
                if offset < 0 or peptide.amino_acids.find(sequence, offset + 1) >= 0:
                    target_slice = None
                    excl_slice = ()
                else:
                    end = offset + len(sequence)
                    target_slice = _clip_intervals(eff['target'], offset, end)
                    excl_slice = _clip_intervals(eff['excluded'], offset, end)
                for index, source in enumerate(antigen.get('sources', ())):
                    if source_matches_candidate(source, peptide.antigen, peptide.amino_acids, prediction_ids, input_ids):
                        masks[antigen['name']][index].append(target_slice)
                        excluded_by_name[antigen['name']].extend(excl_slice)
                        genes_by_name[antigen['name']].update(eff['excluded_genes'])
                        kinds_by_name[antigen['name']].add(peptide.antigen.kind)

    result = {}
    for name, sources in masks.items():
        values = [value for source in sources for value in source]
        excluded = _merge_intervals(excluded_by_name[name])
        agreed = (
            values[0] if sources and all(sources) and values[0] is not None
            and all(v == values[0] for v in values) else None
        )
        if agreed is not None and excluded:
            agreed = _subtract_intervals(agreed, excluded)
        kinds = kinds_by_name[name]
        result[name] = dict(
            target=agreed,
            excluded=excluded,
            excluded_genes=tuple(sorted(genes_by_name[name])),
            is_cta=bool(kinds) and kinds == {'CTA'},
        )
    return result


def product_target_intervals(product, ranked, *, protein_contexts=None):
    """Map candidate masks to products only when all matching placements agree."""
    highlights = product_sequence_highlights(product, ranked, protein_contexts=protein_contexts)
    return {name: info['target'] for name, info in highlights.items()}
