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


def marked_sequence(sequence, intervals=(), css_class="mutant"):
    """Render exact intervals, including deletion junctions, without filling gaps."""
    starts = {start for start, end in intervals if start == end}
    pieces = []
    for index, residue in enumerate(sequence):
        if index in starts:
            pieces.append(Markup('<span class="junction" title="Target junction">|</span>'))
        if index and index % 20 == 0:
            pieces.append(Markup('<wbr>'))
        value = escape(residue)
        if any(start <= index < end for start, end in intervals):
            value = Markup('<span class="{}">{}</span>').format(css_class, value)
        pieces.append(value)
    if len(sequence) in starts:
        pieces.append(Markup('<span class="junction" title="Target junction">|</span>'))
    return Markup('').join(pieces)


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
    wt = epitope.wt.sequence if epitope.wt is not None else ''
    differences = [(i, i + 1) for i, (target, reference) in enumerate(zip(epitope.sequence, wt))
                   if target != reference] if len(wt) == len(epitope.sequence) else []
    return dict(
        values=row, version=prediction.predictor_version,
        target=marked_sequence(epitope.sequence, intervals),
        reference=marked_sequence(wt, differences, "reference-change"),
        location_known=located,
        comparison_note=('Different lengths; no residue alignment inferred.'
                         if wt and len(wt) != len(epitope.sequence) else ''),
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


def product_target_intervals(product, ranked):
    """Map candidate masks to products only when all matching placements agree."""
    from .assembly_antigens import source_matches_candidate

    masks = {antigen['name']: [[] for _ in antigen.get('sources', ())] for antigen in product.antigens}
    for _, peptides in ranked:
        for peptide in peptides:
            prediction_ids = {e.prediction_group_source for e in peptide.epitopes}
            input_ids = {e.input_provenance.input_id for e in peptide.epitopes if e.input_provenance}
            for antigen in product.antigens:
                sequence = antigen['aa']
                offset = peptide.amino_acids.find(sequence)
                intervals = None
                if offset < 0 or peptide.amino_acids.find(sequence, offset + 1) >= 0:
                    intervals = None
                else:
                    mask = peptide.antigen.targetable_mask.sliced(offset, offset + len(sequence))
                    intervals = tuple((i.start, i.end) for i in mask.intervals)
                for index, source in enumerate(antigen.get('sources', ())):
                    if source_matches_candidate(source, peptide.antigen, peptide.amino_acids, prediction_ids, input_ids):
                        masks[antigen['name']][index].append(intervals)
    result = {}
    for name, sources in masks.items():
        values = [value for source in sources for value in source]
        result[name] = (values[0] if sources and all(sources) and values[0] is not None
                        and all(v == values[0] for v in values) else None)
    return result
