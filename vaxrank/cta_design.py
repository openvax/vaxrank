"""CTA design decisions and unfiltered measurements of actual final products."""

from pathlib import Path

import pandas as pd
from topiary import TopiaryResult

from .cta_expression import CTAExpressionResult
from .native_serialization import from_native_json, to_native_json


def selected_regions(ranked, dataset):
    """Recover exact selected coordinates from persisted native occurrences."""
    originals = {e.prediction_id: e for e in dataset.epitopes}
    from .vaccine_library import iter_named_antigens
    names = {id(p): name for name, _, p in iter_named_antigens(
        ranked, candidates_per_slot=max((len(peptides) for _, peptides in ranked), default=1))}
    regions = []
    for source, peptides in ranked:
        for peptide in peptides:
            coordinates = {originals[e.prediction_id].offset - e.offset
                           for e in peptide.epitopes if e.prediction_id in originals}
            if len(coordinates) != 1:
                raise ValueError('CTA selected window requires unambiguous native occurrence coordinates')
            start = coordinates.pop()
            end = start + len(peptide.amino_acids)
            if source.amino_acids[start:end] != peptide.amino_acids:
                raise ValueError('CTA selected window differs from the recorded native protein')
            regions.append(dict(name=names[id(peptide)], antigen=source, start=start, end=end,
                                sequence=peptide.amino_acids, target_epitope_score=peptide.target_epitope_score,
                                epitopes=tuple(peptide.target_epitopes), window_selection=peptide.window_selection_audit))
    return regions


def write_cta_design(args, ranked, dataset, products, predictor=None):
    """Persist every target decision and actual product context for replay.

    Final-context MHC inventories are separate from native ranking. Neither
    predicts antigen processing, mature termini or clinical safety. No missing
    product/native mapping is filled with an arbitrary sequence alignment.
    """
    admission_payload = dataset.selection.get('cta_expression_admission')
    if admission_payload is None:
        return
    admission = from_native_json(admission_payload, CTAExpressionResult)
    regions = selected_regions(ranked, dataset)
    products = products or {}
    sequences = {modality + ':' + p.name: p.cds_aa if modality == 'mrna' else p.sequence
                 for modality, records in products.items() for p in records or ()}
    unsupported = {modality + ':' + p.name: 'unsupported_chemical_modifications'
                   for modality, records in products.items() for p in records or ()
                   if modality == 'peptide' and (p.components.get('n_terminal_acetylation')
                                                 or p.components.get('c_terminal_amidation'))}
    # A replayed audit is usable only for the identical products and native
    # source/region graph; a different policy cannot inherit an old audit.
    from .construct_replay import design_digest
    assembly_identities = {modality: args._construct_assembly[modality]['identity']
                           for modality in products}
    identity = design_digest(args, (products, regions, args._construct_policy, assembly_identities))
    saved = dataset.selection.get('cta_product_audit')
    if saved is not None and saved['identity'] == identity:
        audit = from_native_json(saved['payload'], dict)
    else:
        audit = dict(products=products, sequences=sequences, selected_regions=regions,
                     mhc_status='unassessed', reason='mhc_predictor_not_requested',
                     predictions=(), unassessed_products=unsupported,
                     unassessed_processing_compartments=('proteasomal', 'extracellular_serum', 'endolysosomal'))
        from .construct_sequence import ConstructEvidence, ConstructSequence, ConstructChemicalModification
        from .construct_audit import audit_construct_sequence
        from .version import __version__
        evidence = ConstructEvidence('vaxrank CTA design', 'documented',
                                     'Actual emitted product', vaccine_version=__version__)
        by_name = {r['name']: r for r in regions}
        product_audits = []
        for modality, records in products.items():
            for product in records:
                name = modality + ':' + product.name
                sequence = sequences[name]
                region = by_name.get(product.antigen_names[0]) if len(product.antigen_names) == 1 else None
                mapping = (dict(native_antigen=region['antigen'], native_start=region['start'], native_end=region['end'])
                           if region and sequence == region['sequence'] else
                           dict(mapping_status='unresolved', mapping_reason='composite_or_changed_native_region_mapping'))
                modifications = []
                if modality == 'peptide':
                    for key, label, offset in (
                            ('n_terminal_acetylation', 'N-terminal acetylation', 0),
                            ('c_terminal_amidation', 'C-terminal amidation', len(sequence))):
                        if product.components.get(key):
                            modifications.append(ConstructChemicalModification(offset, offset, label, evidence))
                record = ConstructSequence(name, sequence, modality, evidence,
                                           chemical_modifications=tuple(modifications), **mapping)
                product_audits.append(audit_construct_sequence(record, predictor, genome=args.genome))
        audit['product_audits'] = tuple(product_audits)
        emitted_names = {name for records in products.values() for p in records for name in p.antigen_names}
        audit['assembly_coverage'] = dict(
            selected_regions=len(regions), emitted_regions=sum(r['name'] in emitted_names for r in regions),
            unassembled_region_names=tuple(r['name'] for r in regions if r['name'] not in emitted_names),
            unresolved_product_antigen_names=tuple(sorted(emitted_names - set(by_name))))
        if predictor is not None and sequences:
            evidence = ConstructEvidence('vaxrank CTA design', 'documented',
                                         'Selected native protein region', vaccine_version=__version__)
            audit['native_region_audits'] = tuple(audit_construct_sequence(
                ConstructSequence(r['name'], r['sequence'], 'peptide', evidence,
                                  native_antigen=r['antigen'], native_start=r['start'], native_end=r['end']),
                predictor, genome=args.genome) for r in regions)
            requested = {name: sequence for name, sequence in sequences.items() if name not in unsupported}
            if requested:
                frame = predictor.predict_from_named_sequences(requested)
                from .allele_validation import validate_prediction_frame
                validate_prediction_frame(frame, 'CTA final products')
                for row in frame.itertuples():
                    if sequences[row.source_sequence_name][int(row.peptide_offset):int(row.peptide_offset) + len(row.peptide)] != row.peptide:
                        raise ValueError('Final product prediction differs from its exact sequence context')
                audit.update(mhc_status='predictions_returned' if not frame.empty else 'no_predictions',
                             reason='processing_unassessed', predictions=frame.to_dict('records'))
            else:
                audit['reason'] = 'unsupported_chemical_modifications'
        dataset.selection['cta_product_audit'] = dict(identity=identity, payload=to_native_json(audit))
    if args.output_dir:
        directory = Path(args.output_dir)
        directory.mkdir(parents=True, exist_ok=True)
        admitted_ids = {r['antigen'].source_identifier for r in regions}
        contributed_ids = {r['antigen'].source_identifier for r in regions
                           if r['name'] in {n for records in products.values() for p in records for n in p.antigen_names}}
        decisions = [dict(input_identifier=d.input_identifier, gene_id=d.gene_id,
                          canonical_gene_id=d.gene_identity.get('canonical_gene_id') if d.gene_identity else None,
                          identity_status=d.gene_identity.get('status') if d.gene_identity else 'legacy_unverified',
                          gene_name=d.gene_name, transcript_id=d.transcript_id,
                          expression_value=d.value, expression_unit=admission.input_contract.expression_unit,
                          measurement_level=admission.input_contract.measurement_level,
                          admission_status=d.status,
                          product_status='contributed' if d.assessment and d.assessment.antigen.source_identifier in contributed_ids
                          else 'not_assembled_or_mapping_unresolved' if d.assessment and d.assessment.antigen.source_identifier in admitted_ids
                          else 'not_selected',
                          design_status='selected_region' if d.assessment and d.assessment.antigen.source_identifier in admitted_ids
                          else 'no_eligible_region' if d.status == 'admitted' else d.status)
                     for d in admission.decisions]
        pd.DataFrame(decisions).to_csv(directory / 'cta_target_decisions.csv', index=False)
        (directory / 'cta_design.json').write_text(to_native_json(dict(
            admission=admission, final_product_audit=audit,
            interpretation='Public MS is not a patient measurement. MHC inventories and target masks do not establish processing, immunogenicity or clinical safety.')) + '\n')
        if audit['predictions']:
            TopiaryResult(pd.DataFrame(audit['predictions'])).to_tsv(directory / 'cta_product_predictions.tsv')
    if args.output_epitopes:
        dataset.save(args.output_epitopes)
