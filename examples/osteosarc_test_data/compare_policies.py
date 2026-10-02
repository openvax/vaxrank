"""Reproducible Sid comparison with cached NetMHCpan and optional fresh models.

Run: python -m examples.osteosarc_test_data.compare_policies --output NEW_DIR
Add --processing to explicitly run MHCflurry presentation and Pepsickle. This
experiment uses the small pinned reference, not a complete normal-tissue atlas.
"""

import argparse
from collections import Counter
from dataclasses import asdict, replace
from hashlib import sha256
from importlib.metadata import version
import json
import logging
import os
from pathlib import Path
import shutil
import tempfile

import msgspec
import pandas as pd
from topiary import CachedPredictor, TopiaryPredictor, TopiaryResult

from tests.osteosarc_selection_helpers import DATA, load_selection_inputs, reconstruct_selection
from vaxrank.prediction_input import finite_prediction_value
from vaxrank.cleavage_inference import audit_pepsickle_inputs
from vaxrank.config.loader import load_vaxrank_config, extract_epitope_config_kwargs
from vaxrank.core_logic import vaccine_peptides_for_variant, vaccine_peptides_from_epitopes
from vaxrank.epitope_config import EpitopeConfig
from vaxrank.epitope_dataset import EpitopeDataset
from vaxrank.epitope_dsl import attach_per_allele_scores, epitopes_for_ranking, epitopes_to_topiary_df
from vaxrank.epitope_logic import predict_epitopes, _prediction_from_row
from vaxrank.mutant_protein_fragment import MutantProteinFragment
from vaxrank.input_scope import InputScope, InputProvenance
from vaxrank.vaccine_antigen import VaccineAntigen
from vaxrank.peptide import PeptideConstructConfig, assemble_peptide_constructs
from vaxrank.selection_policy import save_policy_evaluations


def policy_config(name):
    paths = ['builtin:openvax-v1']
    if name != 'openvax-v1':
        paths.append('builtin:' + name)
    return msgspec.convert(extract_epitope_config_kwargs(
        load_vaxrank_config(config_path=paths)), EpitopeConfig)


def add_models(epitopes, mhc_frame, profile):
    """Join predictions by exact source offset, preserving original observations."""
    additional = {}
    for _, row in mhc_frame.iterrows():
        pred = _prediction_from_row(row, peptide=row.peptide, score=row.score,
                                    percentile_rank=finite_prediction_value(row.get('percentile_rank')))
        pred = replace(pred, n_flank=row.get('n_flank'), c_flank=row.get('c_flank'))
        additional.setdefault((row.peptide, int(row.peptide_offset)), []).append(pred)
    updated = [replace(e, predictions=tuple(e.predictions_flat()) + tuple(
        additional.get((e.sequence, e.offset), ()))) for e in epitopes]
    frame = epitopes_to_topiary_df(updated)
    site_scores = {s.bond: s.score for s in profile.sites
                   if s.bond not in profile.padded_bonds}
    # Do not reinterpret the source endpoint or padded fragment context as a
    # well-observed cleavage site. Unknown features remain missing in Topiary.
    frame['pepsickle_cterm_score'] = [site_scores.get(int(row.peptide_offset) + len(row.peptide))
                                     for row in frame.itertuples()]
    return updated, frame


def selected_record(peptides):
    return [dict(sequence=p.amino_acids, target_score=p.target_epitope_score,
                 self_score=p.self_epitope_score, combined_score=p.combined_score,
                 targets=sorted({e.sequence for e in p.target_epitopes})) for p in peptides]


def compact_comparison(output):
    """Summarize measured results without duplicating full prediction tables."""
    source = output / 'comparison.json'
    manifest = json.loads(source.read_text())
    summary = {key: value for key, value in manifest.items() if key != 'records'}
    summary['comparison_sha256'] = sha256(source.read_bytes()).hexdigest()
    summary['records'] = []
    for record in manifest['records']:
        entry = {key: value for key, value in record.items() if key != 'alternatives'}
        entry['alternatives'] = {}
        for name, alternative in record.get('alternatives', {}).items():
            serum = {}
            for weight, constructs in alternative['serum_sensitivity'].items():
                rejected = [row for audit in alternative['serum_decisions'][weight]
                            for row in audit['alternatives'] if not row['eligible']]
                serum[weight] = dict(constructs=[dict(sequence=c['sequence'], **{
                    key: c['components']['window_selection'][key] for key in (
                        'target_score', 'non_cta_self_score', 'at_risk_target_score',
                        'utility', 'unassessed_models')}) for c in constructs],
                    rejected_counts=dict(Counter(row['ineligibility_reason'] for row in rejected)),
                    unassessed_models=sorted({model for row in rejected for model in row['unassessed_models']}))
            entry['alternatives'][name] = dict(windows=alternative['windows'],
                self_trim=alternative['self_trim'], serum_sensitivity=serum)
        summary['records'].append(entry)
    (output / 'comparison-summary.json').write_text(json.dumps(summary, indent=2) + '\n')


def main(output, processing=False):
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    logging.getLogger().setLevel(logging.ERROR)
    documented = json.loads((DATA / 'documented.json').read_text())
    records = documented['records']
    predictor = TopiaryPredictor(models=[CachedPredictor.from_topiary_output(DATA / 'netmhcpan42.tsv')])
    packages = ('vaxrank', 'topiary', 'mhctools', 'isovar', 'varcode')
    if processing:
        packages += ('mhcflurry', 'pepsickle')
    manifest = dict(packages={name: version(name) for name in packages},
        netmhcpan_sha256=sha256((DATA / 'netmhcpan42.tsv').read_bytes()).hexdigest(),
        reference='Pinned small Sid transcript fixture; not complete human self coverage',
        documented_settings='Historical provider settings unavailable; not a tuning target', records=[])
    manifest['hla'] = {key: documented[key] for key in (
        'clinical_class_i', 'prediction_alleles', 'unassessed_hla', 'typing_discrepancy')}
    with tempfile.TemporaryDirectory(prefix='vaxrank-policy-sid-') as directory:
        root = Path(directory)
        os.environ['VAXRANK_REF_PEPTIDES_DIR'] = str(root / 'kmers')
        # Native mutation records retain pyensembl reference handles. Keep their
        # local source and index alive after this process (not in test tempdirs).
        reference = output / 'input-reference'
        shutil.copytree(DATA / 'isovar' / 'reference', reference)
        inputs = load_selection_inputs(output, reference_directory=reference)
        cache = {}
        for record in records:
            key = record['variant_id'], len(record['native_sequence'])
            result, config = reconstruct_selection(inputs, *key)
            evaluations = []
            legacy = vaccine_peptides_for_variant(result, predictor, vaccine_config=config)
            frozen = vaccine_peptides_for_variant(result, predictor, vaccine_config=config,
                                                  epitope_config=policy_config('openvax-v1'),
                                                  policy_evaluations=evaluations)
            assert selected_record(legacy) == selected_record(frozen), record['id']
            case_dir = output / record['id']
            save_policy_evaluations(evaluations, case_dir / 'openvax-v1')
            entry = dict(id=record['id'], variant_id=record['variant_id'],
                         passes_rna_filters=result.passes_all_filters,
                         documented_native=record['native_sequence'],
                         legacy=selected_record(legacy), openvax_v1=selected_record(frozen),
                         baseline_parity=True,
                         documented_match=bool(frozen and frozen[0].amino_acids == record['native_sequence']))
            manifest['records'].append(entry)
            if processing and result.passes_all_filters:
                cache[key] = result, config
        (output / 'comparison.json').write_text(json.dumps(manifest, indent=2) + '\n')
        if processing:
            from mhctools import MHCflurry
            from mhctools.cleavage import CleavageInput
            model = MHCflurry(alleles=documented['prediction_alleles'],
                             default_peptide_lengths=documented['prediction_peptide_lengths'],
                             presentation_allele_mode='per_allele')
            fresh = TopiaryPredictor(models=[model])
            contexts = {str(i): r.top_protein_sequence.amino_acids for i, (r, _) in enumerate(cache.values())}
            fresh_frame = fresh.predict_from_named_sequences(contexts)
            TopiaryResult(fresh_frame).to_tsv(output / 'mhcflurry-evidence.tsv')
            profiles = audit_pepsickle_inputs([CleavageInput(seq, source_id=name)
                                               for name, seq in contexts.items()])
            (output / 'pepsickle-evidence.json').write_text(json.dumps([asdict(p) for p in profiles], indent=2))
            for i, (key, (result, config)) in enumerate(cache.items()):
                fragment = MutantProteinFragment.from_isovar_result(result)
                original = predict_epitopes(predictor, fragment,
                                            EpitopeConfig(score_expr='1.0', min_epitope_score=0.),
                                            genome=result.variant.ensembl)
                provenance = InputProvenance(
                    input_id='sid-policy-' + str(i), source_format='vcf_bam',
                    path=str(DATA / 'netmhcpan42.tsv'), content_sha256=manifest['netmhcpan_sha256'],
                    scope=InputScope(patient_id='Sid-fixture', reference_assembly='GRCh38',
                                     sample_id=key[0], annotation='pinned Sid transcript subset'),
                    observed_mhc_alleles=tuple(documented['prediction_alleles']))
                original = [replace(e, input_provenance=provenance) for e in original]
                source_frame = fresh_frame.loc[fresh_frame.source_sequence_name == str(i)]
                enriched, frame = add_models(original, source_frame, profiles[i])
                entry = next(e for e in manifest['records'] if (e['variant_id'], len(e['documented_native'])) == key)
                entry['processing_status'] = profiles[i].status
                entry['cterm_coverage'] = float(frame.pepsickle_cterm_score.notna().mean())
                entry['alternatives'] = {}
                for name in ('presentation-v1', 'presentation-ba-v1', 'presentation-ba-cterm-v1'):
                    cfg = policy_config(name)
                    evaluations = []
                    scored = attach_per_allele_scores(enriched, cfg, topiary_df=frame,
                                                     policy_evaluations=evaluations)
                    selected = vaccine_peptides_from_epitopes(
                        result.variant, fragment, epitopes_for_ranking(scored, cfg),
                        vaccine_config=config, epitope_config=cfg)
                    case_dir = output / entry['id'] / name
                    save_policy_evaluations(evaluations, case_dir)
                    dataset = EpitopeDataset.from_predictions(scored, result=TopiaryResult(frame),
                                                               config=cfg, policy_evaluations=evaluations)
                    dataset.antigens = {e.prediction_group_source: VaccineAntigen.from_mutant_protein_fragment(fragment)
                                        for e in scored}
                    dataset.mutation_fragments = {e.prediction_group_source: fragment for e in scored}
                    dataset.save(case_dir / 'epitopes.tsv')
                    alternative = dict(windows=selected_record(selected), serum_sensitivity={})
                    self_config = msgspec.structs.replace(config,
                        min_peptide_length=min(15, key[1]),
                        combined_score_expr='sqrt(n_rna_alt) * window_epitope_score',
                        window_selection=dict(self_weight=.25, min_target_fraction=.95))
                    trimmed = vaccine_peptides_from_epitopes(
                        result.variant, fragment, epitopes_for_ranking(scored, cfg),
                        vaccine_config=self_config, epitope_config=cfg)
                    alternative['self_trim'] = selected_record(trimmed)
                    alternative['serum_decisions'] = {}
                    for weight in (0., .25, .5):
                        audit = []
                        constructs = assemble_peptide_constructs([(result.variant, selected)],
                            PeptideConstructConfig(min_antigen_length_aa=min(15, key[1]),
                                max_antigen_length_aa=key[1], window_selection=dict(
                                    self_weight=.25, serum_weight=weight, min_target_fraction=.95)),
                            window_audit=audit)
                        alternative['serum_sensitivity'][str(weight)] = [asdict(c) for c in constructs]
                        alternative['serum_decisions'][str(weight)] = audit
                    entry['alternatives'][name] = alternative
    (output / 'comparison.json').write_text(json.dumps(manifest, indent=2, allow_nan=False) + '\n')
    rows = [dict(id=e['id'], passes_rna_filters=e['passes_rna_filters'],
                 baseline_parity=e['baseline_parity'], documented_match=e['documented_match'],
                 selected=';'.join(p['sequence'] for p in e['openvax_v1'])) for e in manifest['records']]
    pd.DataFrame(rows).to_csv(output / 'comparison.tsv', sep='\t', index=False)
    compact_comparison(output)
    hashes = {str(path.relative_to(output)): sha256(path.read_bytes()).hexdigest()
              for path in output.rglob('*') if path.is_file()}
    (output / 'evidence-sha256.json').write_text(json.dumps(hashes, indent=2) + '\n')
    return manifest


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--processing', action='store_true')
    args = parser.parse_args()
    main(args.output, args.processing)
