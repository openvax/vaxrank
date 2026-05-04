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



import logging
import logging.config

import sys
from importlib.resources import files


import pandas as pd
import serializable


from varcode import VariantCollection
from varcode.cli import variant_collection_from_args
from isovar import isovar_results_to_dataframe, run_isovar
from isovar.cli import (
    protein_sequence_creator_from_args,
    read_collector_from_args,
)
from isovar.cli.rna_args import alignment_file_from_args
from isovar.cli.filter_args import filter_threshold_dict_from_args
from mhctools.cli import (
    mhc_alleles_from_args,
    mhc_binding_predictor_from_args,
)

from .arg_parser import parse_vaxrank_args, check_args
from .epitope_config_args import epitope_config_from_args
from .vaccine_config_args import vaccine_config_from_args
from ..config import load_vaxrank_config

from ..core_logic import run_vaxrank
from ..epitope_io import (
    save_predictions,
    write_neoepitope_report,
)
from ..gene_pathway_check import GenePathwayCheck
from ..report import (
    make_ascii_report,
    make_html_report,
    make_pdf_report,
    make_csv_report,
    make_minimal_neoepitope_report,
    TemplateDataCreator,
)
from ..patient_info import PatientInfo
from ..mrna import RNAConstructConfig, assemble_mrna_constructs, write_mrna_outputs
from ..peptide import (
    PeptideConstructConfig,
    assemble_peptide_constructs,
    write_peptide_outputs,
)
from ..vaf import extract_dna_vaf_by_variant

logger = logging.getLogger(__name__)


def _filter_unannotatable_variants(variants):
    """
    Filter out variants whose contigs cannot be resolved by varcode/pyensembl.

    Some VCF files include variants on alt/decoy contigs (e.g.
    chr14_GL000194v1_random) that varcode cannot normalize to a contig
    in the pyensembl annotation database. These variants crash during
    gene annotation lookups and cannot contribute vaccine peptides
    since no gene annotations are available for them.

    We use varcode's own contig validation (via variant.gene_names)
    rather than comparing raw contig names, since varcode normalizes
    names (e.g. chr1 -> 1) before checking.

    See: https://github.com/openvax/vaxrank/issues/193
    """
    if len(variants) == 0:
        return variants

    valid = []
    skipped_contigs = set()
    for v in variants:
        try:
            # Triggers varcode's contig validation and normalization.
            # Raises ValueError if the contig can't be resolved.
            v.gene_names
            valid.append(v)
        except ValueError:
            skipped_contigs.add(v.contig)

    if skipped_contigs:
        logger.warning(
            "Skipping %d variant(s) on contigs not recognized by %s: %s. "
            "These cannot be annotated with gene information.",
            len(variants) - len(valid),
            variants[0].ensembl.reference_name,
            ", ".join(sorted(skipped_contigs)),
        )
        return VariantCollection(valid)
    return variants


def _has_external_input(args):
    return bool(
        getattr(args, 'input_pvacseq', None)
        or getattr(args, 'input_lens', None))


def _annotate_predictions_with_processing(ranked_vaccine_peptides,
                                          lens_predictions,
                                          human_only=False,
                                          threshold=0.5):
    """Run pepsickle credibility tagging across all EpitopePrediction
    records. Pulls them from BOTH the ranked-vaccine-peptides
    intermediate (one VP per variant; each VP has its own predictions)
    AND the LENS/pVACseq predictions list (which is the union).

    Annotation is in-place on the EpitopePrediction objects.

    Dedup strategy: identity (``id()``) wins when the two paths share
    the same Python object, but a future external-input loader that
    *copies* a prediction into the VP would break that assumption.
    To stay correct in both cases, we *also* dedup on the content
    tuple ``(peptide_sequence, allele, source_sequence, offset)``
    so a prediction reachable as two distinct objects with the same
    semantic identity still gets annotated only once. Both checks
    are O(1).
    """
    from ..processing import annotate_processing
    all_predictions = []
    seen_ids = set()
    seen_keys = set()

    def _key(p):
        return (
            getattr(p, 'peptide_sequence', '') or '',
            getattr(p, 'allele', '') or '',
            getattr(p, 'source_sequence', '') or '',
            getattr(p, 'offset', 0) or 0,
        )

    def _maybe_add(p):
        pid = id(p)
        if pid in seen_ids:
            return
        seen_ids.add(pid)
        # Empty / missing peptide_sequence: degenerate record. Don't
        # apply content-key dedup (every empty-peptide record would
        # collapse to one bucket); just pass it through and let the
        # annotation skip it on its own.
        if not getattr(p, 'peptide_sequence', None):
            all_predictions.append(p)
            return
        # Content-key dedup: skip if a *different* object with the
        # same (peptide, allele, source, offset) already in flight.
        k = _key(p)
        if k in seen_keys:
            return
        seen_keys.add(k)
        all_predictions.append(p)

    for _, peptides in (ranked_vaccine_peptides or []):
        for vp in peptides or []:
            for p in (
                    getattr(vp, 'mutant_epitope_predictions', None) or []):
                _maybe_add(p)
    for p in (lens_predictions or []):
        _maybe_add(p)
    if not all_predictions:
        return 0
    return annotate_processing(
        all_predictions, human_only=human_only, threshold=threshold)


def _epitope_config_from_args_safe(args):
    from ..epitope_config import EpitopeConfig
    try:
        return epitope_config_from_args(args)
    except Exception:
        return EpitopeConfig()


def _resolve_vaccine_type(args):
    """Return the active vaccine type (single value).

    ``--vaccine-type`` is single-valued (one mode per run; default
    ``peptide``). Argparse already validates the choice, but tolerate
    a missing attribute for callers that bypass the parser in tests.

    Known types come from ``_VACCINE_TYPE_DISPATCH`` so adding a new
    vaccine writer is a single registration there — no second
    set/list to keep in sync.
    """
    vtype = getattr(args, 'vaccine_type', 'peptide') or 'peptide'
    if vtype not in _VACCINE_TYPE_DISPATCH:
        raise ValueError(
            "Unknown vaccine type %r. Known: %s." % (
                vtype, sorted(_VACCINE_TYPE_DISPATCH)))
    return vtype


def _emit_neoepitope_report_external(args, report_df, predictions):
    """LENS / pVACseq specific report path: writes the per-(peptide,
    allele) neoepitope CSV / XLSX. Independent from the modality
    dispatch — these are *report* outputs, not vaccine-design outputs.
    """
    if report_df is None or report_df.empty:
        return
    epitope_config = _epitope_config_from_args_safe(args)
    write_neoepitope_report(
        report_df=report_df,
        predictions=predictions,
        excel_report_path=(
            getattr(args, 'output_neoepitope_report', '') or None),
        csv_report_path=getattr(args, 'output_csv', '') or None,
        epitope_config=epitope_config,
    )
    if getattr(args, 'output_epitopes', ''):
        save_predictions(predictions, args.output_epitopes)


def _resolve_axis(args, per_type_attr, shared_attr, fallback):
    """Resolve a shared antigen-design axis with per-type override.

    Resolution: per-type explicit > shared > fallback. Per-type and
    shared values of None mean "unset" (defer to the next layer).
    """
    per = getattr(args, per_type_attr, None)
    if per is not None:
        return per
    shared = getattr(args, shared_attr, None)
    if shared is not None:
        return shared
    return fallback


def _emit_peptide_constructs(args, ranked):
    """Build + write peptide-vaccine constructs. Modality-specific."""
    # Resolve the orthogonal antigen-design axes:
    # per-type override > shared default > config default. The legacy
    # --peptide-mode is consumed by PeptideConstructConfig.__post_init__
    # when no explicit axis flag is set.
    antigen_content = _resolve_axis(
        args, 'peptide_antigen_content', 'antigen_content', None)
    epitopes_per_antigen = _resolve_axis(
        args, 'peptide_epitopes_per_antigen', 'epitopes_per_antigen', 1)
    config_kwargs = dict(
        mode=args.peptide_mode,
        linker=args.peptide_linker,
        min_antigen_length_aa=args.peptide_min_antigen_length_aa,
        max_antigen_length_aa=args.peptide_max_antigen_length_aa,
        antigens_per_construct=args.peptide_antigens_per_construct,
        max_constructs=args.peptide_max_constructs,
        candidates_per_slot=args.peptide_candidates_per_slot,
        n_terminal_acetylation=args.peptide_n_terminal_acetyl,
        c_terminal_amidation=args.peptide_c_terminal_amide,
        scale_mg=args.peptide_scale_mg,
        purity_percent=args.peptide_purity_percent,
        counterion=args.peptide_counterion,
        epitopes_per_antigen=epitopes_per_antigen,
    )
    if antigen_content is not None:
        config_kwargs['antigen_content'] = antigen_content
    peptide_options = PeptideConstructConfig(**config_kwargs)
    constructs = assemble_peptide_constructs(ranked, options=peptide_options)
    fasta_path = args.vaccine_output
    write_peptide_outputs(
        constructs,
        fasta_path=fasta_path,
        manifest_path=getattr(args, 'vaccine_manifest', '') or None,
        order_form_path=getattr(args, 'vaccine_order_form', '') or None,
        options=peptide_options,
    )
    logger.info(
        "Wrote %d peptide construct(s) to %s",
        len(constructs), fasta_path)


def _emit_mrna_constructs(args, ranked):
    """Build + write mRNA-vaccine constructs. Modality-specific."""
    junction_candidates = tuple(
        s.strip() for s in (args.mrna_junction_candidates or "").split(",")
        if s.strip()
    )
    antigen_content = _resolve_axis(
        args, 'mrna_antigen_content', 'antigen_content', 'mutation_spanning')
    epitopes_per_antigen = _resolve_axis(
        args, 'mrna_epitopes_per_antigen', 'epitopes_per_antigen', 1)
    options = RNAConstructConfig(
        signal_peptide=(args.mrna_signal_peptide or None),
        linker=args.mrna_linker,
        include_mitd=args.mrna_include_mitd,
        mitd=args.mrna_mitd,
        utr_5p=args.mrna_5p_utr,
        utr_3p=args.mrna_3p_utr,
        codon_species=args.mrna_codon_species,
        codon_method=args.mrna_codon_method,
        antigen_content=antigen_content,
        epitopes_per_antigen=epitopes_per_antigen,
        min_antigen_length_aa=args.mrna_min_antigen_length_aa,
        max_antigen_length_aa=args.mrna_max_antigen_length_aa,
        antigens_per_construct=args.mrna_antigens_per_construct,
        max_constructs=args.mrna_max_constructs,
        candidates_per_slot=args.mrna_candidates_per_slot,
        max_length_nt=args.mrna_max_length_nt,
        optimize_linkers=args.mrna_optimize_linkers,
        junction_swap_candidates=junction_candidates,
        junction_rank_strong=args.mrna_junction_rank_strong,
        junction_rank_mild=args.mrna_junction_rank_mild,
        poly_a_length=args.mrna_poly_a_length,
        poly_a_segmented=args.mrna_poly_a_segmented,
        poly_a_first_segment=args.mrna_poly_a_first_segment,
        poly_a_segment_linker=args.mrna_poly_a_segment_linker,
    )
    if options.optimize_linkers:
        try:
            mhc_predictor = mhc_binding_predictor_from_args(args)
            mhc_alleles = mhc_alleles_from_args(args)
        except Exception as e:
            logger.warning(
                "Could not load MHC predictor / alleles for "
                "per-junction linker optimization (%s). The optimizer "
                "will fall back to the shared linker at every junction. "
                "Set --mhc-predictor + --mhc-alleles to enable it.", e)
            logger.debug(
                "Predictor / alleles load traceback:", exc_info=True)
            mhc_predictor = None
            mhc_alleles = None
    else:
        mhc_predictor = None
        mhc_alleles = None
    constructs = assemble_mrna_constructs(
        ranked, options=options,
        mhc_predictor=mhc_predictor, mhc_alleles=mhc_alleles)
    output_dir = args.vaccine_output
    write_mrna_outputs(
        constructs,
        output_dir=output_dir,
        manifest_path=getattr(args, 'vaccine_manifest', '') or None,
        csv_path=getattr(args, 'vaccine_csv', '') or None,
        csv_include_full_rows=getattr(
            args, 'vaccine_csv_full_rows', True),
    )
    logger.info(
        "Wrote %d mRNA construct(s) to %s/{cds,no_polyA,full}.fasta",
        len(constructs), output_dir)


# Single dispatch table: vaccine-type → writer. One mode per run, so
# this is selected (not iterated) by ``_emit_outputs``. Adding a new
# vaccine type means adding the writer and registering it here, plus
# extending the ``--vaccine-type`` choices in arg_parser.
_VACCINE_TYPE_DISPATCH = {
    'peptide': _emit_peptide_constructs,
    'mrna': _emit_mrna_constructs,
}


def _emit_outputs(args, ranked, source):
    """Single fan-out point for everything downstream of ranking.

    ``ranked`` is the shared ``[(Variant, [VaccinePeptide])]``
    intermediate; ``source`` is 'pipeline' (VCF/BAM) or 'external'
    (LENS/pVACseq).

    The active vaccine type is single-valued (``--vaccine-type``);
    its writer fires only when ``--vaccine-output`` is set.
    Ranking-only / report-only runs are valid — the writer just
    no-ops without an output path.
    """
    vaccine_type = _resolve_vaccine_type(args)

    # Rank reports run for both sources; on the external-input path
    # the LENS-native per-(peptide, allele) report already consumed
    # ``--output-csv`` (via _emit_neoepitope_report_external) so the
    # rank-report writer would collide — skip it there.
    if (args.output_csv or args.output_xlsx_report) and source != 'external':
        make_csv_report(
            ranked,
            excel_report_path=args.output_xlsx_report,
            csv_report_path=args.output_csv)
    elif args.output_xlsx_report and source == 'external':
        logger.info(
            "--output-xlsx-report ignored on external-input path "
            "(use --output-neoepitope-report for the LENS / pVACseq "
            "XLSX format).")

    if args.output_neoepitope_report and source == 'pipeline':
        num_epitopes = getattr(args, 'num_epitopes_per_vaccine_peptide', None)
        make_minimal_neoepitope_report(
            ranked,
            num_epitopes_per_peptide=num_epitopes,
            excel_report_path=args.output_neoepitope_report)

    writer = _VACCINE_TYPE_DISPATCH.get(vaccine_type)
    fired = '(none)'
    if writer is None:
        logger.warning(
            "vaccine type %r recognized but has no writer "
            "registered (future-type stub); skipping.", vaccine_type)
    elif getattr(args, 'vaccine_output', ''):
        writer(args, ranked)
        fired = vaccine_type
    # No --vaccine-output → ranking-only for this type. Silent (vs.
    # warning) because users running only reports shouldn't be nagged.

    logger.info(
        "Vaccine-type dispatch [%s]: type=%s wrote=%s",
        source, vaccine_type, fired)


def configure_logging(args):
    logging.config.fileConfig(
        str(files('vaxrank').joinpath('logging.conf')),
        defaults={'logfilename': args.log_path})
    if getattr(args, 'verbose', False):
        for name in ['vaxrank', 'isovar', 'varcode', 'pyensembl', 'mhctools']:
            for handler in logging.getLogger(name).handlers:
                if isinstance(handler, logging.StreamHandler) \
                        and not isinstance(handler, logging.FileHandler):
                    handler.setLevel(logging.DEBUG)


def _resolve_ensembl_release(args):
    """If --ensembl-release was given, override args.genome with a specific
    pyensembl EnsemblRelease so that varcode/isovar use that release."""
    ensembl_release = getattr(args, 'ensembl_release', None)
    if ensembl_release is not None:
        import pyensembl
        args.genome = pyensembl.EnsemblRelease(ensembl_release)
        logger.info("Using Ensembl release %d", ensembl_release)


def main(args_list=None):
    """
    Rank personalized cancer neoantigens from somatic variants + tumor
    RNA (or a pre-computed LENS / pVACseq report) and emit the ranked
    candidates as one or more vaccine types plus optional analysis
    reports.

    Inputs (one of):
      - --vcf + --bam: full pipeline (variant calling → Isovar
        transcript assembly → MHC prediction → ranking)
      - --input-lens / --input-pvacseq: external neoepitope report,
        skipping the predict-from-genome stages

    Vaccine-type dispatch is single-valued (``--vaccine-type``;
    default ``peptide``); the writer fires only when
    ``--vaccine-output`` is also set.

    Example (full pipeline, peptide vaccine):

        vaxrank \\
            --vcf somatic.vcf \\
            --bam rnaseq.bam \\
            --mhc-predictor netmhc \\
            --mhc-alleles HLA-A*02:01 \\
            --vaccine-peptide-length 25 \\
            --vaccine-type peptide \\
            --output-pdf-report report.pdf \\
            --vaccine-output vaccine-peptides.fasta

    Example (LENS-driven mRNA design):

        vaxrank \\
            --input-lens patient.lens.tsv \\
            --vaccine-type mrna \\
            --vaccine-output mrna_dir/ \\
            --vaccine-csv layers.csv
    """
    if args_list is None:
        args_list = sys.argv[1:]

    args = parse_vaxrank_args(args_list)
    configure_logging(args)
    _resolve_ensembl_release(args)
    logger.info(args)

    # Fail fast when no output path is set, *before* loading inputs or
    # running predictions. Otherwise a long --input-lens run silently
    # produces nothing on disk (every output is opt-in via its own
    # --output-* flag, with no default destination). Applies to both
    # the VCF/BAM pipeline path and the external-input path.
    check_args(args)

    # Architecture (post-#252):
    #
    #   INPUT
    #     ├── --vcf + --bam  →  ranked_vaccine_peptides_with_metadata_from_parsed_args
    #     └── --input-lens / --input-pvacseq  →  external_input.load_external_ranked
    #                                            └── synthesize ranked from the report
    #
    #   SHARED:  ranked_variants_with_vaccine_peptides
    #
    #   DOWNSTREAM (modality-agnostic reports + peptide / mRNA construct dispatch)
    #     _emit_outputs(args, ranked, source)
    #
    # The external-input path no longer hard-short-circuits — it
    # produces the same shape as the VCF/BAM path so peptide and mRNA
    # vaccine designs work end-to-end from LENS / pVACseq files.
    from ..external_input import load_external_ranked

    patient_info = None
    args_for_report = args
    report_df = None
    predictions = None
    # ``data`` is the pipeline path's intermediate metadata bundle
    # (variants + patient_info + dna_vaf_by_variant). The external
    # path doesn't produce one, so default to an empty dict so the
    # shared template-report block can ``data.get(...)`` either way.
    data = {}

    if _has_external_input(args):
        loaded = load_external_ranked(args)
        if loaded is None:
            ranked_variants_with_vaccine_peptides = []
        else:
            (ranked_variants_with_vaccine_peptides, report_df,
             predictions, patient_info, external_dna_vaf) = loaded
            # Mirror the pipeline path's ``data['dna_vaf_by_variant']``
            # so the shared template-report block picks it up.
            data['dna_vaf_by_variant'] = external_dna_vaf
        if not ranked_variants_with_vaccine_peptides:
            logger.warning(
                "External input produced no ranked vaccine peptides; "
                "writing only the per-(peptide, allele) neoepitope "
                "report (if requested).")
        # Per-(peptide, allele) CSV / XLSX report is unique to the
        # external-input path; emit it before the shared dispatch.
        _emit_neoepitope_report_external(args, report_df, predictions)
        # Template reports want a dict here, not a Namespace —
        # ``vars(args)`` matches the pipeline path's
        # ``data['args']`` shape.
        args_for_report = vars(args)
        source = 'external'
    else:
        data = ranked_vaccine_peptides_with_metadata_from_parsed_args(args)
        ranked_variants_with_vaccine_peptides = data['variants']
        patient_info = data['patient_info']
        args_for_report = data['args']
        source = 'pipeline'

    # Issue #249: annotate epitope predictions with pepsickle proteasome-
    # cleavage credibility. Mutates EpitopePredictions in place; doesn't
    # affect ranking — purely additional information surfaced in the
    # per-epitope report tables. On by default; opt-out via
    # --no-processing-aware-annotation.
    if getattr(args, 'processing_aware_annotation', True):
        _annotate_predictions_with_processing(
            ranked_variants_with_vaccine_peptides, predictions,
            human_only=getattr(args, 'pepsickle_human_only', False),
            threshold=getattr(args, 'pepsickle_threshold', 0.5))

    _emit_outputs(args, ranked_variants_with_vaccine_peptides, source)

    ########################
    # Template-based reports (PDF / HTML / ASCII)
    ########################
    # Run on both pipeline and external (LENS / pVACseq) paths.
    # ``external_input.load_external_ranked`` resolves transcript IDs
    # to pyensembl ``Transcript`` objects and synthesizes a
    # ``PatientInfo`` with proxy variant counts so the template
    # builder has the same shape inputs either way. Antigens whose
    # transcript IDs can't be resolved (release mismatch, ERV, …)
    # render with "—" placeholders rather than crashing — see
    # ``TemplateDataCreator._effect_data``.
    if not (args.output_ascii_report or args.output_html_report or args.output_pdf_report):
        return

    template_data_creator = TemplateDataCreator(
        ranked_variants_with_vaccine_peptides=ranked_variants_with_vaccine_peptides,
        patient_info=patient_info,
        final_review=getattr(args, 'output_final_review', '') or '',
        reviewers=getattr(args, 'output_reviewed_by', '') or '',
        args_for_report=args_for_report,
        input_json_file=getattr(args, 'input_json_file', None),
        cosmic_vcf_filename=getattr(args, 'cosmic_vcf_filename', ''),
        dna_vaf_by_variant=data.get('dna_vaf_by_variant') or {})

    template_data = template_data_creator.compute_template_data()

    if args.output_ascii_report:
        make_ascii_report(
            template_data=template_data,
            ascii_report_path=args.output_ascii_report)

    if args.output_html_report:
        make_html_report(
            template_data=template_data,
            html_report_path=args.output_html_report)

    if args.output_pdf_report:
        make_pdf_report(
            template_data=template_data,
            pdf_report_path=args.output_pdf_report,
            backend=args.pdf_backend)


def run_vaxrank_from_parsed_args(args):
    merged_config = load_vaxrank_config(args)
    epitope_config = epitope_config_from_args(args, merged_config=merged_config)
    vaccine_config = vaccine_config_from_args(args, merged_config=merged_config)

    # Sync args with merged config values so downstream consumers (Isovar,
    # reports, etc.) see the effective configuration.
    args.vaccine_peptide_length = vaccine_config.preferred_peptide_length
    args.padding_around_mutation = vaccine_config.padding_around_mutation
    args.max_vaccine_peptides_per_variant = vaccine_config.max_vaccine_peptides_per_variant
    # Keep legacy key for backward compatibility in JSON/report args.
    args.max_vaccine_peptides_per_mutation = vaccine_config.max_vaccine_peptides_per_variant
    args.num_epitopes_per_vaccine_peptide = vaccine_config.num_mutant_epitopes_to_keep

    mhc_predictor = mhc_binding_predictor_from_args(args)

    prediction_cache = getattr(args, 'prediction_cache', None)
    if prediction_cache:
        from topiary import CachedPredictor
        cached = CachedPredictor.from_topiary_output(
            prediction_cache, fallback=mhc_predictor)
        logger.info("Loaded prediction cache from %s", prediction_cache)
        mhc_predictor = cached

    args.protein_sequence_length = (
            args.vaccine_peptide_length + 2 * args.padding_around_mutation
    )

    # Vaxrank is going to evaluate multiple vaccine peptides containing
    # the same mutation so need a longer sequence from Isovar.
    # We load variants ourselves (instead of run_isovar_from_parsed_args)
    # so we can filter out alt contigs that would crash downstream.
    variants = variant_collection_from_args(args)
    variants = _filter_unannotatable_variants(variants)
    isovar_results = run_isovar(
        variants=variants,
        alignment_file=alignment_file_from_args(args),
        read_collector=read_collector_from_args(args),
        protein_sequence_creator=protein_sequence_creator_from_args(args),
        filter_thresholds=filter_threshold_dict_from_args(args),
    )

    if args.output_isovar_csv:
        df = isovar_results_to_dataframe(isovar_results)
        df.to_csv(args.output_isovar_csv, index=False)

    vaxrank_results = run_vaxrank(
        isovar_results=isovar_results,
        mhc_predictor=mhc_predictor,
        vaccine_peptide_length=args.vaccine_peptide_length,
        max_vaccine_peptides_per_variant=args.max_vaccine_peptides_per_variant,
        num_mutant_epitopes_to_keep=args.num_epitopes_per_vaccine_peptide,
        epitope_config=epitope_config,
        vaccine_config=vaccine_config,
        allow_dna_only_fallback=getattr(args, 'allow_dna_only_fallback', False),
    )

    if getattr(args, 'output_epitopes', ''):
        # Collect all epitope predictions across all variants
        all_predictions = []
        for _variant, peptides in vaxrank_results.ranked_vaccine_peptides:
            for vp in peptides:
                all_predictions.extend(vp.mutant_epitope_predictions)
        save_predictions(all_predictions, args.output_epitopes)

    return vaxrank_results

def ranked_vaccine_peptides_with_metadata_from_parsed_args(args):
    """
    Computes all the data needed for report generation.

    Parameters
    ----------
    args : Namespace
      Parsed user args from this run

    Returns a dictionary containing 3 items:
    - ranked variant/vaccine peptide list
    - a dictionary of command-line arguments used to generate it
    - patient info object
    """

    if hasattr(args, 'input_json_file'):
        with open(args.input_json_file) as f:

            data = serializable.from_json(f.read())
            # the JSON data from the previous run will have the older args saved, which may need to
            # be overridden with args from this run (which all be output related)
            data['args'].update(vars(args))

            # if we need to truncate the variant list based on max_mutations_in_report, do that here
            if args.max_mutations_in_report is not None and len(data['variants']) > args.max_mutations_in_report:
                data['variants'] = data['variants'][:args.max_mutations_in_report]
            return data
    # get various things from user args
    mhc_alleles = mhc_alleles_from_args(args)
    logger.info("MHC alleles: %s", mhc_alleles)

    all_variants = variant_collection_from_args(args)
    variants = _filter_unannotatable_variants(all_variants)
    logger.info("Variants: %d loaded, %d after filtering invalid contigs",
                len(all_variants), len(variants))

    dna_vaf_by_variant = extract_dna_vaf_by_variant(
        all_variants, tumor_sample_name=getattr(args, 'tumor_sample_name', None))

    vaxrank_results = run_vaxrank_from_parsed_args(args)

    variants_count_dict = vaxrank_results.variant_counts()
    assert len(variants) == variants_count_dict['num_total_variants'], \
        "Len(variants) is %d but variants_count_dict came back with %d" % (
            len(variants), variants_count_dict['num_total_variants'])

    if args.output_passing_variants_csv:
        variant_metadata_dicts = vaxrank_results.variant_properties(
            gene_pathway_check=GenePathwayCheck(),
            dna_vaf_by_variant=dna_vaf_by_variant)
        df = pd.DataFrame(variant_metadata_dicts)
        df.to_csv(args.output_passing_variants_csv, index=False)

    ranked_variants_with_vaccine_peptides = vaxrank_results.ranked_vaccine_peptides
    ranked_variants_with_vaccine_peptides_for_report = \
        ranked_variants_with_vaccine_peptides[:args.max_mutations_in_report]
    patient_info = PatientInfo(
        patient_id=args.output_patient_id,
        vcf_paths=variants.sources,
        bam_path=args.bam,
        mhc_alleles=mhc_alleles,
        num_somatic_variants=variants_count_dict['num_total_variants'],
        num_coding_effect_variants=variants_count_dict['num_coding_effect_variants'],
        num_variants_with_rna_support=variants_count_dict['num_variants_with_rna_support'],
        num_variants_with_vaccine_peptides=variants_count_dict['num_variants_with_vaccine_peptides']
    )
    # return variants, patient info, and command-line args
    data = {
        # TODO:
        #  change this field to 'ranked_variants_with_vaccine_peptides'
        #  but figure out how to do it in a backwards compatible way
        'variants': ranked_variants_with_vaccine_peptides_for_report,
        'patient_info': patient_info,
        'args': vars(args),
        'dna_vaf_by_variant': dna_vaf_by_variant,
    }
    logger.info('About to save args: %s', data['args'])

    # save JSON data if necessary. as of time of writing, vaxrank takes ~25 min to run,
    # most of which is core logic. the formatting is super fast, and it can
    # be useful to save the data to be able to iterate just on the formatting
    if args.output_json_file:
        with open(args.output_json_file, 'w') as f:
            f.write(serializable.to_json(data))
            logger.info('Wrote JSON report data to %s', args.output_json_file)

    return data
