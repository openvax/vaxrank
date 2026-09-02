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

from collections import OrderedDict
from dataclasses import dataclass
from importlib import import_module
import logging
import os
import sys
import tempfile

from astropy.io import ascii as asc
import jinja2
import pandas as pd
import roman
from mhctools.pred import Prediction
from varcode import load_vcf_fast

from .cancer_hotspots import get_hotspot_url
from .manufacturability import ManufacturabilityScores
from .processing import PEPSICKLE_PREDICTOR_NAME
from .varcode_effects import (
    OUTCOME_SELECTION_MULTI_OUTCOME,
    is_multi_outcome_effect,
    summarize_varcode_effect_outcomes,
)

logger = logging.getLogger(__name__)


JINJA_ENVIRONMENT = jinja2.Environment(
    loader=jinja2.FileSystemLoader(os.path.dirname(__file__)),
    extensions=[],
    autoescape=False,
    trim_blocks=True,
    lstrip_blocks=True,
)


@dataclass(frozen=True)
class EpitopeReportRowInput:
    """Prediction evidence and patient allele used for one report row."""

    prediction: Prediction
    allele: str


def epitope_report_row_inputs(epitope):
    """Return explicit row inputs for a candidate's template table.

    Select one anchor per ``(patient allele, predictor, version)``. Affinity
    is preferred, followed by presentation, followed by allele-independent
    processing for combinations not covered by either allele-specific kind.
    This preserves every scored allele and predictor in mixed-evidence inputs
    without duplicating the canonical processing Prediction.
    """
    predictions = epitope.predictions_flat()
    row_inputs_by_key = {}
    for kind in ('pMHC_affinity', 'pMHC_presentation'):
        for prediction in predictions:
            if prediction.kind != kind:
                continue
            if not prediction.allele:
                # An allele-scoped kind that arrived blank is malformed —
                # see vaxrank.allele_evidence for the one definition of
                # which kinds carry an allele, and epitope.allele_attributions
                # for which alleles the allele-free kinds were credited to.
                raise ValueError(
                    f"{kind} report evidence requires a patient allele")
            key = (
                prediction.allele,
                prediction.predictor_name,
                prediction.predictor_version,
            )
            row_inputs_by_key.setdefault(
                key,
                EpitopeReportRowInput(prediction, prediction.allele),
            )

    processing_predictions = tuple(
        prediction for prediction in predictions
        if prediction.kind == 'antigen_processing')
    if processing_predictions and not epitope.patient_alleles:
        raise ValueError(
            "Cannot render allele-independent antigen-processing evidence "
            "without explicit patient alleles")
    for prediction in processing_predictions:
        for allele in epitope.patient_alleles:
            key = (
                allele,
                prediction.predictor_name,
                prediction.predictor_version,
            )
            row_inputs_by_key.setdefault(
                key,
                EpitopeReportRowInput(prediction, allele),
            )
    return tuple(row_inputs_by_key.values())


class TemplateDataCreator(object):
    def __init__(
            self,
            ranked_variants_with_vaccine_peptides,
            patient_info,
            final_review,
            reviewers,
            args_for_report,
            input_json_file,
            cosmic_vcf_filename=None,
            dna_vaf_by_variant=None,
            source_vaf_by_variant=None,
            processing_predictions_by_key=None,
            mrna_ranking_decisions=None,
            vaccine_constructions=None,
            target_alleles=None,
            include_manufacturability=None):
        """
        Construct a TemplateDataCreator object, from the output of the vaxrank pipeline.

        ``processing_predictions_by_key`` is the map returned by
        :func:`vaxrank.processing.annotate_processing` —
        ``(peptide, source, predictor_name) -> ProcessingPrediction``.
        Report writers join this in at render time. ``None`` means
        the processing-aware annotation pass didn't run; the
        per-epitope tables then omit the processing columns.
        """
        self.ranked_variants_with_vaccine_peptides = ranked_variants_with_vaccine_peptides
        self.patient_info = patient_info
        self.dna_vaf_by_variant = dna_vaf_by_variant or {}
        self.source_vaf_by_variant = source_vaf_by_variant or {}
        self.processing_predictions_by_key = (
            processing_predictions_by_key or {})
        # Issue #270: mRNA ranking-decisions section. Legacy field;
        # superseded by ``vaccine_constructions['mrna']`` in the new
        # report layout (#269). Both reach the template so existing
        # tests / consumers keep working.
        self.mrna_ranking_decisions = mrna_ranking_decisions
        # Issue #269: per-modality "Vaccine construction" sections —
        # coverage block + selected antigens with allele tiers
        # + dropped-past-cap list. Keys are modality names; values
        # are the dicts produced by
        # :func:`vaxrank.coverage.summarize_construction_decisions`.
        self.vaccine_constructions = vaccine_constructions or {}
        self.target_alleles = list(target_alleles or [])

        # ``vaccine_type`` is rendered in the patient-info block so
        # the reader can tell at a glance what's being designed.
        # Multi-valued (default ``['peptide']``); render as a
        # comma-list so multi-mode runs (peptide+mrna) read cleanly.
        raw_types = args_for_report.get('vaccine_type') or ['peptide']
        if isinstance(raw_types, str):
            raw_types = [raw_types]
        self.vaccine_type = ', '.join(raw_types)

        # Show the *effective* run parameters, not the raw argparse
        # Namespace. We drop: output-path args (covered by the Inputs /
        # file listing), the argparse ``version`` SUPPRESS sentinel,
        # internal config-plumbing keys, and any value that is unset
        # (None / ''). A None here means "not set / inherits a default"
        # — noise in a record of what actually ran, and the source of
        # the confusing ``manufacturability: None`` / ``version:
        # ==SUPPRESS==`` lines.
        _internal_keys = {
            'version', 'config', 'config_set_overrides',
            'config_expr_overrides',
        }
        args_to_display_in_report = {
            k: v for k, v in args_for_report.items()
            if not k.startswith("output")
            and k not in _internal_keys
            and v is not None
            and v != ''
            and v != '==SUPPRESS=='
        }

        self.template_data = {
            'args': sorted(args_to_display_in_report.items()),
            'reviewers': reviewers.split(',') if reviewers else [],
            'final_review': final_review,
            'input_json_file': input_json_file,
            # these report sections are optional. ``include_manufacturability``
            # can be overridden per report (the split-report layout turns it
            # off for the modality-agnostic core report and the mRNA report,
            # on for the peptide report); when None it follows the run's
            # ``--manufacturability`` resolution.
            'include_manufacturability': (
                args_for_report['manufacturability']
                if include_manufacturability is None
                else include_manufacturability),
            'include_wt_epitopes': args_for_report['wt_epitopes'],
        }

        # map from peptide objects to their COSMIC IDs if they exist
        if cosmic_vcf_filename:
            logger.info('Loading COSMIC data...')
            self.cosmic_variant_collection = load_vcf_fast(
                cosmic_vcf_filename, allow_extended_nucleotides=True, include_info=False)
            logger.info('COSMIC data loaded.')
        else:
            self.cosmic_variant_collection = None

    def _patient_info(self):
        """
        Returns an OrderedDict with patient info.

        Inputs section: prefer the explicit ``inputs`` list when
        populated (LENS / pVACseq external paths set this directly so
        the file is labelled accurately); fall back to the legacy
        ``vcf_paths`` + ``bam_path`` fields for the VCF/BAM pipeline
        path and for backward-compatible loading of older
        serialized PatientInfo JSON.
        """
        patient_info = OrderedDict()
        patient_info['Patient ID'] = self.patient_info.patient_id
        # Modality is part of the report's identity — what the user
        # is actually designing. Generic up here; modality-specific
        # construct details (signal peptide, UTRs, etc. for mRNA;
        # synthesis flags for peptide) belong in a per-construct
        # section rendered next to the assembled FASTA / directory.
        patient_info['Vaccine type'] = self.vaccine_type

        if self.patient_info.inputs:
            for label, path in self.patient_info.inputs:
                patient_info[label] = path
        else:
            if self.patient_info.vcf_paths:
                patient_info['VCF (somatic variants) path(s)'] = (
                    '; '.join(self.patient_info.vcf_paths))
            if self.patient_info.bam_path:
                patient_info['BAM (RNAseq reads) path'] = (
                    self.patient_info.bam_path)

        patient_info['MHC alleles'] = ' '.join(self.patient_info.mhc_alleles)
        # Name the processing-credibility predictor explicitly so the
        # ``Processing: …`` columns in the per-VP epitope tables can
        # use predictor-agnostic header names without losing
        # provenance. Today there's only one (pepsickle, see
        # ``vaxrank.processing``); a future per-position predictor
        # plugged in via the same flat record fields would
        # show its own name here.
        patient_info['Processing predictor'] = PEPSICKLE_PREDICTOR_NAME
        patient_info['Total number of somatic variants'] = (
            self.patient_info.num_somatic_variants)
        patient_info['Somatic variants with predicted coding effects'] = (
            self.patient_info.num_coding_effect_variants)
        patient_info[
            'Somatic variants with predicted coding effects and RNA support'
        ] = self.patient_info.num_variants_with_rna_support
        patient_info[
            'Somatic variants with predicted coding effects, RNA support and '
            'predicted MHC ligands'
        ] = self.patient_info.num_variants_with_vaccine_peptides
        return patient_info

    def _variant_data(self, variant, top_vaccine_peptide):
        """
        Returns an OrderedDict with info used to populate variant info section.
        """
        mutant_protein_fragment = top_vaccine_peptide.mutant_protein_fragment
        top_score = _sanitize(top_vaccine_peptide.combined_score)
        igv_locus = "chr%s:%d" % (variant.contig, variant.start)
        rna_vaf = mutant_protein_fragment.rna_vaf
        dna_vaf = self.dna_vaf_by_variant.get(variant)
        source_vaf = self.source_vaf_by_variant.get(variant)
        variant_data = OrderedDict([
            ('IGV locus', igv_locus),
            ('Gene name', mutant_protein_fragment.gene_name),
            ('Top score', top_score),
            ('DNA VAF', '%.3f' % dna_vaf if dna_vaf is not None else 'n/a'),
            ('RNA VAF', '%.3f' % rna_vaf if rna_vaf is not None else 'n/a'),
            ('RNA reads supporting variant allele', mutant_protein_fragment.n_alt_reads),
            ('RNA reads supporting reference allele', mutant_protein_fragment.n_ref_reads),
        ])
        # Only when the reference count was counted rather than derived.
        # Where a source reports a depth and a fraction, ref is depth minus
        # alt, so this subtraction is structurally zero — and the whole
        # depth when no fraction was reported at all. Printing either as
        # "reads supporting other alleles" describes the locus wrongly.
        if mutant_protein_fragment.n_other_reads is not None:
            variant_data['RNA reads supporting other alleles'] = (
                mutant_protein_fragment.n_other_reads)
        # A variant fraction the source did not qualify. Shown under its
        # own label rather than as DNA VAF, which is what LENS's bare
        # `vaf` used to be printed as — asserting an assay the file never
        # stated.
        if source_vaf is not None:
            variant_data['Variant fraction (assay unstated)'] = (
                '%.3f' % source_vaf)
        # Say how those counts were obtained. "51 reads" from an alignment
        # and "51" from depth x VAF are the same integer and not the same
        # claim, and the combined-score DSL weights them identically.
        if mutant_protein_fragment.rna_evidence_method:
            variant_data['RNA read count method'] = (
                mutant_protein_fragment.rna_evidence_method)
        # What those counts counted. Five fragments and five reads are
        # different bars, and the rows above say neither on their own.
        if mutant_protein_fragment.rna_evidence_subject:
            variant_data['RNA read count subject'] = (
                mutant_protein_fragment.rna_evidence_subject)
        if mutant_protein_fragment.sequence_source:
            variant_data['Protein sequence source'] = (
                mutant_protein_fragment.sequence_source)
        return variant_data

    def effect_data(self, predicted_effect, selected_effect=None):
        """OrderedDict with info about the given varcode effect.

        ``predicted_effect`` may be ``None`` on external-input paths
        when the transcript IDs in the source file (LENS / pVACseq)
        couldn't be resolved against the configured pyensembl release.
        Render those rows with ``"—"`` placeholders rather than
        crash.
        """
        display_effect = selected_effect or predicted_effect
        if display_effect is None:
            return OrderedDict([
                ('Effect type', '—'),
                ('Transcript name', '—'),
                ('Transcript ID', '—'),
                ('Effect description', '—'),
            ])
        effect_data = OrderedDict([
            ('Effect type', display_effect.__class__.__name__),
            ('Transcript name', display_effect.transcript_name or '—'),
            ('Transcript ID', display_effect.transcript_id or '—'),
            ('Effect description', display_effect.short_description),
        ])
        if is_multi_outcome_effect(predicted_effect):
            effect_data.update(summarize_varcode_effect_outcomes(
                predicted_effect))
        return effect_data

    def _peptide_header_display_data(self, vaccine_peptide, rank):
        """
        Returns a dictionary with info used to populate the header section of a peptide table.

        Parameters
        ----------
        vaccine_peptide : VaccinePeptide
          The given peptide to convert to display form

        rank : int
          Rank of vaccine peptide in list
        """
        mutant_protein_fragment = vaccine_peptide.mutant_protein_fragment
        amino_acids = mutant_protein_fragment.amino_acids
        mutation_start = mutant_protein_fragment.mutant_amino_acid_start_offset
        mutation_end = mutant_protein_fragment.mutant_amino_acid_end_offset
        aa_before_mutation = amino_acids[:mutation_start]
        aa_mutant = amino_acids[mutation_start:mutation_end]
        aa_after_mutation = amino_acids[mutation_end:]

        header_display_data = {
            'num': roman.toRoman(rank + 1).lower(),
            'aa_before_mutation': aa_before_mutation,
            'aa_mutant': aa_mutant,
            'aa_after_mutation': aa_after_mutation,
        }
        return header_display_data

    def _peptide_data(self, vaccine_peptide, transcript_name):
        """
        Returns a dictionary with info used to populate peptide table contents.

        Parameters
        ----------
        vaccine_peptide : VaccinePeptide
          The given peptide to convert to display form

        transcript_name : str
            RNA transcript name (should match that displayed in effect section)
        """
        mutant_protein_fragment = vaccine_peptide.mutant_protein_fragment
        amino_acids = mutant_protein_fragment.amino_acids
        peptide_data = OrderedDict([
            ('Transcript name', transcript_name),
            ('Length', len(amino_acids)),
            ('Expression score', _sanitize(vaccine_peptide.expression_score)),
            ('Mutant epitope score', _sanitize(vaccine_peptide.target_epitope_score)),
            ('Combined score', _sanitize(vaccine_peptide.combined_score)),
            ('Max coding sequence coverage',
                mutant_protein_fragment.n_alt_reads_supporting_protein_sequence),
            ('Mutant amino acids', mutant_protein_fragment.n_mutant_amino_acids),
            ('Mutation distance from edge',
                mutant_protein_fragment.mutation_distance_from_edge),
        ])
        return peptide_data

    def _manufacturability_data(self, vaccine_peptide):
        """
        Returns an OrderedDict with manufacturability data for the given peptide.
        """
        scores = vaccine_peptide.manufacturability_scores
        manufacturability_data = OrderedDict([
            ('C-terminal 7mer GRAVY score', _sanitize(scores.cterm_7mer_gravy_score)),
            ('Max 7mer GRAVY score', _sanitize(scores.max_7mer_gravy_score)),
            ('N-terminal Glutamine, Glutamic Acid, or Cysteine',
                int(scores.difficult_n_terminal_residue)),
            ('C-terminal Cysteine', int(scores.c_terminal_cysteine)),
            ('C-terminal Proline', int(scores.c_terminal_proline)),
            ('Total number of Cysteine residues', scores.cysteine_count),
            ('N-terminal Asparagine', int(scores.n_terminal_asparagine)),
            ('Number of Aspartate-Proline bonds', scores.aspartate_proline_bond_count),
        ])
        return manufacturability_data

    def _processing_prediction_for(self, epitope, prediction):
        """Look up the ProcessingPrediction for this epitope+prediction
        by ``(peptide, source, predictor_name)``. Returns ``None`` when
        no record exists (annotation pass didn't run, or this
        prediction had no usable source sequence).

        Predictor name is currently always
        :data:`vaxrank.processing.PEPSICKLE_PREDICTOR_NAME`; when a
        second per-position cleavage predictor (NetChop, PAProC, …)
        lands, this method becomes the swap point — read the
        predictor name from a config knob and fall back as needed.
        """
        key = (
            epitope.sequence or '',
            epitope.source_sequence or '',
            PEPSICKLE_PREDICTOR_NAME,
        )
        return self.processing_predictions_by_key.get(key)

    def _wt_ic50_for_allele(self, epitope, allele, predictor=None):
        """Return the WT pMHC_affinity IC50 for ``allele`` (optionally
        scoped to ``predictor``) from this epitope's WT comparator,
        or ``None`` when no WT record exists (e.g. WT peptide was
        too short, or peptide doesn't overlap the mutation).

        Iterates ``predictions_flat()`` rather than
        ``predictions_for('pMHC_affinity')`` so multi-predictor WT
        comparators don't trip the disambiguation check.
        """
        if epitope.wt is None:
            return None
        for p in epitope.wt.predictions_flat():
            if p.kind != 'pMHC_affinity':
                continue
            if p.allele != allele:
                continue
            if predictor is not None and p.predictor_name != predictor:
                continue
            if p.value is not None:
                return p.value
        return None

    def epitope_data(self, epitope, prediction, include_processing=False,
                     include_additional_prediction_axes=False, allele=None):
        """Returns an OrderedDict with epitope data for one
        (CandidateEpitope, mutant Prediction) row.

        One mutant ``CandidateEpitope`` carries N per-allele × per-predictor
        ``mhctools.Prediction`` records. The report keeps one row per explicit
        :func:`epitope_report_row_inputs` result, using affinity as the
        preferred anchor and presentation or processing when affinity is
        absent.

        ``allele`` defaults to ``prediction.allele``. Report rows anchored by
        canonical allele-independent processing evidence pass an explicitly
        retained patient allele here.

        ``include_processing``: when True, always emit the three
        proteasomal-cleavage credibility columns
        (``Processing: C-term`` / ``Processing: max internal`` /
        ``Processing: combined``), using ``'—'`` as the placeholder
        for predictions that weren't annotated. The caller decides
        per-list (per-VaccinePeptide) whether *any* prediction in
        the list has a ProcessingPrediction; if so, every row gets
        the columns so the rendered table has consistent headers
        and column widths.

        When False (default), the legacy 7-column dict is returned
        unchanged — reports that don't run pepsickle keep their
        original shape.

        Processing scores come from
        ``self.processing_predictions_by_key`` (built by
        :func:`vaxrank.processing.annotate_processing`), looked up
        by ``(peptide, source, predictor_name)``. Pre-2.23 these
        lived as ``pepsickle_*`` attributes on the flat record
        itself; the join is the new contract (#272).

        ``include_additional_prediction_axes`` appends presentation
        score/percentile and integrated antigen-processing score leaves
        from the same predictor/version. Missing axes render as ``'—'``
        so mixed-predictor tables retain consistent columns.
        """
        row_allele = prediction.allele if allele is None else allele
        if not row_allele:
            raise ValueError(
                "Template epitope rows require an explicit patient allele")
        wt_ic50 = self._wt_ic50_for_allele(
            epitope, row_allele, predictor=prediction.predictor_name)
        wt_ic50_str = _format_ic50(wt_ic50)
        ic50_str = _format_ic50(prediction.value)
        wt_peptide_sequence = (
            epitope.wt.sequence if epitope.wt is not None else '')
        epitope_data = OrderedDict([
            ('Sequence', epitope.sequence),
            # ``Predictor`` names which MHC affinity tool produced
            # ``IC50`` / ``Score`` / ``WT IC50`` for this row —
            # mhcflurry / netmhcpan / mhcflurry-presentation / etc.
            # Pipeline runs are usually single-predictor; LENS files
            # are sometimes multi-predictor (e.g. mhcflurry +
            # netmhcpan), so making the source explicit per row is
            # the only honest answer.
            ('Predictor', prediction.predictor_name or '—'),
            ('IC50', ic50_str),
            # ``Score`` is the per-allele DSL score for this epitope
            # as computed at predict time by the configured
            # ``EpitopeConfig.score_expr`` (default: logistic of
            # IC50, range [0, 1]; higher = stronger). Read directly
            # from ``epitope.per_allele_scores`` — the single source
            # of truth — rather than recomputing from the raw IC50.
            ('Score',
                _sanitize(epitope.per_allele_scores[row_allele])
                if row_allele in epitope.per_allele_scores else '—'),
            ('Allele', row_allele.replace('HLA-', '')),
            ('WT sequence', wt_peptide_sequence),
            ('WT IC50', wt_ic50_str),
        ])
        if include_additional_prediction_axes:
            presentation = None
            antigen_processing = None
            for related in epitope.predictions_flat():
                if (
                    related.predictor_name != prediction.predictor_name
                    or related.predictor_version != prediction.predictor_version
                ):
                    continue
                if (
                    related.kind == 'pMHC_presentation'
                    and related.allele == row_allele
                ):
                    presentation = related
                elif (
                    related.kind == 'antigen_processing'
                    and not related.allele
                ):
                    antigen_processing = related
            epitope_data['Presentation score'] = (
                '%.3f' % presentation.score
                if presentation is not None else '—')
            epitope_data['Presentation %ile'] = (
                '%.3f' % presentation.percentile_rank
                if (
                    presentation is not None
                    and presentation.percentile_rank is not None
                ) else '—')
            epitope_data['Integrated processing score'] = (
                '%.3f' % antigen_processing.score
                if antigen_processing is not None else '—')
        if include_processing:
            # Column headers are predictor-agnostic ("Processing: …")
            # so a future per-position predictor (NetChop, PAProC)
            # can render through the same table without renaming.
            # The active predictor is named in the patient-info
            # header (``Processing predictor: pepsickle``) so the
            # reader can trace which model produced the values.
            pp = self._processing_prediction_for(epitope, prediction)
            if pp is None:
                c_term = max_int = proc = None
            else:
                c_term = pp.c_term_cleavage_prob
                max_int = pp.max_internal_cut_prob
                proc = pp.processing_score
            epitope_data['Processing: C-term'] = (
                '%.2f' % c_term if c_term is not None else '—')
            epitope_data['Processing: max internal'] = (
                '%.2f' % max_int if max_int is not None else '—')
            epitope_data['Processing: combined'] = (
                '%.2f' % proc if proc is not None else '—')
        return epitope_data

    def _query_cosmic(self, variant):
        if not self.cosmic_variant_collection:
            return None
        if variant in self.cosmic_variant_collection.metadata:
            # IDs in the DB are of the form 'COSM725245'
            cosmic_id = self.cosmic_variant_collection.metadata[variant]['id']
            link_for_report = \
                "http://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=%s" % cosmic_id[4:]
            logger.info("Link for report: %s", link_for_report)
            return link_for_report

        logger.info("Variant not in COSMIC")
        return None

    def _query_cancer_hotspots(self, predicted_effect, gene_name):
        """
        Check if variant is a known cancer hotspot and return link if so.

        Uses bundled cancer hotspots data from Chang et al. 2016/2017.
        Returns ``None`` when ``predicted_effect`` is None (no
        transcript context — common on external-input paths).
        """
        if predicted_effect is None:
            return None
        amino_acids = predicted_effect.short_description
        hotspot_url = get_hotspot_url(gene_name, amino_acids)
        if hotspot_url:
            logger.info("Cancer hotspot found: %s %s -> %s", gene_name, amino_acids, hotspot_url)
        return hotspot_url

    def _databases(self, variant, predicted_effect, gene_name):
        databases = {}
        hotspot_link = self._query_cancer_hotspots(predicted_effect, gene_name)
        if hotspot_link:
            databases['Cancer Hotspots'] = hotspot_link

        cosmic_link = self._query_cosmic(variant)
        if cosmic_link:
            databases['COSMIC'] = cosmic_link

        return databases

    def compute_template_data(self):
        patient_info = self._patient_info()

        # list ranked variants with their peptides
        variants = []
        num = 0
        for (variant, vaccine_peptides) in self.ranked_variants_with_vaccine_peptides:
            variant_short_description = variant.short_description
            if len(vaccine_peptides) == 0:
                try:
                    gene_names = variant.gene_names
                except ValueError:
                    gene_names = ["(unknown — invalid contig)"]
                logger.info(
                    "Skipping gene(s) %s, variant %s: no vaccine peptides",
                    gene_names, variant_short_description)
                continue
            num += 1

            top_peptide = vaccine_peptides[0]
            variant_data = self._variant_data(variant, top_peptide)
            mutant_protein_fragment = top_peptide.mutant_protein_fragment
            predicted_effect = mutant_protein_fragment.predicted_effect()
            predicted_effect_outcomes = mutant_protein_fragment.predicted_effect(
                outcome_selection=OUTCOME_SELECTION_MULTI_OUTCOME)
            effect_data = self.effect_data(
                predicted_effect_outcomes,
                selected_effect=predicted_effect)

            databases = self._databases(
                variant, predicted_effect, mutant_protein_fragment.gene_name)

            peptides = []
            for j, vaccine_peptide in enumerate(vaccine_peptides):
                if not vaccine_peptide.contains_target_epitopes():
                    logger.info('No epitopes for peptide: %s', vaccine_peptide)
                    continue

                header_display_data = self._peptide_header_display_data(vaccine_peptide, j)
                transcript_name = (
                    predicted_effect.transcript_name
                    if predicted_effect is not None else '—')
                peptide_data = self._peptide_data(vaccine_peptide, transcript_name)
                manufacturability_data = self._manufacturability_data(vaccine_peptide)

                # Issue #249 / #272: processing-credibility columns
                # are added to every row in this VP's per-epitope
                # table iff *any* mutant epitope in the list has a
                # ProcessingPrediction in the map. Keeps the rendered
                # HTML/PDF/ASCII table headers consistent with the
                # rows even when some predictions failed annotation
                # (per-source pepsickle errors are graceful, so
                # mixed annotated/unannotated lists are possible).
                #
                # Pepsickle keys off (peptide, source) — independent
                # of allele — so we look up once per CandidateEpitope using
                # the first available leaf Prediction as a probe.
                def _has_processing(e):
                    for p in e.predictions_flat():
                        if self._processing_prediction_for(e, p) is not None:
                            return True
                    return False

                any_processing = any(
                    _has_processing(e)
                    for e in vaccine_peptide.target_epitopes)
                any_additional_prediction_axes = any(
                    prediction.kind in {
                        'pMHC_presentation', 'antigen_processing'}
                    for epitope in vaccine_peptide.target_epitopes
                    for prediction in epitope.predictions_flat())

                epitopes = []
                wt_epitopes = []
                # One row per explicit report input. Affinity remains the
                # preferred anchor; presentation-only and processing-only
                # candidates still render rather than disappearing.
                for e in vaccine_peptide.target_epitopes:
                    for row_input in epitope_report_row_inputs(e):
                        epitopes.append(self.epitope_data(
                            e, row_input.prediction,
                            include_processing=any_processing,
                            include_additional_prediction_axes=(
                                any_additional_prediction_axes),
                            allele=row_input.allele))

                for e in vaccine_peptide.self_epitopes:
                    for p in e.predictions_flat():
                        if p.kind != 'pMHC_affinity':
                            continue
                        epitope_data = self.epitope_data(e, p)
                        key_list = ['Allele', 'IC50', 'Sequence']
                        wt_epitopes.append(
                            {key: epitope_data[key] for key in key_list})

                # hack: make a nicely-formatted fixed width table for epitopes, used in ASCII report
                with tempfile.TemporaryFile(mode='r+') as temp:
                    asc.write(epitopes, temp, format='fixed_width_two_line', delimiter_pad=' ')
                    temp.seek(0)
                    ascii_epitopes = temp.read()

                ascii_wt_epitopes = None
                if len(wt_epitopes) > 0:
                    with tempfile.TemporaryFile(mode='r+') as temp:
                        asc.write(
                            wt_epitopes, temp, format='fixed_width_two_line', delimiter_pad=' ')
                        temp.seek(0)
                        ascii_wt_epitopes = temp.read()
                peptide_dict = {
                    'header_display_data': header_display_data,
                    'peptide_data': peptide_data,
                    'manufacturability_data': manufacturability_data,
                    'epitopes': epitopes,
                    'ascii_epitopes': ascii_epitopes,
                    'wt_epitopes': wt_epitopes,
                    'ascii_wt_epitopes': ascii_wt_epitopes,
                }
                peptides.append(peptide_dict)

            # if there are no peptides for this variant, exclude from report
            if len(peptides) == 0:
                logger.info('No peptides for variant: %s', variant)
                continue

            variant_dict = {
                'num': num,
                'short_description': variant_short_description,
                'variant_data': variant_data,
                'effect_data': effect_data,
                'peptides': peptides,
                'databases': databases,
            }
            variants.append(variant_dict)

        # add package metadata to the report
        package_versions = {}
        for name in ['vaxrank', 'isovar', 'mhctools', 'varcode', 'pyensembl']:
            module = import_module(name)
            version = getattr(module, '__version__')
            package_versions[name] = version

        self.template_data.update({
            'patient_info': patient_info,
            'variants': variants,
            'package_versions': package_versions,
            # Issue #270: legacy mRNA-only ranking block. Kept for
            # back-compat with older templates / tests; new
            # consumers read ``vaccine_constructions['mrna']``.
            'mrna_ranking': self.mrna_ranking_decisions,
            # Issue #269: per-modality "Vaccine construction"
            # blocks. Keyed by modality name (``'peptide'`` /
            # ``'mrna'`` / future ``'dna'``). Empty dict when no
            # ranked antigens exist; the template renders blocks
            # only for keys present here.
            'vaccine_constructions': self.vaccine_constructions,
            'target_alleles': self.target_alleles,
        })
        return self.template_data


def render_report(
        template_data,
        file_handle,
        template_path):
    """Render one packaged Jinja template into an open text stream."""
    template = JINJA_ENVIRONMENT.get_template(template_path)
    report = template.render(template_data)
    file_handle.write(report)

def make_ascii_report(
        template_data,
        ascii_report_path):
    with open(ascii_report_path, "w") as f:
        render_report(template_data, f, 'templates/template.txt')
    logger.info('Wrote ASCII report to %s', ascii_report_path)

def make_html_report(
        template_data,
        html_report_path):
    with open(html_report_path, "w") as f:
        render_report(template_data, f, 'templates/template.html')
    logger.info('Wrote HTML report to %s', html_report_path)

def _pdf_via_pdfkit(html_path, pdf_report_path):
    import pdfkit

    options = {
        'zoom': 0.55,
        'margin-top': '20mm'
    }

    if sys.platform in ('linux', 'linux2'):
        # pdfkit uses wkhtmltopdf, which doesn't work on headless servers;
        # recommended workaround is to use xvfb, as documented here:
        # https://github.com/wkhtmltopdf/wkhtmltopdf/issues/2037#issuecomment-62019521
        from xvfbwrapper import Xvfb
        logger.info('Running pdfkit inside xvfb wrapper')
        with Xvfb():
            pdfkit.from_file(html_path, pdf_report_path, options=options)
    else:
        pdfkit.from_file(html_path, pdf_report_path, options=options)


def _pdf_via_weasyprint(html_path, pdf_report_path):
    import weasyprint
    doc = weasyprint.HTML(filename=html_path)
    doc.write_pdf(pdf_report_path)


def make_pdf_report(
        template_data,
        pdf_report_path,
        backend='pdfkit'):
    with tempfile.NamedTemporaryFile(mode='w', suffix='.html') as f:
        render_report(template_data, f, 'templates/template.html')
        f.flush()

        if backend == 'pdfkit':
            _pdf_via_pdfkit(f.name, pdf_report_path)
        elif backend == 'weasyprint':
            _pdf_via_weasyprint(f.name, pdf_report_path)
        else:
            raise ValueError("Unknown PDF backend: %s (choose 'pdfkit' or 'weasyprint')" % backend)
    logger.info('Wrote PDF report to %s (backend=%s)', pdf_report_path, backend)

def new_columns():
    columns = OrderedDict([
        ("amino_acids", []),
        ("chr", []),
        ("pos", []),
        ("ref", []),
        ("alt", []),
        ("gene_name", []),
        ("variant_rank", []),
        ("peptide_rank", []),
        ("mutation_start", []),
        ("mutation_end", []),
        ("combined_score", []),
        ("expression_score", []),
        ("target_epitope_score", []),
    ])
    for field in ManufacturabilityScores._fields:
        columns[field] = []
    return columns

def _str_sig_figs(input, n_sig_figs):
    return '{:g}'.format(float('{:.{p}g}'.format(input, p=n_sig_figs)))

def _sanitize(val):
    """
    Converts values into display-friendly
    """
    if isinstance(val, bool):
        return int(val)
    else:
        return _str_sig_figs(val, 5)

def _format_ic50(value):
    try:
        if pd.isna(value):
            return 'No prediction'
    except (TypeError, ValueError):
        pass
    return '%.2f nM' % value

def resize_columns(worksheet, amino_acids_col, pos_col):
    """
    Resizes amino acid and mutant position columns in the Excel sheet so that they don't
    have to be expanded.
    """
    worksheet.column_dimensions[amino_acids_col].width = 40
    worksheet.column_dimensions[pos_col].width = 12


def make_minimal_neoepitope_report(
        ranked_variants_with_vaccine_peptides,
        num_epitopes_per_peptide=None,
        excel_report_path=None):
    """
    Creates a simple Excel spreadsheet containing one neoepitope per row

    Parameters
    ----------
    ranked_variants_with_vaccine_peptides :
      Ranked list of (variant, list of its vaccine peptides)

    num_epitopes_per_peptide : int
      The number of epitopes to include for each vaccine peptide; these are sorted before cutoff.
      If None, all epitopes will be included in the output

    excel_report_path : str
      Path to which to write the output Excel file
    """
    rows = []
    # each row in the spreadsheet is one (peptide, allele) record
    for (variant, vaccine_peptides) in ranked_variants_with_vaccine_peptides:
        for vaccine_peptide in vaccine_peptides:
            epitopes = vaccine_peptide.target_epitopes
            if num_epitopes_per_peptide is not None:
                epitopes = epitopes[:num_epitopes_per_peptide]
            for epitope in epitopes:
                wt_peptide_sequence = (
                    epitope.wt.sequence
                    if epitope.wt is not None else '')
                affinity_leaves = [
                    p for p in epitope.predictions_flat()
                    if p.kind == 'pMHC_affinity']
                for p in affinity_leaves:
                    wt_ic50 = None
                    if epitope.wt is not None:
                        for wt_p in epitope.wt.predictions_flat():
                            if (wt_p.kind == 'pMHC_affinity'
                                    and wt_p.allele == p.allele
                                    and wt_p.predictor_name == p.predictor_name
                                    and wt_p.value is not None):
                                wt_ic50 = wt_p.value
                                break
                    # Build the shared core columns via the same helper
                    # the LENS / pVACseq report builders use, so all
                    # three input paths emit identical core columns and
                    # the same missing-value convention. The RNA-read
                    # count is this path's only source-specific column.
                    from .epitope_io import neoepitope_core_row
                    row = neoepitope_core_row(
                        allele=p.allele,
                        mutant_peptide=epitope.sequence,
                        mutant_affinity=p.value,
                        wt_peptide=wt_peptide_sequence,
                        wt_affinity=wt_ic50,
                        gene_name=(
                            vaccine_peptide.mutant_protein_fragment.gene_name),
                        variant=variant.short_description,
                        score=vaccine_peptide.target_epitope_score)
                    row['Variant allele RNA read count'] = (
                        vaccine_peptide.mutant_protein_fragment.n_alt_reads)
                    rows.append(row)

    if len(rows) > 0:
        df = pd.DataFrame.from_dict(rows)
        writer = pd.ExcelWriter(excel_report_path, engine='openpyxl')
        df.to_excel(writer, sheet_name='Neoepitopes', index=False)

        # resize columns to be not crappy
        worksheet = writer.sheets['Neoepitopes']
        worksheet.column_dimensions['B'].width = 23
        worksheet.column_dimensions['D'].width = 27
        worksheet.column_dimensions['E'].width = 26
        worksheet.column_dimensions['F'].width = 17
        worksheet.column_dimensions['G'].width = 30
        worksheet.column_dimensions['H'].width = 9
        worksheet.column_dimensions['I'].width = 18
        writer.close()
        logger.info('Wrote XLSX neoepitope report file to %s', excel_report_path)


def make_csv_report(
        ranked_variants_with_vaccine_peptides,
        excel_report_path=None,
        csv_report_path=None):
    """
    Writes out CSV/XLSX reports as needed.
    """
    # make a bunch of pd frames, save them in an OrderedDict with keys being descriptive Excel
    # sheet names (will be used later for making the Excel report if needed)
    frames = OrderedDict()
    for i, (variant, vaccine_peptides) in enumerate(ranked_variants_with_vaccine_peptides):
        any_vaccine_peptides = False
        if not vaccine_peptides:
            continue

        sheet_name = '%d_%s_chr%s_%d_%s_%s' % (
            i + 1, vaccine_peptides[0].mutant_protein_fragment.gene_name,
            variant.contig, variant.start, variant.ref, variant.alt)
        columns = new_columns()
        for j, vaccine_peptide in enumerate(vaccine_peptides):

            # if there are no predicted epitopes, exclude this peptide from the report
            if not vaccine_peptide.contains_target_epitopes():
                logger.info('No epitopes for peptide: %s', vaccine_peptide)
                continue

            columns["chr"].append(variant.contig)
            columns["pos"].append(variant.original_start)
            columns["ref"].append(variant.original_ref)
            columns["alt"].append(variant.original_alt)
            columns["gene_name"].append(vaccine_peptide.mutant_protein_fragment.gene_name)
            columns["variant_rank"].append(i + 1)
            columns["peptide_rank"].append(j + 1)
            columns["amino_acids"].append(vaccine_peptide.mutant_protein_fragment.amino_acids)
            columns["mutation_start"].append(
                vaccine_peptide.mutant_protein_fragment.mutant_amino_acid_start_offset)
            columns["mutation_end"].append(
                vaccine_peptide.mutant_protein_fragment.mutant_amino_acid_end_offset)
            columns["combined_score"].append(_sanitize(vaccine_peptide.combined_score))
            columns["expression_score"].append(_sanitize(vaccine_peptide.expression_score))
            columns["target_epitope_score"].append(_sanitize(vaccine_peptide.target_epitope_score))
            for field in ManufacturabilityScores._fields:
                columns[field].append(
                    _sanitize(getattr(vaccine_peptide.manufacturability_scores, field)))
            any_vaccine_peptides = True

        if not any_vaccine_peptides:
            continue

        df = pd.DataFrame(columns, columns=columns.keys())
        frames[sheet_name] = df

    if not frames:
        logger.info('No data for CSV or XLSX report')
        return

    all_dfs = pd.concat(frames.values())
    if csv_report_path:
        all_dfs.to_csv(csv_report_path, index=False)
        logger.info('Wrote CSV report file to %s', csv_report_path)

    if excel_report_path:
        writer = pd.ExcelWriter(excel_report_path, engine='openpyxl')

        # copy the variant rank column to position 0, make first sheet called "All"
        all_dfs[''] = all_dfs['variant_rank']
        colnames = all_dfs.columns.tolist()
        colnames.insert(0, colnames.pop(colnames.index('')))
        all_dfs = all_dfs.reindex(columns=colnames)
        all_dfs.to_excel(writer, sheet_name='All', index=False)
        resize_columns(writer.sheets['All'], 'B', 'D')

        # add one sheet per variant
        for sheet_name, df in frames.items():
            # trim sheet names to 31 characters due to limit in Excel
            # should still be unique since they start with the variant
            # index
            shortened_sheet_name = sheet_name[:31]
            df.to_excel(writer, sheet_name=shortened_sheet_name, index=False)
            resize_columns(writer.sheets[shortened_sheet_name], 'A', 'C')

        writer.close()
        logger.info('Wrote manufacturability XLSX file to %s', excel_report_path)
