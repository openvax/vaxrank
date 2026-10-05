"""Evidence-preserving inputs for the existing epitope scoring pipeline.

Topiary owns the table and its encoding. Native Vaxrank objects supplement it;
they never replace source rows or reconstruct missing biological evidence.
"""

from dataclasses import dataclass, field, replace
from pathlib import Path

import msgspec
import pandas as pd
from mhctools.pred import Prediction
from topiary import TopiaryResult, combine_sources, is_stated, read_csv, read_tsv

from .allele_validation import validate_peptide_alleles
from .candidate_epitope import CandidateEpitope, Peptide, stated_or_blank
from .epitope_config import EpitopeConfig
from .epitope_dsl import epitopes_to_topiary_df
from .input_scope import InputProvenance
from .native_serialization import from_native_json, to_native_json
from .vaccine_antigen import VaccineAntigen


DATASET_SCHEMA = "vaxrank.epitope_dataset.v1"
DATASET_METADATA = "vaxrank_epitope_dataset"
ANTIGEN_COLUMN = "vaxrank_antigen_json"


def read_table(path):
    """Read a normalized table using Topiary's typed CSV/TSV codec."""
    return (read_tsv if str(path).endswith(".tsv") else read_csv)(path)


def _value(row, name, default=None):
    value = row.get(name)
    return value if is_stated(value) else default


@dataclass
class EpitopeDataset:
    """Historical evidence plus the native objects that consume it.

    Candidate IDs refer to ``prediction_id`` in the scoring view. Antigens
    are optional, explicitly supplied construction evidence keyed by that ID.
    A missing antigen leaves the candidate available for table-only ranking.
    """

    result: TopiaryResult
    epitopes: tuple[CandidateEpitope, ...] = ()
    provenance: tuple[InputProvenance, ...] = ()
    config: EpitopeConfig | None = None
    antigens: dict[str, VaccineAntigen] = field(default_factory=dict)
    selection: dict = field(default_factory=dict)
    mutation_fragments: dict = field(default_factory=dict)
    construction_reports: list = field(default_factory=list)
    direct_sources: list = field(default_factory=list)
    policy_evaluations: list = field(default_factory=list)
    native_references: object = field(default=None, repr=False, compare=False)

    def add_mutation(self, fragment, epitopes):
        """Capture direct predictions before vaccine windows discard context."""
        import hashlib
        provenance = self.provenance[0]
        context_id = hashlib.sha256((provenance.input_id + to_native_json(fragment)).encode()).hexdigest()
        antigen = VaccineAntigen.from_mutant_protein_fragment(fragment)
        candidates = []
        for epitope in epitopes:
            identity = 'direct:' + context_id + ':' + str(epitope.offset) + ':' + epitope.sequence
            candidates.append(replace(epitope, prediction_id=identity,
                                      input_provenance=provenance))
            self.mutation_fragments[identity] = fragment
            self.antigens[identity] = antigen
        self.epitopes += tuple(candidates)

    @classmethod
    def from_predictions(cls, epitopes, *, result=None, **kwargs):
        epitopes = tuple(epitopes)
        if result is None:
            frame = epitopes_to_topiary_df(epitopes)
            if frame.empty:
                # Native payloads can exist without prediction leaves. Keep
                # an empty, readable evidence table without inventing rows.
                frame = frame.reindex(columns=[
                    "prediction_id", "peptide", "peptide_offset", "allele", "kind"])
            result = TopiaryResult(frame, form="long")
        provenance = {}
        for epitope in epitopes:
            if epitope.input_provenance is not None:
                p = epitope.input_provenance
                provenance[p.input_id] = p
        return cls(result, epitopes, tuple(provenance.values()), **kwargs)

    @classmethod
    def from_topiary(cls, result, *, label="input", sample_name=None):
        """Adapt reported observations; no inferred intervals or admission."""
        if isinstance(result, pd.DataFrame):
            result = TopiaryResult(result)
        if "source_observation_id" not in result.df:
            result = combine_sources({label: result}, sample_name=sample_name)
        # Pandas may represent missing strings as NaN (notably with its
        # string dtype). Native objects use None, while '' stays distinct.
        # Normalize only this consumer view; retain the original table.
        frame = result.long_df.astype(object)
        frame = frame.where(frame.notna(), None)
        epitopes, antigens = [], {}
        for identity, rows in frame.groupby("source_observation_id", sort=False):
            row = rows.iloc[0]
            peptide = _value(row, "peptide")
            if peptide is None:
                continue  # ORF/RNA-only observations remain in the evidence.
            offset = _value(row, "peptide_offset", 0)
            if int(offset) != offset or offset < 0:
                raise ValueError("peptide_offset must be a nonnegative integer")
            source = _value(row, "source_sequence", "")
            if source and source[int(offset):int(offset) + len(peptide)] != peptide:
                raise ValueError("Peptide does not match its stated source_sequence")
            predictions, wt_predictions = [], []
            for record in rows.to_dict("records"):
                kind = _value(record, "kind")
                if kind is None:
                    continue
                method = stated_or_blank(record.get("prediction_method_name"))
                version = stated_or_blank(record.get("predictor_version"))
                # The typed table is the lossless authority for partial
                # measurements. A percentile-only quantitative prediction
                # cannot be represented by mhctools.Prediction (#566). Keep
                # its occurrence, allele and rank in result; do not fabricate
                # a quantitative value or score for a native leaf.
                if any(_value(record, name) is not None for name in ('value', 'score')):
                    predictions.append(Prediction(
                        kind=kind, peptide=peptide,
                        allele=stated_or_blank(record.get("candidate_allele")),
                        predictor_name=method, predictor_version=version,
                        value=_value(record, "value"), score=_value(record, "score"),
                        percentile_rank=_value(record, "percentile_rank"),
                        n_flank=record.get("n_flank"), c_flank=record.get("c_flank")))
                if any(_value(record, 'wt_' + name) is not None
                       for name in ('value', 'score')):
                    wt_predictions.append(Prediction(
                        kind=kind, peptide=_value(record, 'wt_peptide', ''),
                        allele=stated_or_blank(record.get('candidate_allele')),
                        predictor_name=_value(record, 'wt_prediction_method_name',
                                              method),
                        predictor_version=_value(record, 'wt_predictor_version',
                                                 version),
                        value=_value(record, 'wt_value'), score=_value(record, 'wt_score'),
                        percentile_rank=_value(record, 'wt_percentile_rank'),
                        n_flank=record.get('wt_n_flank'), c_flank=record.get('wt_c_flank')))
            wt_sequence = _value(row, 'wt_peptide', '')
            comparators = ({'wt': Peptide(sequence=wt_sequence, predictions=tuple(wt_predictions),
                                         n_flank=row.get('wt_n_flank'), c_flank=row.get('wt_c_flank'))}
                           if (wt_sequence or wt_predictions) and wt_sequence != peptide else {})
            candidate = CandidateEpitope(
                sequence=peptide, source_sequence=source, offset=int(offset),
                source_name=_value(row, "source_sequence_name", ""),
                n_flank=row.get("n_flank"), c_flank=row.get("c_flank"),
                prediction_id=identity, predictions=tuple(predictions), comparators=comparators,
                patient_alleles=tuple(sorted({stated_or_blank(value)
                    for value in rows.candidate_allele if stated_or_blank(value)})),
                source_class=_value(row, "source_class"),
                overlaps_targetable=False)
            validate_peptide_alleles(candidate, f"Topiary observation {identity}")
            payload = _value(row, ANTIGEN_COLUMN)
            if payload is not None:
                antigen = from_native_json(payload, VaccineAntigen)
                if antigen.amino_acids[int(offset):int(offset) + len(peptide)] != peptide:
                    raise ValueError("Peptide does not match its supplied VaccineAntigen")
                antigens[identity] = antigen
                candidate = replace(
                    candidate, source_sequence=antigen.amino_acids,
                    overlaps_targetable=antigen.targetable_mask.overlaps(
                        int(offset), int(offset) + len(peptide)))
            epitopes.append(candidate)
        return cls(result, tuple(epitopes), antigens=antigens)

    def scoring_frame(self):
        """Make a consumer view without changing original IDs or annotations."""
        frame = self.result.long_df.copy()
        if "source_observation_id" in frame:
            original = frame.get("prediction_id", pd.Series(None, index=frame.index, dtype=object))
            frame["prediction_id"] = frame["source_observation_id"].where(
                frame["source_observation_id"].notna(), original)
        if "peptide_offset" not in frame:
            frame["peptide_offset"] = 0
        else:
            frame["peptide_offset"] = frame["peptide_offset"].fillna(0)
        # The candidate representation uses canonical blank allele identities.
        if "allele" in frame:
            frame["allele"] = frame["allele"].map(stated_or_blank)
            if "candidate_allele" in frame:
                frame["allele"] = frame["candidate_allele"].where(
                    frame["source_observation_id"].notna(), frame["allele"])
        return frame

    def scoring_frames(self):
        """Retain source-local model defaults when several inputs are saved."""
        frame = self.scoring_frame()
        groups = self.selection.get("score_groups")
        if not groups:
            return [frame]
        return [frame.loc[frame.prediction_id.isin(group)].copy() for group in groups]

    def select_representatives(self, duplicates=None):
        """Ask Topiary to select observations from the already scored evidence."""
        from topiary import rank_candidates
        if 'candidate_id' not in self.result.df:
            return None
        if self.config is not None and self.config.selection_policy is not None:
            from topiary import select_policy_representatives
            from .selection_policy import combine_evaluations
            evaluation = combine_evaluations(self.policy_evaluations, self.config.selection_policy,
                                             duplicates=duplicates)
            if evaluation is None:
                raise ValueError('Score named policies before selecting representative observations')
            ranked = select_policy_representatives(evaluation)
            self.policy_evaluations = [evaluation]
            self.selection['duplicates'] = evaluation.policy.duplicates
            self.selection['representatives'] = [dict(
                source_observation_id=row.source_observation_id,
                candidate_id=row.candidate_id, candidate_allele=row.allele,
                candidate_observations=row.alternative_occurrences,
                representative_reason=row.representative_reason)
                for row in ranked.itertuples()]
            return set(zip(ranked.source_observation_id, ranked.allele))
        frame = self.result.long_df
        frame = frame.loc[frame.source_observation_id.notna()].copy()
        scores = {(e.prediction_group_source, allele): score
                  for e in self.epitopes for allele, score in e.per_allele_scores.items()}
        frame['vaxrank_selection_score'] = [scores.get(key) for key in frame[[
            'source_observation_id', 'candidate_allele']].itertuples(index=False, name=None)]
        policy = duplicates or self.selection.get('duplicates', 'error')
        ranked = rank_candidates(
            TopiaryResult(frame, metadata=self.result.metadata),
            'vaxrank_selection_score', duplicates=policy)
        self.selection['duplicates'] = policy
        self.selection['representatives'] = ranked[[
            'source_observation_id', 'candidate_id', 'candidate_allele',
            'candidate_observations']].to_dict('records')
        return set(zip(ranked.source_observation_id, ranked.candidate_allele))

    def save(self, path):
        """Persist evidence with Topiary and objects with the native codec."""
        from .selection_policy import encode_evaluations
        extra = dict(self.result.extra)
        extra[DATASET_METADATA] = {
            "schema": DATASET_SCHEMA,
            "epitopes": [to_native_json(e) for e in self.epitopes],
            "provenance": [to_native_json(p) for p in self.provenance],
            "config": msgspec.to_builtins(self.config) if self.config is not None else None,
            "antigens": {key: to_native_json(a) for key, a in self.antigens.items()},
            "selection": self.selection,
            "mutation_fragments": {key: to_native_json(fragment)
                                   for key, fragment in self.mutation_fragments.items()},
            "construction_reports": self.construction_reports,
            "direct_sources": self.direct_sources,
            "policy_evidence": encode_evaluations(self.policy_evaluations),
        }
        from .native_references import NativeReferences
        extra[DATASET_METADATA] = NativeReferences(path, self.native_references).transform(
            extra[DATASET_METADATA])
        result = TopiaryResult(self.result.df.copy(), metadata=self.result.metadata)
        result.extra = extra
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        (result.to_tsv if str(path).endswith(".tsv") else result.to_csv)(path)

    @classmethod
    def load(cls, path):
        from .selection_policy import decode_evaluations
        from .mutant_protein_fragment import MutantProteinFragment
        result = read_table(path)
        payload = result.extra.pop(DATASET_METADATA, None)
        if not isinstance(payload, dict) or payload.get("schema") != DATASET_SCHEMA:
            raise ValueError("Expected a native Vaxrank epitope dataset")
        from .native_references import NativeReferences
        references = NativeReferences(path)
        payload = references.transform(payload, decode=True)
        epitopes = tuple(from_native_json(e, CandidateEpitope) for e in payload["epitopes"])
        for epitope in epitopes:
            validate_peptide_alleles(epitope, str(path))
        antigens = {key: from_native_json(a, VaccineAntigen)
                    for key, a in payload["antigens"].items()}
        fragments = {key: from_native_json(value, MutantProteinFragment)
                     for key, value in payload.get("mutation_fragments", {}).items()}
        # Older datasets already persisted these IDs on the construction
        # antigen. Reuse that recorded provenance only for the same transcripts;
        # do not query an absent annotation to recreate it during replay.
        for key, fragment in fragments.items():
            antigen = antigens.get(key)
            transcript_ids = tuple(sorted({t.id for t in fragment.supporting_reference_transcripts}))
            if (fragment.supporting_reference_protein_ids is None and antigen is not None
                    and antigen.transcript_ids == transcript_ids):
                fragments[key] = replace(fragment, supporting_reference_protein_ids=antigen.protein_ids)
        return cls(
            result=result, epitopes=epitopes,
            provenance=tuple(from_native_json(p, InputProvenance) for p in payload["provenance"]),
            config=msgspec.convert(payload["config"], EpitopeConfig) if payload["config"] else None,
            antigens=antigens,
            selection=payload["selection"],
            mutation_fragments=fragments,
            construction_reports=payload.get("construction_reports", []),
            direct_sources=payload.get("direct_sources", []),
            policy_evaluations=decode_evaluations(payload.get("policy_evidence", [])),
            native_references=references)

    def index_native_references(self):
        """Prepare bundled local annotation indexes without acquiring data."""
        if self.native_references is not None:
            self.native_references.index()

    def report_frame(self):
        """Report each occurrence/allele once; retain full evidence for scoring."""
        frame = self.scoring_frame()
        if frame.empty:
            return frame
        frame = frame.loc[frame.prediction_id.isin(
            e.prediction_group_source for e in self.epitopes)].copy()
        alleles = {e.prediction_group_source: e.patient_alleles or ('',)
                   for e in self.epitopes}
        # Shared allele-free processing rows are evidence for these candidates,
        # not extra peptide/allele report observations. The full scoring frame
        # and native result keep every measurement kind and predictor row.
        frame = frame.loc[[row.allele in alleles[row.prediction_id]
                           for row in frame.itertuples()]]
        frame = frame.drop_duplicates(['prediction_id', 'peptide', 'peptide_offset', 'allele'])
        frame["Prediction identity"] = frame.prediction_id
        frame["Peptide offset"] = frame.peptide_offset
        frame["Mutant peptide sequence"] = frame.peptide
        frame["Allele"] = frame.allele
        from .external_report import ExternalRecord
        adapter_ids = {from_native_json(r, ExternalRecord).key.identifier
                       for saved in self.construction_reports for r in saved['records']}
        def limitation(key):
            if key in adapter_ids:
                return "source adapter evaluates construction evidence"
            antigen = self.antigens.get(key)
            if antigen is None:
                return "missing antigen construction evidence"
            if not antigen.tumor_specificity.admits_construct:
                return "antigen held out by tumor-specificity evidence"
            return ""
        frame["Construction limitation"] = frame.prediction_id.map(limitation)
        frame.attrs = {"topiary_df": self.scoring_frame()}
        return frame


def read_epitope_report(path, *, source_format, scope_declarations=None, **_):
    """Format adapters produce the same report consumed by existing inputs."""
    from .epitope_io import load_predictions
    from .external_report import ExternalReport

    if source_format == "epitopes":
        with open(path) as stream:
            enriched = stream.readline().startswith("#")
        dataset = (EpitopeDataset.load(path) if enriched else
                   EpitopeDataset.from_predictions(load_predictions(path)))
        if not enriched and any(e.per_allele_scores for e in dataset.epitopes):
            # Older native files retained results but not the expression or
            # annotations that produced them. Expose those results to the
            # same DSL; an explicit new config can still score predictions.
            scores = {(e.prediction_group_source, e.sequence, e.offset, allele): score
                      for e in dataset.epitopes for allele, score in e.per_allele_scores.items()}
            frame = dataset.result.df
            frame['vaxrank_saved_score'] = [scores.get(key) for key in frame[[
                'prediction_id', 'peptide', 'peptide_offset', 'allele'
            ]].itertuples(index=False, name=None)]
            dataset.config = EpitopeConfig(
                score_expr='vaxrank_saved_score',
                filter_expr='vaxrank_saved_score == vaxrank_saved_score')
    else:
        import hashlib
        import json
        from dataclasses import asdict
        from .input_scope import combine_scopes, report_declarations, scope_from_mapping
        identity = "topiary:" + hashlib.sha256(Path(path).read_bytes()).hexdigest()
        result = read_table(path)
        producer_scope, _ = report_declarations(path, source_format, frame=result.df)
        scope = combine_scopes(producer_scope, scope_from_mapping(
            scope_declarations or {}, str(path)), str(path))
        label = identity + ':' + hashlib.sha256(
            json.dumps(asdict(scope), sort_keys=True).encode()).hexdigest()
        dataset = EpitopeDataset.from_topiary(
            result, label=label, sample_name=scope.patient_id or label)
    return ExternalReport(source_format, str(path),
                          report_df=dataset.report_frame(),
                          epitopes=dataset.epitopes, dataset=dataset)


def dataset_ranking_result(report, epitopes, genome=None, options=None):
    """Use explicit antigen evidence with the existing construct consumer."""
    from .external_input import (
        ExternalRankingAccumulator, ExternalVariantEntry,
        source_agnostic_construct_options,
    )
    from .vaccine_peptide import VaccinePeptide

    dataset = report.dataset
    from .core_logic import vaccine_peptides_from_epitopes
    from .external_report import ExternalReport, ExternalRecord
    from .external_input import lens_ranking_result, pvacseq_ranking_result

    source_options = source_agnostic_construct_options(options)
    accumulator = ExternalRankingAccumulator(
        require_target_epitopes=options.require_target_epitopes_in_variant)
    representatives = dataset.selection.get('representatives')
    selected = (set((r['source_observation_id'], r['candidate_allele']) for r in representatives)
                if representatives is not None else dataset.select_representatives())
    # Saved source records re-enter their existing occurrence-selection
    # adapters, retaining the original IDs and evidence without a file reread.
    source_ids = set()
    for saved in dataset.construction_reports:
        records = tuple(from_native_json(r, ExternalRecord) for r in saved['records'])
        ids = {r.key.identifier for r in records}
        source_ids.update(ids)
        report_source = ExternalReport(
            saved['source_format'], saved['path'], records=records,
            rows=tuple(from_native_json(saved['rows'], list)))
        ranker = {'lens': lens_ranking_result, 'pvacseq': pvacseq_ranking_result}[saved['source_format']]
        source_genome = genome
        if source_genome is None and saved.get('genome') is not None:
            from pyensembl import Genome, EnsemblRelease
            source_genome = from_native_json(saved['genome'], (Genome, EnsemblRelease))
        outcome = ranker(report_source, [e for e in epitopes if e.prediction_group_source in ids],
                         genome=source_genome, options=options)
        for entry in outcome.entries:
            accumulator.add(entry)
    from varcode import Variant, StructuralVariant
    for source in dataset.direct_sources:
        properties = source['properties']
        accumulator.add(ExternalVariantEntry(
            variant=from_native_json(source['variant'], (Variant, StructuralVariant)),
            resolved_protein_context=properties['is_coding_nonsynonymous'],
            has_rna_support=properties['rna_support'], dna_vaf=properties['dna_vaf']))
    mutation_groups = {}
    groups = {}
    generalized_ids = (set(dataset.result.df.source_observation_id.dropna())
                       if 'source_observation_id' in dataset.result.df else set())
    for epitope in epitopes:
        key = epitope.prediction_group_source
        if key in source_ids:
            continue
        if selected is not None and key in generalized_ids:
            scores = {a: s for a, s in epitope.per_allele_scores.items() if (key, a) in selected}
            if not scores:
                continue
            epitope = replace(epitope, per_allele_scores=scores)
        fragment = dataset.mutation_fragments.get(key)
        if fragment is not None:
            identity = (epitope.input_provenance.input_id if epitope.input_provenance else '',
                        to_native_json(fragment))
            mutation_groups.setdefault(identity, (fragment, []))[1].append(epitope)
            continue
        antigen = dataset.antigens.get(key)
        if antigen is None or not antigen.tumor_specificity.admits_construct:
            continue
        identity = to_native_json(antigen)
        groups.setdefault(identity, (antigen, []))[1].append(epitope)
    for fragment, candidates in mutation_groups.values():
        vaccines = vaccine_peptides_from_epitopes(
            fragment.variant, fragment, candidates, epitope_config=dataset.config,
            vaccine_config=options.vaccine_config,
            vaccine_peptide_length=options.vaccine_peptide_length,
            num_target_epitopes_to_keep=(options.num_target_epitopes_to_keep
                                        if options.num_target_epitopes_to_keep is not None else 1000),
            manufacturability_config=options.manufacturability_config)
        accumulator.add(ExternalVariantEntry(
            variant=fragment.variant, vaccine_peptide=vaccines[0] if vaccines else None,
            vaccine_peptides=tuple(vaccines),
            resolved_protein_context=True, has_rna_support=bool(fragment.n_rna_alt),
            dna_vaf=fragment.dna_vaf))
    options = source_options
    for antigen, candidates in groups.values():
        vaccine = VaccinePeptide(
            antigen=antigen, epitopes=candidates,
            num_target_epitopes_to_keep=options.num_target_epitopes_to_keep,
            combined_score_expr=options.combined_score_expr,
            ranking_rules=options.ranking_rules,
            manufacturability_thresholds=options.manufacturability_thresholds,
            manufacturability_rules=options.manufacturability_rules)
        accumulator.add(ExternalVariantEntry(
            source=antigen,
            vaccine_peptide=vaccine, resolved_protein_context=True))
    return accumulator.result(report.source_format)
