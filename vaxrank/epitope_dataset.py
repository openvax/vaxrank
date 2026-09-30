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

    @classmethod
    def from_predictions(cls, epitopes, *, result=None, **kwargs):
        epitopes = tuple(epitopes)
        if result is None:
            result = TopiaryResult(epitopes_to_topiary_df(epitopes))
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
        frame = result.long_df
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
                predictions.append(Prediction(
                    kind=kind, peptide=peptide,
                    allele=stated_or_blank(record.get("candidate_allele")),
                    predictor_name=stated_or_blank(record.get("prediction_method_name")),
                    predictor_version=stated_or_blank(record.get("predictor_version")),
                    value=_value(record, "value"), score=_value(record, "score"),
                    percentile_rank=_value(record, "percentile_rank"),
                    n_flank=record.get("n_flank"), c_flank=record.get("c_flank")))
                if any(_value(record, 'wt_' + name) is not None
                       for name in ('value', 'score', 'percentile_rank')):
                    wt_predictions.append(Prediction(
                        kind=kind, peptide=_value(record, 'wt_peptide', ''),
                        allele=stated_or_blank(record.get('candidate_allele')),
                        predictor_name=_value(record, 'wt_prediction_method_name',
                                              predictions[-1].predictor_name),
                        predictor_version=_value(record, 'wt_predictor_version',
                                                 predictions[-1].predictor_version),
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
                patient_alleles=tuple(sorted({p.allele for p in predictions if p.allele})),
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
        extra = dict(self.result.extra)
        extra[DATASET_METADATA] = {
            "schema": DATASET_SCHEMA,
            "epitopes": [to_native_json(e) for e in self.epitopes],
            "provenance": [to_native_json(p) for p in self.provenance],
            "config": msgspec.to_builtins(self.config) if self.config is not None else None,
            "antigens": {key: to_native_json(a) for key, a in self.antigens.items()},
            "selection": self.selection,
        }
        result = TopiaryResult(self.result.df.copy(), metadata=self.result.metadata)
        result.extra = extra
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        (result.to_tsv if str(path).endswith(".tsv") else result.to_csv)(path)

    @classmethod
    def load(cls, path):
        result = read_table(path)
        payload = result.extra.pop(DATASET_METADATA, None)
        if not isinstance(payload, dict) or payload.get("schema") != DATASET_SCHEMA:
            raise ValueError("Expected a native Vaxrank epitope dataset")
        epitopes = tuple(from_native_json(e, CandidateEpitope) for e in payload["epitopes"])
        for epitope in epitopes:
            validate_peptide_alleles(epitope, str(path))
        return cls(
            result=result, epitopes=epitopes,
            provenance=tuple(from_native_json(p, InputProvenance) for p in payload["provenance"]),
            config=msgspec.convert(payload["config"], EpitopeConfig) if payload["config"] else None,
            antigens={key: from_native_json(a, VaccineAntigen)
                      for key, a in payload["antigens"].items()},
            selection=payload["selection"])

    def report_frame(self):
        """Expose original annotations beside the familiar report columns."""
        frame = self.scoring_frame()
        if frame.empty:
            return frame
        frame = frame.loc[frame.prediction_id.isin(
            e.prediction_group_source for e in self.epitopes)].copy()
        frame["Prediction identity"] = frame.prediction_id
        frame["Peptide offset"] = frame.peptide_offset
        frame["Mutant peptide sequence"] = frame.peptide
        frame["Allele"] = frame.allele
        def limitation(key):
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
    options = source_agnostic_construct_options(options)
    accumulator = ExternalRankingAccumulator(
        require_target_epitopes=options.require_target_epitopes_in_variant)
    representatives = dataset.selection.get('representatives')
    selected = (set((r['source_observation_id'], r['candidate_allele']) for r in representatives)
                if representatives is not None else dataset.select_representatives())
    groups = {}
    for epitope in epitopes:
        key = epitope.prediction_group_source
        antigen = dataset.antigens.get(key)
        if antigen is None or not antigen.tumor_specificity.admits_construct:
            continue
        if selected is not None:
            scores = {a: s for a, s in epitope.per_allele_scores.items() if (key, a) in selected}
            if not scores:
                continue
            epitope = replace(epitope, per_allele_scores=scores)
        identity = to_native_json(antigen)
        groups.setdefault(identity, (antigen, []))[1].append(epitope)
    for antigen, candidates in groups.values():
        vaccine = VaccinePeptide(
            antigen=antigen, epitopes=candidates,
            num_target_epitopes_to_keep=options.num_target_epitopes_to_keep,
            combined_score_expr=options.combined_score_expr,
            ranking_rules=options.ranking_rules,
            manufacturability_thresholds=options.manufacturability_thresholds,
            manufacturability_rules=options.manufacturability_rules)
        accumulator.add(ExternalVariantEntry(
            source=antigen.source_identifier or to_native_json(antigen),
            vaccine_peptide=vaccine, resolved_protein_context=True))
    return accumulator.result(report.source_format)
