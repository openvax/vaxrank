"""Vaxrank's consumer boundary for Topiary selection policies."""

from dataclasses import replace
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from topiary import SelectionPolicy, TopiaryResult, evaluate_selection_policy


def canonical_digest(definition):
    """Identify the complete effective configuration, independent of key order."""
    text = json.dumps(definition, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def evaluate_frame(frame, cfg, *, group_keys, alleles=None, kind_support=None,
                   default_methods=None, default_versions=None):
    """Evaluate named criteria with the exact occurrence context of the caller.

    The returned score series contains only eligible, finite scores. Every
    decision, including rejected and unknown observations, stays in Topiary's
    evaluation result. Policy evaluation never executes predictors.
    """
    policy = SelectionPolicy.from_dict(cfg.selection_policy)
    # Inferred defaults fill gaps; a portable policy's explicit choices win.
    methods = {**(default_methods or {}), **(policy.default_methods or {}),
               **(cfg.default_methods or {})}
    versions = {**(default_versions or {}), **(policy.default_versions or {})}
    policy = replace(policy, default_methods=methods or None,
                     default_versions=versions or None)
    evaluation = evaluate_selection_policy(
        TopiaryResult(frame, form="long"), policy,
        group_keys=group_keys, alleles=alleles, kind_support=kind_support)
    selected = evaluation.selected
    if selected.empty:
        scores = pd.Series(dtype=float)
    else:
        selected = selected.loc[np.isfinite(selected.score.astype(float))]
        scores = selected.set_index(list(group_keys)).score.astype(float)
    scores.attrs["policy_evaluation"] = evaluation
    return scores


def combine_evaluations(evaluations, definition, *, duplicates=None):
    """Preserve per-source runtime choices for one representative selection."""
    if not evaluations:
        return None
    policy = SelectionPolicy.from_dict(definition)
    if duplicates is not None:
        policy = replace(policy, duplicates=duplicates)
    frames, contexts = [], {}
    group_keys = None
    partitions = []
    for evaluation in evaluations:
        for saved in evaluation.evidence.extra['policy_evaluation']['partitions']:
            frame = evaluation.evidence.long_df
            if saved['source_label'] is not None:
                frame = frame.loc[frame.source_label == saved['source_label']]
            partitions.append((frame, saved['context']))
    for i, (original, context) in enumerate(partitions):
        if group_keys is not None and group_keys != context['group_keys']:
            raise ValueError('Source policy evaluations use different occurrence identities')
        group_keys = context['group_keys']
        label = 'vaxrank-policy-source-%d' % i
        frame = original.copy()
        if 'source_label' in frame and 'vaxrank_original_source_label' not in frame:
            frame['vaxrank_original_source_label'] = frame.source_label
        frame['source_label'] = label
        frame.attrs = {}
        frames.append(frame)
        alleles = context['alleles']
        if alleles and isinstance(alleles[0], dict):
            peptide_keys = [k for k in group_keys if k != 'allele']
            alleles = {tuple(row['keys'][k] for k in peptide_keys): row['alleles'] for row in alleles}
        contexts[label] = dict(
            alleles=alleles, default_methods=context['default_methods'],
            default_versions={(r['kind'], r['method']): r['version'] for r in context['default_versions']},
            kind_support=context['kind_support'])
    return evaluate_selection_policy(TopiaryResult(pd.concat(frames, ignore_index=True)), policy,
                                     group_keys=group_keys, source_contexts=contexts)


def save_policy_evaluations(evaluations, directory):
    """Write replayable Topiary evidence and human-readable decision tables."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    records = []
    for i, evaluation in enumerate(evaluations, 1):
        stem = "policy-%04d" % i
        evidence_path = directory / (stem + "-evidence.tsv")
        evaluation.evidence.to_tsv(evidence_path)
        evaluation.occurrences.to_csv(directory / (stem + "-decisions.tsv"), sep="\t", index=False)
        evaluation.audit.to_csv(directory / (stem + "-criteria.tsv"), sep="\t", index=False)
        records.append(dict(name=evaluation.policy.name, sha256=evaluation.policy.sha256,
                            evidence=evidence_path.name,
                            evidence_sha256=hashlib.sha256(evidence_path.read_bytes()).hexdigest()))
    (directory / "index.json").write_text(json.dumps(records, indent=2) + "\n")
    return records


def encode_evaluations(evaluations):
    """Embed Topiary's typed format, using its public lossless serializer."""
    from tempfile import TemporaryDirectory
    with TemporaryDirectory() as directory:
        path = Path(directory) / "evidence.tsv"
        encoded = []
        for evaluation in evaluations:
            evaluation.evidence.to_tsv(path)
            encoded.append(path.read_text())
        return encoded


def decode_evaluations(encoded):
    from tempfile import TemporaryDirectory
    from topiary import replay_selection_policy
    from topiary import read_tsv
    with TemporaryDirectory() as directory:
        path = Path(directory) / "evidence.tsv"
        evaluations = []
        for text in encoded:
            path.write_text(text)
            evaluations.append(replay_selection_policy(read_tsv(path)))
        return evaluations


def write_run_policy(args, evaluations, epitope_config, vaccine_config):
    """Record resolved settings, CLI overrides and predictor decision context."""
    import msgspec
    from .version import __version__
    output_dir = getattr(args, "output_dir", None)
    if not output_dir:
        output = getattr(args, "output_epitopes", None)
        if not output:
            return
        output_dir = str(output) + ".policy"
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    from .config.loader import configuration_provenance
    resolved = configuration_provenance(args)
    resolved["effective_epitopes"] = msgspec.to_builtins(epitope_config)
    resolved["effective_vaccine_peptides"] = msgspec.to_builtins(vaccine_config)
    resolved["vaxrank_version"] = __version__
    resolved["effective_sha256"] = configuration_digest(resolved)
    resolved["evaluations"] = save_policy_evaluations(evaluations, directory / "policy_evidence")
    (directory / "selection_policy.json").write_text(json.dumps(resolved, indent=2) + "\n")
    return resolved


def configuration_digest(record):
    return canonical_digest({key: value for key, value in record.items()
                             if key == 'configuration' or (
                                 key.startswith('effective_') and key != 'effective_sha256')})


def record_construct_configuration(args, modality, options):
    """Include actual CLI-resolved construct settings in the run identity."""
    from dataclasses import asdict
    directory = getattr(args, 'output_dir', None)
    if not directory:
        return
    path = Path(directory) / 'selection_policy.json'
    if not path.exists():
        return
    record = json.loads(path.read_text())
    record.setdefault('effective_constructs', {})[modality] = asdict(options)
    record['effective_sha256'] = configuration_digest(record)
    path.write_text(json.dumps(record, indent=2) + '\n')
