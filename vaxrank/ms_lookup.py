"""Capture source-verified Hitlist evidence for the MS scoring profile."""

from functools import lru_cache
import hashlib
from importlib.metadata import version
import json
from pathlib import Path

import pandas as pd

from .file_digest import file_sha256
from .ms_evidence import _texts


RECORDS_COLUMN = 'hitlist_ms_records_json'

_UNPERTURBED_CONDITION_CATEGORIES = frozenset({
    '', 'unperturbed', 'baseline', 'control', 'untreated', 'unknown', 'unannotated',
})
_UNPERTURBED_PERTURBATION_TEXTS = frozenset({
    '', 'none', 'unperturbed', 'untreated', 'control', 'baseline',
})
_APM_PERTURBATION_KEYWORDS = (
    'apm', 'hla-dm', 'hla_dm', 'hladm', 'erap', 'tap1', 'tap2', 'tapbp',
    'tapasin', 'b2m', 'psmb', 'immunoproteasome', 'ifng', 'ifn-g', 'ifn_g',
    'interferon', 'nlrc5', 'ciita', 'jak1', 'jak2', 'stat1',
)
_PERTURBATION_FACTOR_TYPES = frozenset({
    'drug', 'drug_treatment', 'chemical', 'perturbation', 'gene_perturbation',
    'genetic_perturbation', 'knockout', 'knockdown', 'crispr', 'sirna', 'shrna',
    'cytokine', 'cytokine_exposure', 'stimulus', 'treatment', 'apm_perturbation',
})
_INTERVENTION_CONDITION_COLUMNS = (
    'condition_drugs',
    'condition_cytokines',
    'condition_knockout_genes',
    'condition_knockdown_genes',
    'condition_overexpression_genes',
    'condition_genetic_variants',
    'condition_transfection',
    'condition_transduction',
    'condition_infection',
    'condition_stimulation',
    'condition_antigen_exposure',
)
_EXPLICIT_PERTURBATION_COLUMNS = (
    'perturbation_type',
    'perturbation_agent',
    'perturbation_target',
    'drug_treatment',
)


def _is_truthy_flag(value):
    if value is True:
        return True
    if isinstance(value, str):
        return value.strip().lower() == 'true'
    return False


def _infer_species_and_class_from_alleles(allele_values):
    """Derive canonical MHC species and class ('I'/'II') via mhcgnomes without guessing."""
    from mhcgnomes import ParseError, parse
    species_found = set()
    classes_found = set()
    for raw in _texts(allele_values):
        try:
            parsed = parse(raw, infer_class2_pairing=False)
        except (ParseError, TypeError, ValueError):
            continue
        species = getattr(parsed, 'species', None)
        if species is not None and getattr(species, 'name', None):
            species_found.add(species.name)
        if getattr(parsed, 'is_class1', False):
            classes_found.add('I')
        elif getattr(parsed, 'is_class2', False):
            classes_found.add('II')
        else:
            mhc_class = getattr(parsed, 'mhc_class', None)
            if mhc_class in ('I', 'II'):
                classes_found.add(mhc_class)
    species_name = next(iter(species_found)) if len(species_found) == 1 else ''
    mhc_class_name = next(iter(classes_found)) if len(classes_found) == 1 else None
    return species_name, mhc_class_name


def _normalize_lineage_fields(row):
    """Extract reviewed donor/specimen/origin/acquisition identities without inference.

    Only explicit reviewed axes with status == 'resolved' and a single unit ID
    qualify as ``unit_id``. Multi-unit, ambiguous, or missing axes remain
    unresolved (no ``unit_id``). Original/reanalysis connectivity is preserved
    via ``reuse_ids``.
    """
    engineered = row.get('is_engineered_mhc') is True
    donor_status = text_or_empty(row.get('donor_status'))
    specimen_status = text_or_empty(row.get('specimen_status'))
    origin_status = text_or_empty(row.get('experimental_origin_status'))
    acquisition_status = text_or_empty(row.get('acquisition_status'))

    donor_ids = _texts(row.get('donor_ids'))
    specimen_ids = _texts(row.get('specimen_ids'))
    origin_ids = _texts(row.get('experimental_origin_ids'))
    acquisition_ids = _texts(row.get('acquisition_ids'))
    explicit_reuse = _texts(row.get('reuse_ids'))

    unit_id = None
    unit_kind = None
    if not engineered and donor_status == 'resolved' and len(donor_ids) == 1:
        unit_id = donor_ids[0]
        unit_kind = 'native_donor'
    elif origin_status == 'resolved' and len(origin_ids) == 1:
        unit_id = origin_ids[0]
        unit_kind = 'engineered_origin'
    elif engineered and donor_status == 'resolved' and len(donor_ids) == 1:
        unit_id = donor_ids[0]
        unit_kind = 'engineered_origin'

    # Connect original and reanalyzed observations sharing reviewed origin,
    # specimen, acquisition, or explicit reuse identifiers.
    reuse_ids = sorted(set(explicit_reuse) | set(origin_ids) | set(specimen_ids) | set(acquisition_ids))
    return dict(
        unit_id=unit_id,
        unit_kind=unit_kind,
        donor_ids=list(donor_ids),
        specimen_ids=list(specimen_ids),
        origin_ids=list(origin_ids),
        acquisition_ids=list(acquisition_ids),
        lineage_statuses={
            k: v for k, v in (
                ('donor', donor_status),
                ('specimen', specimen_status),
                ('experimental_origin', origin_status),
                ('acquisition', acquisition_status),
            ) if v
        },
        reuse_ids=reuse_ids,
    )


def text_or_empty(value):
    from .cells import text
    return text(value)


def normalize_observation(row):
    """Normalize a Hitlist positive-MS row without inventing source identity."""
    from .cells import text
    reasons = []
    if row.get('has_ptm') is not False:
        reasons.append('modified_or_unknown_chemical_form')
    if row.get('is_non_peptide_ligand') is True:
        reasons.append('non_peptide_ligand')
    if row.get('is_potential_contaminant') is True:
        reasons.append('potential_contaminant')
    sample_attribution = text(row.get('sample_attribution'))
    if sample_attribution == 'elution_conditions_excluded':
        reasons.append('excluded_elution_condition')
    condition_category = text(row.get('condition_category'))
    perturbation = text(row.get('perturbation'))
    gene_perturbation = text(row.get('gene_perturbation'))
    cytokine_exposure = text(row.get('cytokine_exposure'))
    apm_perturbed_raw = row.get('apm_perturbed')
    apm_perturbed = _is_truthy_flag(apm_perturbed_raw)
    is_control_arm = _is_truthy_flag(row.get('is_control_arm'))
    perturbation_texts = [
        p for p in (perturbation, gene_perturbation, cytokine_exposure)
        if p and p.lower() not in _UNPERTURBED_PERTURBATION_TEXTS
    ]
    intervention_fields = {}
    for col in _INTERVENTION_CONDITION_COLUMNS + _EXPLICIT_PERTURBATION_COLUMNS:
        raw_val = text(row.get(col))
        if not raw_val:
            continue
        tokens = [
            tok.strip() for tok in raw_val.split(';')
            if tok.strip() and tok.strip().lower() not in _UNPERTURBED_PERTURBATION_TEXTS
        ]
        if tokens:
            intervention_fields[col] = '; '.join(tokens)
    all_perturbation_texts = perturbation_texts + list(intervention_fields.values())
    if intervention_fields:
        if apm_perturbed or any(
            any(kw in p.lower() for kw in _APM_PERTURBATION_KEYWORDS)
            for p in all_perturbation_texts
        ):
            reasons.append('apm_perturbed')
        reasons.append('perturbed_condition')
    elif not is_control_arm:
        if apm_perturbed or any(
            any(kw in p.lower() for kw in _APM_PERTURBATION_KEYWORDS)
            for p in perturbation_texts
        ):
            reasons.append('apm_perturbed')
        if ((condition_category and condition_category.lower() not in _UNPERTURBED_CONDITION_CATEGORIES)
                or bool(perturbation_texts)):
            reasons.append('perturbed_condition')
    references = []
    pmid = text(row.get('pmid') if row.get('pmid') is not None else row.get('source_pmid'))
    if pmid:
        try:
            pmid = str(int(float(pmid)))
        except ValueError:
            raise ValueError(f'Invalid Hitlist PMID: {pmid!r}') from None
        references.append('pmid:' + pmid)
    ref_iri = text(row.get('reference_iri') or row.get('source_reference_iri') or row.get('source_reference_id'))
    if ref_iri:
        references.append(ref_iri)
    identity = (
        text(row.get('observation_id'))
        or text(row.get('provenance_id'))
        or text(row.get('assay_iri'))
        or text(row.get('source_assay_id'))
    )
    if not identity:
        # An unresolved source row can support presence; its content digest is
        # an audit locator, never a donor, specimen or independent observation.
        identity = 'unresolved:' + hashlib.sha256(
            json.dumps(row, sort_keys=True, allow_nan=False).encode()).hexdigest()

    # Preserve experimental vs prediction-derived allele attribution and typing.
    if row.get('cellular_mhc_genotype') is not None:
        typing = row.get('cellular_mhc_genotype')
    elif row.get('mhc_allele_provenance') in ('sample_allele_match', 'sample_locus_match', 'reviewed_arm_genotype'):
        typing = row.get('mhc_allele_set')
    elif row.get('candidate_mhc_set') is not None and row.get('mhc_allele_provenance') != 'restriction':
        typing = row.get('candidate_mhc_set')
    else:
        typing = None

    restriction = text(
        row.get('peptide_mhc_restriction')
        or row.get('mhc_restriction')
        or row.get('reported_mhc')
        or row.get('tested_mhc')
    )
    mhc_species = text(row.get('mhc_species'))
    mhc_class_explicit = 'mhc_class' in row
    mhc_class = text(row.get('mhc_class')) or None
    if not mhc_species or (not mhc_class_explicit and mhc_class not in ('I', 'II')):
        inferred_species, inferred_class = _infer_species_and_class_from_alleles(
            list(_texts(restriction)) + list(_texts(typing))
        )
        if not mhc_species:
            mhc_species = inferred_species
        if not mhc_class_explicit and mhc_class not in ('I', 'II'):
            mhc_class = inferred_class

    basis = text(row.get('restriction_evidence')) or row.get('restriction_evidence')
    if not basis:
        basis = 'unknown'
    # Unresolved class supports sequence transfer only.
    if mhc_class not in ('I', 'II'):
        typing, restriction, basis = None, None, 'unknown'

    monoallelic_verified = (
        (row.get('is_monoallelic') is True or text(row.get('mhc_claim_role')) == 'monoallelic_host')
        and basis == 'monoallelic'
        and bool(restriction)
    )
    source_prediction_qualified = bool(row.get('source_prediction_qualified') is True)
    peptide_seq = row.get('peptide') if row.get('peptide') is not None else row['peptide_sequence']
    result = dict(
        observation_id=identity,
        peptide=peptide_seq,
        eligible=not reasons,
        exclusion_reasons=reasons,
        mhc_species=mhc_species,
        mhc_class=mhc_class,
        restriction=restriction,
        restriction_basis=basis,
        monoallelic_verified=monoallelic_verified,
        typed_alleles=list(_texts(typing)),
        references=references,
        source_prediction_qualified=source_prediction_qualified,
        quality=1.,
        context_weight=1.,
        source_dataset=text(row.get('source') or row.get('source_dataset')),
        source_assay_id=text(row.get('assay_iri') or row.get('source_assay_id')),
        chemical_form_status='unmodified' if row.get('has_ptm') is False else 'modified_or_unknown',
    )
    if condition_category:
        result['condition_category'] = condition_category
    if perturbation:
        result['perturbation'] = perturbation
    if gene_perturbation:
        result['gene_perturbation'] = gene_perturbation
    if cytokine_exposure:
        result['cytokine_exposure'] = cytokine_exposure
    if apm_perturbed_raw is not None:
        result['apm_perturbed'] = apm_perturbed
    if sample_attribution:
        result['sample_attribution'] = sample_attribution
    if is_control_arm:
        result['is_control_arm'] = True
    for col, val in intervention_fields.items():
        result[col] = val

    # Preserve source context and genotype uncertainty metadata when provided.
    source_context = {
        key: text(row.get(key))
        for key in (
            'arm_id', 'condition_id', 'elution_condition_id',
            'cell_type', 'cell_line_name', 'arm_link_status',
        )
        if text(row.get(key))
    }
    if source_context:
        result['source_context'] = source_context
    genotype_uncertainty = text(row.get('genotype_uncertainty'))
    if genotype_uncertainty:
        result['genotype_uncertainty'] = genotype_uncertainty
    if row.get('genotype_complete_loci') is not None:
        result['genotype_complete_loci'] = list(_texts(row.get('genotype_complete_loci')))

    lineage = _normalize_lineage_fields(row)
    if lineage['unit_id']:
        result['unit_id'] = lineage['unit_id']
        result['unit_kind'] = lineage['unit_kind']
    if lineage['reuse_ids']:
        result['reuse_ids'] = lineage['reuse_ids']
    for key in ('donor_ids', 'specimen_ids', 'origin_ids', 'acquisition_ids', 'lineage_statuses'):
        if lineage[key]:
            result[key] = lineage[key]
    return result


@lru_cache(maxsize=4)
def verified_package(directory, manifest_digest):
    """Verify a complete legacy MHC producer package once per manifest identity."""
    from hitlist.mhc_evidence import verify_mhc_evidence_package
    path = Path(directory)
    verified = verify_mhc_evidence_package(path)
    if file_sha256(path / 'manifest.json') != manifest_digest:
        raise ValueError('Hitlist package changed during verification')
    return verified


@lru_cache(maxsize=4)
def verified_evidence_store(directory, manifest_digest):
    """Open and verify a Hitlist EvidenceStore package once per manifest identity."""
    try:
        from hitlist import EvidenceStore
    except ImportError as exc:
        raise ValueError(
            'Hitlist EvidenceStore packages require hitlist with EvidenceStore support'
        ) from exc
    path = Path(directory)
    store = EvidenceStore.open(path)
    if file_sha256(path / 'manifest.json') != manifest_digest:
        raise ValueError('Hitlist package changed during verification')
    return store


def _query_store_df(store, evidence_type, *, where=None):
    """Query an EvidenceStore relation/observation table as a pandas DataFrame."""
    result = store.query(evidence_type, where=where, limit=None)
    if isinstance(result, pd.DataFrame):
        return result
    if hasattr(result, 'to_pandas'):
        return result.to_pandas()
    if hasattr(result, 'table'):
        return result.table.to_pandas()
    return pd.DataFrame(list(result))


def _summarize_coverage_diagnostics(normalized_rows, *, available_tables=(), unavailable_facets=()):
    """Summarize lineage, attribution, and eligibility coverage for replay & audit."""
    exclusion_counts = {}
    basis_counts = {}
    resolved_native = 0
    resolved_engineered = 0
    unresolved = 0
    with_reuse = 0
    unique_native_units = set()
    unique_engineered_units = set()
    eligible = 0
    for rec in normalized_rows:
        if rec.get('eligible') is True:
            eligible += 1
        else:
            for reason in rec.get('exclusion_reasons') or ():
                exclusion_counts[reason] = exclusion_counts.get(reason, 0) + 1
        basis = rec.get('restriction_basis') or 'unknown'
        basis_counts[basis] = basis_counts.get(basis, 0) + 1
        unit_id = rec.get('unit_id')
        unit_kind = rec.get('unit_kind')
        if unit_id and unit_kind == 'native_donor':
            resolved_native += 1
            unique_native_units.add(unit_id)
        elif unit_id and unit_kind == 'engineered_origin':
            resolved_engineered += 1
            unique_engineered_units.add(unit_id)
        else:
            unresolved += 1
        if rec.get('reuse_ids'):
            with_reuse += 1
    return dict(
        total_observations=len(normalized_rows),
        eligible_observations=eligible,
        excluded_observations=len(normalized_rows) - eligible,
        exclusion_reason_counts=dict(sorted(exclusion_counts.items())),
        restriction_basis_counts=dict(sorted(basis_counts.items())),
        resolved_native_donor_observations=resolved_native,
        resolved_engineered_origin_observations=resolved_engineered,
        unresolved_lineage_observations=unresolved,
        unique_resolved_units=dict(
            native_donor=len(unique_native_units),
            engineered_origin=len(unique_engineered_units),
        ),
        observations_with_reuse_links=with_reuse,
        available_tables=sorted(available_tables),
        unavailable_facets=sorted(unavailable_facets),
    )


def _capture_evidence_store_package(directory, digest, peptides, species):
    """Capture mhc_positive_ms + reviewed typing/lineage from a verified EvidenceStore."""
    store = verified_evidence_store(str(directory), digest)
    bindings = set(getattr(store, 'bindings', {}) or ())
    if 'mhc_positive_ms' not in bindings:
        raise ValueError(
            f'Hitlist EvidenceStore package at {directory} does not materialize mhc_positive_ms'
        )
    available_tables = [
        name for name in (
            'mhc_positive_ms', 'assay_context', 'reviewed_lineage',
            'arm_mhc_claim', 'arm_factor', 'arm_control_link',
        )
        if name in bindings
    ]
    unavailable_facets = []
    if 'assay_context' not in bindings:
        unavailable_facets.append('assay_context')
    if 'reviewed_lineage' not in bindings:
        unavailable_facets.append('reviewed_lineage')
    if 'arm_mhc_claim' not in bindings:
        unavailable_facets.append('arm_mhc_claim')
    if 'arm_factor' not in bindings:
        unavailable_facets.append('arm_factor')

    obs_df = _query_store_df(store, 'mhc_positive_ms')
    if not obs_df.empty and 'peptide_sequence' in obs_df.columns:
        obs_df = obs_df.loc[obs_df['peptide_sequence'].isin(peptides)].copy()
    elif not obs_df.empty and 'peptide' in obs_df.columns:
        obs_df = obs_df.loc[obs_df['peptide'].isin(peptides)].copy()

    ctx_by_obs = {}
    if not obs_df.empty and 'assay_context' in bindings:
        obs_ids = set(obs_df['observation_id'].dropna().astype(str))
        ctx_df = _query_store_df(store, 'assay_context')
        if not ctx_df.empty and 'observation_id' in ctx_df.columns:
            ctx_df = ctx_df.loc[ctx_df['observation_id'].astype(str).isin(obs_ids)]
            for rec in json.loads(ctx_df.to_json(orient='records')):
                ctx_by_obs[str(rec['observation_id'])] = rec

    arm_ids = {
        str(rec['arm_id']) for rec in ctx_by_obs.values()
        if rec.get('arm_id') is not None and str(rec['arm_id']).strip()
    }

    lineage_by_arm = {}
    if arm_ids and 'reviewed_lineage' in bindings:
        lin_df = _query_store_df(store, 'reviewed_lineage')
        if not lin_df.empty and 'arm_id' in lin_df.columns:
            lin_df = lin_df.loc[lin_df['arm_id'].astype(str).isin(arm_ids)]
            for rec in json.loads(lin_df.to_json(orient='records')):
                lineage_by_arm.setdefault(str(rec['arm_id']), []).append(rec)

    mhc_claims_by_arm = {}
    if arm_ids and 'arm_mhc_claim' in bindings:
        claim_df = _query_store_df(store, 'arm_mhc_claim')
        if not claim_df.empty and 'arm_id' in claim_df.columns:
            claim_df = claim_df.loc[claim_df['arm_id'].astype(str).isin(arm_ids)]
            for rec in json.loads(claim_df.to_json(orient='records')):
                mhc_claims_by_arm.setdefault(str(rec['arm_id']), []).append(rec)

    factors_by_arm = {}
    if arm_ids and 'arm_factor' in bindings:
        factor_df = _query_store_df(store, 'arm_factor')
        if not factor_df.empty and 'arm_id' in factor_df.columns:
            factor_df = factor_df.loc[factor_df['arm_id'].astype(str).isin(arm_ids)]
            for rec in json.loads(factor_df.to_json(orient='records')):
                factors_by_arm.setdefault(str(rec['arm_id']), []).append(rec)

    control_arm_ids = set()
    if arm_ids and 'arm_control_link' in bindings:
        ctrl_df = _query_store_df(store, 'arm_control_link')
        if not ctrl_df.empty and 'control_arm_id' in ctrl_df.columns:
            control_arm_ids = {
                str(v) for v in ctrl_df['control_arm_id'].dropna() if str(v).strip()
            }

    species_set = set(species)
    rows_by_species = {sp: [] for sp in species}
    normalized_for_diag = []
    raw_obs_rows = json.loads(obs_df.to_json(orient='records')) if not obs_df.empty else []

    for obs in raw_obs_rows:
        merged = dict(obs)
        merged['peptide'] = merged.get('peptide') or merged.get('peptide_sequence')
        obs_id = text_or_empty(merged.get('observation_id'))
        ctx = ctx_by_obs.get(obs_id)
        arm_id = ''
        if ctx:
            arm_id = text_or_empty(ctx.get('arm_id'))
            for key in (
                'arm_id', 'condition_id', 'elution_condition_id',
                'sample_attribution', 'arm_link_status',
                'cell_type', 'cell_line_name', 'donor_status', 'specimen_status',
                'acquisition_status',
            ) + _INTERVENTION_CONDITION_COLUMNS + _EXPLICIT_PERTURBATION_COLUMNS:
                if ctx.get(key) is not None and not merged.get(key):
                    merged[key] = ctx[key]
            if ctx.get('experimental_origin_ids'):
                merged['experimental_origin_ids'] = ctx['experimental_origin_ids']
            if ctx.get('cellular_mhc_genotype') is not None:
                merged['cellular_mhc_genotype'] = ctx['cellular_mhc_genotype']
            elif ctx.get('candidate_mhc_set') is not None and not merged.get('candidate_mhc_set'):
                merged['candidate_mhc_set'] = ctx['candidate_mhc_set']

        # Apply reviewed arm MHC claims (experimental vs predicted, monoallelic host, genotype uncertainty).
        claims = mhc_claims_by_arm.get(arm_id, ())
        if claims:
            typed_alleles = list(_texts(merged.get('cellular_mhc_genotype')))
            complete_loci = set()
            for claim in claims:
                role = text_or_empty(claim.get('mhc_claim_role'))
                basis = text_or_empty(claim.get('evidence_basis'))
                if role == 'monoallelic_host':
                    merged['mhc_claim_role'] = 'monoallelic_host'
                    merged['is_monoallelic'] = True
                    merged['is_engineered_mhc'] = True
                    if basis:
                        merged['restriction_evidence'] = basis
                    if claim.get('reported_mhc') and not merged.get('peptide_mhc_restriction'):
                        merged['peptide_mhc_restriction'] = claim['reported_mhc']
                elif role == 'restriction' and basis and not text_or_empty(merged.get('restriction_evidence')):
                    merged['restriction_evidence'] = basis
                if claim.get('cellular_mhc_genotype'):
                    typed_alleles.extend(_texts(claim['cellular_mhc_genotype']))
                    merged['mhc_allele_provenance'] = 'reviewed_arm_genotype'
                elif claim.get('candidate_mhc_set'):
                    typed_alleles.extend(_texts(claim['candidate_mhc_set']))
                    merged['mhc_allele_provenance'] = 'reviewed_arm_genotype'
                if claim.get('genotype_complete_loci'):
                    complete_loci.update(_texts(claim['genotype_complete_loci']))
            if typed_alleles:
                merged['cellular_mhc_genotype'] = sorted(set(typed_alleles))
            if complete_loci:
                merged['genotype_complete_loci'] = sorted(complete_loci)
                merged['genotype_uncertainty'] = 'complete' if len(complete_loci) >= 3 else 'partial'
            else:
                merged['genotype_uncertainty'] = 'unknown'
        elif not merged.get('cellular_mhc_genotype'):
            merged['genotype_uncertainty'] = 'unknown'

        # Apply reviewed lineage claims strictly per axis (never infer across axes).
        lin_claims = lineage_by_arm.get(arm_id, ())
        if lin_claims:
            axis_ids = {}
            axis_statuses = {}
            for claim in lin_claims:
                axis = text_or_empty(claim.get('lineage_axis') or claim.get('axis'))
                status = text_or_empty(claim.get('lineage_status') or claim.get('status'))
                if not axis:
                    continue
                ids = set(_texts(claim.get('identity_ids')))
                for col, axis_name in (
                    ('donor_id', 'donor'),
                    ('specimen_id', 'specimen'),
                    ('acquisition_id', 'acquisition'),
                ):
                    if axis == axis_name and claim.get(col):
                        ids.update(_texts(claim[col]))
                if axis in axis_statuses and axis_statuses[axis] != status:
                    axis_statuses[axis] = 'ambiguous'
                else:
                    axis_statuses[axis] = status
                axis_ids.setdefault(axis, set()).update(ids)
            if 'donor' in axis_statuses:
                merged['donor_status'] = axis_statuses['donor']
                merged['donor_ids'] = sorted(axis_ids.get('donor', ()))
            if 'specimen' in axis_statuses:
                merged['specimen_status'] = axis_statuses['specimen']
                merged['specimen_ids'] = sorted(axis_ids.get('specimen', ()))
            if 'experimental_origin' in axis_statuses:
                merged['experimental_origin_status'] = axis_statuses['experimental_origin']
                merged['experimental_origin_ids'] = sorted(axis_ids.get('experimental_origin', ()))
            if 'acquisition' in axis_statuses:
                merged['acquisition_status'] = axis_statuses['acquisition']
                merged['acquisition_ids'] = sorted(axis_ids.get('acquisition', ()))
            if 'reuse' in axis_ids or 'original_study' in axis_ids:
                merged['reuse_ids'] = sorted(
                    set(_texts(merged.get('reuse_ids')))
                    | axis_ids.get('reuse', set())
                    | axis_ids.get('original_study', set())
                )

        # Apply reviewed perturbation factors & control-arm links.
        if arm_id and arm_id in control_arm_ids:
            merged['is_control_arm'] = True
        arm_factors = factors_by_arm.get(arm_id, ())
        if arm_factors:
            pert_vals = []
            for factor in arm_factors:
                ftype = text_or_empty(factor.get('factor_type')).lower()
                fval = text_or_empty(factor.get('factor_value'))
                fstatus = text_or_empty(factor.get('claim_status')).lower()
                if fstatus in ('absent', 'negative', 'control', 'unperturbed'):
                    continue
                if ftype in _PERTURBATION_FACTOR_TYPES and fval and fval.lower() not in _UNPERTURBED_PERTURBATION_TEXTS:
                    pert_vals.append(f'{ftype}:{fval}')
            if pert_vals and not merged.get('is_control_arm'):
                merged['perturbation'] = '; '.join(sorted(set(pert_vals)))
                if not merged.get('condition_category'):
                    merged['condition_category'] = 'perturbed'

        row_species = text_or_empty(merged.get('mhc_species'))
        if not row_species:
            row_species, inferred_class = _infer_species_and_class_from_alleles(
                list(_texts(merged.get('peptide_mhc_restriction')))
                + list(_texts(merged.get('reported_mhc')))
                + list(_texts(merged.get('tested_mhc')))
                + list(_texts(merged.get('cellular_mhc_genotype')))
                + list(_texts(merged.get('candidate_mhc_set')))
            )
            if row_species:
                merged['mhc_species'] = row_species
            if inferred_class and not merged.get('mhc_class'):
                merged['mhc_class'] = inferred_class
        if row_species in species_set:
            rows_by_species[row_species].append(merged)
            normalized_for_diag.append(normalize_observation(merged))

    if file_sha256(directory / 'manifest.json') != digest:
        raise ValueError('Hitlist package changed during capture')

    manifest = dict(store.manifest)
    diagnostics = _summarize_coverage_diagnostics(
        normalized_for_diag,
        available_tables=available_tables,
        unavailable_facets=unavailable_facets,
    )
    return dict(
        schema_version=2,
        hitlist_version=version('hitlist'),
        status='captured',
        query=dict(peptides=peptides, mhc_species=species),
        data_directory=str(directory),
        index_metadata_sha256=digest,
        index_metadata=manifest,
        source_contract='verified_hitlist_evidence_store',
        package_schema_version=manifest.get('schema_version'),
        package_family=manifest.get('package_family') or manifest.get('parent_kind') or 'evidence_store',
        coverage_diagnostics=diagnostics,
        species_observations={
            species_name: dict(
                columns=sorted({k for r in rows for k in r}),
                rows=rows,
            )
            for species_name, rows in rows_by_species.items()
        },
    )


def capture_package(directory, sequences, alleles):
    """Query a verified Hitlist MHC package or EvidenceStore with hash verification."""
    import pyarrow.dataset as ds
    from mhcgnomes import parse
    from .allele_validation import parse_prediction_allele

    directory = Path(directory).resolve()
    manifest_path = directory / 'manifest.json'
    if not manifest_path.exists():
        raise ValueError(f'Hitlist package is missing manifest.json: {directory}')
    digest = file_sha256(manifest_path)
    envelope = json.loads(manifest_path.read_text(encoding='utf-8'))
    species = sorted({parse(parse_prediction_allele(a), infer_class2_pairing=False).species.name
                      for a in alleles if a})
    peptides = sorted(set(sequences))

    legacy_path = directory / 'native' / 'observations.parquet'
    is_legacy_mhc = (
        envelope.get('schema_version') == 2
        and legacy_path.exists()
        and 'package_family' not in envelope
    )
    if not is_legacy_mhc:
        return _capture_evidence_store_package(directory, digest, peptides, species)

    manifest = verified_package(str(directory), digest)
    before = legacy_path.stat()
    frame = ds.dataset(legacy_path, format='parquet').to_table(filter=(
        ds.field('peptide').isin(peptides) & ds.field('mhc_species').isin(species))).to_pandas()
    # Current package verification establishes the positive-ligand partition's
    # modality/outcome contract. Retain every field and the package identity.
    after = legacy_path.stat()
    if (before.st_ino, before.st_mtime_ns, before.st_size) != (
            after.st_ino, after.st_mtime_ns, after.st_size) or (
                file_sha256(manifest_path) != digest):
        raise ValueError('Hitlist package changed during capture')
    species_observations = {
        species_name: dict(columns=list(frame), rows=json.loads(
            frame.loc[frame.mhc_species.eq(species_name)].to_json(orient='records')))
        for species_name in species
    }
    normalized_for_diag = [
        normalize_observation(row)
        for table in species_observations.values()
        for row in table['rows']
    ]
    diagnostics = _summarize_coverage_diagnostics(
        normalized_for_diag,
        available_tables=('native/observations.parquet',),
        unavailable_facets=('assay_context', 'reviewed_lineage', 'arm_mhc_claim', 'arm_factor'),
    )
    return dict(schema_version=1, hitlist_version=version('hitlist'), status='captured',
        query=dict(peptides=peptides, mhc_species=species), data_directory=str(directory),
        index_metadata_sha256=digest, index_metadata=manifest,
        source_contract='verified_hitlist_mhc_package',
        coverage_diagnostics=diagnostics,
        species_observations=species_observations)


class MSLookup:
    """One profile lookup with portable source captures and normalized rows."""

    def __init__(self, config):
        from .ms_evidence import validate_config
        self.config = validate_config(config)
        self.queries = []

    def __call__(self, frame):
        from .hitlist_evidence import capture_local_hitlist_evidence, annotate_hitlist_frame
        missing = (frame[RECORDS_COLUMN].isna() if RECORDS_COLUMN in frame
                   else pd.Series(True, index=frame.index))
        if not missing.any():
            return
        rows = frame.loc[missing].copy()
        if self.config.get('lookup', 'auto') == 'off':
            evidence = dict(status='disabled', query=dict(peptides=sorted(set(rows.peptide)),
                mhc_species=[]), species_observations={})
        elif self.config.get('package'):
            evidence = capture_package(self.config['package'], rows.peptide, rows.allele.fillna(''))
        else:
            evidence = capture_local_hitlist_evidence(rows.peptide, rows.allele.fillna(''))
        if self.queries:
            previous = self.queries[0]
            identity = ('status', 'hitlist_version', 'data_directory', 'index_metadata_sha256')
            if any(evidence.get(key) != previous.get(key) for key in identity):
                raise ValueError('Hitlist data changed between MS queries')
        if self.config.get('lookup') == 'required' and evidence['status'] != 'captured':
            raise ValueError('Required MS evidence is ' + evidence['status'])
        if evidence['status'] != 'disabled':
            annotate_hitlist_frame(rows, evidence)
        self._annotate_records(rows, evidence)
        for name in rows:
            if name.startswith('hitlist_'):
                if name not in frame:
                    frame[name] = rows[name].reindex(frame.index)
                else:
                    frame[name] = frame[name].astype(object)
                    frame.loc[missing, name] = rows[name].to_numpy()
        self.queries.append(evidence)

    @staticmethod
    def _annotate_records(frame, evidence):
        records = {}
        for table in evidence.get('species_observations', {}).values():
            for row in table['rows']:
                records.setdefault(row['peptide'], []).append(normalize_observation(row))
        frame[RECORDS_COLUMN] = [json.dumps(records.get(p, []), sort_keys=True, allow_nan=False)
                                 for p in frame.peptide]
        frame['hitlist_ms_availability'] = evidence['status']

    def annotate_captured(self, frame):
        for evidence in self.queries:
            selected = frame.peptide.isin(evidence['query']['peptides'])
            rows = frame.loc[selected].copy()
            self._annotate_records(rows, evidence)
            for name in (RECORDS_COLUMN, 'hitlist_ms_availability'):
                if name not in frame:
                    frame[name] = rows[name].reindex(frame.index)
                else:
                    frame[name] = frame[name].astype(object)
                    frame.loc[selected, name] = rows[name].to_numpy()


def complete_candidate_frame(frame, sequences, predictor, *, alleles=None):
    """Enumerate declared native windows even when a model emits no result.

    Missing rows carry no prediction kind, method or value. They are candidate
    identities, not synthetic predictions. Existing predictor rows are retained.
    """
    from .allele_validation import parse_prediction_allele
    models = getattr(predictor, 'models', (predictor,))
    requests = set()
    for model in models:
        model_alleles = getattr(model, 'alleles', ()) or alleles or ()
        lengths = getattr(model, 'default_peptide_lengths', ())
        if model_alleles and not lengths:
            raise ValueError('MS candidate enumeration requires declared model peptide lengths')
        for allele in model_alleles:
            for length in lengths:
                if type(length) is not int or length < 1:
                    raise ValueError('MS candidate peptide lengths must be positive integers')
                requests.add((parse_prediction_allele(allele), length))
    if not requests:
        raise ValueError('MS native prediction requires declared patient alleles and peptide lengths')
    keys = ['source_sequence_name', 'peptide_offset', 'peptide', 'allele']
    existing = set(frame[keys].itertuples(index=False, name=None)) if not frame.empty else set()
    missing = []
    for name, sequence in sequences.items():
        for allele, length in sorted(requests):
            for offset in range(len(sequence) - length + 1):
                peptide = sequence[offset:offset + length]
                if (name, offset, peptide, allele) not in existing:
                    missing.append(dict(source_sequence_name=name, peptide_offset=offset,
                        peptide=peptide, peptide_length=length, allele=allele,
                        kind=None, prediction_method_name=None, predictor_version=None,
                        value=None, score=None, percentile_rank=None))
    result = pd.concat([frame, pd.DataFrame(missing)], ignore_index=True) if missing else frame.copy()
    result.attrs = dict(frame.attrs)
    return result
