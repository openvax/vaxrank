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

"""Consolidated per-target longitudinal library x timepoint support matrix (#818).

Surfaces observed below-threshold CTA expression (distinguishing missing from
measured zero), per-library exact-window vs. event/different-window mutation
RNA support, qualifying timepoints, and targets qualified only through earlier
samples without pooling libraries or treating repeated pipelines as independent
samples. Omitted on classic single-BAM VCF runs so vanilla OpenVax 1 reports
remain uncluttered.
"""

from collections import OrderedDict, defaultdict
import json
import re

from varcode import Variant


def _timepoint_sort_key(timepoint: str):
    if not timepoint:
        return (1, (), "")
    tokens = [
        (0, int(part)) if part.isdigit() else (1, part.lower())
        for part in re.findall(r"\d+|\D+", timepoint)
    ]
    return (0, tuple(tokens), timepoint)


def _event_identity(source, antigen):
    if antigen is not None and antigen.kind == "CTA":
        return ("CTA", antigen.gene_id or antigen.gene_name or antigen.source_identifier)
    if antigen is not None and antigen.origin is not None and antigen.origin.variants:
        variants = tuple(
            (v.reference_assembly, v.contig, v.start, v.ref, v.alt, v.sv_type, v.end)
            for v in antigen.origin.variants
        )
        return ("variant", variants)
    if isinstance(source, Variant):
        return (
            "variant",
            ((source.genome.reference_name, source.contig, source.start, source.ref, source.alt,
              getattr(source, "sv_type", ""), getattr(source, "end", None)),),
        )
    if antigen is not None:
        return (antigen.kind, antigen.source_identifier or antigen.gene_id or antigen.gene_name)
    return ("unknown", str(source))


def _format_number(value):
    if value is None:
        return None
    if value == 0:
        return "0"
    return f"{value:.6g}"


def build_longitudinal_support_view(
    ranked,
    products=None,
    selection=None,
    patient_info=None,
):
    """Build consolidated longitudinal/CTA library support matrix for reports.

    Returns ``None`` when no multi-library, timepoint-scoped, or CTA-admission
    evidence is present (preserving uncluttered vanilla single-BAM VCF reports).
    """
    from .assembly_antigens import emitted_antigens_for_candidate
    from .cta_expression import CTAExpressionResult, load_cta_expression_result
    from .species_cta import SpeciesCTAExpressionResult

    def _parse_cta_admission(raw):
        if isinstance(raw, (CTAExpressionResult, SpeciesCTAExpressionResult)):
            return raw
        payload = raw if isinstance(raw, str) else json.dumps(raw)
        return load_cta_expression_result(payload)

    selection = selection or {}
    coverage = getattr(patient_info, "input_coverage", None) or {}
    original_inputs = coverage.get("original_inputs") or ()
    provenances = list(getattr(patient_info, "input_provenance", None) or ())
    source_evidence = selection.get("source_evidence") or {}

    # Collect CTA admission records keyed by input_id (or standalone key).
    cta_admissions = OrderedDict()
    for ev_id, ev_record in source_evidence.items():
        if isinstance(ev_record, dict) and "cta_expression_admission" in ev_record:
            cta_admissions[ev_id] = _parse_cta_admission(
                ev_record["cta_expression_admission"]
            )
    if not cta_admissions and selection.get("cta_expression_admission"):
        standalone_id = provenances[0].input_id if len(provenances) == 1 else "cta_expression"
        cta_admissions[standalone_id] = _parse_cta_admission(
            selection["cta_expression_admission"]
        )

    # Collect library/sample columns from declared input_provenance, CTA contracts,
    # Isovar source_evidence, and attached RNAFeatureEvidence scopes.
    # Keyed by (sample_id, library_id) so repeated pipelines on one library
    # remain dependent runs in a single column rather than separate samples.
    libraries_by_key = OrderedDict()
    input_to_lib_key = {}

    def _register_library(sample_id, library_id=None, timepoint=None, input_id=None, pipeline=None):
        sample = sample_id or library_id or input_id or "unknown"
        library = library_id or sample
        key = (sample, library)
        entry = libraries_by_key.get(key)
        if entry is None:
            entry = dict(
                key=f"{sample}:{library}" if library != sample else sample,
                sample_id=sample,
                library_id=library,
                timepoint=timepoint or "",
                input_ids=[],
                pipelines=[],
            )
            libraries_by_key[key] = entry
        elif timepoint and not entry["timepoint"]:
            entry["timepoint"] = timepoint
        if input_id and input_id not in entry["input_ids"]:
            entry["input_ids"].append(input_id)
            input_to_lib_key[input_id] = key
        if pipeline and pipeline not in entry["pipelines"]:
            entry["pipelines"].append(pipeline)
        return key

    if original_inputs:
        for prov_item in original_inputs:
            get_prov = prov_item.get if isinstance(prov_item, dict) else lambda k: getattr(prov_item, k, None)
            scope_obj = get_prov("scope") or {}
            get_scope = scope_obj.get if isinstance(scope_obj, dict) else lambda k: getattr(scope_obj, k, None)
            input_id = get_prov("input_id")
            sample = get_scope("sample_id") or get_scope("patient_id") or input_id
            _register_library(
                sample,
                library_id=get_scope("library_id") or get_scope("sample_id") or input_id,
                timepoint=get_scope("timepoint"),
                input_id=input_id,
                pipeline=get_prov("source_format"),
            )
    else:
        for prov in provenances:
            scope = prov.scope
            _register_library(
                scope.sample_id or scope.patient_id or prov.input_id,
                library_id=scope.library_id or scope.sample_id or prov.input_id,
                timepoint=scope.timepoint,
                input_id=prov.input_id,
                pipeline=prov.source_format,
            )

    for _, peptides in ranked:
        for vp in peptides:
            for ep in vp.epitopes:
                prov = ep.input_provenance
                if prov is not None and prov.input_id not in input_to_lib_key:
                    scope = prov.scope
                    _register_library(
                        scope.sample_id or scope.patient_id or prov.input_id,
                        library_id=scope.library_id or scope.sample_id or prov.input_id,
                        timepoint=scope.timepoint,
                        input_id=prov.input_id,
                        pipeline=prov.source_format,
                    )

    def _cta_pipeline_label(contract):
        src = str(contract.evidence_source or "").strip()
        ver = str(contract.evidence_version or "").strip()
        if src and "/" not in src and len(src) <= 40:
            return f"{src}:{ver}" if ver and len(ver) <= 24 else src
        assay = str(getattr(contract, "assay", "") or "").strip()
        return assay or "cta_expression"

    for input_id, admission in cta_admissions.items():
        contract = admission.input_contract
        existing_key = input_to_lib_key.get(input_id)
        pipeline_label = _cta_pipeline_label(contract)
        if existing_key is not None:
            lib = libraries_by_key[existing_key]
            if pipeline_label not in lib["pipelines"]:
                lib["pipelines"].append(pipeline_label)
        else:
            _register_library(
                contract.sample_id,
                library_id=contract.sample_id,
                input_id=input_id,
                pipeline=pipeline_label,
            )

    for ev_id, ev_record in source_evidence.items():
        if not isinstance(ev_record, dict) or "isovar_result" not in ev_record:
            continue
        res = ev_record["isovar_result"]
        sample = res.get("sample_id") or res.get("provenance", {}).get("sample_id")
        src = res.get("source") or res.get("provenance", {}).get("source")
        if sample:
            # Match existing sample column if one was registered with a specific library_id
            existing = next((k for k in libraries_by_key if k[0] == sample), None)
            if existing is None:
                _register_library(sample, library_id=sample, pipeline=src or "isovar")

    for _, peptides in ranked:
        for vp in peptides:
            for rna_ev in getattr(vp.antigen, "rna_support", ()):
                sample, src = rna_ev.measurement["evidence_scope"]
                existing = next((k for k in libraries_by_key if k[0] == sample), None)
                if existing is None:
                    _register_library(sample, library_id=src, pipeline="rna_feature_support")

    has_timepoints = any(lib["timepoint"] for lib in libraries_by_key.values())
    if len(libraries_by_key) <= 1 and not has_timepoints and not cta_admissions:
        return None

    ordered_lib_keys = sorted(
        libraries_by_key,
        key=lambda k: (_timepoint_sort_key(libraries_by_key[k]["timepoint"]), k[0], k[1]),
    )
    libraries = [libraries_by_key[k] for k in ordered_lib_keys]

    distinct_timepoints = sorted(
        {lib["timepoint"] for lib in libraries if lib["timepoint"]},
        key=_timepoint_sort_key,
    )
    if len(distinct_timepoints) > 1 or (not distinct_timepoints and len(libraries) > 1):
        design_scope = "longitudinal_union"
        design_scope_label = (
            "Longitudinal multi-sample design (targets admitted in any qualifying sample "
            "remain eligible; library counts are never pooled)"
        )
    else:
        design_scope = "single_timepoint"
        tp_text = distinct_timepoints[0] if distinct_timepoints else (
            libraries[0]["sample_id"] if libraries else "single sample"
        )
        design_scope_label = f"Single-timepoint / single-sample design ({tp_text})"

    latest_timepoint = distinct_timepoints[-1] if len(distinct_timepoints) > 1 else None

    # Index CTA decisions by (lib_key, canonical_or_source_gene_id)
    cta_by_lib_gene = defaultdict(list)
    cta_gene_meta = OrderedDict()
    cta_gene_observations = defaultdict(list)
    small_cta_panel = True
    for input_id, admission in cta_admissions.items():
        lib_key = input_to_lib_key.get(input_id) or (
            admission.input_contract.sample_id,
            admission.input_contract.sample_id,
        )
        cta_decisions = [
            d for d in admission.decisions
            if d.status not in {"not_cta", "unresolved_reference", "unsafe_reference", "not_current_cta"}
        ]
        if len(cta_decisions) > 50:
            small_cta_panel = False
        for decision in cta_decisions:
            canonical_gene = (
                (decision.gene_identity or {}).get("canonical_gene_id")
                or decision.gene_id
            )
            gene_key = canonical_gene or decision.input_identifier
            cta_by_lib_gene[(lib_key, gene_key)].append((admission, decision))
            cta_gene_observations[gene_key].append(decision)
            if gene_key not in cta_gene_meta and (
                decision.gene_id or decision.status in {"admitted", "held_out", "unknown_expression"}
            ):
                cta_gene_meta[gene_key] = dict(
                    gene_id=gene_key,
                    gene_name=decision.gene_name or gene_key,
                    transcript_id=decision.transcript_id,
                )

    all_cta_genes = OrderedDict(
        (gene_key, meta)
        for gene_key, meta in cta_gene_meta.items()
        if small_cta_panel or any(
            d.status in {"admitted", "unknown_expression"}
            or (d.status == "held_out" and d.value is not None and d.value > 0)
            for d in cta_gene_observations[gene_key]
        )
    )

    # Index Isovar source results by (sample_id, variant_coords_or_event)
    isovar_by_sample_event = defaultdict(list)
    for ev_id, ev_record in source_evidence.items():
        if not isinstance(ev_record, dict) or "isovar_result" not in ev_record:
            continue
        res = ev_record["isovar_result"]
        sample = res.get("sample_id") or res.get("provenance", {}).get("sample_id")
        var = res.get("variant")
        if sample and isinstance(var, dict) and var.get("contig") and var.get("start") is not None:
            vkey = (
                var.get("reference_name") or var.get("reference_assembly") or "",
                str(var["contig"]),
                int(var["start"]),
                str(var.get("ref") or ""),
                str(var.get("alt") or ""),
            )
            isovar_by_sample_event[(sample, vkey)].append((ev_id, ev_record))

    all_products = [p for group in (products or {}).values() for p in group]

    # Group ranked candidates by biological target event + exact window sequence
    target_entries = OrderedDict()
    ranked_windows_by_lib_event = defaultdict(list)

    for source, peptides in ranked:
        for vp in peptides:
            antigen = vp.antigen
            ev_key = _event_identity(source, antigen)
            input_ids = {
                e.input_provenance.input_id for e in vp.epitopes if e.input_provenance
            }
            matched_libs = {
                input_to_lib_key[iid] for iid in input_ids if iid in input_to_lib_key
            }
            for rna_ev in getattr(antigen, "rna_support", ()):
                sample = rna_ev.measurement["evidence_scope"][0]
                for lk in ordered_lib_keys:
                    if lk[0] == sample:
                        matched_libs.add(lk)
            for lk in matched_libs:
                ranked_windows_by_lib_event[(lk, ev_key)].append(vp)

            selected_products = sorted({
                prod_name
                for prod_name, names in emitted_antigens_for_candidate(
                    antigen, vp.amino_acids, vp.epitopes, all_products
                ).items()
                if names
            })
            row_key = (
                ev_key if antigen.kind == "CTA" else (ev_key, vp.amino_acids)
            )
            entry = target_entries.get(row_key)
            if entry is None:
                entry = dict(
                    event_key=ev_key,
                    kind=antigen.kind,
                    gene=antigen.display_gene_name,
                    gene_id=antigen.gene_id,
                    transcript_ids=list(antigen.transcript_ids),
                    window=vp.amino_acids,
                    peptides=[],
                    selected_products=[],
                )
                target_entries[row_key] = entry
            entry["peptides"].append((source, vp))
            for prod_name in selected_products:
                if prod_name not in entry["selected_products"]:
                    entry["selected_products"].append(prod_name)
            for tid in antigen.transcript_ids:
                if tid and tid not in entry["transcript_ids"]:
                    entry["transcript_ids"].append(tid)

    # Also include held-out / below-threshold CTA targets that produced 0 ranked windows
    # so longitudinal/single-sample below-threshold CTA evidence is retained in the matrix!
    for gene_key, meta in all_cta_genes.items():
        ev_key = ("CTA", gene_key)
        if ev_key not in target_entries:
            target_entries[ev_key] = dict(
                event_key=ev_key,
                kind="CTA",
                gene=meta["gene_name"],
                gene_id=gene_key,
                transcript_ids=[meta["transcript_id"]] if meta["transcript_id"] else [],
                window="",
                peptides=[],
                selected_products=[],
            )

    rows = []
    for entry in target_entries.values():
        ev_key = entry["event_key"]
        kind = entry["kind"]
        window = entry["window"]
        cells_out = []
        qualifying_timepoints = []
        qualifying_samples = []
        latest_qualified = False
        latest_below_threshold = False

        for lk in ordered_lib_keys:
            lib = libraries_by_key[lk]
            tp = lib["timepoint"]
            sample = lib["sample_id"]

            if kind == "CTA":
                gene_key = ev_key[1]
                records = cta_by_lib_gene.get((lk, gene_key), ())
                # Also match by raw gene_id if canonical differed
                if not records and entry["gene_id"] and entry["gene_id"] != gene_key:
                    records = cta_by_lib_gene.get((lk, entry["gene_id"]), ())

                window_rna = []
                for _, vp in entry["peptides"]:
                    if vp.amino_acids == window or not window:
                        for rna_ev in getattr(vp.antigen, "rna_support", ()):
                            if rna_ev.measurement["evidence_scope"][0] == sample:
                                vals = rna_ev.measurement["values"]
                                window_rna.append(
                                    f"{vals['reads']}r/{vals['fragments']}f ({rna_ev.scope})"
                                )

                if not records:
                    if window_rna:
                        state = "window_rna_only"
                        summary_text = "; ".join(dict.fromkeys(window_rna))
                        qualified = True
                    else:
                        state = "unmeasured"
                        summary_text = "unmeasured"
                        qualified = False
                    cells_out.append(dict(
                        library_key=lib["key"],
                        state=state,
                        qualified=qualified,
                        summary=summary_text,
                        details={},
                    ))
                else:
                    parts = []
                    qualified = False
                    below = False
                    cell_details = []
                    for admission, dec in records:
                        contract = admission.input_contract
                        policy = admission.admission_policy
                        thresh = _format_number(policy.min_tumor_expression)
                        receipt = admission.input_sha256[:12]
                        feat_label = dec.transcript_id or dec.gene_id or dec.input_identifier
                        if dec.value is None:
                            state_i = "missing"
                            parts.append(
                                f"missing ({contract.measurement_level} {feat_label}; "
                                f"threshold >={thresh} {contract.expression_unit}; sha256:{receipt})"
                            )
                        else:
                            val_str = _format_number(dec.value)
                            if dec.status == "admitted":
                                state_i = "qualified"
                                qualified = True
                            elif dec.value < policy.min_tumor_expression:
                                state_i = "below_threshold"
                                below = True
                            else:
                                state_i = dec.status or "held_out"
                            zero_note = " [measured zero]" if dec.value == 0 else ""
                            parts.append(
                                f"{val_str} {contract.expression_unit}{zero_note} "
                                f"({state_i.replace('_', ' ')}; threshold >={thresh} {contract.expression_unit}; "
                                f"{contract.measurement_level} {feat_label}; "
                                f"{contract.evidence_source} {contract.evidence_version}; sha256:{receipt})"
                            )
                        cell_details.append(dict(
                            feature_id=dec.input_identifier,
                            gene_id=dec.gene_id,
                            transcript_id=dec.transcript_id,
                            level=contract.measurement_level,
                            value=dec.value,
                            measured_zero=(dec.value == 0),
                            missing=(dec.value is None),
                            unit=contract.expression_unit,
                            threshold=policy.min_tumor_expression,
                            status=dec.status,
                            evidence_source=contract.evidence_source,
                            evidence_version=contract.evidence_version,
                            input_sha256=admission.input_sha256,
                        ))
                    if window_rna:
                        parts.append("window RNA: " + ", ".join(dict.fromkeys(window_rna)))
                    overall_state = (
                        "qualified" if qualified
                        else "below_threshold" if below
                        else "missing" if all(d["missing"] for d in cell_details)
                        else cell_details[0]["status"] or "held_out"
                    )
                    cells_out.append(dict(
                        library_key=lib["key"],
                        state=overall_state,
                        qualified=qualified,
                        summary="; ".join(parts),
                        decisions=cell_details,
                    ))
                    if below and tp and tp == latest_timepoint:
                        latest_below_threshold = True
                if cells_out[-1]["qualified"]:
                    if tp and tp not in qualifying_timepoints:
                        qualifying_timepoints.append(tp)
                    if sample not in qualifying_samples:
                        qualifying_samples.append(sample)
                    if tp and tp == latest_timepoint:
                        latest_qualified = True
            else:
                # Non-CTA target: separate exact-window vs. different-window vs. event RNA vs. unmeasured
                lib_vps = ranked_windows_by_lib_event.get((lk, ev_key), ())
                exact_vps = [vp for vp in lib_vps if vp.amino_acids == window]
                other_vps = [vp for vp in lib_vps if vp.amino_acids != window]

                # Also check attached RNAFeatureEvidence on any peptide for this sample
                exact_rna_ev = []
                for _, vp in entry["peptides"]:
                    if vp.amino_acids != window:
                        continue
                    for rna_ev in getattr(vp.antigen, "rna_support", ()):
                        if rna_ev.measurement["evidence_scope"][0] == sample:
                            exact_rna_ev.append(rna_ev)

                if exact_vps or exact_rna_ev:
                    parts = []
                    for rna_ev in exact_rna_ev:
                        vals = rna_ev.measurement["values"]
                        scope_tag = "allele/event span" if rna_ev.scope == "event" else "exact window"
                        parts.append(
                            f"{vals['reads']} reads, {vals['fragments']} fragments ({scope_tag})"
                        )
                    if not parts:
                        for vp in exact_vps:
                            frag = vp.mutant_protein_fragment
                            if frag is not None:
                                item = f"{frag.n_alt_reads} alt reads"
                                if frag.n_alt_fragments is not None:
                                    item += f", {frag.n_alt_fragments} alt fragments"
                                parts.append(item + " (variant fragment)")
                            else:
                                cands = [
                                    json.loads(val)
                                    for key_m, val in vp.antigen.source_metadata
                                    if key_m == "isovar_candidate"
                                ]
                                for cand in cands:
                                    sup = cand.get("occurrence_rna_support") or cand.get("rna_support") or {}
                                    if sup.get("fragments") is not None:
                                        r_txt = f"{sup['reads']} reads, " if sup.get("reads") is not None else ""
                                        parts.append(f"{r_txt}{sup['fragments']} fragments (reconstructed ORF)")
                    summary_text = "; ".join(dict.fromkeys(parts)) or "ranked (exact window)"
                    qualified = bool(exact_vps)
                    cells_out.append(dict(
                        library_key=lib["key"],
                        state="exact_window",
                        qualified=qualified,
                        summary=summary_text,
                    ))
                elif other_vps:
                    other_seqs = ", ".join(dict.fromkeys(vp.amino_acids for vp in other_vps))
                    cells_out.append(dict(
                        library_key=lib["key"],
                        state="different_window_only",
                        qualified=True,
                        summary=(
                            f"different ranked window ({other_seqs}); "
                            f"exact window {window} not ranked (does not imply absent RNA)"
                        ),
                    ))
                else:
                    # Check Isovar source_evidence for this sample + variant event
                    iso_matches = []
                    if ev_key[0] == "variant" and len(ev_key[1]) == 1:
                        vtuple = ev_key[1][0][:5]
                        iso_matches = isovar_by_sample_event.get((sample, vtuple), ())
                    if iso_matches:
                        parts = []
                        for _, ev_rec in iso_matches:
                            res = ev_rec["isovar_result"]
                            cands = res.get("protein_candidates") or res.get("candidates") or ()
                            frags = [
                                (c.get("rna_support") or {}).get("fragments")
                                for c in cands
                                if (c.get("rna_support") or {}).get("fragments") is not None
                            ]
                            if frags:
                                parts.append(
                                    f"no ranked window; source RNA observed (max {max(frags)} fragments)"
                                )
                            else:
                                parts.append(
                                    f"no ranked window ({res.get('assembly_state') or res.get('status') or '0 coding candidates'})"
                                )
                        cells_out.append(dict(
                            library_key=lib["key"],
                            state="no_ranked_window",
                            qualified=False,
                            summary="; ".join(dict.fromkeys(parts)),
                        ))
                    else:
                        cells_out.append(dict(
                            library_key=lib["key"],
                            state="unmeasured",
                            qualified=False,
                            summary="unmeasured in library",
                        ))
                if cells_out[-1]["qualified"]:
                    if tp and tp not in qualifying_timepoints:
                        qualifying_timepoints.append(tp)
                    if sample not in qualifying_samples:
                        qualifying_samples.append(sample)
                    if tp and tp == latest_timepoint:
                        latest_qualified = True

        qualified_only_earlier = (
            design_scope == "longitudinal_union"
            and latest_timepoint is not None
            and bool(qualifying_timepoints)
            and not latest_qualified
        )
        from .input_coverage import report_source_key
        source_keys = list(dict.fromkeys(
            report_source_key(src) for src, _ in entry["peptides"]
        ))
        rows.append(dict(
            gene=entry["gene"],
            gene_id=entry["gene_id"],
            kind=kind,
            source_keys=source_keys,
            transcript_ids=entry["transcript_ids"],
            window=window,
            selected=bool(entry["selected_products"]),
            selected_products=entry["selected_products"],
            qualifying_timepoints=qualifying_timepoints,
            qualifying_samples=qualifying_samples,
            qualified_only_earlier=qualified_only_earlier,
            below_threshold_at_latest=latest_below_threshold,
            cells=cells_out,
        ))

    return dict(
        design_scope=design_scope,
        design_scope_label=design_scope_label,
        latest_timepoint=latest_timepoint,
        libraries=libraries,
        rows=rows,
    )
