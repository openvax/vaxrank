# Exploratory SV and neoORF tier

Status: proposal. Not implemented.

## Problem

The Nov 2026 Sid designs (`output/sid-ms-design-2026-11-08`,
`output/sid-rna-design-2026-11-09`) are 80-100% CTA. A 20-product peptide design
has 16-18 CTA, 2-4 mutations, and no SV, fusion, neoORF or viral product.
RNA-supported SV leads exist (TPST1::CRCP, a 30-aa upstream ORF; FOXO3::STRADA
and CCDC47, 21 aa) but never reach a slot.

## What already exists (read before building)

- **A per-type quota already exists.** `--min-<type>-antigens` and
  `--max-<type>-antigens` (`vaxrank/antigen_allocation.py`) serve minimum goals
  first, enforce hard maxima, then follow ranked order. #834 removed the
  *automatic* balancing mode because it could displace a stronger candidate
  merely to represent another kind (`ranked-antigen-limits.md`). #758 is closed.
  This spec adds **no new allocation algorithm**.
- **SV and neoORF candidates already enter** through the Isovar input
  (`vaxrank/isovar_antigens.py`: `_sv`, `_sv_export`, `_variants`). Each is
  admitted or `held_out` by `AdmissionPolicy`: `min_fragments`, `allow_partial`,
  `allow_ambiguous`, plus `no_targetable_amino_acid_change` and
  `rna_fragment_support_unavailable`. A held-out candidate is never selected, and
  allocation never admits one to meet a goal.

So the gap is **admission, not allocation**. Before designing anything, find
out which gate rejects TPST1::CRCP and the FOXO3 leads in the Sid run (partial
ORF? ambiguous frame? fragment threshold? no tumor-specificity attestation?).
The run provenance says only "original admission policies and known holdouts
are retained; no new viral or SV targets inferred", so this is unmeasured.

## Proposal

1. **Measure first.** A report listing, per held-out SV/fusion/neoORF candidate,
   its reasons and fragment support. No policy change; answers which gate to move.
2. **An `exploratory` admission status**, distinct from `admitted` and
   `held_out`. A candidate qualifies only if it is RNA-supported at the stated
   threshold, targetable, and not safety-flagged. Partial or ambiguous
   candidates stay held out unless the existing `allow_*` flags are set. RNA
   support is not proof of translation or tumor specificity, so the status and
   its reasons appear in the report and the native audit and cannot be
   presented as validated.
3. **Opt-in use, through existing controls.** Exploratory candidates are
   eligible only when the user passes `--include-antigen-type sv` (and `fusion`,
   etc.). They then compete in ranked order like any other candidate. To
   guarantee representation, the user sets `--min-sv-antigens N`; unmet goals
   are already reported (`insufficient_eligible_sequences` or
   `construct_capacity`).
4. **Named derived policy, not a change to `openvax-v1`.** A policy such as
   `exploratory-sv-v1` may bundle the include flag and a small default minimum
   (suggest 1-2). `builtin:openvax-v1` stays frozen (golden hash).

## RNA support must be comparable across kinds

#757 shipped (#822, `docs/design/rna-feature-support.md`): a common
`RNAFeatureEvidence` record (reads, fragments, per-million, effective-length
densities) and identical DSL bindings for every kind. An exploratory SV can sit
in a mixed list only if its support is counted the same way as a mutation's.
Known differences in the admission gate today (`isovar_antigens.py`):

- Mutation, SV export and fusion use fragment support, but with different
  scopes: mutation = alt-allele fragments; `_sv` = **minimum** direct-junction
  fragments; `_sv_export` = union of junction-linked witnesses per occurrence.
  `min_fragments` therefore does not mean one thing across kinds.
- CTA admission uses an expression table with its own unit (TPM or similar),
  not fragments.
- The default score `sqrt(n_rna_alt) * target_epitope_score` uses RNA only for
  mutations; other kinds default to `target_epitope_score` alone, so RNA does
  not currently affect their rank against mutations at all.

Required before mixing: (a) gate and rank every kind on the same
`RNAFeatureEvidence` measurement and unit, scoped to the allele-bearing
junction/interval, not the whole gene; (b) name a derived policy whose score
expression uses that binding for all kinds; (c) report the unit and scope next
to each product; (d) never convert TPM to fragments or the reverse. Missing
measurement stays "unknown", not zero. Check whether the Sid RNA design used
such a named expression; I have not verified that.

## Non-goals

- No automatic category balancing or "underrepresented kind" priority.
- No change to scores, the Topiary DSL, or RNA normalization (#757).
- No inferred viral or SV targets. Candidates come from supplied Isovar results.

## Open questions

- Is the Sid blocker admission, or are SV candidates absent from the input
  files? Step 1 answers it; the answer decides whether steps 2-4 are needed.
- Should an exploratory product carry a wet-lab/validation flag in the report?

## Verification

Synthetic Isovar fixtures: SV and neoORF each admitted, held out for each
reason, and exploratory. Check that a weak exploratory candidate cannot displace
a stronger one without `--min-*`, and that an unmet goal is reported and not
filled from held-out candidates. Peptide and mRNA modalities, CLI/YAML
equivalence, native replay of a saved design, frozen `openvax-v1` hash.
