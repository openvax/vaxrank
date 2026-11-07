# MS evidence for every antigen category

Status: proposed algorithm, 2026-11-06. This document changes no runtime behavior.
Vaxrank baseline: 3.51.2, including global Hitlist annotation in PR #809.
Evidence/scoring contracts inspected with Hitlist 2.1.9, Topiary 5.109.0 and
MHCtools 4.8.1.
Configuration and feature names below are proposals.

The current `--hitlist-local-evidence` option already annotates supported input
paths before filtering/scoring and saves evidence for replay; CTA expression
inputs query automatically. This plan extends that foundation with automatic
lookup across categories, recurrence-aware scoring and direct-MS admission.

## Scope and decisions

Enrich exact candidate peptides from mutations, fusions, CTAs, viruses, supplied
antigens and imported predictions through one Hitlist path before filtering.
Use installed data automatically, without implicit downloads. A named
`ms-evidence-v1` selection policy applies the algorithm below; annotation alone
leaves the existing ranking unchanged.

There is one scoring algorithm: qualify positive ligand-MS observations, assign
patient-allele support weights, deduplicate source units, apply diminishing
returns, then combine the result with a bounded prediction score. All antigen
categories use the same rules. Weights and thresholds are explicit policy
parameters, not separate algorithms for different input types.

The result is a ranking heuristic, not a calibrated probability of presentation,
immune response or tumor specificity. Source observations do not establish
presentation in the current patient. Antigen admission, expression, exact-self
and tissue-exclusion rules remain separate requirements.

## Hitlist evidence contract

Hitlist is the acquisition and curation layer for IEDB, CEDAR and supplementary
sources. Vaxrank must not maintain a second IEDB importer. Target coverage is the
union of eligible source evidence, with overlapping records deduplicated.

| Available contract | Use in this algorithm |
|---|---|
| `load_ms_observations`, or verified `EvidenceStore` relation `mhc_positive_ms` | Positive MHC-ligand MS observations eligible for scoring |
| `mhc_reported_measurement` | Original qualitative/numeric measurement, endpoint, units, inequality and reporting statuses |
| `mhc_allele_provenance`, `restriction_evidence`, reported restriction and candidate sets | Distinguish candidate-set construction from experimental, monoallelic, predicted or unknown restriction |
| Source contributors, assay/material links and `arm_mhc_claim` where materialized | Original record identities, source sample attribution, independent cellular typing and capture context |
| Reviewed lineage axes/statuses | Experimental origin, specimen, donor and acquisition; resolved, unknown, pooled or ambiguous |
| Source-specific quantitative MS relations | Additional ligand-MS measurements only when modality, identity and context satisfy the same contract |
| Binding assays, T-cell assays, bulk proteomics, synthetic standards and processing assays | Retain their own modality; they do not contribute ligand-MS observations |
| `mhc_other` / `load_other_assays` | Retain negative or unestablished MS outcomes for audit; no positive contribution |

The [Hitlist evidence-store contract](https://github.com/pirl-unc/hitlist/blob/main/docs/evidence-store.md)
specifies mandatory positive-MS identifiers and conditional measurement fields.
A declared catalog field is not proof that an installed generation contains it:
inspect `EvidenceStore.catalog()` and capture availability. The existing flat
reader supports basic observation annotation; richer claims require the relevant
materialized relations. Both readers must emit one normalized contract and
converge on the same source identities. Merge overlapping package/index records
by producer identities rather than adding copies.

Retain source/assay/observation IDs and contributors; exact peptide and chemical
form; modality/outcome; raw measurements and quality scope; publication IDs;
material/condition and native/engineered context; reported restriction, candidate
sets, typing provenance and restriction basis; lineage IDs/statuses; source
protein mappings and species roles. Missing fields remain missing. Supplemental
positive admission does not manufacture a reported qualitative measurement.

Before claiming a superset of an earlier IEDB extraction, reconcile eligible
source assay IDs and exact peptide/restriction records against a pinned export.
Every omission needs a disposition: duplicate, changed sequence/restriction,
unsupported modality/outcome, source not installed, or ingestion gap.
Protein-only MS and negative assays are not positive ligand evidence merely
because an earlier extraction matched "mass spectrometry". Verified ingestion
gaps belong upstream in Hitlist.

## Algorithm

### 1. Capture and qualify observations

Batch unique exact peptide queries by compatible MHC species, retaining all
source restrictions rather than querying only patient-matching restrictions.
Keep MHC species, host organism and peptide-source species distinct: a viral
peptide presented by human HLA is relevant to a human patient.

Positive MHC-ligand MS with compatible known MHC species is the scoring input.
Known incompatible MHC class cannot support an allele of another class; unknown
class permits only sequence-level transfer until an association resolves it.
Retain other captured observations with exclusion reasons.

Do not transfer mutant evidence from WT, nested peptides, or modified forms
without an explicit chemical-equivalence rule. Preserve unknown modification
status separately from explicitly unmodified sequence. Retain multiple
protein/gene/transcript occurrences. A shared peptide gets the same sequence
evidence in each occurrence, with distinct provenance and admission rules.

Set an observation quality multiplier `q[o]` in [0, 1]. Initially q=1 for a
producer-established positive passing the configured quality gates, and q=0 for
an excluded observation. This is a policy weight, not a probability. Missing
numerical quality is allowed and reported by default; a strict gate must name
the required field, scope and missing-value behavior. Source-reviewed mappings
may downweight weak identifications. Qualitative strength, intensity, spectra,
FDR and PEP are not interchangeable confidence scales; study-level FDR is not
per-peptide confidence.

Intensity, abundance, spectra and repeated acquisitions remain reportable
measurements, without automatic cross-study abundance bonuses or extra
biological observations. Explicit negatives are not subtracted: assay absence,
failed identification and untested peptides differ. Disputed positives are
reviewed or excluded with an explicit reason.

### 2. Assign support for each patient allele

Use mhcgnomes identity/resolution, including class-II pairs. Evaluate all patient
alleles; do not discard an observed association because another allele has the
best prediction. Retain all applicable facets but choose only the highest weight:

| Facet | Required support for patient allele a | Initial weight |
|---|---|---:|
| Sequence observed | Positive exact peptide; no patient-allele association established | 0.25 |
| Sequence plus patient prediction | Above, plus qualifying peptide–a prediction; source typing may be unknown or nonmatching | 0.5 |
| Typed overlap | a occurs in typing attributed to source material; restriction remains ambiguous | 1 |
| Recorded association, basis unknown | Exact reported peptide–a association with unestablished restriction basis | 1 |
| Predicted source-compatible association | Qualifying peptide–a prediction with source-attributed compatible typing, or source-reported prediction with sufficient method/threshold provenance | 3 |
| Experimental restriction | Source-backed experimental restriction of this MS-observed peptide to a | 6 |
| Monoallelic MS | Verified monoallelic source context and expressed molecule matching a | 8 |

Sequence-only credit is weak transfer support, not an assertion that a presents
the peptide. It cannot activate direct-MS admission. Setting its weight to zero
disables this transfer without losing the observation.

A precise allele name is not proof of experimental attribution:
`mhc_allele_provenance=exact` does not imply experimental restriction.
A singleton inferred set is not monoallelic; a publication-wide allele pool is
not a source sample's genotype. Source predictions lacking sufficient
qualification metadata fall back to recorded-unknown support, with their
predicted origin retained.

Use exact molecular identity under a named mhcgnomes equivalence rule for
allele-specific weights. Coarse typing, serotypes and unresolved alternatives
retain compatible sets without exact-match credit. Independently typed exact
source alleles can still support typed overlap despite coarse peptide
restriction. Do not invent phase or narrow a class-II pair or class-I alternative
to one favored molecule.

Prediction qualification names method/version, binding versus presentation
endpoint, metric, units and threshold by MHC class. Reuse patient predictions
and compatible captured source predictions. The initial affinity-based gate is
predicted IC50 <= 50 nM for the run's explicitly selected class-I affinity
predictor. It is a heuristic cutoff, not a universal definition of presentation.
Percentile ranks, presentation scores and other MHC classes require their own
named qualification rule; absent or unsupported rules give no predicted-facet
bonus. Source-reported assignments must satisfy a configured rule in their
native metric, not a guessed conversion between methods.
Extra source prediction is opt-in
through MHCtools/Topiary and its cache. Use source context when a model needs
flanks or expression. Absolute gates may qualify multiple alleles; retain all
without converting set size into a probability. A best-source-allele or margin
rule would require the full supported source genotype; this algorithm uses
absolute gates.

Let `t[o]` in [0, 1] be a material-context multiplier, initially 1 for eligible
native and engineered material. Report those contexts separately. Supplied or
private evidence is not stronger merely because of where it was stored.

    c[o,a] = q[o] * t[o] * max(applicable facet weights)

No qualifying observation gives no contribution. Prediction qualification does
not add an experiment; its bonus is model-dependent and can be zeroed for
sensitivity analysis.

### 3. Deduplicate and aggregate support

Use producer-backed identity/lineage. Repeated exports, alternate searches,
technical acquisitions and multiple allele rows do not create biological
recurrence. Preserve their contributing records.

Choose resolved donor for native material and resolved experimental origin for
engineered material. These are separate axes, not a pooled count of independent
samples. Reconcile material represented through both contexts before it can
contribute twice. A specimen-based configuration must report donor clustering.
Do not use Hitlist's descriptive `count_samples` fallbacks as donor identities.

    unit_support[u,a] = max(c[o,a] for observations attributed to resolved unit u)
    R[a] = sum(unit_support[u,a] over eligible resolved units)
    P[a] = max(c[o,a] over all accepted observations, default=0)

Retain native-donor and engineered-origin subtotals. P preserves strong evidence
when lineage is unresolved and is not added to R. Pooled evidence contributes
once through P, never once per possible donor. Twenty unresolved sample labels
do not become twenty resolved units.

Publications form a separate corroboration channel. Start from canonical
reference IDs; link references when Hitlist records source reuse or they share
a resolved recurrence unit for this peptide. Each connected group g contributes
only its strongest observation. Construct these identity groups before quality
or scoring filters, using all captured contributor links; filtering out a bridge
record must not split one known reused experiment into several groups:

    reference_support[g,a] = max(c[o,a] in reference group g)
    J[a] = reference_discount * sum(reference_support[g,a])
    M[a] = max(R[a], P[a], J[a])
    S[a] = min(ms_cap, log2(1 + M[a]))

Initially `reference_discount=0.25`, allowing limited bibliographic
corroboration when biological identity is incomplete. This does not label papers
independent experiments. Missing reference IDs add nothing to J; zero discount
disables this channel. Taking the maximum avoids adding the same support as
samples, papers and presence. Authoritative reuse/lineage corrections can lower
a previous score and must be recorded as evidence revisions.

Report matching, known-nonmatching and unresolved reference sets using unions
and differences, never sums of per-allele counts. One reference naming two patient
alleles remains one reference. Keep original bibliography and group membership.

Initially `ms_cap=8`. Combining weights before the logarithm prevents
facet-splitting bonuses. At fixed unit/reference grouping, qualification and
policy, additional distinct units or stronger facets cannot reduce support.

### 4. Combine with prediction

B[p,a] is a named, bounded [0, 1] binding/presentation score before eligibility
gates zero it out. The initial profile can reuse Vaxrank's normalized affinity
transform. Other metrics need an explicit bounded criterion; do not silently
clamp or reinterpret arbitrary user scores.

    H[a] = 1 - 2 ** (-S[a] / ms_scale)
    C[p,a] = B[p,a] + (1 - B[p,a]) * H[a]

Initially `ms_scale=4`; cap=8 limits the boost to 75% of remaining distance to 1.
No support yields C=B. Increasing support improves C with diminishing returns.
B=0 can receive positive support; B=1 stays 1. This is score arithmetic, not an
independent-probability calculation.

Keep missing B missing. Only an admitted direct-MS candidate (step 5) may use
B_effective=0 here, with explicit `score_basis=ms_only`; otherwise it has no
rankable combined score.

Compute each C from B and S for the same allele. Keep the allele vector; do not
combine the best prediction for one allele with another allele's MS support.
The peptide/window aggregation below applies the MS increment once while
preserving the existing prediction-only aggregation.

### 5. Apply eligibility before final ranking

Separate hard antigen requirements from prediction eligibility:

    admitted = hard_antigen_requirements AND (
        existing_prediction_gate OR qualified_direct_MS
    )

The proposed profile's direct-MS route requires at least one accepted exact
patient-allele experimental or verified monoallelic observation passing the
quality/context gates, with positive effective direct-facet weight and S>0.
Zeroing the MS contribution therefore cannot leave a hidden admission bypass.
Unknown donor identity does not erase verified restriction.
Predicted restriction, typed overlap and generic sequence presence
cannot activate this route. Expose a switch to disable it.

Enrich before the first prediction filter. When a direct-MS candidate has no
numerical prediction, create an evidence-backed peptide/allele candidate row,
never a fabricated affinity. All hard antigen requirements still apply.

### 6. Add one MS increment per peptide

Vaxrank currently sums per-allele scores into epitope scores and sums epitopes
into vaccine-window scores (`vaccine_peptide.py`). Preserve this prediction
baseline; summing C across alleles would award the same MS observation repeatedly.
For each allele admitted by step 5:

    b0[p,a] = B[p,a] if the existing prediction gate passes, otherwise 0
    D[p,a] = C[p,a] - b0[p,a]
    D[p] = max(D[p,a] over admitted patient alleles, default=0)
    peptide_score[p] = sum_a b0[p,a] + D[p]

The no-prediction/direct-MS case has b0=0 and D=C. If a direct observation admits
an otherwise prediction-filtered allele, C is its evidence-backed contribution.
Alleles failing both routes do not participate in the maximum. Record the allele
providing D, with deterministic ties by canonical allele identity. This can
differ from the best predicted allele; neither attribution is hidden.

At the window level, retain the existing prediction-only sum W0 and add D once
per unique eligible minimal-peptide sequence in that window. Repeated sequence
positions keep their source mappings but do not receive extra MS increments.
Compose the resulting presentation contribution with existing expression,
clonality and context terms. No MS leaves the original baseline unchanged.
Window/peptide sums need not lie in [0, 1]; only B, H and C have that domain.
Keep assembly and separate allele-coverage objectives explicit.

## Worked example

For B=0.20, q=t=1, cap=8, scale=4 and the weights above:

| Evidence | M | S | C, for one admitted allele |
|---|---:|---:|---:|
| No qualifying MS | 0 | 0 | 0.200 |
| Exact peptide, source typing unknown, patient prediction qualifies | 0.5 | 0.585 | 0.277 |
| One typed sample sharing the patient allele | 1 | 1.000 | 0.327 |
| Same sample, qualifying source-compatible prediction | 3 | 2.000 | 0.434 |
| That predicted support in 20 resolved native donors | 60 | 5.931 | 0.714 |
| Twenty technical acquisitions from the first donor | 3 | 2.000 | 0.434 |
| Verified monoallelic source, one origin or unresolved origin | 8 | 3.170 | 0.538 |
| Twenty unresolved sample labels within one reference | 3 | 2.000 | 0.434 |
| Twenty reference groups with predicted support, lineage unresolved | 15 | 4.000 | 0.600 |

These weights let extensive predicted corroboration outweigh one monoallelic
observation. Weights, reference discount and cap control that tradeoff; the
numbers are uncalibrated. The example assumes prediction eligibility passes
where there is no direct experimental restriction.

## Configuration, audit and replay

| Setting | Initial behavior |
|---|---|
| Lookup | `auto`, `off`, `required`; auto reads current installed data without acquiring it |
| Inputs | Installed Hitlist data/package and optional supplied MS under the same contract |
| Quality/context | Explicit gates and source-reviewed multipliers; missing quality allowed and reported |
| Allele matching | Named exact molecular equivalence; uncertain compatibility retained separately |
| Prediction qualification | Selected class-I IC50 <= 50 nM initially; named rules for other metrics/classes; extra source inference opt-in |
| Facet weights | 0.25 / 0.5 / 1 / 1 / 3 / 6 / 8 as listed |
| Recurrence | Native donors and engineered experimental origins; presence floor for unresolved identity |
| Reference discount | 0.25; zero disables bibliographic corroboration |
| Transform | cap=8, scale=4 |
| Base criterion | Explicit bounded binding/presentation score; prediction gate separately named |
| Direct-MS admission | Enabled in this opt-in profile; exact experimental/monoallelic support required |

Acquisition and evidence definitions belong under `evidence.ms`. Numeric
criteria, weights, gates and arithmetic compile into existing Topiary selection
policies, with no new expression language or ad-hoc ranking engine. Recompute
grouping/identity features from captured observations when policy changes,
rather than reweighting overlapping raw counts.

Reject non-finite/negative weights or caps, nonpositive scale, and multipliers
or reference discounts outside [0, 1]. A present B must be finite and in [0, 1].
Missing B follows the explicit admission rule, not a global missing-value fill.

Reports show availability, strongest facet and basis, all compatible alleles,
prediction provenance, counted unit/reference IDs, unresolved coverage,
R/P/J/M/S, B/H/C/D, admission route and final contribution. Show best overall,
best evidence-supported and increment-providing alleles, including disagreements.
Spectra, detected samples, specimens, donors, origins and references retain their
own units.

Distinguish unavailable, stale/unverifiable, query failure, successful zero
matches and positive evidence with incomplete metadata. A zero match is not an
experimental negative. Auto mode without usable data preserves prediction-only
ranking; required mode fails with the precise reason.

Capture source/package hashes, schema/curation versions, query scope, normalized
observations and excluded records, predictions, resolved parameters and feature
definitions. Ordinary replay uses this snapshot without installed Hitlist data
or predictors. Explicit refresh produces a new traceable snapshot.

## Implementation sequence

1. Verify public Hitlist field/identity contracts and source coverage with
   pinned real-format fixtures. Extend missing fields upstream rather than
   recreating curation joins. Resolve the allele boundary work tracked in
   Vaxrank [#789](https://github.com/openvax/vaxrank/issues/789). Preserve the
   explicit-filter fix for [#804](https://github.com/openvax/vaxrank/issues/804),
   which landed with [#809](https://github.com/openvax/vaxrank/pull/809).
2. Extend the shared `hitlist_evidence.py` path introduced in #809 with automatic
   installed-data lookup and the richer evidence contract. Keep direct VCF/BAM,
   generic antigen and imported prediction annotation before selection, with
   bundle precedence, query scopes and occurrence identities preserved.
3. Implement normalized observations, lineage/reference grouping, per-allele
   facets and R/P/J. Use Hitlist identities, mhcgnomes matching and
   MHCtools/Topiary predictions; retain missingness and contributor links.
4. Add this single named profile and direct-MS admission, including candidates
   without numerical predictions. Use Topiary `join_annotations` and named
   criteria. Add report explanations and replay fields.
5. Compare with the unchanged baseline. Test sensitivity to weights, reference
   discount, cap/scale, quality, context and predicted association. Use
   source/lineage-separated validation where available; record training overlap
   and do not claim independence for unresolved reused evidence.

Hitlist owns modality, curation, typing provenance and lineage. Relevant coverage
work includes [#855](https://github.com/pirl-unc/hitlist/issues/855),
[#680](https://github.com/pirl-unc/hitlist/issues/680),
[#705](https://github.com/pirl-unc/hitlist/issues/705) and
[#706](https://github.com/pirl-unc/hitlist/issues/706).
Incomplete lineage need not block presence-level annotation.

## Acceptance criteria

- Every category/input path gets evidence before selection; identical sequences
  and MHC contexts get identical MS features.
- IEDB/CEDAR/supplement overlap, duplicate exports, alternate searches, repeated
  reviews and multiple allele rows do not multiply observations.
- Typed, predicted-match, twenty-donor and monoallelic cases reproduce the
  example; twenty technical runs do not mimic twenty donors.
- A PMID naming two patient alleles counts once; reference remainders stay
  nonnegative. Reused experiments and known shared units merge reference groups.
  Filtering an observation cannot split a known identity group.
- Unknown/pooled lineage preserves presence without manufacturing units.
  Authoritative reuse corrections revise scores with an audit.
- Exact/coarse/missing/conflicting typing remain distinct. Inferred singleton
  sets never become monoallelic; all qualifying alleles remain available.
- No support yields B; zeroed weights/multipliers yield no bonus. Support is
  monotone at fixed identities/policy and C stays in [0, 1]. Missing B is
  rankable only through direct-MS admission.
- Non-MS, bulk-protein, synthetic-standard and negative/unestablished records
  cannot enter positive ligand-MS. Preserve measurement statuses/quality scope.
- Patient prediction plus sequence MS never becomes experimental restriction;
  unsupported source alleles are not failed predictions.
- Aggregation cannot mix B/S from different alleles, count a peptide twice in
  the window's MS increment, or override hard antigen requirements. With no MS,
  the prior per-allele, peptide and window prediction sums remain unchanged.
- Replay reproduces annotations, admission, scores and constructs with data
  removed and predictors disabled; refresh changes provenance.
- Implementation PRs include regression tests, lint/full tests, CLI/report
  smoke, dependency-contract checks, version bump, green CI and deployment
  under AGENTS.md.

## Scientific basis

Monoallelic elution supports allele attribution while retaining engineered
context ([Abelin et al., 2017](https://pubmed.ncbi.nlm.nih.gov/28228285/)).
Multiallelic MS requires attribution/deconvolution; genotype alone does not
establish peptide restriction
([Alvarez et al., 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6885703/)).
The [HLA Ligand Atlas study](https://jitc.bmj.com/content/9/4/e002071.long)
distinguishes observed ligands from predicted allotype associations. These
support the evidence distinctions; none validates the numeric weights.
