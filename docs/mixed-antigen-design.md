# Mutation, CTA and viral candidates in one vaccine design

Vaxrank can assemble admitted native mutation, CTA and viral antigens together.
Admission, ranking and allocating product slots answer different questions.
A source label or a binding prediction alone does not establish tumor specificity.

## Information needed at each stage

| Decision | Shared requirements | CTA-specific evidence | Viral-specific evidence |
| --- | --- | --- | --- |
| Retain a hypothesis for consideration | Patient/sample identity; exact sequence or an explicit unresolved sequence hypothesis; source identifiers and evidence provenance. Incomplete hypotheses may remain report-only. | Gene/transcript and annotation identity; CTA reference membership with its version. | Virus species/type/strain, accession/version, ORF/protein and patient sequence differences where known. |
| Admit an antigen to construct selection | Actual translated amino-acid sequence, occurrence coordinates, explicit targetable intervals and an evidence-bearing admission decision. | Patient tumor expression, declared measurement level/units/assay and threshold policy; OncoRef identity and admission result. Gene abundance is not isoform or protein abundance. | Evidence linking the expressed viral ORF to this tumor, with sample/assay and decision provenance. A reference FASTA, infection history or sequence match alone does not supply this decision. |
| Co-rank with mutation candidates | Same patient and declared host genotype/reference; comparable model, version, endpoint, units and sequence context; explicit missing-data, duplicate and combined-score policies. | Keep expression and normal-tissue restriction separate from binding/presentation scores. | Keep tumor association and viral expression separate from binding/presentation scores. |
| Select across source strata | Explicit allocation of slots, per-source eligibility/ranking, product modality and size limits; an approach to shared sequences and allele coverage. | State whether CTA slots are reserved, capped or compete freely. | State whether viral slots are reserved, capped or compete freely. |

For mutation candidates, retain somatic call and reference provenance, mutant
sequence/context and available RNA evidence. Missing RNA support is not a measured
zero. For every source, normal-tissue/self sharing, presentation-pathway integrity,
processing and T-cell recognition are separate assessments, not consequences of
passing the admission or rank gate.

The expression-first [CTA workflow](design/expression-cta-inputs.md) implements
OncoRef-backed human CTA admission. Its expression threshold and gene exclusions
are explicit design policies. Other species require their own evidence contract;
removing the human reference guard would not establish one (see
[#596](https://github.com/openvax/vaxrank/issues/596)).

Viral inputs currently enter through native `VaccineAntigen` evidence or the
`vaxrank_antigen_json` column of a Topiary table. Supply `kind="viral"`, the
sequence, `TargetableMask`, source species/identifiers/metadata and
`TumorSpecificityAttestation` with its `TumorSpecificityEvidence` records.
Use `held_out` while evidence is missing; an `admitted` or explicitly justified
`overridden` decision is an assertion by the producing workflow, not a biological
conclusion Vaxrank derives from the label. Keep strain/accession, sample, assay,
evidence version and rationale in the record. There is no automatic viral FASTA
admission or complete oncovirus catalogue in this release
([#258](https://github.com/openvax/vaxrank/issues/258)).

## A shared ranking

First save each design's candidate evidence using `--output-epitopes`; preserve
adjacent `.references` directories. Use all candidate observations, not only
previously selected products. Declare common scope in a manifest:

```yaml
schema: vaxrank.input_manifest.v1
patient_id: patient-001
reference_assembly: GRCh38
mhc_alleles: [HLA-A*02:01, HLA-B*07:02]
inputs:
  - {format: epitopes, path: mutation.tsv, sample_id: tumor-rna}
  - {format: epitopes, path: cta.tsv, sample_id: tumor-rna}
  - {format: epitopes, path: viral.tsv, sample_id: tumor-rna}
```

`reference_assembly` is the shared **host** coordinate reference, not the viral
accession. Keep the viral reference with the antigen's own provenance. Known
mutation annotation can come from native evidence; annotation-free CTA/viral
reports do not require `--ensembl-release`. Mutation reports without retained
annotation still require a suitable configured reference. Custom retained
mutation references are validated against their actual saved identity.

When every input already carries compatible MHCflurry affinity measurements:

```sh
vaxrank --input-manifest inputs.yaml \
  --config-text 'epitopes.score_expr=affinity["mhcflurry"].value.logistic_normalized(350, 150)' \
  --config-text 'vaccine_peptides.combined_score_expr=target_epitope_score' \
  --vaccine-type peptide --output-dir mixed \
  --output-epitopes mixed/native.tsv --no-processing-aware-annotation
```

Select a specific recorded predictor version as well when inputs differ; the
DSL accepts `affinity["mhcflurry", "VERSION"]`. Binding nM, presentation scores,
percentiles and observed MS presentation remain different measurements. The
formula above is an explicit binding-only design score, not a probability of
immune response. It removes the mutation default's RNA-count multiplier from
this comparison; it does not normalize different candidate universes or infer
missing expression. Shared sequences and competing source occurrences use the
existing Topiary representative/duplicate policy.

For incompatible or missing model results, use the documented
[additive prediction path](input-workflows.md#saved-scoring-policies)
with a uniquely named run and select that run's columns in the DSL. MHCtools
supplies the models and Topiary retains their evidence. Additive prediction
scores existing candidate/allele occurrences; it does not discover new ORFs,
scan absent peptides, or create missing HLA assignments. Ensure those candidates
exist upstream. Inspect missingness before interpreting a pooled rank.

The manifest's genotype is a declaration. Observed prediction alleles can be a
subset; they do not establish a complete genotype. The assembler may reorder
ranked candidates for allele coverage, so inspect actual products as well as
the epitope ranking table. Peptide and mRNA constructs have independent limits.

## Source strata without inventing cross-source weights

Native scoring and audit tables expose `vaxrank_antigen_kind` and
`vaxrank_admission_status`, derived from the saved typed antigen objects.
These consumer fields do not trust a producer's similarly named columns.
Missing construction evidence stays unknown. Admission is still enforced after
DSL scoring: a high-scoring held-out candidate cannot enter a construct.
Legacy report-only rows without typed antigens can have unknown source fields;
use their documented adapter fields or enrich construction evidence first.

If the scores are not comparable, select separate panels with the same saved
pool, explicit source filters and product caps. For example, reserve at most two
single-antigen peptide products for CTAs:

```sh
vaxrank --input-epitopes mixed/native.tsv \
  --config-text 'epitopes.filter_expr=vaxrank_antigen_kind == "CTA"' \
  --vaccine-type peptide --peptide-antigens-per-construct 1 \
  --peptide-max-constructs 2 --output-dir cta-stratum \
  --output-epitopes cta-stratum/native.tsv --no-processing-aware-annotation
```

Use `"mutation"` or `"viral"` in separate runs with independently chosen caps.
Compose any additional eligibility predicates into the same filter; setting
`filter_expr` replaces an inherited scalar filter. A named selection policy
uses its own `filter_by` instead; do not also supply a scalar filter.
The cap is a maximum, not a guarantee: insufficient admitted candidates leave
unfilled slots. Review the resulting panels together for repeated sequences,
shared allele coverage and total capacity. These are independently selected
panels; Vaxrank does not yet optimize source minimum/maximum quotas jointly
inside one concatenated mRNA construct. Free competition in a pooled run and
separate source allocations are distinct design policies.

## Scientific interpretation and audit limits

CTA reference membership cannot substitute for patient tumor evidence. The
original NY-ESO-1 discovery found restricted normal-tissue transcription
([Chen et al., 1997](https://pubmed.ncbi.nlm.nih.gov/9050879/)); a subsequent
breast-tumor study found that transcript detection often lacked detectable
protein ([Sugita et al., 2004](https://aacrjournals.org/cancerres/article/64/6/2199/512367/NY-ESO-1-Expression-and-Immunogenicity-in)).
This motivates retaining measurement level and assay instead of promoting RNA
measurements into protein or presentation evidence.

HPV16 E6/E7 peptide vaccination was studied in patients with HPV16-positive
lesions ([Kenter et al., 2009](https://doi.org/10.1056/NEJMoa0810097)). Requiring
an explicit link between a viral target and the patient's tumor is a workflow
policy informed by that setting, not a claim that any detected viral sequence
is an eligible cancer target.

MHCflurry's presentation model combines binding and processing information;
its validation on ligand datasets does not establish T-cell recognition
([O'Donnell et al., 2020](https://doi.org/10.1016/j.cels.2020.06.010)). Preserve
these distinctions in both ranking and reports. Saved admission and construct
replay provide traceability; automatic final-product safety enforcement remains
separate work ([#311](https://github.com/openvax/vaxrank/issues/311)).
