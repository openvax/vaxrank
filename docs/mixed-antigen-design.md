# Experimental mixed-antigen design

**Experimental:** CTA/viral integration and mixed-source selection are research
workflows. Admission and ranking do not establish clinical safety or efficacy.

## What information is needed?

| Stage | Minimum information |
| --- | --- |
| Consider a candidate | Patient/sample, sequence or unresolved sequence hypothesis, source identifiers and evidence provenance. Incomplete candidates remain report-only. |
| Admit to a construct | Exact amino-acid sequence, source coordinates, targetable intervals and a documented tumor-specificity decision. Held-out candidates cannot enter constructs. |
| Rank sources together | Same patient and declared host genotype/reference; compatible predictor versions, endpoints and units; explicit scoring, missing-data and duplicate policies. |
| Allocate by source | Eligibility rules and slot limits for mutation, CTA and viral candidates, plus product size and allele-coverage requirements. |

- **Mutation:** somatic-call/reference provenance, mutant sequence and available
  RNA evidence. Missing RNA is not zero expression.
- **CTA:** gene/transcript identity, versioned OncoRef reference and patient tumor
  expression with assay, units and threshold. Gene RNA does not establish protein
  abundance. See the [CTA input guide](design/expression-cta-inputs.md).
- **Viral:** virus/strain, accession/version, ORF sequence and evidence linking its
  expression to this tumor. Native `VaccineAntigen` evidence must carry the
  targetable mask and admission decision. A FASTA or infection history alone is
  insufficient; automatic viral admission is not available.

## Rank together or select separate panels

For resolved sequences without predictions, use the experimental
[unified antigen input](antigen-inputs.md). It accepts all supported kinds
in one dataset and predicts them through the shared pipeline.

Save candidates with `--output-epitopes`, keep their `.references` directories,
and combine them using an [input manifest](input-workflows.md). Preserve viral
accessions separately from the manifest's host reference.

For inputs with compatible saved MHCflurry affinity measurements, this explicitly
uses binding alone for ranking:

```sh
vaxrank --input-manifest inputs.yaml \
  --config-text 'epitopes.score_expr=affinity["mhcflurry"].value.logistic_normalized(350, 150)' \
  --config-text 'vaccine_peptides.combined_score_expr=target_epitope_score' \
  --vaccine-type peptide --output-dir mixed \
  --output-epitopes mixed/native.tsv --no-processing-aware-annotation
```

Select a recorded version with `affinity["mhcflurry", "VERSION"]` when needed.
For incompatible scores, use
[common-model prediction](input-workflows.md#generate-common-predictions-and-re-rank).
Binding, presentation and immune response remain distinct; this score does not
infer missing expression or tumor specificity.

For separate source allocations, filter the saved pool and cap each panel:

```sh
vaxrank --input-epitopes mixed/native.tsv \
  --config-text 'epitopes.filter_expr=vaxrank_antigen_kind == "CTA"' \
  --vaccine-type peptide --peptide-antigens-per-construct 1 \
  --peptide-max-constructs 2 --output-dir cta-stratum \
  --output-epitopes cta-stratum/native.tsv --no-processing-aware-annotation
```

Use `"mutation"` or `"viral"` for the other panels. These source fields come from
native antigen evidence. Compose additional eligibility rules into the same
filter; it replaces the inherited scalar filter. Named policies use `filter_by`.
Caps are maxima, not guaranteed slots. Review duplicates, allele coverage and
total capacity across panels. **Joint source quotas inside one mRNA construct
are not implemented.**

Saved CTA/viral reports need no Ensembl release; mutation reports still need
annotation when it was not retained. Optional [window exposure](core-window-exposure.md)
is also experimental and does not change selection.
