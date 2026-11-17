# One vaccine run for mutations, ncORFs, SVs, fusions and CTAs

All supplied antigen sources enter the same prediction, window selection,
ranking, peptide/mRNA assembly and report workflow. Use `--antigen-types` to
choose the types eligible for vaccine products:

```sh
# Ordinary SNVs and indels only
vaxrank ... --antigen-types mutation

# Add the other sources to the same design
vaxrank ... --antigen-types mutation ncORF fusion SV CTA
```

The default includes all supplied types. `ncORF` selects mutations translated
through an ORFdata extra-ORF model, including uORFs and dORFs; it is separate from
ordinary `mutation` allocation. These switches select supplied evidence, not
discoveries to fabricate. `--exclude-antigen-type CTA` and the existing per-type
minimum/maximum options remain available. Identical product sequences share a
slot and retain all source mappings; allocation uses the first ranked source.

## Combine inputs

Use a patient-scoped manifest when combining direct VCF/BAM with other sources.
The reference annotation and patient genotype must agree across the inputs.

```yaml
schema: vaxrank.input_manifest.v1
patient_id: patient-001
sample_id: tumor-rna
reference_assembly: GRCh38
annotation: ensembl:115
mhc_alleles: [HLA-A*02:01]
inputs:
  - format: isovar
    path: ncorfs.json
  - format: isovar
    path: structural-variants.json
  - format: isovar
    path: fusions.json
  - format: cta_expression
    path: expression.tsv
    cta_expression:
      measurement_level: gene
      id_column: gene_id
      value_column: tumor_tpm
      expression_unit: TPM
      evidence_source: Salmon
      evidence_version: "1.10.3"
      assay: bulk RNA-seq
      min_expression: 2
```

An [Isovar input](isovar-inputs.md) can reconstruct directly from BAM/reference
models or load a saved producer result. For a source of called somatic events,
its envelope can simply declare:

```json
{
  "schema": "vaxrank.isovar_input.v1",
  "somatic": true,
  "result": "isovar-result.json"
}
```

The somatic declaration replaces a handwritten source-wide attestation. It
does not establish translation or bypass RNA support, sequence consistency,
targetable-change or unresolved-frame checks. Supported local coding windows
are accepted by default; use `admission.allow_partial: false` to require complete
coding sequences. Distinct supplied ORF models retain separate candidates even
when they translate the same RNA in different frames. Each candidate still
requires a resolved carrier strand unless ambiguity is explicitly enabled.

For ncORFs, a single request applies a pinned ORFdata catalogue to every
overlapping variant in the VCF, so there is no per-variant input to prepare:

```json
{
  "schema": "vaxrank.isovar_input.v1",
  "somatic": true,
  "reconstruct": {
    "family": "variant_orfs",
    "inputs": {
      "vcf": "somatic.vcf",
      "bam": "tumor.bam",
      "catalog": "orfs.json.gz",
      "expected_sha256": "<SHA-256 of the catalogue file>"
    }
  }
}
```

The manifest supplies the RNA sample and annotation. ORFdata validates each
catalogue placement against that annotation; incompatible placements raise an
error rather than transferring translation evidence to a different protein.

## Current CTA sequence policy

Acquire the owner-pinned sequence bundle once, explicitly:

```sh
python -c 'from oncoref import sequence_bundle; print(sequence_bundle.ensure(download=True))'
```

Then select it for the run:

```sh
vaxrank --vcf somatic.vcf --bam tumor.bam \
  --input-manifest inputs.yaml --ensembl-release 115 \
  --mhc-predictor mhcflurry --mhc-alleles HLA-A*02:01 \
  --cta-sequence-bundle current \
  --antigen-types mutation ncORF fusion SV CTA \
  --vaccine-type peptide mrna --output-dir vaccine \
  --output-html-report vaccine/report.html --output-pdf-report vaccine/report.pdf
```

`--cta-sequence-bundle PATH` selects a verified local copy. Loading never downloads
implicitly. The bundle determines the current qualified CTA representatives and
the exact protein-reference identity; the installed annotation must match it.
This policy replaces the historical gene-pattern exclusions for that run.
Historical CTA admission options cannot be mixed into the current bundle policy.
Omitting the bundle option retains the existing historical admission API.

For every predicted CTA peptide, all contiguous reference occurrences are checked.
OncoRef's `exclude` role (Unsafe) takes precedence, followed by `review`, then
`allow`. An allowed CTA occurrence cannot cancel a match in an Unsafe protein.
Missing matches remain unassessed. Matching never crosses ambiguous residues,
stops or record boundaries. Only eligible peptides contribute target score;
the ordinary expression, MHC prediction and ranking requirements still apply.

## What the reports show

The default HTML, PDF and text reports include a compact biological view:

- Fusion/SV partner loci, junction position or split codon, and RNA support.
- ncORF type and source, catalogued translation evidence separately from patient
  RNA support, and the recorded tumor/reference sequence comparison.
- CTA peptide decisions with matching Unsafe/review genes and protein counts
  (exact accessions are retained in the workbook and native evidence), plus
  ranges of tested peptide starts for every assessed source, including sources
  with no ranked window.

Comparisons label an annotation reference, a focal-allele reversion on the same
RNA background, or a retained patient background explicitly. They do not claim
that patient-normal protein was measured. Unavailable comparisons stay missing.
CSV/Excel include the same biological fields; Excel also includes a CTA targeting
sheet. Detailed producer evidence remains in native outputs.

Saved native predictions retain the source evidence, reference identity,
containment decisions and selection configuration. Replaying them does not load
the current OncoRef bundle or reconstruct Isovar sequences again. Changing the
sequence-reference policy requires new CTA inputs:

```sh
vaxrank --input-epitopes vaccine/candidate_predictions.tsv --output-dir replay
```

Reference translation catalogues and exact Unsafe containment do not establish
patient presentation or clinical safety. OncoRef's recorded coverage gaps remain
attached to the native source evidence.
