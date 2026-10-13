# Experimental unified antigen inputs

Mutation, fusion, CTA, ERV, splice and viral antigens can use the same
`AntigenDataset` API and `--input-antigens` CLI. Each record carries a
`VaccineAntigen`: sequence, targetable intervals, source identity and an
evidence-bearing admission decision. These are experimental research interfaces.

Source adapters or callers supply the admission evidence. A source label or
sequence alone does not establish tumor specificity. Held-out records remain
in the evidence and cannot enter constructs. This input requires resolved
sequences; unresolved hypotheses remain in upstream evidence tables.

## Prepare and predict

Given `mutation` and `viral` objects from your evidence adapters:

```python
from vaxrank import AntigenDataset
from vaxrank.input_scope import InputScope

inputs = AntigenDataset(
    antigens={"mutation-1": mutation, "viral-1": viral},
    scope=InputScope(
        patient_id="patient-001", reference_assembly="GRCh38",
        annotation="ensembl:93", mhc_alleles=("HLA-A*02:01",),
    ),
)
inputs.save("antigens.json")
```

```sh
vaxrank --input-antigens antigens.json \
  --mhc-predictor mhcflurry --ensembl-release 93 \
  --vaccine-type peptide mrna --output-dir design \
  --output-antigens design/antigens.json
```

The annotation describes the host; each antigen keeps its own source species
and accession. New predictions require an installed host annotation for
exact-self assessment. Conflicting reference or genotype declarations fail.
The Python equivalent is `inputs.predict(predictor, genome=genome)`, using a
Topiary predictor. It returns the shared `EpitopeDataset`.

## Combine new sequences and saved evidence

An experimental manifest can mix new antigen sequences with saved mutation,
CTA or viral predictions. Paths are relative to the manifest:

```yaml
schema: vaxrank.input_manifest.v1
patient_id: patient-001
reference_assembly: GRCh38
annotation: ensembl:93
mhc_alleles: [HLA-A*02:01]
inputs:
  - format: epitopes
    path: mutations-and-cta.tsv
  - format: antigens
    path: viral-antigens.json
```

```sh
vaxrank --input-manifest inputs.yaml \
  --mhc-predictor mhcflurry --ensembl-release 93 \
  --vaccine-type peptide mrna --output-dir mixed-design
```

Only the `antigens` entries receive new predictions. Saved measurements stay
unchanged. All declarations are checked before loading models; a manifest may
fill missing scope but cannot override conflicting evidence. Duplicate files
and loading an antigen input alongside its native prediction export are rejected.

When each file already declares its scope, repeat `--external-input
antigens=PATH` or combine `--input-antigens` with `--input-epitopes` instead.
The Python `vaxrank.external_rescoring.prepare_reports` API accepts the same
`("antigens", path)` inputs and an `antigen_predictor_factory(scope)` returning
a Topiary predictor. The factory runs once after validation.

Direct VCF/BAM still runs its source adapter first; add
its saved candidate file to this manifest. Fresh/additive table rescoring
is a separate step. Pooling does not calibrate different models or evidence
types: use an explicit [selection policy](mixed-antigen-design.md).

## Add CTA expression

An experimental `cta_expression` entry admits targets from patient expression
before prediction. Add it beside `antigens` or `epitopes` in the manifest above:

```yaml
  - format: cta_expression
    path: expression.tsv
    sample_id: tumor-rna
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

Use `transcript` for transcript measurements. Each entry declares its own
columns, units, sample and threshold; one matrix may supply several sample
columns. Existing human Ensembl/OncoRef admission rules apply. Expression is
not converted into protein abundance or mutation read support.

Optional `exclusion_policy` lists `gene_patterns` and `gene_exceptions`;
defaults remain `MAGE*` and `MAGEA4`. An optional `hitlist_evidence_bundle` path
is relative to the manifest and is verified against that input's contract.
All input scopes are checked before loading one predictor for new sources.

`--output-dir` writes `cta_admission.json` for one CTA source, or separate files
under `cta_admissions/` for multiple sources. `--output-cta-admission` selects a
single-source file explicitly. Held-out and unresolved decisions remain in
the native evidence, even when a CTA source contributes no admitted proteins
to an otherwise usable mixed run. A wholly unadmitted CTA-only run stops after
writing requested admission records.

## Shared outputs

The ordinary ranking, peptide/mRNA assembly and reports apply to every kind.
Reports distinguish resolved, admitted and held-out records by source kind.
`design/candidate_predictions.tsv` retains predictions, policies and actual
products for replay:

```sh
vaxrank --input-epitopes design/candidate_predictions.tsv --output-dir replay
```

`--output-antigens` exports resolved native antigen records in the input JSON
format. It preserves admission evidence and input history, including held-out
and unmeasured records. It excludes prediction values and design policies;
keep the native candidate file for replay. Unresolved table rows do not become
invented antigen sequences. This export also works with native mutation and
CTA expression workflows.

Pooled native files retain each input's CTA admission decisions and captured
Hitlist evidence in `dataset.selection["source_evidence"]`, including unresolved
features. Prior product audits do not establish an audit for a new mixed design.

See [mixed-source selection](mixed-antigen-design.md) for evidence requirements
and current source-allocation limits.
