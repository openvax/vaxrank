# Experimental unified antigen inputs

Mutation, fusion, CTA, ERV, splice and viral antigens can use the same
`AntigenDataset` API and `--input-antigens` CLI. Each record carries a
`VaccineAntigen`: sequence, targetable intervals, source identity and an
admission decision and its evidence.

Callers supply resolved sequences and admission evidence. Held-out records stay
in the evidence but cannot enter constructs; unresolved sequences stay upstream.

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

A manifest can mix new sequences with saved predictions. Paths are relative
to the manifest:

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

Add a `cta_expression` entry beside `antigens` or `epitopes` to select CTA targets
from patient expression before prediction:

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

`format: cta_expression` selects the CTA expression adapter. The file can be CSV
or TSV; the nested `cta_expression` settings describe its columns and units.
Use `transcript` for transcript measurements. Each entry declares its own
columns, units, sample and threshold; one matrix may supply several sample
columns. Existing human Ensembl/OncoRef admission rules apply. Expression is
not converted into protein abundance or mutation read support.

Optional `exclusion_policy` lists `gene_patterns` and `gene_exceptions`;
defaults remain `MAGE*` and `MAGEA4`. An optional `hitlist_evidence_bundle` path
is relative to the manifest and is verified against that input's contract.

`--output-dir` writes `cta_admission.json` for one CTA source, or separate files
under `cta_admissions/` for multiple sources. `--output-cta-admission` selects a
single-source file explicitly. Held-out and unresolved decisions remain in
the native evidence, even when a CTA source contributes no admitted proteins
to an otherwise usable mixed run. A wholly unadmitted CTA-only run stops after
writing requested admission records.

## Declare the host species (experimental draft)

Declare a host by name. For an ordinary mouse, this is enough:

```yaml
host_species: mouse
```

For an HLA mouse, add the intended MHC species:

```yaml
host_species: mouse
mhc_species: human
```

Only use a list when targeting multiple MHC species: `mhc_species: [human, mouse]`.
Both fields accept common names, scientific names or NCBI taxonomy IDs;
`host_taxon_id` remains an alias. Equivalent declarations normalize to the same
stored identity. The CLI equivalent is `--host-species mouse --mhc-species human`.
The API accepts the same declarations through `InputScope.from_dict(...)`.

The host must agree with annotation; declared and observed alleles must match
the intended MHC species. Antigen source species remain separate, allowing
viral sequences in a human or mouse design. Declarations survive native replay.
For fresh mouse predictions, select mouse annotation with
`--genome GRCm38 --ensembl-release 93` (or another matching build/release).

Producer tables may leave host cells blank. Integer-valued numeric taxon cells
are normalized when read; fractional IDs and booleans are rejected. Manifest
and API taxon IDs must be positive integers.

Donor/transgene sequences and extra ORFs still need explicit self-reference
integration, tracked in [#625](https://github.com/openvax/vaxrank/issues/625).
Current CTA expression admission remains human-only. Inputs without new scope
fields keep their historical identities; known host annotation still enables
species checks.

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
keep the native candidate file for replay. This export also supports native mutation and CTA expression workflows.

Pooled native files retain each input's CTA admission decisions and captured
Hitlist evidence in `dataset.selection["source_evidence"]`, including unresolved
features. A new mixed design requires new product audits.

See [mixed-source selection](mixed-antigen-design.md) for evidence requirements
and current source-allocation limits.
