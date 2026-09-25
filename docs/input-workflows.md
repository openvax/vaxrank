# Choosing and combining inputs

Vaxrank supports direct variant analysis and imported prediction tables as
separate CLI modes. LENS and pVACseq tables can be combined in one run. A run
containing direct, LENS, pVACseq and native Exacto inputs together is not yet
supported.

| Input | Entry point | Current boundary |
| --- | --- | --- |
| Somatic variants + tumor RNA | `--vcf`, `--bam`, `--mhc-predictor`, `--mhc-alleles` | Reconstructs protein context with Isovar and predicts epitopes. Cannot be mixed with report inputs in the same CLI run. |
| LENS report | `--external-input lens=PATH` or `--input-lens PATH` | Imports reported peptide occurrences, available context and evidence. |
| pVACseq report | `--external-input pvacseq=PATH` or `--input-pvacseq PATH` | Imports all-epitope or aggregated TSVs. An aggregated table already selected its best epitopes upstream; importing it cannot recover omitted candidates. |
| Exacto output | No native importer | Requires upstream normalization work; treating an Exacto protein as generic FASTA does not preserve its full provenance. |
| Saved Vaxrank candidate predictions | Python `load_predictions(path)` | Native CSV/TSV reload is a library API; there is no `--input-epitopes` CLI mode. |
| Already constructed `VaccineAntigen` objects | Python `predict_epitopes(..., antigen=...)` and `vaccine_peptides_for_antigen(...)` | Library integration; callers provide context, targetable intervals and admission evidence. |

## Declare the scope of combined reports

The ordinary VCF + BAM command is unchanged. It needs no manifest:

```sh
vaxrank --vcf somatic.vcf --bam tumor.bam \
  --mhc-predictor netmhcpan --mhc-alleles 'HLA-A*02:01,HLA-B*07:02' \
  --output-csv ranked.csv
```

A single `--input-lens` or `--input-pvacseq` report also needs no manifest.
Unstated patient, reference, sample and genotype metadata remain unknown.

To combine reports, declare their **shared patient, reference assembly and
patient MHC allele set**. For example, save this as `inputs.yaml`:

```yaml
schema: vaxrank.input_manifest.v1
patient_id: patient-001
reference_assembly: GRCh38
mhc_alleles: ['HLA-A*02:01', 'HLA-B*07:02']
inputs:
  - format: lens
    path: patient.lens.tsv
    sample_id: tumor-biopsy
    library_id: tumor-rna-1
  - format: pvacseq
    path: patient.all_epitopes.tsv
    sample_id: tumor-biopsy
    library_id: tumor-rna-1
```

Paths resolve relative to the manifest. JSON works too. Shared declarations
can appear at the top level; each input can also declare `patient_id`,
`reference_assembly`, `mhc_alleles`, `annotation`, `sample_id`, `library_id`
and `timepoint`. Sample/library/timepoint defaults can differ per input.
Known patient, reference, annotation and genotype declarations must agree.
Use the exact reference name; Vaxrank does not equate different assembly
labels or lift coordinates between builds. For Ensembl annotations, use
`annotation: ensembl:110` (with the actual release used by the producer).
An explicitly configured `--ensembl-release` must match the declared assembly
and annotation. Without one, a declared assembly still determines variant
identity; it does not invent a producer annotation release.

Unknown fields may be omitted or set to `null`. Combining inputs requires
known patient, assembly and genotype for every input. Other missing fields
remain unknown in saved provenance; their compatibility has not been
established. A manifest is a user assertion, not independent verification of
the patient's identity or typing assay. `--output-patient-id` remains a display
label and cannot authorize pooling.

Producer tables may carry declarations in columns with those exact names
(`mhc_alleles` uses a JSON list in a cell). Vaxrank checks every row, including
LENS `Hsap37.`/`Hsap38.` origin markers. Conflicting producer declarations
cannot be overridden by a manifest. Reports that already declare all required
scope can still use repeated `--external-input`; otherwise use the manifest.
List every input in the manifest rather than mixing it with individual input
flags. Conflicts fail before scoring, live model initialization or exports.

Report allele coverage may be a **subset** of the declared genotype. A union
of observed prediction alleles is never promoted into a declared genotype.
Alleles are normalized with MHCgnomes; differing resolutions are not guessed
to be equivalent. Reported alleles outside the declared set are rejected.

## Rank using original predictions

No predictor installation or raw sequencing inputs are needed to reuse the
table values. Use the manifest above for multiple files:

```sh
vaxrank \
  --input-manifest inputs.yaml \
  --external-predictions input \
  --output-csv original-ranked.csv \
  --output-input-predictions originals.tsv
```

`--output-csv` writes a per-peptide/allele report in external mode;
`--output-neoepitope-report report.xlsx` writes its Excel counterpart.

`input` is the default. It preserves historical prediction values and applies
Vaxrank's configured filters/scores; it does **not** reproduce each producer's
original ranking algorithm. Model selection is resolved within each source
table, so the pooled `rank` is conditional on that scoring policy. The report
also provides `source_rank`. Different model scores are not automatically
calibrated onto one scale.

Do not pass `--mhc-predictor`, `--mhc-alleles`, or `--mhc-alleles-file` in this
mode: predictions retain their recorded alleles. Declare the patient genotype
in the manifest. To compute predictions for additional patient alleles, use
`fresh` mode below.

Use the Topiary DSL through `epitopes.filter_expr` and `epitopes.score_expr`
to choose which original signals drive ranking. For example, this adds an
expression threshold to a LENS run whose normalized table contains `gene_tpm`:

```sh
vaxrank --input-lens patient.lens.tsv \
  --config-text 'epitopes.filter_expr=(affinity.value < 500) & (gene_tpm > 1)' \
  --config-text 'epitopes.score_expr=affinity.logistic_normalized(350, 150)' \
  --output-csv expression-ranked.csv
```

Unqualified `affinity` uses the selected default method. Use a qualified
reference such as `affinity[mhcflurry]` when that method is present in every
table being scored. Check feature availability and missing values in each
source before applying a shared policy; source-specific signals are not
guaranteed to exist in every format.

### Construct selection and final order

`vaccine_peptides.combined_score_expr` determines final construct order for
direct, LENS, pVACseq and mixed external inputs. The first selected vaccine
peptide represents each variant or antigen. Constructs rank by descending
combined score; this order feeds the template reports and the peptide/mRNA
assemblers before optional HLA-coverage selection.

Occurrence/window selection happens first and remains source-specific:

- Direct inputs generate windows, apply `score_fraction_of_best`, and choose
  among the retained windows using `vaccine_peptides.ranking_rules`.
- Within each external report, the strongest eligible epitope selects its
  reported occurrence/context. Only compatible epitopes from that context
  enter the construct. A combined-score expression does not rescan windows
  or recover candidates omitted by the producer.
- When multiple reports supply a construct for the same variant or antigen,
  the shared final ranking policy chooses one. Counts and predictions stay
  with that source-derived construct; observations are not summed or blended.

The existing `require_target_epitopes_in_variant` setting applies to external
inputs too. Its default excludes constructs with no target epitopes before
final ranking, including known-self-only constructs. Their input observations
remain in the audit report. Set the option to `false` only when intentionally
retaining constructs without target epitopes.

Exact combined-score ties use descending RNA support only when **all tied
representatives have counts with the same stated unit and derivation**.
Otherwise that tie skips RNA and uses descending target-epitope score.
Complete ties preserve input order, including the supplied file order for
repeated-source alternatives. Missing RNA is not assigned a count by this
ranking step, and reads are not converted to fragments.

Common ordering does not make evidence from different sources comparable.
The existing mutation default is `sqrt(n_rna_alt) * target_epitope_score`;
source-agnostic antigens without mutation counts default to
`target_epitope_score`. The mutation DSL retains its legacy numeric-zero
binding for unavailable counts; provenance still distinguishes missing from
measured zero. For a shared score that does not weight RNA, use:

```sh
--config-text 'vaccine_peptides.combined_score_expr=target_epitope_score'
```

Custom expressions using unavailable mutation fields fail for source-agnostic
antigens. The external per-epitope/allele CSV and Excel `rank`/`source_rank`
columns remain epitope ranks, separate from construct order.

## Generate common predictions and re-rank

Request fresh predictions explicitly, with installed models and the desired
patient HLA set:

```sh
vaxrank \
  --input-manifest inputs.yaml \
  --external-predictions fresh \
  --mhc-predictor mhcflurry \
  --mhc-alleles 'HLA-A*02:01,HLA-B*07:02' \
  --config-text 'epitopes.filter_expr=affinity[mhcflurry].value < 500' \
  --config-text 'epitopes.score_expr=affinity[mhcflurry].logistic_normalized(350, 150)' \
  --output-csv fresh-ranked.csv \
  --output-input-predictions originals.tsv \
  --output-epitopes fresh.tsv
```

This predicts the **reported peptides**, with their available flanks, on the
requested HLA set, which must fit the declared genotype when one is available.
It can add peptide-HLA pairs, but does not scan new peptide
windows, recover omitted pVACseq candidates, or extend protein context.
Haplotype-scoped predictors are not supported by this external rescoring path.

The active `CandidateEpitope` objects contain fresh predictions. Original
values appear separately under `Input ...` report columns for matching
alleles, and under `input_...` columns in the in-memory scoring frame. Keep
`originals.tsv` as well: original alleles outside the fresh HLA set do not get
rows in the fresh report. This is not yet a single lossless file containing
all original evidence plus all new predictions and ranking policies.

## Construct outputs and provenance

Create the destination with `mkdir -p vaccines`, then add
`--vaccine-type peptide mrna --output-dir vaccines --ensembl-release N` to
construct admitted antigens. Replace `N` with the release matching the input
annotation: `--output-dir` automatically requests ASCII/PDF reports, which
require it. Explicit ASCII/HTML/PDF output flags have the same requirement.

External-report construction uses the shared linker without new junction
predictions by default. To optimize junctions while retaining every historical
candidate score, add `--mrna-junction-predictor mhcflurry` (or another
mhctools model that supplies percentile ranks). The manifest's declared
genotype supplies the query alleles; `--mrna-junction-alleles 'HLA-A*02:01'`
can select a subset. For a single report with unknown genotype, an explicit
junction allele set is required. Report prediction coverage is not used as a
genotype declaration.

Cached JSON rerendering (`--input-json-file`) can also request new junction
predictions by supplying both `--mrna-junction-predictor` and
`--mrna-junction-alleles`; its default performs no new junction prediction.

`--mrna-optimize-linkers` explicitly requires a configured junction model;
`--mrna-no-optimize-linkers` disables junction queries even if a model is
configured. Candidate `--mhc-predictor` settings still apply only to explicit
`--external-predictions fresh` requests. Executable and weights paths can be
set independently with `--mrna-junction-predictor-path` and
`--mrna-junction-predictor-models-path`. The corresponding YAML keys live under
`mrna:` (`junction_predictor`, `junction_alleles`, and the two path keys);
`optimize_linkers: null` selects auto mode.

VCF/BAM auto mode reuses a single candidate predictor. With multiple candidate
predictors, it keeps the shared linker and logs that junction prediction is
disabled. Select `--mrna-junction-predictor` to optimize those junctions;
an explicit `--mrna-optimize-linkers` request requires that choice.

The mRNA manifest's `elements.junction_swap` records the selected policy and
whether optimization ran. Its `prediction` record preserves observed model
names/versions/kinds, query alleles and lengths, and the peptide queries for
each candidate linker at each junction. Unknown model versions remain null.
This metadata also survives native RNAConstruct serialization. Optional
cleavage annotation can be disabled with `--no-processing-aware-annotation`.
Creating the directory first
also avoids the direct pipeline's pending [output-path fix](https://github.com/openvax/vaxrank/pull/456).

LENS supplies reported context windows. pVACseq tables generally supply only
the epitope; Vaxrank does not retrieve the separate mutant/wild-type FASTA
automatically or invent longer flanks. See the
[pVACseq output reference](https://pvactools.readthedocs.io/en/latest/pvacseq/output_files.html)
for the distinction between these artifacts.

Both readers produce `CandidateEpitope` objects containing sequence/context,
prediction identity, allele/model/version predictions and known comparators.
`ExternalReport` retains normalized source rows for scoring and construction.
Reports carry input path, format, content hash and resolved scope. Original
producer declarations and manifest assertions are retained separately.
Candidate occurrence identities include scope; sequence hashes link exact
peptide/window matches without claiming full ORF equivalence.

Native `--output-epitopes` and `--output-input-predictions` files preserve each
candidate's `input_provenance` and `input_evidence`, including Topiary's
normalized RNA measurements, stated units/derivation and source-specific
evidence. `load_predictions(path)` restores those fields. The per-epitope CSV
also carries scope columns and `input_provenance_json`. External output
directories include `candidate_predictions.tsv` and `input_provenance.json`;
the run summary and template reports distinguish declared genotype from
reported prediction coverage. These provenance files do not replace the
original input reports or provide the generalized table-reload CLI in #346.

Different samples, libraries and timepoints stay separate observations. One
representative construct may be selected for their shared patient/reference
event, but its counts come from that observation. RNA support is never summed
between reports, even when libraries or read sets may overlap. All imported
candidate observations remain available in the native prediction export.

Construction still uses format-specific adapters and selects one
source-derived window per grouped source. It does not merge RNA counts or
assemble a consensus ORF across reports. LENS splice, CTA/self and ERV
construction also requires explicit category opt-ins; see
[antigen inputs](https://github.com/openvax/vaxrank/blob/main/README.md#upstream-inputs).

## Remaining integration work

- [Vaxrank #497](https://github.com/openvax/vaxrank/issues/497) and
  [Topiary #366](https://github.com/openvax/topiary/issues/366): common input
  normalization, additive predictions and evidence-preserving reload.
- [Topiary #365](https://github.com/openvax/topiary/issues/365): native Exacto.
- [Topiary #370](https://github.com/openvax/topiary/issues/370): ORF/occurrence
  identity and reconciliation across source observations.

The [unified evidence design](unified_evidence.md) describes the intended
model, not an additional supported CLI workflow. For current configuration
keys and defaults, run `vaxrank --print-default-config`.
