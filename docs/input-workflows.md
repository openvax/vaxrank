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

Use reports from the **same patient and compatible reference/annotation**.
Vaxrank does not currently validate this across files. Report allele sets can
be partial prediction coverage; their union is not an independently verified
patient genotype. `--output-patient-id` labels the run, rather than checking
input identity.

## Rank using original predictions

No predictor installation or raw sequencing inputs are needed to reuse the
table values. Repeat `--external-input` for multiple files of either format:

```sh
vaxrank \
  --external-input lens=patient.lens.tsv \
  --external-input pvacseq=patient.all_epitopes.tsv \
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
mode: predictions retain their recorded alleles. To compute predictions for a
different HLA set, use `fresh` mode below.

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
guaranteed to exist in every format. The separate
`vaccine_peptides.combined_score_expr` computes construct scores. Currently,
external final ordering still uses target-epitope score, while the direct
pipeline orders by combined score. A custom combined-score expression alone
will therefore not reorder imported constructs; this inconsistency is tracked
in [#506](https://github.com/openvax/vaxrank/issues/506).

## Generate common predictions and re-rank

Request fresh predictions explicitly, with installed models and the desired
patient HLA set:

```sh
vaxrank \
  --external-input lens=patient.lens.tsv \
  --external-input pvacseq=patient.all_epitopes.tsv \
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
requested HLA set. It can add peptide-HLA pairs, but does not scan new peptide
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

For original-only mRNA construction, also pass `--mrna-no-optimize-linkers`.
Junction optimization is a separate prediction step and otherwise attempts to
use a live MHC model, even in `input` mode. Optional cleavage annotation can be
disabled with `--no-processing-aware-annotation`. Choosing an independent
junction model while preserving historical candidate scores needs
[#507](https://github.com/openvax/vaxrank/issues/507). Creating the directory first
also avoids the direct pipeline's pending [output-path fix](https://github.com/openvax/vaxrank/pull/456).

LENS supplies reported context windows. pVACseq tables generally supply only
the epitope; Vaxrank does not retrieve the separate mutant/wild-type FASTA
automatically or invent longer flanks. See the
[pVACseq output reference](https://pvactools.readthedocs.io/en/latest/pvacseq/output_files.html)
for the distinction between these artifacts.

Both readers produce `CandidateEpitope` objects containing sequence/context,
prediction identity, allele/model/version predictions and known comparators.
`ExternalReport` retains normalized source rows for scoring and construction.
Reports carry input path, format and content hash; sequence hashes link exact
peptide/window matches without claiming full ORF equivalence. Gene,
transcript, RNA and antigen evidence remain source-dependent. They are not
all stored on the candidate object or guaranteed to survive native reload.

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
- [Vaxrank #505](https://github.com/openvax/vaxrank/issues/505): validate input
  patient/reference/genotype scope before pooling.

The [unified evidence design](unified_evidence.md) describes the intended
model, not an additional supported CLI workflow. For current configuration
keys and defaults, run `vaxrank --print-default-config`.
