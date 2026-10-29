# Choosing and combining inputs

Vaxrank accepts variants, scored reports and saved evidence. Experimental antigen
and CTA expression inputs use the same scoring, selection and construct outputs.

| Input | Entry point | Supplies |
| --- | --- | --- |
| Somatic variants + tumor RNA | `--vcf`, `--bam`, `--mhc-predictor`, `--mhc-alleles` | Isovar protein context and new predictions |
| LENS report | `--input-lens PATH` | Reported peptides, context and evidence |
| pVACseq report | `--input-pvacseq PATH` | All-epitope or aggregated TSVs; omitted candidates cannot be recovered |
| Topiary table | `--input-topiary FILE` | Normalized predictions, metadata and features |
| Native Exacto table | Manifest `format: exacto` | Peptide, ORF and RNA evidence; requires an explicit sample |
| Saved Vaxrank candidates | `--input-epitopes FILE` | Native candidates; enriched files also retain evidence and policies |
| Resolved antigens (experimental) | `--input-antigens JSON` or manifest `format: antigens` | Sequences, targetable intervals and admission evidence |
| CTA expression (experimental) | `--input-cta-expression FILE` or manifest `format: cta_expression` | Patient expression with [declared columns, units and threshold](antigen-inputs.md#add-cta-expression) |

Table-only candidates can be ranked. Construction also requires sequence context
and admission evidence. Unresolved sequence hypotheses stay in upstream tables.

## Flag reference

| Flag | Purpose |
|---|---|
| `--input-manifest PATH` | YAML/JSON input list with shared and per-input scope |
| `--external-input FORMAT=PATH` | Repeat for `lens`, `pvacseq`, `topiary`, `epitopes` or `antigens` files that already declare compatible scope |
| `--duplicate-candidates error\|best\|worst` | Choose a representative among duplicate observations; default requires scores to agree |
| `--external-predictions input\|fresh\|additive` | Reuse predictions (default), replace LENS/pVACseq predictions, or add named features |
| `--output-input-predictions PATH` | Save original predictions separately |

A manifest replaces individual input flags. Direct `--vcf`/`--bam` inputs can
still be added alongside it.

## Declare the scope of combined reports

A single VCF/BAM run or imported report needs no manifest:

```sh
vaxrank --vcf somatic.vcf --bam tumor.bam \
  --mhc-predictor netmhcpan --mhc-alleles 'HLA-A*02:01,HLA-B*07:02' \
  --output-csv ranked.csv
```

Combining reports requires a shared **patient, reference assembly and MHC allele
set**. Save this as `inputs.yaml`:

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

`format` selects the input adapter, not the file extension. Paths are relative
to the manifest. Shared fields can appear at the top level or on each input:

- `patient_id`, `reference_assembly`, `annotation` and `mhc_alleles` must agree
  wherever declared. Use the exact assembly name and, if known,
  `annotation: ensembl:110` with the producer's release.
- `sample_id`, `library_id` and `timepoint` may differ per input.
- `rna_evidence_extent` optionally declares `whole_library`, `regional_subset`
  or `targeted` evidence per input. This is a user declaration, not an inference
  from the filename or read count. Omit it when the extent is unknown.
- Optional unknown fields may be omitted or `null`. `--output-patient-id` is
  only a display label, not a declaration for combining evidence.

An explicit `--ensembl-release` must match the declared assembly and annotation.
Assembly alone does not imply the producer's annotation release. Producer
columns with these names are checked on every row, including LENS assembly
markers; manifest values cannot override conflicts. An `mhc_alleles` table cell
uses a JSON list. Conflicts fail before models or exports run.

Report alleles may cover a subset of the declared genotype. MHCgnomes normalizes
alleles; Vaxrank rejects reported alleles outside that set and does not infer a
genotype from report coverage.

## Rank using original predictions

```sh
vaxrank --input-manifest inputs.yaml --external-predictions input \
  --output-csv original-ranked.csv --output-input-predictions originals.tsv
```

`input` is the default. It applies Vaxrank's filters and scoring to saved
measurements, retaining both pooled `rank` and within-source `source_rank`.
Different models are not automatically calibrated. CSV and
`--output-neoepitope-report report.xlsx` contain per-peptide/allele results.

Saved-report-only runs need no predictor or prediction allele flags. New
`antigens` or `cta_expression` entries still require a predictor and host
annotation; only those entries receive new predictions.

Choose ranking signals with the Topiary DSL:

```sh
vaxrank --input-lens patient.lens.tsv \
  --config-text 'epitopes.filter_expr=(affinity.value < 500) & (gene_tpm > 1)' \
  --config-text 'epitopes.score_expr=affinity.logistic_normalized(350, 150)' \
  --output-csv expression-ranked.csv
```

Unqualified `affinity` uses each source's selected method. A qualified reference
such as `affinity[mhcflurry]` requires that method in every scored table. Check
that the chosen features exist across inputs.

### Construct selection and final order

`vaccine_peptides.combined_score_expr` orders constructs across direct and
imported inputs, before optional HLA-coverage selection:

- Direct inputs generate windows, apply `score_fraction_of_best`, then select
  with `vaccine_peptides.ranking_rules`.
- External inputs select the reported context containing their strongest
  eligible epitope. Only compatible epitopes from that context enter a construct.
- Where sources offer the same variant or antigen, final ranking selects one
  source's construct. Counts and predictions remain with that source.

By default, `require_target_epitopes_in_variant` excludes constructs with no
target epitopes, including known-self-only constructs. Their observations remain
in the reports.

The mutation default is `sqrt(n_rna_alt) * target_epitope_score`. Antigens without
mutation counts default to `target_epitope_score`. To use binding alone across
sources:

```sh
--config-text 'vaccine_peptides.combined_score_expr=target_epitope_score'
```

The mutation DSL retains its legacy numeric zero for unavailable counts;
provenance distinguishes missing from measured zero. Custom expressions using
unavailable mutation fields fail for other antigen kinds.

Exact score ties use RNA support only when all tied counts share stated units
and derivation, then target-epitope score, then input order. Counts are never
summed across discoveries or converted from reads to fragments. CSV/Excel
`rank` columns are epitope ranks, separate from construct order.

## Generate common predictions and re-rank

```sh
vaxrank --input-manifest inputs.yaml --external-predictions fresh \
  --mhc-predictor mhcflurry --mhc-alleles 'HLA-A*02:01,HLA-B*07:02' \
  --config-text 'epitopes.filter_expr=affinity[mhcflurry].value < 500' \
  --config-text 'epitopes.score_expr=affinity[mhcflurry].logistic_normalized(350, 150)' \
  --output-csv fresh-ranked.csv --output-input-predictions originals.tsv \
  --output-epitopes fresh.tsv
```

`fresh` is supported for LENS/pVACseq. It predicts reported peptides with their
available flanks on the requested allele set, within any declared genotype. It
can add peptide/allele pairs but cannot recover omitted peptides or extend
context. Haplotype models retain their genotype and deconvolved presenter;
allele-free models retain their scope.

Fresh predictions populate the active candidates. Matching original values
appear in `Input ...` report columns and `input_...` scoring features. Keep the
separate original export: alleles outside the fresh set have no fresh report row.

## Combine direct VCF/BAM with tables

```sh
vaxrank --vcf somatic.vcf --bam tumor.bam \
  --mhc-predictor netmhcpan --mhc-alleles 'HLA-A*02:01,HLA-B*07:02' \
  --input-manifest inputs.yaml --duplicate-candidates best \
  --output-epitopes combined.tsv --output-csv ranked.csv
```

Shared manifest declarations also apply to direct inputs. An optional `direct:`
object sets their sample/library/timepoint. Reference and allele compatibility
are checked before prediction.

Predictor, Isovar, DNA fallback and cache settings apply to direct candidates;
imported tables retain their measurements unless additive scoring is requested.
Direct fragments, counts, comparators and reference classifications are saved
before window selection. Reload `combined.tsv` with `--input-epitopes` to repeat
scoring and construction without VCF/BAM or models. Use the same vaccine
configuration to reproduce window selection.

## Construct outputs and provenance

Add `--vaccine-type peptide mrna --output-dir vaccines` for constructs and
ASCII/PDF reports. Native CTA/viral reports use retained provenance. Mutation
reports need `--ensembl-release N` when transcript annotation is missing; use
the producer's release. New sequence prediction requires host annotation.
These rules are the same for single files and manifests.

External inputs use the shared mRNA linker by default. Explicit junction
prediction has independent model and allele settings; see
[junction linker optimization](vaccine-designs.md#junction-linker-optimization).

LENS supplies reported context windows. pVACseq generally supplies the epitope
alone; Vaxrank does not fetch its separate mutant/wild-type FASTAs. LENS splice,
CTA/self and ERV construction require [source opt-ins](how-it-works.md#lens-antigen-source-opt-ins).

Native exports retain source rows, scope, content hashes, producer declarations,
manifest assertions, comparators and RNA units/derivation. Output directories
include `candidate_predictions.tsv` for replay and `ranked_candidates.csv` for
reviewing source windows. `--output-neoepitope-report predictions.csv` exports
prediction evidence with scope columns and `input_provenance_json`. See the
[evidence model](unified_evidence.md) for identities and source relationships.

## Reload or rank a normalized table

```sh
vaxrank --input-topiary candidates.tsv --output-epitopes saved.tsv --output-csv ranked.csv
vaxrank --input-epitopes saved.tsv --output-csv reloaded.csv
```

Topiary tables use `peptide`, `allele`, `kind`, `value` and optional prediction
and annotation columns. Its CSV/TSV writers preserve typed metadata and unknown
versus empty flanks. Enriched Vaxrank exports add native candidates, the evidence
frame, scoring policy and antigen references.

`load_predictions(path)` returns candidates; `EpitopeDataset.load(path)` exposes
the richer evidence. Older candidate-only files remain readable, but cannot
recover construction evidence they did not save.

## Replay a saved vaccine design

```sh
vaxrank --input-topiary candidates.tsv --output-dir design \
  --output-epitopes design/native.tsv --peptide-max-constructs 1 \
  --mrna-poly-a-length 73 --mrna-no-optimize-linkers
vaxrank --input-epitopes design/native.tsv --output-dir replay
```

A native file saved with construct outputs retains actual products, assembly
settings, source mappings and audits. Unchanged replay reuses those products
without prediction, codon optimization or annotation downloads, including
optimized linkers and their provenance.

A successful native run with zero candidates is also replayable. It retains its
input scope and settings, and writes empty product manifests when construct
outputs are requested. It does not acquire predictions or create placeholder
evidence. Empty external tables and legacy prediction-only files do not carry
this dataset contract and are rejected.

Current YAML and explicit CLI options override saved defaults. Changed sources,
policies or assembly settings invalidate retained products and audits. Combined
designs need compatible settings or explicit overrides, and are assembled anew.
Checksums detect altered product payloads. Older exports without saved products
can reconstruct from available evidence but cannot recover exact historical
products or omitted settings.

The original `--max-mutations-in-report` limit also applies to assembly; other
evidence stays in the native file. Move sibling reference bundles with the table
as described in [Moving native mutation datasets](native-references.md).

## Saved scoring policies

Native reload restores saved epitope policies unless explicitly overridden.
Combining different policies requires one explicit policy. Source-local model
preferences survive export, as do all observations after duplicate selection.

Table-only candidates need a `VaccineAntigen` for construction: sequence,
targetable mask and admission evidence. Supply `vaxrank_antigen_json` with a
matching `peptide_offset`, or use `EpitopeDataset.antigens` in Python. Missing or
held-out evidence appears in `Construction limitation`.

Add named prediction features without replacing originals:

```sh
vaxrank --input-topiary candidates.tsv \
  --external-predictions additive --external-prediction-prefix comparison \
  --mhc-predictor mhcflurry --mhc-alleles 'HLA-A*02:01,HLA-B*07:02' \
  --config-value 'epitopes.score_expr=comparison__mhcflurry__pMHC_affinity__score' \
  --output-epitopes comparison.tsv --output-csv ranked.csv
```

The prefix must be unique. Features affect ranking only when the DSL selects
them; omit the config override to retain the current formula. Model, version,
units and invocation metadata survive export. Additive inference scores existing
peptide/allele candidates; it creates neither new windows nor allele assignments
for allele-free Exacto observations.

Fresh and additive inference retain distinct occurrence flanks and comparator
contexts. `--external-peptide-only` explicitly discards flanks for inference.
When a native candidate has a source window, Vaxrank derives current inference
flanks from that exact sequence and local peptide offset. Its saved
`n_flank`/`c_flank` fields remain historical prediction metadata: replay never
rewrites them, and a sliced source window is interpreted at its rebased offset.
Without a source window, `None` still means unknown context and `""` a known
sequence terminus.
`--prediction-cache predictions.tsv` accepts contextual Topiary output with one
configured fallback model. The same cache policy applies to direct variants,
native antigen and CTA proteins, and fresh/additive external occurrences. Each
stage first validates its patient genotype and source scope, then queries exact
peptide, endpoint, genotype and flanking context. A reusable cache may contain
other genotypes; those rows never expand the current patient repertoire.
Distinct contexts stay distinct, and unknown flanks do not stand in for known
sequence termini. Missing coverage invokes the configured model, and Topiary
requires matching model/version/endpoint evidence before adding its results.
Cached measurements and input source identities remain unchanged. The explicit
`--external-peptide-only` mode instead queries without flanks; it does not reuse
context-dependent measurements as peptide-only results.

Original-only native replay needs neither a model nor a cache. Optional
processing annotation requires `--processing-aware-annotation`.

## Native Exacto inputs

Exacto peptide-variant, primary-structure and translation tables, including gzip,
use Topiary's native reader:

```yaml
schema: vaxrank.input_manifest.v1
patient_id: patient-001
reference_assembly: GRCh38
mhc_alleles: [HLA-A*02:01, HLA-B*07:02]
inputs:
  - format: exacto
    path: peptide-variants.tsv
    sample_id: tumor-rna
    library_id: tumor-long-reads
    exacto:
      primary_structures: primary-structures.tsv
      transcript_read_support: transcript-read-support.tsv
      read_set_id: tumor-bam-read-names
      tag: exacto-run-001
  - format: pvacseq
    path: pvacseq.tsv
```

Companion paths are relative to the manifest. Optional `schema` selects the
native schema; extra fields remain in metadata. Read-support companions require
library and read-set identities. Counts remain transcript-level evidence with
exact read membership, not variant support or TPM.

Without primary structures, context and geometry remain unknown. Full ORFs and
partial translations stay distinct. Exacto supplies no tumor-specificity
admission, reference comparator or pMHC prediction; construction needs those in
enriched evidence. ORF-only records can be exported without candidate peptides.

Inspect event, ORF, occurrence and RNA relationships through
`dataset.evidence_views()`. Shared sequence alone does not merge ORFs or transfer
another source's prediction, admission or support. Legacy report IDs remain in
`vaxrank_prediction_id` when additive inference normalizes evidence.

## Evidence model

See [Evidence and provenance](unified_evidence.md) for the data model.
Run `vaxrank --print-default-config` for current configuration defaults.

## Experimental window exposure

`vaxrank-window-exposure` explores delivery of admitted SLP windows under explicit
kinetic assumptions, without changing selection. See the
[window exposure guide](core-window-exposure.md).

## Coverage in design reports

HTML, PDF and text reports lead with input coverage, then emitted products,
ranking, source outcomes, optional candidate details, and methods/provenance.
The source ledger retains held-out antigens, records without retained predictions,
producer-unresolved ORF models, and sources that did not reach the ranked list.
LENS source-policy exclusions and unparseable external sources appear explicitly.
CTA expression decisions retain their measurement level, units, sample and assay.
Reasons can overlap; they must not be added as disjoint losses.

The family table counts recorded source identities, deduplicated by the shared
ranker's identity rules. These are not counts of independent events, reference
models, RNA reconstructions or unique vaccine products. The producer table and
actual product section report those different units separately. A family with
no supplied evidence recorded has not thereby been screened and excluded.
Retained epitope records do not establish complete prediction coverage.

New direct VCF/BAM runs retain parsed variant counts and per-variant outcomes
in the native source-evidence catalog, before shared ranking and product
assembly. Replay preserves those original observations. Different observations
of the same input remain separate histories; do not sum them. Historical files
without this metadata show unknown original denominators. The report does not
reconstruct upstream VCF counts from surviving antigens.

Original receipts and current replay-container declarations are displayed
separately. Neither later manifests nor unknown scope fields overwrite the
original sample, library, timepoint or reference record. All of this uses the
existing native dataset; no additional evidence sidecar is required.

## Inspect a saved design without copying its prediction table

Use `--no-output-epitopes` to keep the report/product bundle while disabling its
default native evidence export:

```sh
vaxrank --input-epitopes saved/candidate_predictions.tsv \
  --output-dir review --no-output-epitopes \
  --output-html-report review/report.html
```

This writes the usual product files, text/PDF reports and ranked CSV, plus the
requested HTML report. It does not rewrite `candidate_predictions.tsv`.
Original native inputs and their reference bundles remain the evidence source;
keep them together. Effective policy, saved product mappings and replay checks
remain active. Changing design controls can still request a new assembly;
this flag only controls the native output file. It cannot be combined with an
explicit `--output-epitopes PATH`. Without it, the default output bundle is
unchanged.

The Python report API can also render existing products directly: pass the
retained product objects to `TemplateDataCreator(products=...)` and then call
`make_html_report`, `make_pdf_report` or `make_ascii_report`. Report rendering
itself does not export native evidence.

## Reading event kind and coding context

An additional variant ORF can describe a mutation in an annotated CDS or in a
noncanonical ORF. Reports keep these dimensions separate: for example,
`mutation / noncanonical dORF`. The Isovar producer's reference kind, ORF type,
ORF ID, catalog, reference identity and source transcript travel with the
candidate and its product source mappings. Catalog context does not establish
translation or tumor specificity in this sample.

Candidate CSV/XLSX exports retain `source_kind` and add `coding_context` plus
the structured `coding_references` JSON column. A shared product can have
multiple coding contexts; neither the first display name nor a shared peptide
sequence erases the other origins. Historical exports recover context only
from an unambiguous linked reference model; otherwise the report says
`coding context not recorded`. Rendering historical products does not rewrite
their saved names, sequences or manifests.
