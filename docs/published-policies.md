# Published comparison policies

Inspect the bundled definitions without providing biological inputs:

```sh
vaxrank --list-policies
vaxrank --show-policy lens-v1
vaxrank --show-policy pvacseq-aggregate-v1
vaxrank --show-policy tesla-presentation-v1
vaxrank --show-policy tesla-recognition-v1
```

Select a bundle with `--config builtin:NAME`. `--show-policy` prints its resolved
configuration, expanded Topiary expressions, required evidence, source version
and citation. The frozen `openvax-v1` bundle remains available. Overrides create
a derived configuration; the saved effective definition and digest identify it.

| Bundle | Evidence and behavior |
| --- | --- |
| `lens-v1` | LENS 1.9 documented score: `max(0, 1 − IC50/1000) × log2(coding-support reads)/max(log2(coding-support reads)) × CCF`. Missing optional CCF contributes 1. |
| `pvacseq-aggregate-v1` | Retains upstream best-peptide/Tier decisions and reproduces pVACtools 7.1.4 default aggregate ordering: Tier, summed dense ranks of expression/IC50/percentile, IC50 rank, Gene, AA Change. |
| `tesla-presentation-v1` | Experimental Wells et al. presentation gates: affinity `<34 nM`, tumor epitope abundance `>33 TPM`, and pMHC binding stability `>1.4 hours`. Passing epitopes order by affinity. |
| `tesla-recognition-v1` | Adds recognition gates: mutant/WT affinity `<0.1` with positive WT affinity, or published foreignness `>10⁻¹⁶`. A known passing OR branch is sufficient; otherwise missing evidence excludes the occurrence. |

These bundles compare epitope selection. Vaccine construction, multiple-epitope
aggregation and manufacturability remain Vaxrank procedures. Their combined
score uses `target_epitope_score`, preventing an additional RNA factor. Their
scores are not calibrated immunogenicity probabilities.

## LENS evidence and normalization

Use a LENS report or normalized Topiary evidence with the reader's actual
`lens_rna_reads_covering_genomic_origin_with_peptide_cds` column. This is the
reported coding-sequence support count, not Isovar partial support, gene TPM or
arbitrary coverage. Affinity is in nM and optional `ccf` is a fraction.

Normalization covers all occurrences in each input report/source partition.
Evaluate patients/reports separately; report-wide operations reject mixed named
samples. Repeated predictor rows do not add support. No pseudocount is added:
zero support or a zero log-support maximum stays unscorable. Upstream discovery
and filtering are retained in the input, not reconstructed by this score.

[LENS 1.9 scoring specification](https://uselens.io/en/lens-v1.9.0/faq.html#how-is-the-prioritization-score-calculated)

## pVACseq aggregate scope

Use the aggregate report columns `Tier`, `Allele Expr`, `IC50 MT`, `%ile MT`,
`Gene`, and `AA Change`. The Topiary reader exposes these as `pvacseq_tier`,
`rna_alt_expression`, `affinity.value`, `affinity.rank`, `gene`, and `aa_change`.
The preset requires the default median IC50/combined percentile metric choices.
Nondefault source metrics need a derived policy.

The upstream Tier carries presentation, binding, expression, clonality,
transcript, reference-match and mutation-anchor decisions. This bundle consumes
those decisions; it does not assign tiers to raw all-epitopes rows. Missing
numeric sort measurements rank at the bottom. Unknown tiers are excluded.
Expression ranks descending, matching the pinned implementation despite the
documentation's ascending wording. The score `1/ordinal_rank(...)` represents
the ordering and exact sort-key ties remain tied. Conflicting scores for multiple
occurrences of one candidate require an explicit duplicate-resolution choice.

Tests compare a synthetic aggregate report covering every supported Tier,
missing measurements and ties against the actual v7.1.4 sorting function. The
fixture records the source commit and hash of that function's file.

[pVACseq aggregate criteria](https://pvactools.readthedocs.io/en/stable/pvacseq/output_files.html#the-pvacseq-aggregate-report-tiers),
[pVACtools 7.1.4 sorting implementation](https://github.com/griffithlab/pVACtools/blob/v7.1.4/pvactools/lib/sort.py)

## TESLA evidence

Provide normalized Topiary rows with affinity in nM, `pMHC_stability` values in
hours, and `tumor_abundance_tpm` containing tumor epitope abundance. Gene TPM and
read counts are not automatically substituted. Recognition also requires matched
WT affinity and/or the published IEDB similarity-derived `foreignness` feature.
The CLI does not calculate foreignness or obtain missing experimental measures.

All thresholds are strict and missing required evidence remains unknown. These
are explicit experimental reproductions of cohort-derived published gates,
not newly fitted or externally validated Vaxrank models. Ordering survivors by
affinity is a disclosed comparison convention.

[Wells et al., Cell 2020](https://doi.org/10.1016/j.cell.2020.09.015)

## Save and replay

For a report, select its corresponding bundle and save native evidence:

```sh
vaxrank --input-lens report.tsv --config builtin:lens-v1 \
  --output-epitopes evidence.tsv --output-csv lens-comparison.csv
vaxrank --input-pvacseq aggregate.tsv --config builtin:pvacseq-aggregate-v1 \
  --output-epitopes evidence.tsv --output-csv pvacseq-comparison.csv
vaxrank --input-topiary measured.tsv --config builtin:tesla-recognition-v1 \
  --output-epitopes evidence.tsv --output-csv tesla-comparison.csv
vaxrank --input-epitopes evidence.tsv --output-csv replay.csv
```

Native evidence retains the effective policy, all source measurements and
per-occurrence decisions. The adjacent `.policy` directory records configuration
derivation, source citations, model choices and replayable Topiary audit tables.
Evaluation and replay do not run predictors. Required missing evidence is not
filled with plausible measurements; percentile-only source rows remain available
even when they cannot create a quantitative native prediction object.

These table-based comparisons do not require lifting the Isovar dependency cap
(#562), relocating mutation annotation references (#563), or completing the
remaining exact-window prediction consumer integration (#497).
