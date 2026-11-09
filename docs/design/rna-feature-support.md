# RNA support across antigen categories

Implementation of https://github.com/openvax/vaxrank/issues/757. Preserve the
frozen openvax-v1 configuration, including its RNA and pMHC settings. New read,
fragment, per-million and effective-length-normalized measurements are optional
inputs to the existing DSL and retain the selected feature, sample and library.

Implementation sequence:

1. Consume Isovar's common feature-support record; preserve it on antigen and
   mutation objects, native datasets and per-window scoring views.
2. Expose identical DSL bindings for every source kind. Keep old expressions
   unchanged; missing measurements fail explicitly when requested. Record
   normalization and score policies with reports and replay.
3. Exercise API/CLI input, mixed-source selection, native replay and reports;
   run lint, full tests and an end-to-end CLI/report smoke before release.
4. Rerun the six-library longitudinal Sid analysis and separate T2 comparison
   using acquired local RNA evidence for CTAs as well as mutations. Preserve
   baseline artifacts and use explicitly named RNA expressions; compare raw
   and normalized measurements, source/timepoint support, final selections and
   construct changes. Retain source-specific limitations and known holdouts.

No fixed scoring method, category multiplier or temporal aggregation is implied
by the measurement contract. Read and fragment effective lengths are separate;
SE counterparts must agree. Counts, unit identities, denominators, geometry and
feature scope must refer to the same observation rule.

The producer feature declares `fragment_mode="union"` for joint sequenced-mate
support or `"any_read"` for templates with an individually supporting read.
Both counting and effective length honor this choice. Neither rule is selected
implicitly by antigen category.

## Interface

Isovar 1.85 provides `RNAFeature`, `RNASamplingProfile`,
`effective_rna_feature_lengths`, `measure_rna_feature` and
`measure_rna_path_features`, plus the `isovar rna-features` CLI. Its
[feature contract](https://github.com/openvax/isovar/blob/master/RNA_FEATURE_SUPPORT.md)
defines the counting predicate and `L_effective` for base substitutions,
indels, long insertions, rearrangements, complex SVs, splice changes, altered
ORFs, CTAs, non-mutant isoforms, ERVs, viral proteins and ncORFs.

Attach a validated producer record with the public `RNAFeatureEvidence` type:

```python
from dataclasses import replace
from vaxrank import RNAFeatureEvidence

# Exact 25-aa window within the parent protein. The producer must measure
# the corresponding nucleotide interval on its declared RNA path.
antigen = replace(antigen, rna_support=(
    RNAFeatureEvidence(measurement, scope="window", protein_interval=(40, 65)),
))
```

`scope="event"` with no protein interval means evidence for the declared causal
RNA feature, applicable to this source's targetable windows. It does not assert
that downstream coding bases were observed. `scope="window"` binds only at
that exact protein window. Slicing rebases contained measurements and removes
those that cannot apply; it never prorates counts. `VaccineAntigen` and legacy
`MutantProteinFragment` both carry the same `rna_support` tuple.

The existing `AntigenDataset.save/load` and `--input-antigens` accept these typed
attachments. Native prediction exports preserve them for offline replay.
Measurements with multiple features/libraries remain separate. A scalar score
requires one applicable measurement; callers choose the feature/library or
compare separately scored samples under an explicit temporal policy. No
implicit sum, average, maximum, independence assumption or missing-to-zero
conversion is introduced.

The same columns are available to Topiary's scoring view and the existing
combined-score DSL for every antigen category:

| Quantity | Read binding | Fragment binding |
| --- | --- | --- |
| Observed support | `rna_reads` | `rna_fragments` |
| Whole-library denominator | `rna_library_reads` | `rna_library_fragments` |
| Effective length in bp | `rna_effective_length_reads` | `rna_effective_length_fragments` |
| Per million | `rna_reads_per_million` | `rna_fragments_per_million` |
| Per effective kb | `rna_reads_per_effective_kb` | `rna_fragments_per_effective_kb` |
| Per effective kb per million | `rna_reads_per_effective_kb_per_million` | `rna_fragments_per_effective_kb_per_million` |

For example, an explicitly selected new policy can use:

```sh
vaxrank --input-antigens measured-antigens.json --config-text \
  'vaccine_peptides.combined_score_expr=sqrt(rna_fragments_per_effective_kb_per_million) * target_epitope_score' \
  --report-candidate-details --output-dir design
```

Predictions and MS evidence still determine `target_epitope_score` through the
selected pMHC configuration; RNA supplies an independently configurable factor.
In Topiary expressions, read these columns with `column(rna_fragments)`; for
example `epitopes.filter_expr=column(rna_fragments) >= 2`. Math transforms use
Topiary syntax such as `column(rna_fragments).sqrt()`. Combined-score expressions
use the flat binding names and Python-style functions shown above.

Raw support and FPM are alternative expressions using the other bindings, not
category-dependent defaults. The frozen `openvax-v1` profile retains its complete
original RNA and pMHC configuration. No existing score changes just because a
new measurement is attached.

Detailed reports show sample/library, feature intervals and rule, both counts,
both effective lengths, all available normalizations, missing reasons, geometry
assumptions and measurement provenance. These are observed local RNA counts and
an experimental geometric density, not gene TPM or evidence of translation.

## Acquisition requirements

The common interface is available to every category, but new counts require
original RNA evidence and an explicit feature/path. Existing TPM-only CTA data
or historical assembly counts do not become these measurements automatically.
Missing geometry leaves length-normalized values unavailable; missing original
whole-library counts leaves RPM/FPM unavailable. A cropped regional BAM never
supplies its own whole-library denominator. Sequence/locus/isoform ambiguity
stays in the producer provenance and does not become unique attribution.
