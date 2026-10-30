# Patient-background comparisons

`annotation_reference` is the default comparator policy. The opt-in
`patient_background` policy uses Isovar's exact RNA-derived coding hypothesis,
reverting specified somatic edits while retaining nearby germline and
unexplained sequence changes. It changes the comparator sequence and its
predictions; the existing mutant-evidence scoring formula is unchanged.

The baseline is explicit:

| Value | Reverted edits |
| --- | --- |
| `focal` | Only the nominated somatic variant |
| `supplied_somatic_haplotype` | The nominated variant and other supplied somatic edits attributed to that RNA hypothesis |

Neither choice implies that every somatic variant in the patient was supplied.
The source evidence records exactly which cDNA edits were removed and retained.
A computed baseline, or separate spanning REF-allele RNA, does not establish
observed patient-normal tissue or translated protein expression.

## Live VCF/RNA input

Add these controls to a normal Vaxrank VCF/BAM command:

```bash
--comparison-reference patient_background \
--patient-background-baseline focal \
--rna-sample-id sample-1 \
--rna-source tumor-rna-library-1
```

Pairs are exported from the Isovar results already reconstructed by that run;
there is no second BAM pass. In a pooled input manifest, `direct.sample_id` and
`direct.library_id` supply the same identities. CLI identities must agree with
manifest declarations. Supply germline calls through the existing Isovar
`--germline-vcf` controls when available.

## Saved Isovar hypotheses

For an Isovar protein-hypothesis input, supply its independent matching pair
export:

```bash
vaxrank --input-isovar hypotheses-input.json --mhc-predictor mhcflurry \
  --comparison-reference patient_background \
  --patient-background-baseline focal \
  --patient-background-pairs patient-background-pairs.json.gz \
  --output-dir design --report-candidate-details
```

`hypotheses-input.json` follows the existing [Isovar input contract](isovar-inputs.md).
Repeat `--patient-background-pairs` for distinct sample/read sets. Each export
must contain the full retained set for that sample/read set. Sample, read-set,
event, hypothesis, nucleotide translation, transcript and RNA evidence identities
must agree. An amino-acid sequence match alone never selects a comparator.
Other antigen families do not acquire a patient comparator merely because they
are labelled mutation-derived.

The equivalent YAML controls are:

```yaml
epitopes:
  comparison_reference: patient_background
  patient_background_baseline: focal
```

## Python API

```python
from isovar import export_patient_background_pairs
from vaxrank import EpitopeConfig, PatientBackgroundPairs, run_vaxrank

pairs = PatientBackgroundPairs(export_patient_background_pairs(
    results, sample_id="sample-1", source="tumor-rna-library-1", baseline="focal"))
design = run_vaxrank(
    results, predictor, patient_background=pairs,
    epitope_config=EpitopeConfig(comparison_reference="patient_background"))
```

`AntigenDataset.predict(..., patient_backgrounds=[pairs], config=...)` uses the
same policy for generic Isovar protein hypotheses. For a single explicit
hypothesis, `pairs.match_hypothesis(hypothesis_export, event_id, hypothesis_id)`
returns an exact context whose `compare(start, end)` accepts a zero-based,
half-open amino-acid interval in the original selected hypothesis.

`CandidateEpitope.comparators['patient_background']` contains the available
peptide and recorded Topiary predictions. `comparison_evidence` preserves the
selected baseline, states, original coordinates and all alternative pair IDs.
`reference_comparator` resolves the declared comparison; it returns `None` for
an unavailable patient background. The older `.wt` accessor remains the
annotation/supplied WT role and is not relabelled.

## Uncertainty, reports and replay

Peptide boundaries are projected through the recorded cDNA edits. No protein
alignment or equal-offset guess bridges a reverted indel. Unequal peptide
lengths are retained; predictions outside the requested peptide lengths remain
unmeasured. Boundaries crossing a length-changing edit, incompatible codon
frames, premature stops and missing contexts have explicit reasons.
Comparator predictions use exact peptide-only requests through Topiary, without
invented upstream or downstream flanks. The recorded prediction context makes
that limitation explicit for processing-aware models.

All compatible nucleotide/transcript alternatives participate. Conflicting
sequences, or a mixture of available and unavailable alternatives, remain
ambiguous. Identical target/comparator sequences remain explicit evidence.
Missing sequence and missing prediction are separate states; the presence of
measurements does not assert complete predictor/allele coverage.
Comparator affinity is displayed only for the target row's exact allele,
predictor and version. Conflicting measurements stay in the evidence graph and
are labeled in reports; spreadsheet affinity cells remain blank with an explicit
`conflicting` state.

Full reports label the comparison policy, peptide comparison and bounded
source context, including removed/retained edits, phase and separate original
REF-allele RNA support. Native exports retain the complete source once in their
existing evidence store. Each comparison links back to its exact records.
Peptide coordinates retain their original hypothesis interval when vaccine
windows are sliced.

Replay needs neither the original pair file nor RNA reconstruction or new
predictions. Saved comparisons retain their baseline policy when rescoring;
changing that policy requires new comparisons from the original Isovar input.
When supplying an explicit scoring configuration for a saved patient-background
run, retain its `comparison_reference` and `patient_background_baseline` values.
The native reader validates source identities and saved peptide coordinates.

Contracts: [Isovar patient-background pairs](https://github.com/openvax/isovar/blob/v1.64.0/docs/patient-background-pairs.md)
and [Topiary's reader](https://github.com/openvax/topiary/blob/v5.108.4/topiary/io_isovar_hypotheses.py).
