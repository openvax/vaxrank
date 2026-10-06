# Expression-first CTA admission (#584)

This is the first implementation slice of [#584](https://github.com/openvax/vaxrank/issues/584).
The library adapter admits reference CTA proteins from patient expression
tables without VCF/BAM. CLI, report and Hitlist integration remain to be wired.

```python
from pyensembl import EnsemblRelease
from vaxrank.cta_admission import CTAAdmissionPolicy
from vaxrank.cta_expression import (
    CTAExpressionInput, CTATargetExclusionPolicy, admit_cta_expression,
)

result = admit_cta_expression(
    "patient-expression.tsv",
    input_contract=CTAExpressionInput(
        measurement_level="gene",  # or "transcript"
        id_column="gene_id", value_column="TPM",
        expression_unit="TPM", sample_id="patient-001",
        evidence_source="patient RNA quantification workflow",
        evidence_version="workflow-2026.10", assay="bulk RNA-seq",
    ),
    admission_policy=CTAAdmissionPolicy(min_tumor_expression=2, expression_unit="TPM"),
    genome=EnsemblRelease(93),  # caller's installed, declared annotation
    exclusion_policy=CTATargetExclusionPolicy(
        gene_patterns=("MAGE*",), gene_exceptions=("MAGEA4",)),
)
result.save("cta-admission.json")
antigens = result.admitted_antigens
```

The example threshold is a caller choice, not a universal biological cutoff.
Units, sample, measurement level, column names, assay and producer version must
be declared. The adapter uses Topiary's expression loader with the declared
columns and rejects malformed values or identifiers before coercion can lose
them. Blank/NA values stay unknown; zero is a measured value. Negative,
infinite, nonnumeric, missing-ID and duplicate normalized-ID rows fail.
Ensembl gene and transcript IDs are supported, including version suffixes;
the original feature identifier and source SHA-256 stay in the result.

Gene inputs select OncoRef's canonical transcript as a reference sequence.
That selection does not assign the gene's abundance to an individual isoform
or to a protein. Transcript inputs select only their measured transcript and
hold out a version mismatch. No values are summed across isoforms, genes or
identical proteins. Identical proteins retain independent gene/transcript
occurrences and expression evidence. Neither path fabricates alternate-read
counts or variant evidence.

OncoRef owns CTA membership and canonical admission. Each CTA decision carries
its versioned reference resolution, including noncanonical, expression-unknown
and excluded targets. Non-CTA and unresolved features stay in the decision
list. Canonical CTAs above the threshold produce admitted `VaccineAntigen`
objects through the existing admission model; noncanonical candidates remain
held out. This first adapter does not introduce an override route.

The default gene-symbol policy excludes `MAGE*` except exactly `MAGEA4`.
Empty patterns and exceptions allow all otherwise admissible targets. Matching
uses case-sensitive canonical OncoRef symbols, and exceptions are exact
symbols. This changes target eligibility only: the full OncoRef CTA candidate
universe remains the negative-self-reference exclusion set.

`CTAExpressionResult.load("cta-admission.json")` uses allowlisted native
serialization and reproduces policies, decisions, expression evidence, genome
identity and admitted source sequences without the expression file, new
OncoRef queries or prediction. It replays admission; prediction and construct
selection have their own native-dataset persistence paths.

Scientific basis: [Salmon](https://doi.org/10.1038/nmeth.4197) estimates
transcript abundance. [Ensembl canonical selection](https://www.ensembl.org/info/genome/genebuild/canonical.html)
chooses a representative transcript, rather than measuring patient isoform or
protein abundance.

## Draft readiness and subsequent work

- [Topiary #486](https://github.com/openvax/topiary/issues/486) blocks CSV tables
  with fewer than two data rows. A strict expected-failure regression records
  that upstream bug; adopt its released fix before marking this PR ready.
- Wire expression input into the shared CLI prediction, filtering/ranking,
  native dataset, report, peptide and mRNA construct paths. Use Topiary's DSL
  for configurable expression-aware selection and scores.
- Integrate versioned Hitlist observations with positive MS modality and full
  occurrence, assay, sample and HLA assignment provenance. Binding or
  fluorescence evidence must not count as MS presentation evidence
  ([Hitlist #644](https://github.com/pirl-unc/hitlist/issues/644)).
- Run end-to-end gene/transcript CLI fixtures and a real-data smoke, and retain
  the final construct and safety audit findings.

The OncoRef dependency pin remains unchanged; coordinated data adoption is
tracked separately in [#579](https://github.com/openvax/vaxrank/issues/579).
