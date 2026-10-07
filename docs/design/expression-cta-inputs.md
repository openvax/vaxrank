# Expression-first CTA admission (#584)

Vaxrank designs CTA vaccines from patient expression tables and HLA without
VCF/BAM. OncoRef supplies the CTA reference; Topiary/mhctools predict the
admitted proteins, and the shared DSL selects/ranks peptide or mRNA designs.

```bash
vaxrank --input-cta-expression patient-expression.tsv \
  --cta-expression-level gene --cta-id-column gene_id --cta-value-column TPM \
  --cta-expression-unit TPM --cta-min-expression 2 \
  --cta-sample-id patient-001 --cta-expression-source 'patient RNA workflow' \
  --cta-expression-version workflow-2026.10 --cta-expression-assay 'bulk RNA-seq' \
  --ensembl-release 93 --mhc-predictor mhcflurry \
  --mhc-alleles 'HLA-A*02:01,HLA-B*07:02' \
  --vaccine-type peptide mrna --output-dir cta-design \
  --output-epitopes cta-design/native.tsv
```

The annotation and predictor assets must be available. Change the level to
`transcript` and supply Ensembl transcript IDs for transcript measurements.
The minimum is explicitly required in the declared units. Target exclusions
default to `MAGE*` with an exact `MAGEA4` exception; repeat
`--cta-exclude-gene-pattern` and `--cta-allow-gene` to replace these lists.
`--cta-no-gene-exclusions` explicitly allows all otherwise admitted CTA genes.

CTA designs select bounded windows using the shared window selector and vaccine
length settings. They retain all native source occurrences, including repeated
peptides and identical proteins from different genes. The initial score is the
ordinary binding policy. Expression-aware rules are explicit DSL choices, e.g.
`--config-text 'epitopes.filter_expr=cta_gene_tpm > 3'` or
`--config-text 'epitopes.score_expr=cta_expression_value'`. Gene and transcript
TPM columns are distinct; arbitrary declared units never populate a TPM column.

The output includes all admission decisions (`cta_admission.json`), target
selection/rejection reasons (`cta_target_decisions.csv`), exact selected native
coordinates and final product inventories (`cta_design.json`), unfiltered final
product measurements (`cta_product_predictions.tsv`), the standard reports and
peptide/mRNA files. Final product measurements do not change native ranking.
Native-region and actual-product audits reuse the existing construct audit API.
Composite or changed native mappings and chemically modified peptides retain
explicit unassessed findings. Assembly coverage records contributing regions,
unassembled regions and unresolved product names. Processing
compartments remain unassessed; neither RNA nor predicted binding demonstrates
presentation, immunogenicity or clinical safety.

Reload with `--input-epitopes cta-design/native.tsv --output-dir replay`.
Admission, predictions, effective policies, public evidence and unchanged
product audits and actual assembly (including selected optimized linkers) replay
without the expression/evidence inputs or upstream
queries. Changed products cannot inherit an old final-context audit. Explicit
junction optimization uses the ordinary separate junction-model options; CTA
input does not implicitly optimize linkers with the candidate model.

## Optional Hitlist evidence

Install `vaxrank[hitlist]` and pass `--hitlist-evidence-bundle DIRECTORY` from
Hitlist's `export cta-evidence` command. Vaxrank verifies the public portable
contract (Hitlist >=1.65.1), every artifact hash and captured relationships. It
requires the same expression input SHA-256, columns, units, measurement level,
Ensembl release and OncoRef package version. It does not rebuild evidence indexes
or silently change a reference. The bundle's CTA definition is recorded separately
from Vaxrank's full candidate negative-self-reference universe.

All mappings, contributors, raw presentation and excluded observations, identities,
lineage and tissue evidence are captured in the native file. Positive MS must be
established by the producer's verified modality/measurement contract; fluorescence,
binding, negative and unknown-modality observations do not count. This consumes the
corrected evidence export contract while raw scanner correction remains tracked in
[Hitlist #644](https://github.com/pirl-unc/hitlist/issues/644).

Scoring columns include `hitlist_ms_observations`, `hitlist_query_status`,
`hitlist_cta_specific`, `hitlist_avoid_sequence`, and tissue-donor/review fields.
Per-allele facets distinguish `hitlist_ms_monoallelic`,
`hitlist_ms_experimental_restriction`, `hitlist_ms_predicted_restriction`,
`hitlist_ms_donor_hla_candidate`, and coarse/untyped or unknown restriction
evidence. Donor typing is a candidate set, not a measured restriction. These
facets can overlap and are not a scalar confidence score. Public observations
are not measurements from this patient.

Only exact observed peptide strings acquire MS support. A longer observation
does not support unobserved nested epitopes. Unqueried sequences have null support,
including targets outside a bundle's selection. Zero in a captured indexed
sequence means no positive observations in that captured source, not an
experimentally measured negative. MS/tissue-based selection is an explicit DSL
choice. The benign-tissue interpretation follows the primary
[HLA Ligand Atlas study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8054196/).

## Library admission API

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

## Released prerequisites

- Requires Topiary >=5.94.3, which fixes delimiter inference for empty and
  single-feature expression tables ([Topiary #486](https://github.com/openvax/topiary/issues/486)).
  Consumer regressions cover CSV/TSV gene and transcript inputs with zero,
  one and multiple rows, with no expected failures.
- Requires OncoRef 1.8.207, validated in
  [#579](https://github.com/openvax/vaxrank/issues/579). The extended CTA panel
  remains a separate upstream alternative; default admission is unchanged.
- Contextual prediction caches remain tracked in
  [Topiary #468](https://github.com/openvax/topiary/issues/468); CTA input currently
  rejects `--prediction-cache` rather than discarding its native flank context.
