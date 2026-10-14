<a id="experimental-expression-first-cta-admission"></a>

# Experimental CTA expression inputs

Vaxrank accepts patient expression tables without VCF/BAM. OncoRef supplies
human CTA membership; Topiary/MHCtools predict admitted proteins for the shared
ranking and peptide/mRNA design workflows.

To combine expression with viral sequences or saved mutations, use a
[`cta_expression` manifest entry](../antigen-inputs.md#add-cta-expression).
For one table:

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

Install the annotation and predictor assets first. Use `transcript` with
Ensembl transcript IDs for transcript measurements. The threshold is required
in the declared units; `2 TPM` here is an example choice.

Default exclusions are `MAGE*` with an exact `MAGEA4` exception. Repeated
`--cta-exclude-gene-pattern` and `--cta-allow-gene` options replace those lists;
`--cta-no-gene-exclusions` removes them. Matching uses canonical, case-sensitive
OncoRef symbols.

Windows use the shared selector and vaccine length settings. Ranking defaults
to binding. Expression rules are explicit, for example
`--config-text 'epitopes.filter_expr=cta_gene_tpm > 3'` or
`--config-text 'epitopes.score_expr=cta_expression_value'`. Gene and transcript
TPM columns are distinct; other units do not populate TPM columns.

Alongside standard reports and constructs, outputs include:

| File | Contents |
|---|---|
| `cta_admission.json` | Every admission decision, including held-out targets |
| `cta_target_decisions.csv` | Selection and rejection reasons |
| `cta_design.json` | Selected coordinates and final product inventories |
| `cta_product_predictions.tsv` | Final product measurements, without changing native ranking |

Reload with `--input-epitopes cta-design/native.tsv --output-dir replay`.
Unchanged products, settings and audits replay without original inputs or new
queries. Changed products invalidate prior audits. See
[junction optimization](../vaccine-designs.md#junction-linker-optimization) for
separate junction-model settings.

## Optional Hitlist evidence

Install `vaxrank[hitlist]` and supply `--hitlist-evidence-bundle DIRECTORY`,
created by Hitlist's `export cta-evidence` command (>=1.65.1). Vaxrank verifies
artifact hashes, relationships, expression input, columns, units, measurement
level, Ensembl release and OncoRef version. It retains the captured evidence
without rebuilding indexes.

Scoring columns include `hitlist_ms_observations`, `hitlist_query_status`,
`hitlist_cta_specific` and `hitlist_avoid_sequence`. Per-allele facets distinguish
monoallelic MS, experimental restriction, predicted restriction and donor-HLA
candidates. These overlapping facets are not summed into a confidence score.
Donor typing gives candidate alleles, not measured restriction.

Only exact observed peptides acquire MS support; longer observations do not
support unobserved nested epitopes. Unqueried sequences have null support. Zero
means no positive observations in the captured source, not an experimental
negative. Only verified MS measurements count; fluorescence, binding and unknown
modalities do not. Public observations remain separate from patient measurements.
MS/tissue filtering is an explicit DSL choice. See the primary
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
    genome=EnsemblRelease(93),
    exclusion_policy=CTATargetExclusionPolicy(
        gene_patterns=("MAGE*",), gene_exceptions=("MAGEA4",)),
)
result.save("cta-admission.json")
antigens = result.admitted_antigens
```

Declare units, sample, measurement level, columns, assay and producer version.
Blank/NA expression stays unknown; zero is measured. Negative, nonfinite,
nonnumeric, missing-ID and duplicate normalized-ID rows fail. Original IDs,
including version suffixes, and the input SHA-256 are retained.

Canonical gene inputs use OncoRef's canonical transcript. Verified alternate
loci use their own longest annotated protein, with transcript ID breaking ties.
Transcript inputs use only the measured transcript; version mismatches are held
out. Gene expression is not assigned to an isoform or protein. Counts are not
summed across genes, transcripts or identical proteins; source occurrences stay
separate. See [Salmon](https://doi.org/10.1038/nmeth.4197) and
[Ensembl canonical selection](https://www.ensembl.org/info/genome/genebuild/canonical.html).

Only canonical CTAs above the threshold are admitted. Noncanonical, unresolved,
unmeasured and excluded targets remain in the decision list. Target exclusions
are separate from the self-reference background: the full verified CTA candidate
set defines source exclusions. Ambiguous or unverified loci and non-CTA sources
with identical sequences remain in that background.

Audits retain annotation identities, mapping decisions and reference hashes.
Earlier native files keep their original identity contract and exclusions on
replay. `CTAExpressionResult.load("cta-admission.json")` restores admission
without the expression file or new OncoRef queries; save the native candidate
dataset to replay predictions and constructs too.

## Released prerequisites

Core version requirements are installed with Vaxrank. Hitlist is optional.
CTA inputs currently reject `--prediction-cache`; contextual cache support is
tracked in [#624](https://github.com/openvax/vaxrank/issues/624).
