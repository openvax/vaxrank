# How Vaxrank works

Choose `--vaccine-type peptide`, `mrna`, or both before setting antigen lengths,
linkers and construct limits. Vaxrank ranks antigen windows, assembles constructs
and writes reports. Mutation/CTA/viral integration and expanded vaccine designs
are experimental.

## Upstream inputs

The VCF/BAM workflow starts after three upstream steps:

1. Call somatic variants with a tool such as
   [Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/5257897695387-Mutect2)
   or [Strelka2](https://github.com/Illumina/strelka).
2. Align tumor RNA-seq reads into a BAM.
3. Obtain the patient's MHC alleles from typing or clinical records.

## Pipeline stages

| Step | What runs |
| --- | --- |
| Prepare sequences | Varcode annotates variants; Isovar reconstructs RNA-supported protein context. Antigen and CTA inputs use their own adapters. |
| Predict and score | MHCtools runs predictors; Topiary applies filters and scoring expressions. Saved predictions skip model calls. |
| Select windows | Vaxrank chooses antigen windows using epitope scores, reference matches and the configured ranking rules. |
| Assemble and export | Vaxrank applies construct limits and linkers, then writes sequences, synthesis files and reports. |

See [library responsibilities](library-responsibilities.md),
[input workflows](input-workflows.md) and [vaccine designs](vaccine-designs.md).

## Previously scored evidence

LENS, pVACseq, Topiary and native Vaxrank files reuse supplied predictions by
default. A manifest can combine them with new antigen sequences or CTA
expression. [Input workflows](input-workflows.md) covers prediction modes and
required metadata.

### LENS antigen source opt-ins

LENS `FUSION` rows can enter constructs using their supplied `pep_context`.
`SPLICE`, `CTA/SELF` and `ERV` rows remain in reports but need an opt-in for
construct selection:

```sh
vaxrank --input-lens patient.lens.tsv \
  --include-antigen-source SPLICE \
  --include-antigen-source CTA/SELF \
  --include-antigen-source ERV
```

`--exclude-antigen-source` removes a category. The YAML equivalent is
`vaccine_peptides.included_antigen_sources`. Source IDs, transcripts, RNA and
normal-tissue evidence remain attached to each record. Targetable intervals
come from the report; CTA/self and ERV opt-ins are recorded as overrides.

### Fusion antigens from supplied assemblies

`VaccineAntigen.from_fusion_sequence(...)` takes a translated sequence and exact
amino-acid junction. Targetable peptides must span both partners.

For results from `isovar.reconstruct_fusion(...)`, use:

```python
from vaxrank import fusion_antigens_from_isovar, predict_epitopes, VaccinePeptide

adapted = fusion_antigens_from_isovar(result, tumor_specificity=attestation)
for antigen in adapted.admitted_antigens:
    epitopes = predict_epitopes(mhc_predictor=predictor, antigen=antigen, genome=genome)
    candidate = VaccinePeptide(
        antigen=antigen, epitopes=epitopes,
        combined_score_expr="target_epitope_score",
        ranking_rules=("target_epitope_score", "manufacturability", "self_epitope_score"),
    )
```

The adapter retains all coding hypotheses in `.antigens` and the full Isovar
result in `.reconstruction`. A separate tumor-specificity attestation is
required for admission. Ambiguous, unsupported and partial-CDS hypotheses stay
held out. Mixed junction codons are targetable; peptides wholly inside an
inserted sequence are not junction-spanning. This API consumes Isovar results;
fusion discovery from BAM files remains upstream.

## Mutant transcript assembly (Isovar)

Isovar assembles reads around each variant, retaining nearby phased variants and
supported splice context. The ordinary Vaxrank path selects
`IsovarResult.top_protein_sequence`. The separate fusion adapter retains multiple
hypotheses. [RNA context settings](configuration.md#rna-reconstruction-context)
control requested length and support thresholds.

## Candidate epitope scoring

MHCtools predicts overlapping peptides at the requested lengths. The default
binding policy converts affinity into a logistic score; custom filters and
scores use the [Topiary DSL](configuration.md#custom-filtering-and-scoring-with-the-topiary-dsl).

## Vaccine peptide ranking

The default mutation workflow ranks windows using target epitope scores and RNA
support, with manufacturability criteria for window selection. Reference
matching includes every available annotated protein sequence, including NMD,
IG/TR and other translated biotypes. Exact self matches are excluded from
mutation target epitopes. Cancer-hotspot annotations flag recurrent mutations.

[Selection policies](selection-policies.md) describes alternative scores and
window rules. [Mixed-antigen selection](mixed-antigen-design.md) explains how to
rank sources that carry different kinds of evidence.

## Data model

| Object | Contents |
| --- | --- |
| `VaccineAntigen` | Sequence, targetable intervals, source identity and admission evidence |
| `MutantProteinFragment` | Mutation positions, variant and RNA support for mutation inputs |
| `CandidateEpitope` | Peptide occurrence, context, predictions, comparators and scores |
| `VaccinePeptide` | Selected antigen window and its epitopes |

Ranked outputs are `(source, [VaccinePeptide, ...])` pairs. Direct mutation runs
can retain multiple windows with `max_vaccine_peptides_per_variant`; external
adapters select one source-derived window. `--peptide-candidates-per-slot` and
`--mrna-candidates-per-slot` control how many alternatives are exported.

See [evidence and provenance](unified_evidence.md) for shared-sequence mappings.

## Key modules

| Area | Modules |
| --- | --- |
| Selection and scoring | `core_logic.py`, `epitope_logic.py`, `vaccine_peptide.py` |
| Inputs and persistence | `epitope_io.py`, `external_input.py`, `antigen_input.py` |
| Reference matching | `reference_proteome.py` |
| Assembly and linkers | `peptide.py`, `mrna.py`, `vaccine_library.py`, `junction_swap.py` |
| Reports | `report.py` |

## Libraries Vaxrank builds on

[PyEnsembl](https://github.com/openvax/pyensembl) supplies annotation,
[Varcode](https://github.com/openvax/varcode) coding consequences,
[Isovar](https://github.com/openvax/isovar) RNA reconstruction,
[MHCtools](https://github.com/openvax/mhctools) prediction, and
[Topiary](https://github.com/openvax/topiary) filtering and ranking.
[DnaChisel](https://github.com/Edinburgh-Genome-Foundry/DnaChisel) optimizes mRNA
coding sequences. Jinja templates produce reports through WeasyPrint or pdfkit.
