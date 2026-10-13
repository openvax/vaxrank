# How Vaxrank works

Vaxrank turns antigen evidence into ranked vaccine peptides and finished
constructs. Variant calling, RNA alignment and HLA typing happen upstream. Vaxrank
handles antigen ranking, selection, construct assembly and reports.

## Pipeline stages

Choose the vaccine type first (`--vaccine-type peptide`, `mrna`, or both). That
choice sets the defaults and constraints for antigen selection, construct
assembly and reporting, including peptide synthesis checks and mRNA sequence
optimization. Vaxrank then:

1. **Prepares candidates.** [Varcode](https://github.com/openvax/varcode)
   annotates coding consequences for mutation inputs;
   [Isovar](https://github.com/openvax/isovar) assembles RNA-supported
   sequences and predicts protein context where the frame is justified.
   Experimental [unified inputs](antigen-inputs.md) admit CTA expression or
   load supplied antigen sequences through source-specific evidence adapters.
2. **Scores epitopes.** Predictors run through
   [MHCtools](https://github.com/openvax/mhctools); their evidence is scored and
   filtered with [Topiary](https://github.com/openvax/topiary). Saved
   predictions can be reused instead.
3. **Selects and assembles.** Candidates are ranked, antigen windows selected,
   and the chosen vaccine type's limits, linkers and assembly settings applied.
4. **Exports designs.** Peptide or mRNA sequences, synthesis planning files and
   the requested analysis reports are written out.

See [Varcode, Isovar, and Vaxrank](library-responsibilities.md) for the
division of responsibility between the libraries, and
[vaccine designs](vaccine-designs.md) for the design space of step 3.

## Upstream inputs

For the VCF + BAM workflow, three things happen before Vaxrank runs, for
example in
[neoantigen-vaccine-pipeline](https://github.com/openvax/neoantigen-vaccine-pipeline):

1. Call somatic variants from tumor and matched-normal DNA with
   [Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/5257897695387-Mutect2)
   or [Strelka2](https://github.com/Illumina/strelka), producing a VCF.
2. Align tumor RNA-seq reads to produce a BAM.
3. Supply the patient's HLA alleles from typing or clinical records.

## Previously scored evidence

Candidates that already carry predictions can enter through saved native
datasets, [LENS](https://github.com/openvax/lens),
[pVACseq](https://github.com/griffithlab/pVACtools), or an input manifest,
skipping reconstruction and prediction. These paths reuse the supplied
sequences and predictions by default; new predictions can be requested
explicitly. They feed the same selection, assembly and report writers.
pVACseq parsing is delegated to Topiary and accepts both `all_epitopes.tsv`
and `all_epitopes.aggregated.tsv`. See
[input workflows](input-workflows.md) for the commands and their limits.

LENS `FUSION` rows are eligible for construct ranking. Vaxrank treats their
caller-supplied `pep_context` as an assembled fusion antigen rather than
fabricating a single-locus genomic variant. Fusion ID, both partner genes,
breakpoints, transcript IDs, fusion type, RNA evidence fields and sequence
source are preserved in the report, and ranked fusion antigens can be included
in both peptide and mRNA constructs. This path imports LENS's peptide
predictions; it does not infer a coding frame or translate raw breakpoints when
the upstream report supplies no established protein context.

### LENS antigen source opt-ins

LENS `SPLICE`, `CTA/SELF` and `ERV` rows are preserved in the neoepitope report
but require an explicit opt-in before they enter construct ranking:

```sh
vaxrank --input-lens patient.lens.tsv \
  --include-antigen-source SPLICE \
  --include-antigen-source CTA/SELF \
  --include-antigen-source ERV
```

`--exclude-antigen-source` removes a default or configured category (for
example `--exclude-antigen-source FUSION`). The same selection can be set with
`vaccine_peptides.included_antigen_sources` in YAML.

Evidence is retained per category: splice targets keep their genomic junction,
description, coding sequence and RNA evidence; CTA/self targets keep gene,
transcript, normal-tissue, mTEC and expression evidence; ERV targets keep the
ORF identifier and all `erv_*` evidence. Only caller-designated peptide
intervals are targetable, because these LENS rows do not state an exact
amino-acid junction offset. CTA/self and ERV opt-ins are recorded as
review-required tumor-specificity overrides, not as proof of sequence novelty or
absence from normal tissue.

### Fusion antigens from supplied assemblies

When an upstream assembler supplies the translated sequence and the exact
amino-acid junction, `VaccineAntigen.from_fusion_sequence(...)` records the
junction as a zero-width target mask. `predict_epitopes(..., antigen=antigen)`
then enumerates the predictor's peptide lengths and marks only peptides
containing residues from both fusion partners as targetable; no fictitious
wild-type comparator is generated.

For validated results from `isovar.reconstruct_fusion(...)` or `isovar fusion`,
use `fusion_antigens_from_isovar(result, tumor_specificity=attestation,
gene_name="DONOR::ACCEPTOR", species="Homo sapiens")`. It returns an
`IsovarFusionAntigens` containing every coding hypothesis in `.antigens`, the
complete JSON-normalized result in `.reconstruction`, and the eligible subset in
`.admitted_antigens`. No best transcript is selected. Nucleotide junctions are
mapped exactly: a mixed codon is targetable, while peptides lying entirely
inside an inserted sequence are not junction-spanning.

RNA validation is not tumor-specificity evidence. Absent a separate
`TumorSpecificityAttestation`, every hypothesis is held out. Ambiguous,
unsupported and partial-CDS reconstructions remain held out even with an
admission attestation; their reasons, alternatives, RNA support and provenance
are preserved. Exploratory predictions can use any hypothesis, but only
admitted antigens can enter a `VaccinePeptide`. This adapter expects trusted
Isovar output: it does not validate arbitrary caller declarations against
reference annotation, and it does not discover fusions from BAM files.

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

Score fusion antigens source-agnostically as above, rather than with the
small-variant default's `n_rna_alt` binding. Full reconstruction provenance
travels in each antigen's `source_metadata`; no wild-type protein,
self-proteome exclusion or clinical neoantigen status is invented. The existing
small-variant `run_vaxrank` API and the LENS CLI path are unchanged.

## Mutant transcript assembly (Isovar)

For each somatic variant, [Isovar](https://github.com/openvax/isovar) extracts
RNA-seq reads overlapping the mutant locus, assembles them into RNA sequence,
then predicts protein context where the frame is justified. Compared with an
isolated DNA edit, this can:

- Preserve nearby variants linked by compatible read/fragment evidence.
- Represent observed splice differences where reconstruction and frame matching
  support them; not every possible splice path is resolved automatically.
- Provide RNA support for the reconstructed sequence. Missing support fails the
  configured RNA evidence criteria; it does not prove biological absence.

The ordinary path selects Isovar's top protein sequence. The separate
supplied-fusion adapter retains coding hypotheses; neither implies that the
ordinary pipeline evaluates every RNA/SV alternative. See
[current limits and the DNA fallback](library-responsibilities.md#available-today-and-remaining-work).

## Candidate epitope scoring

Each mutant protein fragment is sliced into overlapping subsequences of epitope
length (typically 8–15 amino acids). These candidate epitopes are scored for
predicted MHC binding affinity by the selected predictor, and the predictions
are converted to a score between 0 and 1 through a logistic function
parameterized by the `EpitopeConfig` settings (see
[configuration](configuration.md)).

## Vaccine peptide ranking

Candidate antigen windows, typically 25 amino acids, are constructed
around each mutation and scored from the epitope predictions they contain. Candidates are then filtered and ranked by:

1. **Epitope content:** sum of configured epitope scores. Binding predictions
   alone do not establish immunogenicity.
2. **Reference proteome filtering:** exact matches in the selected genome's
   annotated protein sequences are removed from mutation-derived candidates,
   regardless of the source transcript biotype.
3. **Cancer hotspot annotation:** variants at known recurrently mutated
   positions are flagged, using bundled data from
   [cancerhotspots.org](https://www.cancerhotspots.org/) (~2,700 mutations
   across cancer types).
4. **Peptide manufacturability:** tie-breaking by hydropathy-based synthesis difficulty
   (C-terminal and 7-mer window GRAVY scores).

Reference matching includes available annotated NMD, non-stop-decay, IG/TR and
loss-of-function protein sequences. Source-attributed matches and construct
reports retain each transcript's biotype; membership does not establish
expression or presentation. Older saved sources without a biotype remain
explicitly unknown. The translation-selection policy is part of the cache key,
so older indexes are rebuilt automatically. See the
[GENCODE/Ensembl biotype definitions](https://www.gencodegenes.org/pages/biotypes.html).

## Data model

Vaxrank's central data unit is the VaccinePeptide (VP), a ranked antigen window
considered for vaccine construction. A VP bundles:

- a source-agnostic `VaccineAntigen`: amino-acid sequence, targetable
  intervals, source provenance and tumor-specificity attestation;
- for ordinary SNV/indel antigens, a `MutantProteinFragment` with mutation
  positions, gene name, source variant and RNA evidence metrics;
- `CandidateEpitope` objects: peptide occurrences carrying context, provenance,
  per-allele/model/version `mhctools.Prediction` values, known comparators and
  DSL scores. Target epitopes drive vaccine ranking.

The shared output is a list of `(source, [VaccinePeptide, ...])` tuples.
Ordinary variant sources are `varcode.Variant` objects; non-variant antigen
sources retain their own identity. The direct pipeline can retain multiple
windows per variant via `max_vaccine_peptides_per_variant`; external adapters
currently select one source-derived window:

```
ranked_variants_with_vaccine_peptides = [
    (Variant_A, [VP_A1, VP_A2, ...]),    # multiple windows around variant A's mutation
    (Variant_B, [VP_B1]),                 # single SLP for variant B
    ...
]
```

Three settings control how many alternates a variant contributes:

- `--vaccine-peptide-length`: SLP size and the default RNA context request;
  `--protein-sequence-length` or legacy `--padding-around-mutation` can
  override the reconstruction context explicitly.
- `max_vaccine_peptides_per_variant` (config): how many alternate windows per
  variant reach the ranked output.
- `--peptide-candidates-per-slot` / `--mrna-candidates-per-slot` (CLI): how
  many VP alternates per variant slot the construct assembler renders into
  FASTAs.

Reports render one section per variant; within a section, each VP gets its own
per-epitope sub-table. Column counts can differ per VP. For example, when
pepsickle credibility tagging succeeded for one VP and failed for another, only
the successful VP's table shows the processing columns.

## Key modules

Shared evidence and selection:

- `core_logic.py`: main vaccine peptide selection algorithm
- `epitope_logic.py`: candidate epitope scoring and filtering
- `epitope_io.py`: LENS / pVACseq / vaxrank-native I/O for epitope predictions
- `external_input.py`: synthesizes the canonical ranked-vaccine-peptides shape
  from a LENS / pVACseq report, so external-input runs reach the same dispatch
  as VCF + BAM
- `reference_proteome.py`: set-based kmer index for reference proteome
  filtering (O(1) lookup, built once and cached)
- `cancer_hotspots.py`: cancer mutation hotspot annotation
- `vaccine_peptide.py`: vaccine peptide scoring and manufacturability
- `vaccine_library.py`: shared linker vocabulary and compositional grammar with
  primary-source citations

Assembly for each vaccine type:

- `peptide.py`: peptide construct assembly plus FASTA / JSON manifest / vendor
  order-form CSV writers; sub-modes `slp` / `minimal_epitope` / `multi_epitope`
- `mrna.py`: mRNA construct assembly plus three-FASTA / structured manifest /
  long-format CSV writers, DnaChisel codon optimization, 2A frozen-codon
  handling, configurable polyA tail (default A120, optional segmented
  BNT162b2 pattern), per-junction MHC-aware linker swap (issue #247)
- `mrna_library.py`: mRNA-specific elements (5'/3' UTRs including the tandem
  2× HBB FI; signal peptides HLA-A / HLA-B / tPA / IgK / CD8A / CD28; MITD
  HLA-A / HLA-B)
- `junction_swap.py`: per-junction linker optimizer that minimizes predicted
  MHC presentation of chimeric k-mers spanning antigen junctions

Reports:

- `report.py`: analysis-report generation (ASCII, HTML, PDF, XLSX, CSV, JSON)

## Libraries Vaxrank builds on

From the [OpenVax](https://github.com/openvax) ecosystem:

- [pyensembl](https://github.com/openvax/pyensembl): reference genome annotation
- [varcode](https://github.com/openvax/varcode): variant effect prediction from DNA
- [isovar](https://github.com/openvax/isovar): RNA-based mutant transcript
  assembly and variant phasing
- [MHCtools](https://github.com/openvax/mhctools): unified interface to MHC
  binding predictors
- [topiary](https://github.com/openvax/topiary): epitope filtering, scoring and
  ranking DSL

Other key dependencies: `msgspec` for configuration serialization (YAML/JSON),
`pandas` and `numpy` for data processing, `jinja2` with `pdfkit`/`weasyprint`
for report generation, and `dnachisel` for mRNA codon optimization.
