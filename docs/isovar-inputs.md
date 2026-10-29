# Isovar ORFs as Vaxrank inputs

Vaxrank 3.47 adds `--input-isovar` and `format: isovar` manifest entries. They
accept saved Isovar evidence or invoke a released Isovar reconstruction API,
then use the same Topiary prediction, selection, peptide/mRNA design, reports,
and native replay as [other antigen inputs](antigen-inputs.md).

Requires Isovar **1.61.0+**, ORFdata **0.39.1+**, and their compatible dependencies.
Results are interpreted by their exact schema, not by guessing from columns.
Unknown schemas and contradictory sequences, coordinates or declarations fail.
The existing `--vcf`/`--bam` workflow remains available; it does not automatically
search ncORF catalogues or reconstruct structural events.
For mRNA linker optimization, select `--mrna-junction-predictor` explicitly,
as with other native sequence inputs; the candidate predictor alone does not
request new junction predictions.

## Supported sources

| Family | Producer / accepted result | Candidate interpretation |
| --- | --- | --- |
| `variant_orfs` | `reconstruct_variant_orfs`; `isovar.variant_rna_orfs.v1` | Every supplied ordinary CDS/ncORF model and RNA alternative. Kind `mutation`; targetable amino-acid change comes from `mutation_interval`. |
| `sv` | `reconstruct_sv_rna`; `isovar.sv_rna_candidates.v5` | Both annotated-frame translations and exploratory ATG ORFs, kept separately. Kind `SV`; targetable masks cross observed altered junction boundaries. |
| `sv` | `export_sv_rna_orfs`; `isovar.sv_rna_orfs.v4` | Every exported nucleotide candidate and placement. A placement without independently resolvable support stays held out. |
| `fusion` | `reconstruct_fusion`; `isovar.fusion_rna.v3` | All fusion coding hypotheses, with split-codon/junction masks. Kind `fusion`. |
| `viral` | `reconstruct_viral_sequences`; `isovar.viral_sequences.v1` | Observed, reference-anchored coding predictions. Kind `viral`; source species and protein IDs remain separate from the host reference. |
| `isoforms` | `reconstruct_antigen_isoforms`; `isovar.antigen_isoforms.v1` | All annotated competitors and observed alternatives. Kind `isoform`; a reference model with unobserved coding bases cannot be admitted. |
| `protein_hypotheses` | `export_protein_hypotheses`; `isovar.protein_hypotheses.v2` | All ordinary variant protein hypotheses, translations, filter outcomes, phasing and support evidence. Kind `mutation`. Live export accepts typed `IsovarResult` objects in Python; the CLI replays an exported result. |

No adapter picks the longest ORF or highest-ranked Isovar protein in advance.
Synonymous ncORF results and stop-only changes without a targetable residue or
internal junction stay in the evidence with an explicit held-out reason. An
internal zero-width deletion boundary requires residues on both sides.

SV targeting is currently **junction-spanning**. It does not label an entire
host ORF, inserted segment or downstream tail as tumor-specific merely because
a path crosses an SV. Isovar's supplied-reference peptide comparisons remain
available in the source evidence, separately from Vaxrank's host-self assessment.

## Saved result input

Use a JSON envelope to supply patient scope and independent admission evidence:

```json
{
  "schema": "vaxrank.isovar_input.v1",
  "scope": {
    "patient_id": "patient-1",
    "sample_id": "tumor-1",
    "reference_assembly": "GRCh38",
    "annotation": "ensembl:110",
    "mhc_alleles": ["HLA-A*02:01"]
  },
  "result": "isovar-result.json.gz",
  "options": {
    "admission": {
      "allow_partial": false,
      "allow_ambiguous": false,
      "min_fragments": 2
    }
  }
}
```

```bash
vaxrank --input-isovar source.json --ensembl-release 110 \
  --mhc-predictor random --no-processing-aware-annotation \
  --output-dir review --output-antigens review/antigens.json
```

`random` is a software smoke-test predictor. Choose a configured scientific
predictor for an analysis. This example supplies no tumor-specificity evidence,
so its candidates remain held out. RNA assembly alone does not admit a target.

`result` may be an inline object or a relative JSON/JSON.gz path. A raw producer
file is also accepted directly; declare the patient/reference/genotype in an
input manifest. All relative paths resolve beside the file declaring them.
A manifest can combine Isovar, saved predictions and other native antigens:

```yaml
schema: vaxrank.input_manifest.v1
patient_id: patient-1
reference_assembly: GRCh38
annotation: ensembl:110
mhc_alleles: [HLA-A*02:01]
inputs:
  - format: isovar
    path: source.json
    isovar:
      admission:
        allow_partial: true
  - format: epitopes
    path: earlier-predictions.tsv
```

Declare each source's sample/library correctly. Do not load a result alongside
its own native prediction export: that duplicates candidate identities.

## Admission and uncertainty

`IsovarAdmission` is shared by Python, source files, YAML and CLI overrides:

| Control | Default | Meaning |
| --- | --- | --- |
| `allow_partial` | `false` | Permit an observed local coding window to pass the completeness gate. Missing termini remain missing. |
| `allow_ambiguous` | `false` | Permit competing reference frames or unresolved/exploratory SV frames to pass the ambiguity gate. All alternatives remain separate. |
| `min_fragments` | `2` | Minimum producer-scoped RNA fragments; positive integer. Does not mean molecules, TPM, independent libraries, or protein expression. |

These controls **also require an independent tumor-specificity attestation**.
They do not bypass absent targetable content, unobserved isoform coding bases,
unavailable support, failed producer filters, or the minimum support gate.
Relaxing partialness does not relax frame ambiguity.
For variant ORFs, frames are compared within each RNA candidate using the
resolved carrier transcript's strand, the focal allele's codon phase (0, 1 or 2 in
transcript orientation), and the translation table. Trimming a noncoding
prefix cannot change this comparison. Candidate evidence records
`reference_strand` and `variant_codon_phase` so the decision can be inspected.
ORFdata 0.42.1 exports `resolved_transcript.strand` for every ORF kind, including
annotated CDSs and `TranscriptORF.at` placements. This carrier field is authoritative,
including an explicit `null` for unknown orientation. Historical snapshots without
the field retain their original source-model strand fallback.
For unknown orientation (`reference_strand: null`), phase comparisons can establish a shared frame only within the same
identified transcript and reference; alternatives across different transcripts
remain ambiguous because their relative orientation is unavailable.

The envelope accepts the existing `TumorSpecificityAttestation` record as
`tumor_specificity`, applying to every candidate in that source, or a
`candidate_attestations` object keyed by exact candidate ID. Prefer the latter
when a review supports only particular products. An attestation carries its
status, evidence kind/source/version, rationale, patient scope and evidence
records. A held-out or missing attestation never becomes admitted automatically.
Use a source-wide assertion only when its evidence covers **every** candidate.
The original requested assertion survives even when a reconstruction gate holds
that candidate out.

For viral/isoform inputs, an admitted source assertion makes the observed coding
region targetable. It does not add host-gene exclusions to the self reference.
The same exact-self checks and downstream Topiary DSL policies still apply.
An overlapping viral CDS is not automatically proof of the sample's translated
product. Reference annotation, RNA coverage, initiation and protein translation
are distinct claims.

Alternative Isovar candidate IDs do not imply distinct vaccine products.
Identical emitted sequences share a design slot, with all contributing source
mappings retained; related but different windows remain alternatives from the
same recorded event. See [distinct sequences and source mappings](vaccine-designs.md#distinct-sequences-and-source-mappings).

Support retains its original units and scope:

- Variant ORFs: assembled RNA fragments; coding tiling and single-read spans are
  retained separately. Tiling does not establish end-to-end phase.
- Annotated SV products: the minimum direct fragment count across their altered
  junctions; this is not a whole-protein spanning-read count.
- Exploratory SV products: full-interval, junction-linked witnesses for the
  particular occurrence, using the measured occurrence field exported by
  Isovar 1.62.2 or the raw observation ledger.
- Viral/isoform products: sequence-supporting fragments, with coding-base
  coverage and whole-CDS segments retained separately.
- Ordinary protein exports: raw assembly support. Optional fractional catalog
  support and abundance estimates are retained without replacing raw counts.

New v4 exports carry `occurrences[].rna_support` separately from the candidate
union. Vaxrank validates its scope, identities and counts and uses only that
placement's fragments. A measured empty set is zero. Overlapping sets are not
independent evidence, and an unsupported alternative cannot borrow the union.

Historical v4 files without this field still retain multiple unmeasured placements
with `rna_fragment_support_unavailable`; a single occurrence can use its candidate
union. Use the raw v5 reconstruction or re-export with Isovar 1.62.2 to recover
placement support. Explicit inconsistent fields fail instead of falling back.
All ambiguity, partial-sequence and tumor-specificity gates still apply.
Admission opt-ins do not bypass evidence validation: contradictory source scopes
or measured counts are input errors even when ambiguous products are allowed.

## Reconstruction requests and controls

Replace `result` with exactly one `reconstruct` object:

```json
{
  "schema": "vaxrank.isovar_input.v1",
  "reconstruct": {
    "family": "variant_orfs",
    "inputs": {
      "bam": "tumor.bam",
      "sample_id": "tumor-1",
      "source": "RNA-library-1",
      "variant": {"contig": "1", "start": 12345, "ref": "C", "alt": "T"},
      "reference_orfs": [
        {"transcript_id": "TRANSCRIPT_ID", "catalog": "orfs.json.gz", "expected_sha256": "CATALOG_SHA256"}
      ]
    }
  },
  "options": {
    "reads": {"min_mapping_quality": 10},
    "protein": {"min_variant_sequence_coverage": 2, "protein_sequence_preference": "balanced"},
    "admission": {"allow_partial": true}
  }
}
```

This is a contract template: replace the locus/transcript/catalog placeholders
and add the same patient scope used above. Reference ORFs are validated against
the selected annotation. Catalogue discovery/downloads are not implicit.
Supply all competing models you intend to assess; an empty model list retains
an unresolved result instead of guessing a frame.

`inputs` may be an inline object or an existing JSON/JSON.gz file. Python accepts
the corresponding typed objects. CLI alignment requests use Isovar's indexed
BAM/CRAM opener; provide `num_rna_decompression_threads` in `inputs` to control
reading. Python callers can supply already-opened alignment objects, including
an explicit CRAM reference. Scratch/reference cache directories belong in
`inputs`, not in reconstruction controls.

| Family | CLI `inputs` contract |
| --- | --- |
| `variant_orfs` | `bam`, `sample_id`, `source`, `variant`, `reference_orfs`. Each reference entry names `transcript_id` plus either `catalog`/`expected_sha256`, `annotated: true`, or the source-defined `TranscriptORF.at` arguments (`start`, `orf_id`, `source`, optional expected `protein`, `translation_table`, etc.). The host annotation is required. |
| `sv` | `bam`, `source`, `event`. `event` is an inline object or relative file following Isovar's `sv_rna_input_from_dict` contract: event/sample/reference identity, oriented donor/acceptor, regions, references and event provenance. |
| `fusion` | Isovar's existing `fusion_from_dict` input (`fusion`, `references`, `reads`, optional `reference_names`). No new BAM assembly is implied. |
| `viral` | `bam`, `genome_fasta`, `sample_id`, `source`, optional `contigs`, `reference_cache`, `scratch_dir`. The FASTA must match the alignment reference; it never fills missing RNA. For a custom reference set `reconstruction.reference_profile: null` and explicit contigs. |
| `isoforms` | `bam`, `sample_id`, `source`, `targets`, optional `reference_identity`, `scratch_dir`. Targets are serialized `IsoformTarget` records with `FusionReference` objects; include annotated competitors. `isoform_targets_from_genome` can prepare these in Python. |

`options` has four sections. `admission` is the typed policy above;
`reconstruction` uses the selected public producer's optional JSON-valued
controls with their **actual installed defaults**. `reads` configures
`ReadCollector` for `variant_orfs`, `sv`, `viral` and `isoforms`; `protein` configures
`ProteinSequenceCreator` for `variant_orfs`. Required data objects cannot be
passed as controls. Unknown/inapplicable names and incompatible types fail.
Saved results reject reconstruction/read/protein controls rather than pretending
to rerun with them.

Common controls are listed below; see the corresponding Isovar public function's
docstring for its full definition. Vaxrank retains the resolved values and
producer package versions in every live result's evidence.

Viral and targeted-isoform reconstruction inherit
`ReadCollector.for_rna_paths()` defaults from Isovar: MAPQ-zero and secondary
alignments, and reads without base qualities, remain eligible; duplicate-marked
reads and soft-clipped bases are excluded. Override these through `options.reads`
to apply the same policy to discovery and complete-input support. The resolved
policy is recorded even when no override is supplied. A Python `read_collector`
object is not accepted in `inputs`; for callbacks or custom read-end profiles,
run Isovar directly and replay the resulting evidence.

| Section / producer | Useful controls |
| --- | --- |
| `reads` / `ReadCollector` | `min_mapping_quality`, `use_secondary_alignments`, `use_duplicate_reads`, `use_soft_clipped_bases`, `use_reads_without_base_qualities`, `merge_overlapping_fragments`, `infer_read_ends`, `trim_poly_a` |
| `protein` / `ProteinSequenceCreator` | `protein_sequence_length`, `protein_context_peptide_length`, `protein_sequence_preference`, `min_variant_sequence_coverage`, `min_protein_sequence_support_fraction`, `min_transcript_prefix_length`, `max_transcript_mismatches`, `count_mismatches_after_variant`, `variant_sequence_assembly`, `min_assembly_overlap_size` |
| `reconstruction` / SV | `min_anchor_bases`, `min_overlap`, `min_alternative_fragments`, `min_alternative_fraction`, `min_local_variant_fraction`, `max_records`, `max_queries`, `max_paths`, `max_extension_segments`, `assemble`, `include_regional_candidates`, `exon_anchored_clips`, `breakpoint_window`, `max_breakpoint_shift`, `annotated_junction_tolerance`, `min_orf_amino_acids`, `max_orf_candidates`, `dense_support`, `max_evidence_groups`, `inclusion_min_*`, `inclusion_mapq_255_is_unique`, `peptide_lengths` |
| `reconstruction` / fusion | `min_fragments`, `peptide_lengths` |
| `reconstruction` / viral | `reference_profile`, `max_discovery_reads`, `max_paths`, `min_overlap` |
| `reconstruction` / isoforms | `max_discovery_reads`, `max_paths`, `min_overlap` |
| `reconstruction` / protein export (Python) | `partial_read_support`, `support_max_edits`, `support_min_overlap` |

`variant_orfs` keeps every supplied reference model; a nonzero
`max_protein_sequences_per_variant` is rejected. Source-supported genetic codes
are validated where the producer supplies them. Fusion/SV/targeted-isoform
contracts currently use table 1. Python callbacks, custom adapter profiles and
external Salmon execution are not serialized into a source request. Call the
owning Isovar API directly for those operations, then adapt its result. Existing
`add_isoform_abundances` results retain their quantification and provenance.

Source-specific YAML and CLI overrides use the same names:

```yaml
isovar:
  sources:
    variant_orfs:
      admission:
        allow_partial: true
        allow_ambiguous: false
      protein:
        protein_sequence_preference: support
    sv:
      reconstruction:
        max_paths: 128
```

```bash
vaxrank --input-isovar source.json --config analysis.yaml \
  --set isovar.sources.variant_orfs.admission.min_fragments=3 \
  --ensembl-release 110 --mhc-predictor random --output-dir review
```

`--set` aliases the existing `--config-value`. Effective precedence, from low
to high, is producer/policy defaults, source-file options, family options in
YAML with CLI overrides applied, then per-entry manifest options. This allows a
specific manifest entry to declare an intentional exception to a family policy.
The resolved policy is always retained. Family options apply only to inputs of
that family; ordinary VCF/BAM RNA-context controls remain under their existing
`isovar` keys.

## Python and replay

```python
from vaxrank import (
    IsovarAdmission, IsovarSourceOptions,
    antigens_from_isovar, reconstruct_isovar_antigens,
)

# Existing producer result; no annotation or reconstruction is performed.
antigens = antigens_from_isovar(
    result, scope=patient_scope,
    candidate_attestations=reviewed_candidates,
    options=IsovarSourceOptions(admission=IsovarAdmission(allow_partial=True)),
)
antigens.save("antigens.json")

# Typed reference models/read evidence; Isovar owns reconstruction.
antigens = reconstruct_isovar_antigens(
    "variant_orfs",
    inputs=dict(variant=variant, read_evidence=read_evidence,
                reference_orfs=reference_orfs, sample_id="tumor-1", source="RNA-1"),
    options={"protein": {"protein_sequence_preference": "support"}},
    scope=patient_scope, candidate_attestations=reviewed_candidates,
)
predictions = antigens.predict(topiary_predictor, genome=host_genome)
```

Use `load_isovar_input(path, genome=...)` to load the exact CLI source contract
from Python. Predicting new sequences requires the same host annotation or
frozen `HostReference` and MHC repertoire as other antigen inputs.

`AntigenDataset.source_evidence` stores each complete producer result once,
keyed by its content identity, plus admission decisions and live reconstruction
settings/versions. Antigen metadata links to that identity and candidate/model/
path/translation identifiers. Native prediction exports retain this evidence,
including unmeasured proteins, reference-only models and inputs with no proteins.
RNA witnesses, reference comparisons, producer limitations, genetic codes and
phase/abundance evidence are not replaced by summary counts.

`--output-antigens` exports the resolved catalog and source evidence;
`--output-epitopes` retains predictions and design state for offline replay.
Reports and run summaries include Isovar source outcomes and held-out reasons.
Replaying native predictions needs neither the original BAM nor Isovar reference
catalogues. Re-predicting an antigen export still needs a compatible host-self
reference. A no-protein result produces an empty candidate catalog with evidence,
not a fabricated negative-expression or viral-absence claim.

## Evidence interpretation and validation

RNA-backed coding predictions do not establish translation or tumor specificity.
This separation is consistent with studies distinguishing noncanonical ORF
prediction from translation/functional evidence ([Prensner et al., 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8195866/)).
Genetic-code validation follows the producer's declared translation table and
[NCBI's CDS genetic-code conventions](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi).

Offline regression fixtures pin Isovar 1.61.0 outputs, including the complete
COLO829 DYNLL2 result with 300 alternatives, synthetic ncORFs, SV paths, viral
coding and competing isoforms. The viral fixture is synthetic RNA from a public
reference; it is **not** the original SiHa sample. A bounded real SiHa fixture
remains tracked in [Isovar #531](https://github.com/openvax/isovar/issues/531) and
[Vaxrank #656](https://github.com/openvax/vaxrank/issues/656).
