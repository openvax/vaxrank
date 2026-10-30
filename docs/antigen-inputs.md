# Experimental unified antigen inputs

Mutation, fusion, SV, isoform, CTA, ERV, splice and viral antigens can use the same
`AntigenDataset` API and `--input-antigens` CLI. Each record carries a
`VaccineAntigen`: sequence, targetable intervals, source identity and an
admission decision and its evidence.

Callers supply resolved sequences and admission evidence. Held-out records stay
in the evidence but cannot enter constructs. [Isovar source inputs](isovar-inputs.md)
can also reconstruct ncORFs, SV ORFs, viral coding and targeted isoforms directly,
retaining the producer evidence and unresolved/no-protein outcomes in
`AntigenDataset.source_evidence`.

## Prepare and predict

Given `mutation` and `viral` objects from your evidence adapters:

```python
from vaxrank import AntigenDataset
from vaxrank.input_scope import InputScope

inputs = AntigenDataset(
    antigens={"mutation-1": mutation, "viral-1": viral},
    scope=InputScope(
        patient_id="patient-001", reference_assembly="GRCh38",
        annotation="ensembl:93", mhc_alleles=("HLA-A*02:01",),
    ),
)
inputs.save("antigens.json")
```

```sh
vaxrank --input-antigens antigens.json \
  --mhc-predictor mhcflurry --ensembl-release 93 \
  --vaccine-type peptide mrna --output-dir design \
  --output-antigens design/antigens.json
```

The annotation describes the host; each antigen keeps its own source species
and accession. New predictions require an installed host annotation for
exact-self assessment. Conflicting reference or genotype declarations fail.
The Python equivalent is `inputs.predict(predictor, genome=genome)`, using a
Topiary predictor. It returns the shared `EpitopeDataset`.

Reuse contextual Topiary predictions while retaining live coverage for cache
misses:

```sh
vaxrank --input-antigens antigens.json \
  --mhc-predictor mhcflurry --prediction-cache prior-predictions.tsv \
  --ensembl-release 93 --output-dir design
```

The cache must describe the same single model/version as the configured live
predictor. Lookup uses each native protein occurrence's exact neighbours,
endpoint and validated patient genotype. A peptide measured beside different
neighbours is a miss, not a reusable score; the live model fills that occurrence
and Topiary validates its returned allele, endpoint and version before merging
it. Native output retains the completed evidence and replays without either the
cache or predictor.

Exact-self lookup shares repeated peptide queries within the same antigen kind
and gene-exclusion policy. Overlapping ORF windows therefore do not each rescan
the complete host proteome. Every source retains its own candidate identity,
admission evidence and full matching gene/transcript/protein provenance; shared
queries do not pool RNA support or merge candidates. Explicit host-reference
snapshots and custom matching callbacks retain their own matching behavior.

## Direct candidate API

For an admitted antigen and its scored `CandidateEpitope` records:

```python
from vaxrank.vaccine_peptide import VaccinePeptide

candidate = VaccinePeptide(antigen=antigen, epitopes=epitopes)
```

With no mutation fragment, the default combined score is `target_epitope_score`;
ranking uses target epitope score, manufacturability, then self epitope score.
These are the same defaults used for antigen-only candidates through the CLI and
`vaccine_peptides_for_antigen`. No fabricated RNA counts are required. The source
kind alone does not imply that single-variant RNA support is available.

Candidates with a `mutant_protein_fragment` retain the mutation default
`sqrt(n_rna_alt) * target_epitope_score` and its RNA/position ranking rules.
Explicit `combined_score_expr` and `ranking_rules` are preserved. Referencing
mutation-only fields without a fragment raises an error. Resolved defaults and
explicit overrides survive JSON serialization. Admission and self-reference
requirements still apply to direct construction.

## Combine new sequences and saved evidence

A manifest can mix new sequences with saved predictions. Paths are relative
to the manifest:

```yaml
schema: vaxrank.input_manifest.v1
patient_id: patient-001
reference_assembly: GRCh38
annotation: ensembl:93
mhc_alleles: [HLA-A*02:01]
inputs:
  - format: epitopes
    path: mutations-and-cta.tsv
  - format: antigens
    path: viral-antigens.json
```

```sh
vaxrank --input-manifest inputs.yaml \
  --mhc-predictor mhcflurry --ensembl-release 93 \
  --vaccine-type peptide mrna --output-dir mixed-design
```

Only the `antigens` entries receive new predictions. Saved measurements stay
unchanged. All declarations are checked before loading models; a manifest may
fill missing scope but cannot override conflicting evidence. Duplicate files
and loading an antigen input alongside its native prediction export are rejected.

When each file already declares its scope, repeat `--external-input
antigens=PATH` or combine `--input-antigens` with `--input-epitopes` instead.
The Python `vaxrank.external_rescoring.prepare_reports` API accepts the same
`("antigens", path)` inputs and an `antigen_predictor_factory(scope)` returning
a Topiary predictor. The factory runs once after validation.

Direct VCF/BAM still runs its source adapter first; add
its saved candidate file to this manifest. Fresh/additive table rescoring
is a separate step. Pooling does not calibrate different models or evidence
types: use an explicit [selection policy](mixed-antigen-design.md).

## Add CTA expression

Add a `cta_expression` entry beside `antigens` or `epitopes` to select CTA targets
from patient expression before prediction:

```yaml
  - format: cta_expression
    path: expression.tsv
    sample_id: tumor-rna
    cta_expression:
      measurement_level: gene
      id_column: gene_id
      value_column: tumor_tpm
      expression_unit: TPM
      evidence_source: Salmon
      evidence_version: "1.10.3"
      assay: bulk RNA-seq
      min_expression: 2
```

`format: cta_expression` selects the CTA expression adapter. The file can be CSV
or TSV; the nested `cta_expression` settings describe its columns and units.
Use `transcript` for transcript measurements. Each entry declares its own
columns, units, sample and threshold; one matrix may supply several sample
columns. Existing human Ensembl/OncoRef admission rules apply. Expression is
not converted into protein abundance or mutation read support.

Optional `exclusion_policy` lists `gene_patterns` and `gene_exceptions`;
defaults remain `MAGE*` and `MAGEA4`. An optional `hitlist_evidence_bundle` path
is relative to the manifest and is verified against that input's contract.

`--output-dir` writes `cta_admission.json` for one CTA source, or separate files
under `cta_admissions/` for multiple sources. `--output-cta-admission` selects a
single-source file explicitly. Held-out and unresolved decisions remain in
the native evidence, even when a CTA source contributes no admitted proteins
to an otherwise usable mixed run. A wholly unadmitted CTA-only run stops after
writing requested admission records.

### Species-specific CTA expression (experimental)

For dogs, declare `host_species: dog` and a `host_reference` snapshot. The same
`cta_expression` entry accepts NCBI Gene IDs or versioned RefSeq transcript IDs.
It uses OncoRef's canine evidence mapped to RefSeq annotation 106 /
GCF_011100685.1. From a source checkout or extracted source distribution, build
the complete annotated reference from its pinned protein FASTA and feature table
(the transfer budget below permits up to 100 MB):

```sh
python examples/build_cta_host_reference.py dog dog.proteome.json.gz --max-download-bytes 100000000
```

The example uses `Proteome.from_refseq()` and OncoRef's verified public assets.
Vaxrank checks both source hashes, the annotation, identifier namespaces and the
complete protein inventory before prediction. No Ensembl release is needed.

Canine candidates stay held out until an explicit, gene-specific evidence review
is supplied beside `min_expression`:

```yaml
      overrides:
        "477552":
          evidence_source: project target review
          evidence_version: "1"
          evidence_id: review-PRAME
          rationale: Recorded justification for this research target
```

An override cannot bypass the expression threshold or an unresolved mapping.
Proteins with unknown residues stay unresolved, with their accession and sequence
checksum recorded; other targets continue. Gene measurements do not switch to a
different isoform to conceal an unresolved sequence.
Human `MAGE*` exclusions do not apply to dogs. Self checks omit only the reviewed
target gene; other genes and extra ORFs remain searchable. OncoRef's partial
coverage, family ambiguity and normal-tissue observations stay in native outputs.
Its canine evidence includes [kidney expression for CTCFL and DDX53](https://doi.org/10.1016/j.vetimm.2025.110946).
Mouse P1A evidence is available, including normal-tissue observations and strain
provenance. Its single-gene sequence seed is not a complete mouse host inventory;
using it as that inventory fails before acquisition or prediction.

An optional `evidence_snapshot` path selects an exported OncoRef snapshot,
relative to the manifest. Python uses the same adapter:
`admit_cta_expression(..., host_reference=host, overrides={gene_id: evidence})`,
with `CTAOverrideEvidence` records. Saved results replay without current OncoRef
queries or the original expression table.

### Common RNA restriction policy (experimental)

An explicit `restriction_policy` selects the same admission contract for human,
dog and mouse. Install `vaxrank[cta]` to resolve the required OncoRef 1.23.1 and
ORFdata 0.42.6 together; see the [dependency constraints](installation.md#common-cta-evidence). It
uses the snapshot's identifier namespaces: Ensembl gene IDs and exact versioned
transcripts, or NCBI Gene IDs and exact versioned RefSeq transcripts. It does not
use gene symbols, human orthology, or human default CTA exclusions to admit targets.

Build the annotated inventory with the same API for each supported reference:

```python
from vaxrank.cta_reference_adapter import build_cta_host_reference

host = build_cta_host_reference("human", max_bytes=300_000_000)
host.save("human.proteome.json.gz")
```

The equivalent example is
`python examples/build_cta_host_reference.py human human.proteome.json.gz --max-download-bytes 300000000`.
The pinned human Ensembl source set is about 288 MB; the resulting snapshot is
about 45 MB. The builder streams source files without annotation indexes and
checks the total transfer budget before downloading. `offline=True, max_bytes=0`
reuses verified local files; `local_overrides` accepts OncoRef's species-qualified
asset paths. This produces an annotated-only inventory. Compose extra ORFs with
ORFdata when available, and retain their selection and coverage in that snapshot.
Ensembl inventories use PyEnsembl's assembly/annotation labels, checked against
the owner's exact assembly accession and release, so a matching variant input
can share the reference without relabeling its provenance.

Python supplies `CTARestrictionPolicy` and `CTARestrictionEvidence` from
`vaxrank.cta_restriction` to the existing entry point:

```python
result = admit_cta_expression(
    expression_path, input_contract=input_contract, admission_policy=patient_policy,
    host_reference=host, evidence_snapshot=snapshot_path,
    restriction_policy=normal_rna_policy, restriction_evidence=normal_rna_evidence,
)
```

In a manifest, `cta_expression.restriction_policy` contains the policy's JSON
fields and `cta_expression.restriction_evidence` names a JSON file relative to
the manifest. Export a typed bundle with `json.dumps(dataclasses.asdict(bundle))`;
read it with `CTARestrictionEvidence.read(path)`. An `evidence_snapshot` path
selects the exact OncoRef snapshot to which that bundle is bound.

The caller must declare a policy ID/version, native RNA unit/scale, somatic upper
and reproductive lower bounds, required tissues, independent donor/sample minima,
and a rationale. Vaxrank supplies no validated biological thresholds. Observations
retain exact protein/transcript occurrences, quantitative bounds, allocation,
tissue role, donor/sample IDs, pooling, source checksums and locators, study lineage,
and exclusion reasons. Donor and sample IDs must be reconciled across studies;
repeated publications and technical replicates cannot be counted as new donors.

Every full-protein equivalent occurrence in the host must have independently
allocated transcript RNA coverage. Missing or ambiguous RNA, pooled negatives,
unknown transcript mappings, and absent required tissues remain unresolved.
Gene-level patient RNA remains gene-level RNA; it selects a reference protein
without asserting isoform abundance or malignant-cell expression. Native OncoRef
facts remain attached. Compatible, verified numeric transcript RNA can be compared
with the policy, but never replaces independently allocated coverage. Incompatible
or qualitative somatic facts require review. Exclusions cannot erase observed
somatic counterevidence at the target or an identical protein at another locus.

A target is admitted only if both patient expression and the restriction policy
pass. `overrides` cannot be combined with this policy. Only the target's annotated
gene is excluded from self matching; other loci and all extra ORFs remain. Native
outputs freeze the policy, full evidence and findings for offline replay. A changed
reference, omitted equivalent occurrence or altered decision is rejected.

Admission remains a draft research selection requiring review. It does not
establish binding, presentation, processing or immunogenicity. Existing allele
capability checks and audits of every final peptide/mRNA product window still
apply. The released dog and mouse candidate facts alone do not constitute a
validated normal-tissue restriction panel. Existing human defaults and saved
historical decisions retain their recorded semantics.

## Declare the host species (experimental)

Declare a host by name. For an ordinary mouse, this is enough:

```yaml
host_species: mouse
```

For an HLA mouse, add the intended MHC species:

```yaml
host_species: mouse
mhc_species: human
```

Only use a list when targeting multiple MHC species: `mhc_species: [human, mouse]`.
Both fields accept common names, scientific names or NCBI taxonomy IDs;
`host_taxon_id` remains an alias. Equivalent declarations normalize to the same
stored identity. The CLI equivalent is `--host-species mouse --mhc-species human`.
The API accepts the same declarations through `InputScope.from_dict(...)`.

The host must agree with annotation; declared and observed alleles must match
the intended MHC species. Antigen source species remain separate, allowing
viral sequences in a human or mouse design. Declarations survive native replay.
For fresh mouse predictions, select mouse annotation with
`--genome GRCm38 --ensembl-release 93` (or another matching build/release).

Producer tables may leave host cells blank. Integer-valued numeric taxon cells
are normalized when read; fractional IDs and booleans are rejected. Manifest
and API taxon IDs must be positive integers.

Donor/transgene self layers remain tracked in
[#625](https://github.com/openvax/vaxrank/issues/625).
Inputs without new scope fields keep their historical identities; known host
annotation still enables species checks.

## Host proteome snapshots (experimental)

Use an ORFdata snapshot to include selected uORFs and other noncanonical ORFs:

```yaml
host_reference: host.proteome.json.gz
inputs:
  - format: antigens
    path: antigens.json
```

Add this to a manifest with `schema: vaxrank.input_manifest.v1`. Paths are
relative to the manifest. Standalone inputs accept `--host-reference FILE`.
Build the snapshot with ORFdata's `Proteome.build(...).save(...)`; its selected
evidence, source hashes, exclusions and every mapped occurrence are retained.
An optional `{path: FILE, identity: SHA256}` declaration pins its content identity.

The snapshot must match the declared host and annotation. It supplies exact-self
assessment for new native antigen predictions without a live Ensembl annotation.
The historical human CTA path uses annotation for transcript/sequence selection.
Mapped nonhuman inputs and the opt-in common restriction policy select sequences
directly from the snapshot using its declared namespaces.
Custom assembly names and versioned assembly accessions use the snapshot's
reference metadata directly.
CTA gene exclusions remove annotated proteins, leaving additional ORFs at those
genes searchable. Viral entries never become host self.

Native outputs bundle the snapshot in their `.references/` directory. Move that
directory with the table for offline replay. Changed or missing snapshots fail.
Saved predictions require their original reference; select a different reference
by predicting the source antigens again. Legacy producer tables and direct
VCF/BAM runs retain their existing annotated reference path.

Python accepts a `HostReference`, ORFdata `Proteome`, or snapshot path through
`AntigenDataset.predict(..., host_reference=...)` and
`audit_construct_sequence(..., host_reference=...)`.
For an HLA mouse, also pass `mhc_species="human"` to the audit, whether using
a snapshot or a live genome.

## Shared outputs

The ordinary ranking, peptide/mRNA assembly and reports apply to every kind.
Reports distinguish resolved, admitted and held-out records by source kind.
`design/candidate_predictions.tsv` retains predictions, policies and actual
products for replay:

```sh
vaxrank --input-epitopes design/candidate_predictions.tsv --output-dir replay
```

`--output-antigens` exports resolved native antigen records in the input JSON
format. It preserves admission evidence and input history, including held-out
and unmeasured records. It excludes prediction values and design policies;
keep the native candidate file for replay. This export also supports native mutation and CTA expression workflows.

Pooled native files retain each input's CTA admission decisions and captured
Hitlist evidence in `dataset.selection["source_evidence"]`, including unresolved
features. A new mixed design requires new product audits.

See [mixed-source selection](mixed-antigen-design.md) for evidence requirements
and current source-allocation limits.
