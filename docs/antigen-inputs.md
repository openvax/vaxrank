# Experimental unified antigen inputs

Mutation, fusion, CTA, ERV, splice and viral antigens can use the same
`AntigenDataset` API and `--input-antigens` CLI. Each record carries a
`VaccineAntigen`: sequence, targetable intervals, source identity and an
admission decision and its evidence.

Callers supply resolved sequences and admission evidence. Held-out records stay
in the evidence but cannot enter constructs; unresolved sequences stay upstream.

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
(about 14 MB of downloads):

```sh
python examples/build_canine_host_reference.py dog.proteome.json.gz
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
Mouse CTA evidence is currently unavailable; requesting it gives a clear error.

An optional `evidence_snapshot` path selects an exported OncoRef snapshot,
relative to the manifest. Python uses the same adapter:
`admit_cta_expression(..., host_reference=host, overrides={gene_id: evidence})`,
with `CTAOverrideEvidence` records. Saved results replay without current OncoRef
queries or the original expression table.

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
Human CTA expression still uses annotation for transcript/sequence selection;
mapped canine CTA inputs select sequences directly from the snapshot.
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
