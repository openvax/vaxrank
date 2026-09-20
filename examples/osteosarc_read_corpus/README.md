# Expanded osteosarc regional read corpus

This maintainer workflow downloads an **analysis dataset**, not a read-capped
test fixture. The original 49 source–variant cases cover 44 unique vaccine loci.
Those remain included alongside additional coding SNVs/indels, structural
hypotheses and RNA-linked companion alleles. A source-reported candidate is not
automatically a proven somatic variant or a usable antigen.

The current local run is
`output/osteosarc-read-corpus/runs/2026-09-18T112041Z/`.
Its consolidated deliverable is `dataset/REPORT.md`, with machine-readable
events, source inventory, coverage and acquisition receipts beside it. Inspect
the acquisition statuses before using a run: an interrupted download is not
zero coverage.

## Scope

- Tumor bulk/single-cell Illumina RNA, ONT and PacBio across all clinical
  timepoints identified in the inspected inventory; tumor/normal WGS/WES
  controls for germline/somatic attribution. Pooled blood scRNA is excluded
  because patient-specific donor assignment is unresolved.
- Selected representatives and tagged long-read evidence are queried around every eligible event:
  ±2 kb, both SV endpoints, full stored read sequences and paired mates.
- No ALT, quality, duplicate, secondary/supplementary or read-count filtering.
  Unpaired long-read files use a single pass, with paired-flag verification and
  mate retrieval if paired records unexpectedly occur.
- Native GRCh37 coordinates are uniquely chain-mapped and reference-verified;
  unresolved mappings remain explicit. A chromosome-name prefix is not an
  assembly identifier.
- PacBio T1's **entire** published genomic BAM and its complete deduplicated
  unaligned input are cached. Exact molecule accounting and whole-genome
  remapping of omitted molecules audit loss at the alignment stage. This does
  not audit losses in upstream segmentation, primer filtering or deduplication.
- Read-processing products are not biological replicates. Do not pool tagged,
  deduplicated and reprocessed copies to inflate evidence.
- User-selected policy: prioritize library/timepoint coverage, retain distinct
  providers/preparations and DNA controls, and keep processing alternatives
  inventoried. `library_groups.json` records reviewed acquisition families;
  `source-selection.json` pins the decisions for a run. Uncertain relationships
  are not proof of identical libraries. Already cached alternatives are retained.

Regional subsets do not recover hard-clipped bases, unqueried supplementary
partners, all unmapped reads or entire transcript loci. Keep the full original
URLs to expand regions or test alternative alignments when necessary. The
symbolic MUC3A/GAPVD1 entries are anchor queries, not invented literal alleles.

## Reproduce or resume

Requirements: Vaxrank's dependencies, including **osteosarc==0.1.1**;
samtools with HTTPS and `--fetch-pairs` support (CI pins **1.21**); curl; minimap2 for
the optional PacBio audit; local indexed full GRCh38 and GRCh37 references; an
Isovar checkout containing `tests/data/osteosarc/expansion/` and the pinned
osteosarc bucket listing. Pass explicit local paths below. Large files stay
under `output/`, outside version control.

Downloads, header inspection and regional extraction use osteosarc's public
`Cache`, `inspect_alignment` and `extract_reads` APIs. `--cache-root` selects an
explicit shared cache; otherwise osteosarc uses `OSTEOSARC_CACHE`, then
`OPENVAX_DATA_CACHE`, then the platform's OpenVax cache. Downloaded objects use
`objects/sha256/`; extraction receipts and immutable derivatives live under
`osteosarc/`. Verified older whole-file downloads can be imported with their
source URL and checksum instead of downloading them again. `--offline` permits
only cached remote objects/derivatives and explicit local inputs.

Each run still owns its reviewed event panel and source-selection policy.
This acquisition migration does not load or apply the package's allele
corrections implicitly. The published MAP2 fixture and the separately corrected
allele regression are documented in the
[selection pilot](../../tests/data/osteosarc/selection_validation/README.md).

```sh
CORPUS_RUN="output/osteosarc-read-corpus/runs/$(date -u +%Y-%m-%dT%H%M%SZ)"
python3 examples/osteosarc_read_corpus/build.py prepare --run "$CORPUS_RUN" \
  --isovar-repo /path/to/isovar --listing /path/to/bucket-listing-compressed.json \
  --reference /path/to/GRCh38.fa
python3 examples/osteosarc_read_corpus/build.py lift --run "$CORPUS_RUN" \
  --isovar-repo /path/to/isovar --reference /path/to/GRCh37.fa
python3 examples/osteosarc_read_corpus/build.py hg19-mito --run "$CORPUS_RUN" \
  --isovar-repo /path/to/isovar
python3 examples/osteosarc_read_corpus/build.py supplement --run "$CORPUS_RUN" \
  --reference /path/to/GRCh38.fa
python3 examples/osteosarc_read_corpus/build.py lift --run "$CORPUS_RUN/phase-and-sv-context" \
  --isovar-repo /path/to/isovar --reference /path/to/GRCh37.fa
python3 examples/osteosarc_read_corpus/select_sources.py --run "$CORPUS_RUN"
python3 examples/osteosarc_read_corpus/build.py acquire --run "$CORPUS_RUN" --workers 6
python3 examples/osteosarc_read_corpus/build.py acquire --run "$CORPUS_RUN/phase-and-sv-context" --workers 4
python3 examples/osteosarc_read_corpus/build.py pacbio-full --run "$CORPUS_RUN"
python3 examples/osteosarc_read_corpus/audit_pacbio.py --run "$CORPUS_RUN"
python3 examples/osteosarc_read_corpus/audit_pacbio.py --run "$CORPUS_RUN" \
  --realign --reference /path/to/GRCh38.fa
python3 examples/osteosarc_read_corpus/finalize.py --run "$CORPUS_RUN"
python3 examples/osteosarc_read_corpus/quality_audit.py --run "$CORPUS_RUN"
```

Rerun `acquire` to retry failures: verified completed BAMs are reused, and
earlier failed receipts/incomplete files are retained. `--source-id` can target
specific retries within the pinned selection; it does not override deferred
products. Deliberately revise the pinned policy to acquire an alternative.
Do not run two acquisitions of the same source/part at once.
Do not change an existing run's event definitions; start a timestamped new run.

Start a **new run directory** when moving an older direct-samtools acquisition
to this implementation. New receipts pin the osteosarc version, source/index,
event panel, coordinate mapping and snapshot identity. Legacy acquisition
receipts do not establish all of those identities and are not silently reused.
Resuming a new acquisition validates both the exported BAM and index and the
osteosarc derivative; modified evidence fails explicitly. Failed-attempt
receipts are retained. An `integrity_error` requires a new run; retries preserve
the damaged evidence and its failure receipt instead of repairing it silently.
Interrupted whole-file transfers restart through
osteosarc's cache; this is not byte-range resumption.

Offline integration tests compare complete original read records (including
floating-point tag bits), recover off-panel mates, preserve native hg19
mitochondrial coordinates and reject damaged cached evidence. The five real
selection BAMs also pass through this acquisition path before reconstruction
and real cached MHC prediction: all seven documented comparisons retain the
same selections and RNA gates. A tiny selected fixture does not establish
complete RNA coverage or historical ranking-protocol reproduction.

`finalize.py` combines the two acquisition parts per source using the maximum
multiplicity of each identical SAM record across parts. This preserves true
duplicates in the source while avoiding duplicate counting from overlapping
queries. Use the resulting `dataset/alignments/`, not a concatenation of the
parts. Missing parts remain explicitly incomplete.

Coverage only establishes read availability. Isovar assembly/translation,
Varcode comparison, matched-normal genotype assessment, direct phase evidence
and Vaxrank ranking are downstream analyses with their own evidence gates.
Event-directed RNA reconciliation of large/ambiguous variants is tracked in
[Isovar #306](https://github.com/openvax/isovar/issues/306), including biological
soft-clip evidence, competing transcript paths and unresolved coding frames.
When RNA is insufficient, preserve DNA-only candidates explicitly; see
[Vaxrank #482](https://github.com/openvax/vaxrank/issues/482).

The run-local `software/isovar-1.18.1/` installation is isolated from the shared
environment. Set `PYTHONPATH` to that directory and run `quality_audit.py` with
`--label isovar-1.18.1` to reproduce the missing-QUAL comparison. Rerun
`finalize.py` to include both audits in the report. Do not replace unknown
quality values with fabricated scores.

Sources: [osteosarc data guide](https://osteosarc.com/data/),
[BAM catalogue](https://osteosarc.com/bams/),
[variant catalogue](https://osteosarc.com/variants/),
[fusion catalogue](https://osteosarc.com/fusions/),
[pbmm2 alignment behavior](https://github.com/PacificBiosciences/pbmm2),
[minimap2 splice alignment](https://github.com/lh3/minimap2).
