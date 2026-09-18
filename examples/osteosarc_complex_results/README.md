# Osteosarcoma complex-variant Vaxrank results

This result set follows complex variants through four separate claims: the DNA
event, sample-specific RNA evidence, translated protein, and Vaxrank peptide
selection. It includes positive, negative, ambiguous, and unresolved outcomes;
an RNA-supported event is never promoted to a peptide unless its coding frame
and translation are defensible.

The expanded comparison contains:

- assembly-dependent frameshifts: PIP5K1A and TECPR1;
- cross-platform and validated complex-indel controls: GTF3C5, GLIS3,
  RNF213, and H1-2;
- both DYNC1H1 loci and the multi-transcript NAV2 result;
- explicit RNA compound haplotypes for MAP2 and NTF3;
- long-read phasing/context results for CD109 and ZNF436;
- a matched ONT/Illumina KTN1 reconstruction where both platforms encode the
  same frameshift tail but only ONT supplies a standard 25-aa construct; and
- fusion/SV audits for ATP5MG::KMT2A, TPST1::CRCP,
  FOXO3::STRADA/CCDC47, PARD3B::CDKN2B, the AMPH internal deletion, and a
  long-read-rich unnamed chr21 junction, GABBR1::SLC29A1,
  OTUD7A::FMN1, DLG5, AFF3, and KEAP1.

`source/assembled_antigens.json` is the machine-readable bridge between RNA
evidence and Vaxrank. It records exact translated sequence, targetable amino
acid intervals, all assembled transcript identifiers, evidence counts,
platform/source provenance, and the sequence published by osteosarc.com when
available. The site currently publishes protein or vaccine sequence for
TECPR1, MAP2, DYNC1H1 p.Val314Ile, NAV2, GLIS3, RNF213, and ZNF436. The
DYNC1H1 p.Gln3267His, NTF3, and CD109 pages publish the variant annotation but
not a protein context, so those displayed sequences remain attributed to RNA
evidence rather than to the website.

`source/sv_audits.json` retains junction nucleotide sequences and reasons for
withholding translation. DLG5, AFF3, and KEAP1 now have dedicated
transcript-structure panels. It also inventories additional calls that are not
ready for panels: GAPVD1, MUC3A, FAM157A, MYO15B, SPRED1, ITM2B::RB1, and
internal PTPRD, EYS, and EDA structural junctions.
MYO15B includes the 17-aa vaccine peptide published by osteosarc.com, but it
is not reranked because the page does not establish sample-specific RNA
translation provenance.

`source/additional_candidate_audit.json` preserves the exact GTF3C5, RNF213,
GLIS3, and KTN1 alleles and per-BAM fragment counts, the four source BAM URLs,
validated tagged-ONT paths for GABBR1::SLC29A1 and OTUD7A::FMN1, and the
original DLG5, AFF3, and KEAP1 lead inventory. Its checksum is included in every
generated flattened result. PAVE annotations and the transcript-specific
annotations used by the assembled panels remain separate where they differ.

`source/platform_comparison_audit.json` records matched regional ONT and
Illumina reconstruction for KTN1, DLG5, AFF3, and KEAP1. It separates four
possible platform effects: a longer construct-eligible protein context, longer
single-molecule exon phasing, agreement on an unchanged coding splice, and
platform-specific evidence for a noncoding isoform. Counts are never treated as
a direct sensitivity comparison because the ONT single-cell and Illumina bulk
libraries differ in preparation, depth, and deduplication.

`source/orf_platform_audit.json` and `orf-platform-report.md` summarize the
pinned 44-variant × 164-product Isovar audit by ILMN, ONT, and PacBio. They
distinguish a validated local coding window from a full-length ORF, compare
the same 49 checksummed original-read subsets with assembly on and off, and
retain every non-exact RNA/isolated-effect sequence. In the full audit, ILMN,
ONT, and the single genomic-coordinate PacBio product establish local protein
windows for 39, 30, and 3 of 44 loci, respectively. The paired corpus has the
same ORF-availability count with assembly on and off, while assembly extends
the returned ILMN protein context in five cases. These denominators are
reported separately because products are not replicates and the paired corpus
is deliberately enriched.

The audit also records the present attribution boundary: Isovar can phase
multiple supplied somatic variants through shared fragment names, but the 44
nominated loci contain no nearby pair. The public `run_isovar` path does not
accept a matched germline variant set, so additional transcript edits remain
unexplained rather than being guessed germline. NTF3's adjacent RNA base is
sequence-resolved across ILMN and ONT, but its germline/somatic origin remains
unresolved.

`rank_assembled_antigens.py` admits only inputs whose evidence gate is `pass`,
constructs explicit `VaccineAntigen` objects, and runs NetMHCpan through the
normal Vaxrank/Topiary prediction path. These assembled-antigen windows use a
source-agnostic target epitope score; RNA evidence is displayed separately and
is not fabricated as a single-Varcode-variant read count. Original Isovar
mutation panels retain their historical RNA-weighted score.

The explicit assembled-antigen rankings are exploratory binding ranks: no full
human reference proteome was supplied, so they do not claim complete exact-self
screening. That limitation is repeated in the combined figure and recorded in
the generated JSON.

The five assessed patient class-I alleles are recorded in the source JSON. The
clinical null allele HLA-A*01:11N is deliberately not predicted.

Rebuild predictions, the flattened result source, and a timestamped combined
PDF plus per-page SVG and 3600x2280 PNG files:

```bash
PATH=/path/to/netmhc-bundle/bin:$PATH \
  python examples/osteosarc_complex_results/rank_assembled_antigens.py \
  --output examples/osteosarc_complex_results/source/assembled_rankings.json
python examples/osteosarc_complex_results/audit_platform_orfs.py \
  --isovar-repo /path/to/isovar \
  --cache-directory /path/to/local/reference-cache \
  --output examples/osteosarc_complex_results/source/orf_platform_audit.json \
  --markdown-output examples/osteosarc_complex_results/orf-platform-report.md
python examples/osteosarc_complex_results/build_results.py
python examples/osteosarc_complex_results/generate.py \
  --combined-output output/pdf/vaxrank-all-figures.pdf
```

The figure command performs no prediction and is reproducible without the
proprietary predictor. Every variant directory contains the exact record shown
on that page. The run manifest stores source checksums, software version, page
order, raster dimensions, and combined-PDF path.
