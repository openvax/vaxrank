# Osteosarcoma complex-variant Vaxrank results

This result set follows complex variants through four separate claims: the DNA
event, sample-specific RNA evidence, translated protein, and Vaxrank peptide
selection. It includes positive, negative, ambiguous, and unresolved outcomes;
an RNA-supported event is never promoted to a peptide unless its coding frame
and translation are defensible.

The expanded comparison contains:

- assembly-dependent frameshifts: PIP5K1A and TECPR1;
- validated complex-indel controls: GLIS3, RNF213, and H1-2;
- both DYNC1H1 loci and the multi-transcript NAV2 result;
- explicit RNA compound haplotypes for MAP2 and NTF3;
- long-read phasing/context results for CD109 and ZNF436; and
- fusion/SV audits for ATP5MG::KMT2A, TPST1::CRCP,
  FOXO3::STRADA/CCDC47, PARD3B::CDKN2B, the AMPH internal deletion, and a
  long-read-rich unnamed chr21 junction.

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
withholding translation. It also inventories additional calls that are not
ready for panels: GAPVD1, MUC3A, FAM157A, MYO15B, SPRED1, ITM2B::RB1,
GABBR1::SLC29A1, and internal PTPRD, EYS, and EDA structural junctions.
MYO15B includes the 17-aa vaccine peptide published by osteosarc.com, but it
is not reranked because the page does not establish sample-specific RNA
translation provenance.

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
python examples/osteosarc_complex_results/build_results.py
python examples/osteosarc_complex_results/generate.py \
  --combined-output output/pdf/vaxrank-all-figures.pdf
```

The figure command performs no prediction and is reproducible without the
proprietary predictor. Every variant directory contains the exact record shown
on that page. The run manifest stores source checksums, software version, page
order, raster dimensions, and combined-PDF path.
