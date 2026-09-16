# Osteosarcoma complex-variant Vaxrank results

This compact result set follows complex variants through four distinct claims:
the DNA event, sample-specific RNA reconstruction, translated mutant protein,
and Vaxrank vaccine selection. It intentionally includes positive, negative,
ambiguous, and unresolved outcomes.

- **PIP5K1A p.G474fs:** long-read RNA assembly establishes a 49-aa frameshift
  protein. Vaxrank and NetMHCpan select the 25-aa long peptide
  `SSFSRRAAPVATPALLTSHRSLGNT`.
- **H1-2 15-nt deletion:** assembly succeeds and passes the RNA gates, but no
  target epitope passes default filtering. No vaccine peptide is selected.
- **ATP5MG::KMT2A:** a coding fusion hypothesis is retained for audit but held
  out because the Isovar result is ambiguous.
- **TPST1::CRCP T2:** 22 long reads rescue expression of the DNA-supported
  junction, but no CDS/frame can be justified, so translation and ranking are
  withheld.

The five assessed patient class-I alleles are recorded in `source/results.json`.
The clinical null allele HLA-A*01:11N is deliberately not predicted. Scores in
the source record were generated with the local NetMHCpan 4.2c installation;
the figure command performs no prediction and therefore remains reproducible
without proprietary predictor software.

Rebuild a timestamped six-page PDF plus SVG and 3600x2280 PNG pages:

```bash
python examples/osteosarc_complex_results/generate.py \
  --combined-output output/pdf/vaxrank-complex-variant-results.pdf
```

Or use the installed command with another result record:

```bash
vaxrank-complex-result-figure results.json \
  --output-root result-figures \
  --combined-output vaxrank-complex-variant-results.pdf
```

Every per-variant directory includes the exact JSON record rendered on that
page. The run manifest stores the input checksum, software version, page order,
raster dimensions, and combined-PDF path.
