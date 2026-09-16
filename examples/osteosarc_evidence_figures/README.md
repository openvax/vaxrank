# Osteosarcoma cross-platform evidence figures

These examples separate four claims that are often conflated: a DNA
rearrangement exists, short-read RNA detects its expression, local assembly
resolves sequence around the event, and a long read spans the transcript
junction. A detected junction is not automatically a validated coding fusion.

The curated cases are:

- **DSCAML1 :: ESRRG:** a recurrent DNA structural variant with no cataloged
  RNA junction. The transcript and protein claims are withheld.
- **TPST1 :: CRCP at T1:** DNA, short-read RNA, and long-read RNA independently
  support the same junction, but the source records do not establish a full
  coding transcript.
- **TPST1 :: CRCP at T2:** the strongest matched-sample long-read rescue. DNA
  supports the rearrangement and CTAT-LR-Fusion reports 22 ONT junction reads,
  while the matched T2 short-read fusion call is absent.
- **MAF :: ENSG00000261722 at T2:** a long-read-only counterexample. Its 31 ONT
  reads are visible, but absent DNA and short-read support keep the event and
  any protein product unresolved.

Rebuild a timestamped SVG/PDF/high-resolution PNG run:

```bash
python examples/osteosarc_evidence_figures/generate.py
```

Or use the installed command on the JSON evidence schema:

```bash
vaxrank-evidence-figure evidence.json --output-root figures
```

Every protein panel is deliberately withheld because the osteosarc.com source
tables establish breakpoint or splice-junction sequence, not a phased,
transcript-resolved coding sequence. The small-variant examples in
`examples/osteosarc_mutation_figures/` show cases where both transcript
nucleotide sequence and protein translation can be supported.
