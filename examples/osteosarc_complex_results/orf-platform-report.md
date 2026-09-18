# Osteosarc RNA platform and ORF audit

A protein result below is a validated local RNA-derived coding window, not a reconstructed full-length ORF. Products are not biological replicates, and counts are not platform-sensitivity estimates.

## Full 44-variant audit (assembly enabled)

| Platform | Products | Completed rows | ALT loci | Protein loci | Exact in ≥1 product | Non-exact in ≥1 product |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| ILMN | 149 | 5894 | 41 | 39 | 37 | 6 |
| ONT | 6 | 259 | 36 | 30 | 29 | 1 |
| PacBio | 9 | 44 | 3 | 3 | 3 | 0 |

Varcode produces an exact reference-plus-isolated-edit protein for all 44 loci. That sequence is an annotation result, not an observed sample-specific RNA haplotype.

## Paired original-read corpus

The same bounded BAM subsets were rerun with Isovar assembly on and off. This 49-case corpus is enriched for informative loci and is not an unbiased cohort sample.

| Platform | Assembly | Cases | Protein windows | Exact | Non-exact | No window |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| ILMN | assembly_on | 32 | 27 | 26 | 1 | 5 |
| ILMN | assembly_off | 32 | 27 | 26 | 1 | 5 |
| ONT | assembly_on | 16 | 14 | 13 | 1 | 2 |
| ONT | assembly_off | 16 | 14 | 13 | 1 | 2 |
| PacBio | assembly_on | 1 | 0 | 0 | 0 | 1 |
| PacBio | assembly_off | 1 | 0 | 0 | 0 | 1 |

Assembly changed the returned local protein length in five of 32 ILMN cases, but did not create or remove a protein window in this selected corpus. No paired ONT or PacBio case changed sequence; the PacBio denominator is only one case.

## Non-exact RNA protein windows

| Locus | Platform | Product rows | Interpretation |
| --- | --- | ---: | --- |
| NME1 | ILMN | 2 | The non-exact rows end one residue earlier than the isolated-edit window; other products recover the exact local sequence. |
| NR2F2 | ILMN | 1 | The RNA window contains an unexplained three-nucleotide transcript deletion; no matched DNA variant has yet been linked to it. |
| NTF3 | ILMN + ONT | 9 | RNA establishes the adjacent base change and compound AG>GT allele, changing Varcode's isolated Lys-to-Arg call to serine; germline versus somatic origin remains unresolved. |
| ROBO2 | ILMN | 1 | One ILMN product contains an additional amino-acid substitution; other products recover the isolated-edit sequence. |
| ZNF674 | ILMN | 1 | One ILMN product contains an additional amino-acid substitution; other products recover the isolated-edit sequence. |
| ZNF764 | ILMN | 2 | Two ILMN products contain a divergent downstream protein context; other products recover the isolated-edit sequence. |

## Attribution

The current public Isovar path identifies the focal somatic edit and retains additional transcript edits as unexplained. It can phase multiple supplied somatic variants by shared fragment names, but this 44-variant set contains no nearby nominated pair. A matched germline variant input is not wired into `run_isovar`; consequently no observed mismatch can yet be promoted to known germline by this audit.

See the adjacent JSON for every source row, sequence, transcript edit, checksum, and denominator.
