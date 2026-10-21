# WLOS_0001 canine DNA/RNA replay

This offline example carries one published **somatic tumor/normal DNA call**
and **original tumor RNA alignments** from dog WLOS_0001. It runs the VCF and
BAM through Isovar and Vaxrank, then writes `report.txt` and `receipt.json`.
The fixture is one TP53 locus, not a whole-sample benchmark.

From a Vaxrank source checkout with its dependencies installed:

```sh
python -m examples.wlos0001_tp53.run --output-dir /tmp/wlos0001-result
```

The expected receipt has one coding, RNA-supported variant, 107 reference and
297 alternate RNA reads (102 and 289 fragments), and one ranked 25-mer. The
report declares CanFam3.1, Ensembl 98, and each tumor sample's identity.
`DLA-88*001:01` is an example allele passed to a **random test predictor**;
the dog's MHC genotype is unknown, and the reported affinity numbers have no
biological meaning. The local self-reference contains only three TP53 protein
sequences. This run does not validate whole-dog MHC binding or self similarity.

## Source and scope

[Gardner et al. (2019)](https://doi.org/10.1038/s42003-019-0487-2),
Supplementary Data 1, maps WLOS_0001 to normal WGS SRR8691711, tumor WGS
SRR8691710, and tumor RNA [SRR8691780](https://www.ncbi.nlm.nih.gov/sra/SRR8691780).
Supplementary Data 6 reports the matched-normal TP53 `5:32563056 G>A` somatic
call with tumor DNA VAF 0.462 on CanFam3.1. The DNA alignments are **not** in
this small fixture; the VCF transcribes the published call. The BAM retains
1,087 original aligned tumor RNA records extracted from the SRA run around
that locus. Isovar's default read filters determine the counts above.

The [upstream Isovar fixture](https://github.com/openvax/isovar/tree/main/tests/data/wlos0001_tp53)
documents the read extraction, masked chromosome construction, Ensembl release
98 source files, and checksums of the full source downloads. The masked FASTA
retains the original chromosome coordinates and exact TP53 bases but masks
everything else to `N`. It must not be used for analysis outside TP53.
The paper's `p.Arg263Cys` residue number depends on its annotation; Ensembl
98 transcript numbering differs. The stable keys are the genomic substitution
and RNA alternate allele.

These files were copied byte-for-byte from the Apache-2.0-licensed Isovar
fixture. Verify the local copies with `shasum -a 256 examples/wlos0001_tp53/data/*`:

| File | SHA-256 |
|---|---|
| `canfam3-chr5-tp53-masked.fa.gz` | `043b584f03a0828915b055e85ac0e89f25c8369e193566a7917a66b0f5ead9ef` |
| `tp53-ensembl-98.gtf` | `eb647153620db5780f1c1917df6ea55c6998535815c46dea9528566c65d61011` |
| `tp53-proteins-ensembl-98.fa` | `204f9ff2b9e84cfa24cdb84dfa316522924c628310be07243b4b07b91a759782` |
| `tp53-transcripts-ensembl-98.fa` | `61322648f241629e62f311bbb74f9682c693331f6d56f46fbda904aaff4791dd` |
| `wlos0001-tp53.bam` | `f0943702c838bc8dcc2db67410f0dd34d31ab65654b63a87c1b6f7e303ac3e00` |
| `wlos0001-tp53.bam.bai` | `a5809ec5ef9481c80fa934d3e823beba85cfbf6a91c82500dbb599e08f74c971` |
| `wlos0001-tp53.vcf` | `6ba88eaa472e9e238a206e6e31cbe13b99058ee50d1074187ed586cc92e8d3b6` |

The random predictor and TP53-only self-reference make this a workflow test.
Real canine interpretation needs the dog's DLA genotype, a validated canine
binding model, and an assembly-matched complete self proteome; those gaps are
tracked in [#596](https://github.com/openvax/vaxrank/issues/596) and
[#625](https://github.com/openvax/vaxrank/issues/625).
