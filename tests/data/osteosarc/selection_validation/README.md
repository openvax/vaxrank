# Sid final-selection pilot

This is a software-validation comparison, not a vaccine recommendation or a
reproduction of an undocumented historical ranking protocol.

## Inputs

- Five original-read cases from Isovar **1.8.1**, immutable commit
  `9297bb7bf19def6c33dd8295030de2b040d7dee3`. `import_assets.py` reads Git
  objects, checks upstream hashes, and preserves the original case records,
  source receipts/headers and GRCh38 Ensembl-87 reference manifest/models.
- DYNC1H1 V314I and EXOC4 S34I use the source-labelled ONT T1 product;
  H1-2 and MAP2 use source-labelled ONT T2 products; the second DYNC1H1 locus
  uses the 2022-12-16 bulk-RNA product. These are selected-read subsets, not
  unbiased VAF samples. Vendor reprocessing is not biological replication.
- The seven independently documented sequences are transcribed in
  `documented.json`, with explicit modality, provider/version and source URLs.
  Both mRNA entries remain mRNA regardless of their provider/display label
  (draft correction #441); neither establishes a complete historical translated
  product. Peptide boundaries
  and JLF's added KK are explicit. Target intervals are zero-based,
  half-open RNA-affected amino-acid coordinates; an in-frame deletion crossing
  codons can retain an affected codon rather than a zero-width protein interval.
- The clinical HLA table is preserved verbatim in structured form. The null
  A*01:11N allele is not replaced or scored as an expressed allele. Computational
  A*01:01 homozygosity is recorded as a discrepant source claim, not adjudicated.
- `netmhcpan42.tsv` contains **6,140 actual model rows for 614 peptides**, from
  NetMHCpan **4.2c**, `-BA`, five class-I alleles, lengths 8–11. Both binding
  and elution-score leaves are retained. Every requested RNA-context k-mer and
  position-aligned WT comparator is included before filtering. The manifest
  records exact inputs, settings, package versions and checksums. CI performs
  offline cache lookups only, with network calls forbidden and no live fallback.

Maintainer-only regeneration uses installed local models and writes to a new
directory, preserving the previous evidence for comparison:

```sh
python -m tests.data.osteosarc.selection_validation.generate_predictions --output /tmp/new-sid-selection-cache
```

The independent source records are never generated from prediction outputs.

Regenerated on 2026-09-16 with Isovar 1.17.0 and Vaxrank 3.18.1 after the
read-identity correction and coordinated dependency update. All 11 RNA
contexts, 614 requested peptides and 6,140 prediction values are unchanged
from the preceding cache after matching peptide/allele/kind/occurrence keys.
The regenerated cache has a different row order, not changed model scores.
The H1-2 ONT T2 fixture here is distinct from the ambiguous bulk T0 fixture
held out by the RNA-context tests. The manifest records
the actual new runtime versions, generation time and output checksum.

## Comparison configuration

Use each documented **native** peptide length, Vaxrank's adaptive RNA context,
balanced selection, 0.85 compatible read-name support, the independent two-read
per-base floor, assembly enabled, and no DNA fallback. All remaining epitope,
window-scoring and tie-break settings are recorded defaults, not fitted to the
historical choices. The native RNA-derived sequence is compared before any
documented manufacturing addition is represented.

## Initial live-model diagnostic, 2026-09-09

The native window is reconstructable in **all seven** cases. Final selection
matches **four of seven native selections** under the stated configuration:

| Documented selection | Exact native selection? | Explanation |
| --- | --- | --- |
| DYNC1H1 mRNA minimal | Yes | `KRFHATISF` |
| DYNC1H1 CeGaT | Yes | `KHGKRFHATISFDTDTGL` |
| EXOC4 CeGaT | Yes | `SVIRTLSTIDDVEDREN` |
| EXOC4 JLF V2/V3 | Yes | `SVIRTLSTIDDVEDRENEKGR`; leading KK remains a manufacturing addition |
| DYNC1H1 long mRNA | No | Current selection is `LKHGKRFHATISFDTDTGLKQALETVNDYN`, the same 30-mer length placed 7 residues toward the C-terminus of the documented window; both windows contain V314I and share 23 residues |
| DYNC1H1 JLF V2/V3 | No | Current selection is `LKHGKRFHATISFDTDTGLKQA`, the same 22-mer length placed 3 residues toward the N-terminus of the documented window; both windows contain V314I and share 19 residues; the documented entry adds a terminal KK |
| H1-2 CeGaT | No selection | Fixture fails `min_ratio_alt_to_other_fragments` before MHC prediction; not evidence of absent binding |

The two differing DYNC1H1 windows have the same requested lengths as their
documented native comparators. Both disagreements are window placement over one
reconstructed protein at the single DYNC1H1 V314I locus, not a disagreement
about the reconstructed sequence: the documented string is present in the
reconstructed protein in all seven cases, and each competing window contains the
mutation. Because no historical ranking settings are published, neither
placement can be called the correct one here. No sequence or score was changed
to force agreement. H1-2 retains its deletion sequence and the actual filter failure.
MAP2 is a separate no-DNA-fallback regression, outside these seven documented
comparisons. The imported Isovar 1.8.1 manifest deliberately retains the
published GRCh38 `chr2:209694768 CCTGGGCTACTGTGTGTTCAATA>C` allele (22-bp
deletion), preserving the original fixture and its checksums.

Osteosarc 0.1.0's `allele-MAP2-chr2-209694768` correction identifies the complex
replacement `CCTGGGCTACTGTGTGTTCAATAAGTACACAGT>CAGGG` at the same GRCh38
position (net -28 bp). The primary [Tempus Pindel VCF](https://sid-sijbrandij-osteosarc-dataset.s3.us-west-2.amazonaws.com/vendor/tempus/TL-24-ALMY2X4KMV/DNA/TL-24-ALMY2X4KMV.soma.pindel.vcf)
records the equivalent GRCh37 `2:210559493 CTGGGCTACTGTGTGTTCAATAAGTACACAGT>AGGG`,
with `c.2599_2630delinsAGGG` and its anchored representation in `OLD_VARIANT`.
The [CeGaT somatic table](https://sid-sijbrandij-osteosarc-dataset.s3.us-west-2.amazonaws.com/vendor/cegat/P116686_2_S000048/P116686_2_somatic.tsv)
lists consistent component substitutions at GRCh37 chr2:210559493/210559494
and a deletion anchored at chr2:210559496; these are not independent targets.

An offline comparison on 2026-09-19 with Vaxrank 3.19.4, Isovar 1.21.0 and
Varcode 7.0.0 applied both alleles separately to the same pinned BAM/reference
with unchanged defaults. **Both** yielded one alternate read/fragment, zero
reference/other reads, no protein meeting the two-read floor, and no MHC
prediction or vaccine selection. The regression now checks both alleles.
This selected-read subset does not establish absence of mutant expression in
the full RNA source, nor does it settle the corrected protein or vaccine
peptide sequence. The four-of-seven comparison above is unchanged because it
contains no MAP2 record. See [Vaxrank #487](https://github.com/openvax/vaxrank/issues/487).

**Upstream blocker, resolved 2026-09-11:** Topiary #296 duplicated cache
occurrence coordinates, and six end-to-end cache tests failed with that
diagnosis. It is fixed upstream by Topiary #299, commit `2c9f867`, which rebinds
a cached protein scan to the occurrence it was asked about, and released as
Topiary 5.55.1 on PyPI. `requirements.txt` floors topiary at `>=5.55.1`
accordingly. All 18 pilot tests pass against the real published release,
installed non-editable and verified directly against PyPI's index rather than
the topiary source tree, and the offline path composes end to end.

A second upstream gap surfaced while chasing this one: Topiary #300,
`serum_half_life`/`blood_half_life` declared as kinds but unreachable from the
ranking DSL, also fixed in 5.55.1. Reaching that fix required raising
vaxrank's own `requirements.txt` mhctools floor, which now reads `>=3.44.0`
and stands on two separate reasons. From 3.39.0 mhctools' `Kind` class carries
both half-life kinds, without which `KIND_ALIASES` never sees them whatever
topiary version is installed. From 3.44.0 `CleavageModel` requires
`scored_endpoint` for quantitative evidence and a strictness grade for motif
rules, which vaxrank now supplies. `tests/test_declared_dependencies.py`
holds both reasons so relaxing the pin for one cannot silently drop the other.

## Wild-type comparators

A comparator is emitted only where the mutant window is genuinely index
aligned with the reference protein, verified residue by residue outside the
mutated interval rather than assumed. An insertion or deletion shifts every
downstream reference position, so windows straddling or following the indel
have no length-matched reference counterpart; those record no comparator
instead of a shifted slice. The earlier generator sliced the reference at the
mutant's own offsets unconditionally, which cached `AAKPKVVKP` — reference
residues spanning H1-2's deleted `AAKPK` — as a wild-type peptide for a
window that covers none of them. `unaligned_wt_comparators` in the manifest
counts the skipped windows per context so the omission is visible rather than
silent.

## Limits and follow-up

Historical model versions, filters, candidate pools, manual decisions and order
are unknown. Four matching strings do not establish protocol reproduction.
Class-II predictions are unassessed, including for documented CeGaT selections.
The compact reference is not a full-human self-proteome safety audit. PacBio,
single-cell/UMI attribution, additional variants, assemblies and timepoints are
not silently represented by these five cases; they remain in the upstream
Isovar cohort inventory. Full manufactured-context cleavage and non-CTA-self
audits are Vaxrank #422/#423, not implied by this comparison.
