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
- `netmhcpan42.tsv` contains **6,560 actual model rows for 656 peptides**, from
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
Expanded MAP2 has one alternate read object and no protein meeting the two-read
floor; it is never converted into a DNA-backed vaccine selection silently.

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
ranking DSL, also fixed in 5.55.1. Reaching that fix additionally required
raising vaxrank's own `requirements.txt` mhctools floor from `>=3.13.3` to
`>=3.39.0`, matching topiary's real requirement; below that floor mhctools'
`Kind` class predates both kinds, so `KIND_ALIASES` never sees them regardless
of the topiary version installed.

## Limits and follow-up

Historical model versions, filters, candidate pools, manual decisions and order
are unknown. Four matching strings do not establish protocol reproduction.
Class-II predictions are unassessed, including for documented CeGaT selections.
The compact reference is not a full-human self-proteome safety audit. PacBio,
single-cell/UMI attribution, additional variants, assemblies and timepoints are
not silently represented by these five cases; they remain in the upstream
Isovar cohort inventory. Full manufactured-context cleavage and non-CTA-self
audits are Vaxrank #422/#423, not implied by this comparison.
