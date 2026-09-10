# #423: independently documented Sid construct pilot

## Reproduce the current pilot

The independent definitions, selected original-read regressions, full-source
real-model cache and report renderer live in
`tests/data/osteosarc/construct_audit/`. They are maintainer/research assets in
the source repository, not an automatic vaccine-design command.

Render the comparison and the linked complete per-bond/MHC inventory offline:

```sh
python -m tests.data.osteosarc.construct_audit.render_report --output /tmp/sid-report
```

On macOS arm64, use the same `DYLD_FALLBACK_LIBRARY_PATH=/opt/homebrew/lib`
environment as the repository scripts. Run the regression groups through
`./test.sh tests/test_osteosarc_construct_rna.py tests/test_osteosarc_construct_predictions.py`.

Maintainer regeneration requires the installed, licensed local NetMHCpan 4.2c
model, Pepsickle 0.1.3, the recorded published dependencies, full human Ensembl
114 and Isovar's complete regional acquisition cache. The generator checks the
original BAM bytes against the released Isovar receipts before reconstruction:

```sh
python -m tests.data.osteosarc.construct_audit.generate_predictions \
  --regional-cache /path/to/isovar-regional-cache --output /tmp/new-sid-cache
```

Use a new output directory; existing real caches are not overwritten. See the
Isovar 1.8.1 expansion README/acquisition scripts for bounded indexed regional
acquisition. Whole BAM downloads are neither required nor used. The separate
`import_rna.py ISOVAR_REPO ISOVAR_REGIONAL_CACHE --output NEW_DIRECTORY` recreates
selected offline fixtures; it reads pinned Git objects and never consults the
historical vaccine strings when selecting read names.

The current report includes a real model flag: the non-target reference peptide
`DTGLKQAL` in JLF V2/V3 and the long-mRNA antigen segment has NetMHCpan 4.2c EL
percentile rank 0.47 for HLA-B*08:01. Its exact Ensembl matches are attributed to
DYNC1H1 transcripts outside the pinned oncoref CTA-candidate universe. This is
not evidence of tissue presentation, TCR cross-reactivity or toxicity.

The independently intended `KRFHATISF` ligand is overlaid separately from the
single-residue mutation mask. Internal cuts, N/C boundaries, missing context
padding and terminal sentinels remain separate. Complete historical mRNA
products and manufactured terminal chemistry are not established; bare-peptide
model inputs are explicitly conditional. Other variants/constructs, remaining
RNA products and unassessed models are enumerated in the report.

## Scope and implementation contract

This work follows #421 and the complete-context API in #422 / PR #435.
It also fixes #436: standard EnsemblRelease objects must retain a full
content-based reference identity and reuse source snapshots safely.
Review also found #438: the self-reference policy omitted available annotated
translations outside the literal protein_coding biotype. All reference paths now
retain every annotated protein sequence and preserve source transcript biotype;
older disk indexes are invalidated by policy. The pilot self flag has eight
protein_coding and 22 NMD-labelled DYNC1H1 source transcripts. Annotation is not
proof of normal-tissue expression, tolerance or antigen presentation; reference
LoF and IG/TR records also retain their distinct labels. See Ensembl's
[biotype definitions](https://mart.ensembl.org/info/genome/genebuild/biotypes.html)
and [primary NMD-targeted transcript/MHC-I experiments](https://doi.org/10.1073/pnas.1309956110).
The latter support not treating predicted NMD as an absolute absence criterion;
they do not establish processing or presentation of the Sid sequences.

CI also exposed #439, an owned Xvfb child surviving graceful PDF cleanup. Bounded
forced cleanup is tested without dropping the real Linux PDF test; the underlying
wrapper gap is tracked at cgoldberg/xvfbwrapper#76.
It does not replace historical selection validation in #414 or invent the
missing historical ranking settings. Start with DYNC1H1 p.Val314Ile as expressly
required by #423; report the pilot's limits instead of whole-cohort coverage.

## Independent inputs and scientific distinctions

- Pin the published DYNC1H1 sequences and target annotations from
  https://osteosarc.com/variant/DYNC1H1-chr14-101980529/, retrieved 2026-09-09.
  Four documented entries: KRFHATISF (mRNA minimal epitope),
  VLLTLDILKHGKRFHATISFDTDTGLKQAL (mRNA antigen segment),
  GKRFHATISFDTDTGLKQALETKK (JLF V2/V3 final peptide), and
  KHGKRFHATISFDTDTGL (CeGaT peptide). Explicitly map KRFHATISF occurrences;
  do not infer independently intended ligand identity from prediction scores.
- The page separately lists ISFDTDTGL and RFHATISF as tested ELISPOT peptides.
  They are assay comparators, not automatically documented vaccine selections.
  Assay results do not calibrate cleavage or establish a mechanistic explanation.
- JLF's terminal KK is documented; the solubility rationale is attributed to the
  user's 2026-09-09 report. The page mentions CSBio peptides in immune-monitoring
  notes but does not establish that CSBio requested these additions. Do not guess
  another residue or undocumented chemistry. Keep conditional free-terminus
  calculations distinct from manufactured chemistry, which remains unresolved.
- A published mRNA antigen segment is not a full translated product. Do not
  fabricate a CDS, signal peptide, linkers, antigen order or mature termini for
  Sid. Report complete historical mRNA-product processing as unassessed unless
  an independent complete sequence is found and validated.
- Pin clinical HLA typing from https://osteosarc.com/data/. Analyze the five
  documented non-null class-I alleles with explicit model/kind/length requests.
  Preserve A*01:11N as a null call, not an expressed presenter, and retain the
  computational typing discrepancy at https://osteosarc.com/dragen/hla/.
  Class II remains explicitly unassessed in the initial class-I pilot.

## RNA and source provenance

Adopt pinned Isovar 1.8.1 data, not newly synthesized reads. Reuse the prepared
#414 asset importer, checksums and dataset-specific reference identity without
copying a second independent fixture definition. The original #414 worktree and
its blocked cached-selection tests must be preserved. Decide the shared-fixture
commit boundary before moving assets: both PRs should consume the same immutable
data and no known-failing test is silently omitted from its owning issue.

For the DYNC1H1 pilot reconstruct the native context through Vaxrank/Isovar with
balanced selection, 85% compatible read-name support and the independent
two-read-object/base floor. Preserve source/timepoint/read provenance and the
actual source gene/transcripts/species. Additions never acquire RNA support.
Reconstructability, chosen historical window and predicted processing are
separate observations. Do not change thresholds to force a historical match.

The released matrix contains eligible DYNC1H1 PacBio T1, ONT T1/T2/T3,
short-read single-cell T1/T2/T3 and bulk T0/T2 products. Import bounded original
records from these nine products after validating the regional BAM bytes against
the released Isovar acquisition receipts. The existing selected ONT case alone
is not sufficient modality coverage. Keep full-region upstream counts separate
from selected fixture counts, and preserve all 164 upstream RNA-product outcomes
in the coverage inventory. PacBio has two default-quality alternate reads;
many additional aligned alternate bases lack qualities and cannot be silently
promoted to quality-filtered RNA support.

The real-model comparison report must use complete original regional RNA, not
selected regression subsets. A full-data check showed that subsampling changes
some balanced contexts: full ONT T1/T2 lack the complete long-mRNA native window
at default settings, although their selected fixtures contain it. Keep both
observations with explicit input scope. Anchor the four documented final
constructs to the earliest bulk T0 context, which contains all four without
changing thresholds; retain every other source's exact window limitations.

Enumerate absent/unassessed PacBio, ONT, bulk and short-read single-cell cases
against the released upstream manifest. Vendor reprocessing is not a biological
replicate; read-name/fragment/cell/UMI counts are not interchangeable. Selected
fixture subsets are not unbiased VAF estimators. Broader coverage remains a
named limitation; new upstream data gaps should be filed on Isovar if untracked.

## Real inference and self assessment

Generate reproducible real NetMHCpan class-I and Pepsickle caches on complete
native/final contexts. Pin input and output hashes, exact package/backend/model
versions, flags and actual model-asset hashes. Never fabricate prediction scores
or derive expected historical choices from the implementation being tested.
CI consumes the small offline caches and forbids model downloads/network.

Use released canonical enzyme models only inside their stated domain. Label
qualitative CPN terminal recognition as such, not serum kinetics or total
peptide stability. Native/embedded segment termini are not exposed substrates;
unsupported chemistry, other proteases and mature-product localization remain
unassessed. Any iterative trimming scenario needs its own explicit assumptions.

Overlay exact internal/boundary observations on independently documented target
ligands and all returned predicted ligands. Inventory exact self matches against
a full, versioned human reference, preserving every non-CTA gene/transcript and
species. A small Isovar reconstruction reference cannot establish proteome-wide
self absence. Reuse existing attributed self-policy/Topiary DSL primitives, not
ad-hoc ranking filters. Any strong-ligand threshold must be named by kind, native
metric and model, with its primary justification; it is not a clinical verdict.

## Deliverables and gates

An independently sourced construct manifest, reproducible maintainer generator,
small hash-pinned real caches, explicit unassessed-coverage inventory, readable
comparison report, and offline regressions for exact sequences, modifications,
target/ligand cut positions, source-attributed self results and native round trips.
List non-pilot constructs and missing models instead of implying all are covered.
Require lint, full tests/coverage, report smoke, review and green CI before the
next version-bumped PR merge and immediate clean-main/PyPI-verified deployment.
