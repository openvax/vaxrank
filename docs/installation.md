# Installation

```sh
pip install vaxrank
```

Vaxrank requires Python 3.11 or newer.

## Common CTA evidence

The ordinary installation supports OncoRef 1.23.1 through 1.24.0 and requires
ORFdata 0.42.6 or newer,
which support the [common RNA restriction policy](antigen-inputs.md#common-rna-restriction-policy-experimental).
The older `vaxrank[cta]` spelling remains a compatible alias:

```sh
pip install 'vaxrank[cta]'
pip check
```

Common CTA evidence and patient-background comparisons use one environment with
Topiary 5.108.4 or newer. PirlyGenes also constrains the supported OncoRef range,
so resolve
related upgrades together:

```sh
pip install --upgrade 'vaxrank[cta]' topiary pirlygenes oncoref
pip check
```

The resolver may retain an older supported OncoRef release until a compatible
PirlyGenes release is available. Do not bypass these constraints with `--no-deps`.

The common policy remains opt-in: supply `restriction_policy`, its RNA evidence
and the complete host reference. Installing a newer reference package does not
admit a target, change historical CTA defaults or reinterpret saved decisions.
OncoRef 1.24's [TCR-T evidence](https://pirl-unc.github.io/oncoref/cta-release-1.24.0/)
retains its peptide, HLA and product limits; it is not an automatic admission
override for vaccine windows.

## Reference annotation

Vaxrank uses [PyEnsembl](https://github.com/openvax/pyensembl) for reference
genome annotation. Install an Ensembl release matching the reference genome your
variants were called against:

```sh
# GRCh38
pyensembl install --release 113 --species human
# GRCh37 (legacy)
pyensembl install --release 75 --species human
```

## PDF reports

PDF generation uses [wkhtmltopdf](https://wkhtmltopdf.org/downloads.html) by
default. Install a compatible executable and make it available on `PATH`.

Alternatively, pass `--pdf-backend=weasyprint` to use
[WeasyPrint](https://doc.courtbouillon.org/weasyprint/stable/first_steps.html)
(experimental). It uses system libraries such as Pango:

```sh
pip install weasyprint
# macOS also needs: brew install pango
```

On Apple Silicon, WeasyPrint loads Pango through dyld, which does not search
Homebrew's `/opt/homebrew/lib` by default. Add this to your shell profile:

```sh
export DYLD_FALLBACK_LIBRARY_PATH="/opt/homebrew/lib:$DYLD_FALLBACK_LIBRARY_PATH"
```

Intel macOS already includes Homebrew's `/usr/local/lib` in dyld's default
fallback path.

## Isovar version window

RNA reconstruction requires a compatible Isovar release; see
[Isovar support compatibility](isovar-support.md) for the tested range and the
measurements Vaxrank retains from it.
