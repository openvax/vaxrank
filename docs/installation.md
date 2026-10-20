# Installation

```sh
pip install vaxrank
```

Vaxrank requires Python 3.11 or newer.

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
