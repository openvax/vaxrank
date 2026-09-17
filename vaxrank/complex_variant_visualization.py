"""Publication-ready Vaxrank results for assembly-dependent variants."""

from datetime import datetime, timezone
from hashlib import sha256
from html import escape
import json
from pathlib import Path
import re
import shutil
import tempfile

from .mutation_visualization import _configure_native_library_path, _slug, _write_png
from .version import __version__


RUN_TIMESTAMP_PATTERN = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{6}Z$")
OUTCOME_STYLE = {
    "selected": ("SELECTED", "#176b4d", "#dff4e9"),
    "no_target_binder": ("NO TARGET BINDER", "#755400", "#fff2c2"),
    "held_out": ("HELD OUT", "#8a4b16", "#fde9d5"),
    "no_translation": ("NO TRANSLATION", "#7a2938", "#f8e3e7"),
}


def _lines(text, width=72, limit=3):
    words = str(text).split()
    lines = []
    current = []
    for word in words:
        candidate = " ".join([*current, word])
        if current and len(candidate) > width:
            lines.append(" ".join(current))
            current = [word]
        else:
            current.append(word)
    if current:
        lines.append(" ".join(current))
    if len(lines) > limit:
        lines = lines[:limit]
        lines[-1] = lines[-1].rstrip(".") + "..."
    return lines


def _text_lines(text, x, y, css_class="body", width=72, limit=3, step=18):
    return "".join(
        '<text x="%s" y="%s" class="%s">%s</text>' %
        (x, y + index * step, css_class, escape(line))
        for index, line in enumerate(_lines(text, width=width, limit=limit)))


def _sequence_lines(sequence, x, y, width=58, step=19):
    if not sequence:
        return '<text x="%s" y="%s" class="muted">Not established</text>' % (x, y)
    return "".join(
        '<text x="%s" y="%s" class="sequence">%s</text>' %
        (x, y + index * step, escape(sequence[start:start + width]))
        for index, start in enumerate(range(0, len(sequence), width)))


def _base_style():
    return """
text { font-family: Inter, Helvetica, Arial, sans-serif; }
.title { font-size: 27px; font-weight: 700; fill: #17212b; }
.subtitle { font-size: 13px; fill: #445364; }
.kicker, .panel-title { font-size: 11px; font-weight: 700; letter-spacing: 1px; fill: #526273; }
.outcome { font-size: 15px; font-weight: 700; letter-spacing: 1px; }
.body { font-size: 12px; fill: #334252; }
.strong { font-size: 13px; font-weight: 700; fill: #17212b; }
.muted { font-size: 12px; fill: #687786; }
.sequence { font-family: "SFMono-Regular", Consolas, monospace; font-size: 13px; font-weight: 600; fill: #17212b; }
.table-head { font-size: 10px; font-weight: 700; letter-spacing: 0.7px; fill: #526273; }
.table { font-size: 11px; fill: #263544; }
.table-mono { font-family: "SFMono-Regular", Consolas, monospace; font-size: 11px; font-weight: 600; fill: #17212b; }
.footer { font-size: 9px; fill: #687786; }
"""


def render_complex_result_svg(record):
    """Render one complex-variant decision and ranking result as SVG."""
    outcome = record["outcome"]
    if outcome not in OUTCOME_STYLE:
        raise ValueError("Unknown complex-variant outcome %r" % outcome)
    outcome_label, outcome_color, outcome_background = OUTCOME_STYLE[outcome]
    rna = record["rna_evidence"]
    counts = rna.get("counts_text") or (
        "ALT %s | REF %s | OTHER %s | top protein %s" % (
            rna.get("alt_fragments", "n/a"),
            rna.get("ref_fragments", "n/a"),
            rna.get("other_fragments", "n/a"),
            rna.get("top_protein_fragments", "n/a")))
    filters = rna.get("gate_status") or (
        "PASS" if rna.get("passes_vaxrank_filters") else "DOES NOT PASS")
    filter_color = (
        "#176b4d" if rna.get("passes_vaxrank_filters") else "#7a2938")

    epitope_rows = []
    # Keep the publication panel legible while preserving every ranked row in
    # the adjacent record.json.
    for index, epitope in enumerate(record.get("top_epitopes", ())[:5]):
        y = 531 + index * 22
        fill = "#f5f8fa" if index % 2 else "#ffffff"
        epitope_rows.append(
            f'<rect x="628" y="{y - 15}" width="540" height="22" fill="{fill}"/>'
            f'<text x="642" y="{y}" class="table-mono">{escape(epitope["sequence"])}</text>'
            f'<text x="786" y="{y}" class="table">{escape(epitope["allele"])}</text>'
            f'<text x="930" y="{y}" class="table">{epitope["ic50_nm"]:.2f}</text>'
            f'<text x="1022" y="{y}" class="table">{epitope["percentile_rank"]:.3f}</text>'
            f'<text x="1115" y="{y}" text-anchor="end" class="table">{epitope["epitope_score"]:.3f}</text>')
    if not epitope_rows:
        epitope_rows.append(
            '<rect x="628" y="516" width="540" height="58" rx="6" fill="#f5f8fa"/>'
            '<text x="648" y="548" class="muted">%s</text>' % escape(
                record.get(
                    "epitope_message",
                    "No target epitope passed the configured filters.")))

    peptide = record.get("selected_long_peptide")
    peptide_panel = (
        '<text x="52" y="526" class="sequence">%s</text>' % escape(peptide)
        if peptide else
        '<text x="52" y="526" class="muted">No peptide selected for a vaccine construct.</text>')
    scores = ""
    if peptide:
        score_basis = record.get(
            "score_basis", "RNA support × target epitope score")
        scores = (
            '<text x="52" y="551" class="body">Rank score %.3f | target epitope score %.3f</text>'
            '<text x="52" y="570" class="muted">%s</text>' %
            (record["combined_score"], record["target_epitope_score"],
             escape(score_basis)))

    reason = record.get("decision_reason", "")
    sources = " | ".join(record.get("source_labels", ()))
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="760" viewBox="0 0 1200 760">
<rect width="1200" height="760" fill="#ffffff"/>
<style>{_base_style()}</style>
<text x="32" y="45" class="kicker">COMPLEX-VARIANT VAXRANK RESULT</text>
<text x="32" y="78" class="title">{escape(record["title"])}</text>
<text x="32" y="103" class="subtitle">{escape(record["variant"])} | {escape(record["sample"])} | {escape(record["variant_class"])}</text>
<rect x="32" y="122" width="1136" height="55" rx="8" fill="{outcome_background}"/>
<text x="52" y="146" class="outcome" fill="{outcome_color}">{outcome_label}</text>
{_text_lines(reason, 238, 145, width=105, limit=2, step=17)}

<text x="32" y="208" class="panel-title">1  DNA-ONLY INTERPRETATION</text>
<rect x="32" y="222" width="552" height="101" rx="8" fill="#eef3f7"/>
{_text_lines(record["dna_only"], 52, 248, width=76, limit=4)}

<text x="616" y="208" class="panel-title">2  SAMPLE-SPECIFIC RNA EVIDENCE</text>
<rect x="616" y="222" width="552" height="101" rx="8" fill="#eef3f7"/>
<text x="636" y="248" class="strong">{escape(rna["platform"])} | {escape(rna["source"])}</text>
<text x="636" y="272" class="body">{escape(counts)}</text>
<text x="636" y="298" class="body">Vaxrank RNA gates: </text>
<text x="771" y="298" class="strong" fill="{filter_color}">{filters}</text>

<text x="32" y="353" class="panel-title">3  ASSEMBLED OR TRANSLATED PROTEIN</text>
<rect x="32" y="367" width="1136" height="73" rx="8" fill="#f8fafb"/>
{_sequence_lines(record.get("protein_sequence"), 52, 396)}

<rect x="32" y="479" width="552" height="186" rx="8" fill="#f5f8fa"/>
<text x="52" y="505" class="panel-title">4  VACCINE-CONSTRUCT DECISION</text>
{peptide_panel}{scores}
{_text_lines(record["interpretation"], 52, 585, width=72, limit=3, step=18)}

<text x="616" y="493" class="panel-title">5  TARGET EPITOPES PASSING DEFAULT FILTERS</text>
<text x="642" y="511" class="table-head">PEPTIDE</text>
<text x="786" y="511" class="table-head">ALLELE</text>
<text x="930" y="511" class="table-head">IC50 NM</text>
<text x="1022" y="511" class="table-head">RANK %</text>
<text x="1115" y="511" text-anchor="end" class="table-head">SCORE</text>
{''.join(epitope_rows)}

<text x="52" y="646" class="muted">{escape(record.get("limitation", ""))}</text>

<line x1="32" y1="690" x2="1168" y2="690" stroke="#d8e0e7"/>
<text x="32" y="714" class="footer">{escape(sources)}</text>
<text x="1168" y="737" text-anchor="end" class="footer">Vaxrank {escape(__version__)} | research use only</text>
</svg>'''


def render_summary_svg(records, metadata, page_index=1, page_count=1):
    """Render the decision-matrix cover page."""
    if len(records) > 7:
        raise ValueError("Summary pages support at most seven result rows")
    rows = []
    for index, record in enumerate(records):
        label, color, background = OUTCOME_STYLE[record["outcome"]]
        y = 241 + index * 52
        rows.append(
            f'<rect x="32" y="{y - 27}" width="1136" height="46" rx="7" fill="{background}"/>'
            f'<text x="52" y="{y - 6}" class="strong">{escape(record["title"])}</text>'
            f'<text x="52" y="{y + 11}" class="footer">{escape(record["variant_class"])}</text>'
            f'<text x="420" y="{y + 3}" class="body">{escape(record["summary_rna"])} </text>'
            f'<text x="1138" y="{y + 3}" text-anchor="end" class="outcome" fill="{color}">{label}</text>')
    alleles = ", ".join(metadata["hla_alleles"])
    page_label = (
        "" if page_count == 1 else " | decision matrix %d of %d" % (
            page_index, page_count))
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="760" viewBox="0 0 1200 760">
<rect width="1200" height="760" fill="#ffffff"/>
<style>{_base_style()}</style>
<text x="32" y="47" class="kicker">OSTEOSARCOMA COMPLEX-VARIANT RESULTS</text>
<text x="32" y="84" class="title">Assembly changes what Vaxrank can claim</text>
<text x="32" y="112" class="subtitle">Actual RNA reconstruction, translation, epitope prediction, and vaccine-selection outcomes{page_label}</text>
<rect x="32" y="139" width="1136" height="48" rx="8" fill="#dceefa"/>
<text x="52" y="168" class="strong">A DNA event is not an antigen. Vaxrank requires a defensible mutant protein before ranking epitopes.</text>
<text x="52" y="207" class="panel-title">VARIANT</text>
<text x="420" y="207" class="panel-title">RNA / TRANSLATION EVIDENCE</text>
<text x="1115" y="207" text-anchor="end" class="panel-title">VAXRANK OUTCOME</text>
{''.join(rows)}
<rect x="32" y="595" width="1136" height="71" rx="8" fill="#eef3f7"/>
<text x="52" y="620" class="panel-title">PREDICTION CONTEXT</text>
{_text_lines('NetMHCpan %s; patient class-I alleles: %s' % (metadata['predictor_version'], alleles), 52, 646, width=150, limit=1)}
<text x="52" y="662" class="footer">Exploratory binding ranks: full-proteome exact-self screening was not run for explicit assembled-antigen inputs.</text>
<line x1="32" y1="690" x2="1168" y2="690" stroke="#d8e0e7"/>
<text x="32" y="714" class="footer">Generated {escape(metadata["analysis_date"])} from committed evidence records</text>
<text x="1168" y="737" text-anchor="end" class="footer">Vaxrank {escape(__version__)} | research use only</text>
</svg>'''


def render_provenance_svg(records, metadata):
    """Render methods, provenance, and interpretation limits."""
    sources = []
    for record in records:
        sources.extend(record.get("source_urls", ()))
    sources = list(dict.fromkeys(sources))
    source_lines = "".join(
        f'<text x="52" y="{532 + index * 19}" class="footer">{escape(source)}</text>'
        for index, source in enumerate(sources[:7]))
    versions = " | ".join(
        "%s %s" % item for item in metadata["software_versions"].items())
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="760" viewBox="0 0 1200 760">
<rect width="1200" height="760" fill="#ffffff"/>
<style>{_base_style()}</style>
<text x="32" y="47" class="kicker">METHODS AND PROVENANCE</text>
<text x="32" y="84" class="title">What is measured, inferred, and withheld</text>

<rect x="32" y="121" width="355" height="330" rx="8" fill="#eef3f7"/>
<text x="52" y="151" class="panel-title">1  RECONSTRUCT</text>
{_text_lines('Isovar groups RNA fragments by allele, assembles supported sequence, and translates only when the coding frame is defensible.', 52, 183, width=43, limit=6, step=20)}
<text x="52" y="322" class="panel-title">EVIDENCE RETAINED</text>
{_text_lines('ALT, REF, and OTHER fragment counts; supporting-fragment count; platform; sample; transcript or junction provenance.', 52, 351, width=43, limit=4, step=20)}

<rect x="423" y="121" width="355" height="330" rx="8" fill="#e8f3fa"/>
<text x="443" y="151" class="panel-title">2  PREDICT AND RANK</text>
{_text_lines('Vaxrank enumerates mutant-protein windows and predicts class-I binding through mhctools and Topiary. Explicit assembled-antigen panels use target epitope score and do not claim a full-proteome exact-self screen.', 443, 183, width=43, limit=8, step=20)}
<text x="443" y="362" class="panel-title">ALLELE POLICY</text>
{_text_lines('The clinical null allele HLA-A*01:11N is not assessed. The five expressed class-I alleles are evaluated.', 443, 391, width=43, limit=3, step=20)}

<rect x="814" y="121" width="354" height="330" rx="8" fill="#f8e9e5"/>
<text x="834" y="151" class="panel-title">3  WITHHOLD WHEN NEEDED</text>
{_text_lines('Junction support alone does not establish a CDS, reading frame, translated product, or junction-spanning peptide. Ambiguous and unresolved fusion hypotheses remain outside vaccine selection.', 834, 183, width=43, limit=8, step=20)}
<text x="834" y="362" class="panel-title">NO EMPTY-CATEGORY FUDGING</text>
{_text_lines('Assembly success may still produce no target binder. That is a negative Vaxrank result, not a missing result.', 834, 391, width=43, limit=3, step=20)}

<text x="32" y="486" class="panel-title">SOFTWARE AND MODEL</text>
<text x="52" y="507" class="body">{escape(versions)}</text>
{source_lines}
<line x1="32" y1="690" x2="1168" y2="690" stroke="#d8e0e7"/>
<text x="32" y="714" class="footer">Exact evidence and scores are stored beside every figure in record.json.</text>
<text x="1168" y="737" text-anchor="end" class="footer">Vaxrank {escape(__version__)} | research use only</text>
</svg>'''


def _write_multipage_pdf(svgs, path):
    _configure_native_library_path()
    from weasyprint import HTML

    pages = "".join('<div class="page">%s</div>' % svg for svg in svgs)
    html = (
        '<html><head><style>@page { size: 12.5in 7.9167in; margin: 0; } '
        'body { margin: 0; background: white; } '
        '.page { width: 12.5in; height: 7.9167in; break-after: page; } '
        '.page:last-child { break-after: auto; } '
        'svg { width: 100%; height: 100%; display: block; }'
        f'</style></head><body>{pages}</body></html>')
    HTML(string=html).write_pdf(path)


def generate_complex_variant_results(
        input_json, output_root, timestamp=None, png_scale=3,
        combined_output=None):
    """Create timestamped result pages and a combined PDF."""
    input_json = Path(input_json)
    output_root = Path(output_root)
    timestamp = timestamp or datetime.now(timezone.utc).strftime("%Y-%m-%dT%H%M%SZ")
    if not RUN_TIMESTAMP_PATTERN.fullmatch(timestamp):
        raise ValueError("Timestamp must use UTC form YYYY-MM-DDTHHMMSSZ")
    if not 0 < png_scale <= 10:
        raise ValueError("PNG scale must be greater than 0 and at most 10")
    run_directory = output_root / timestamp
    if run_directory.exists():
        raise FileExistsError("Result run already exists: %s" % run_directory)
    payload = json.loads(input_json.read_text())
    records = payload["records"]
    metadata = payload["metadata"]
    if not records:
        raise ValueError("Complex-variant JSON must contain records")

    output_root.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".%s-" % timestamp, dir=output_root))
    try:
        summary_chunks = [records[i:i + 7] for i in range(0, len(records), 7)]
        summary_svgs = [
            render_summary_svg(
                chunk, metadata, page_index=index,
                page_count=len(summary_chunks))
            for index, chunk in enumerate(summary_chunks, start=1)
        ]
        provenance_svg = render_provenance_svg(records, metadata)
        page_entries = [
            ("summary" if index == 1 else "summary-%d" % index, svg, None)
            for index, svg in enumerate(summary_svgs, start=1)
        ]
        page_entries.extend(
            (_slug(record["id"]), render_complex_result_svg(record), record)
            for record in records)
        page_entries.append(("provenance", provenance_svg, None))
        pages = [svg for _, svg, _ in page_entries]
        page_files = []
        for index, (name, svg, record) in enumerate(page_entries, start=1):
            directory = staging / ("%02d-%s" % (index, name))
            directory.mkdir()
            svg_path = directory / "result.svg"
            png_path = directory / "result.png"
            svg_path.write_text(svg)
            _write_png(svg, png_path, png_scale)
            page_files.append({
                "page": index,
                "name": name,
                "svg": str(svg_path.relative_to(staging)),
                "png": str(png_path.relative_to(staging)),
            })
            if record is not None:
                record_path = directory / "record.json"
                record_path.write_text(
                    json.dumps(record, indent=2, sort_keys=True) + "\n")
        pdf_path = staging / "vaxrank-all-figures.pdf"
        _write_multipage_pdf(pages, pdf_path)
        manifest = {
            "schema_version": 1,
            "created_utc": timestamp,
            "vaxrank_version": __version__,
            "source_json": str(input_json),
            "source_json_sha256": sha256(input_json.read_bytes()).hexdigest(),
            "page_count": len(pages),
            "png_scale": png_scale,
            "png_dimensions": {"width": round(1200 * png_scale),
                               "height": round(760 * png_scale)},
            "pages": page_files,
            "combined_pdf": pdf_path.name,
        }
        (staging / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        staging.rename(run_directory)
        if combined_output:
            combined_output = Path(combined_output)
            combined_output.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(run_directory / pdf_path.name, combined_output)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return run_directory
