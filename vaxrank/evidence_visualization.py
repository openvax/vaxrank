"""Publication-ready comparison of DNA, short-read, and long-read evidence."""

from datetime import datetime, timezone
from hashlib import sha256
from html import escape
import json
from pathlib import Path
import re
import shutil
import tempfile

from .mutation_visualization import _slug, _write_pdf, _write_png
from .version import __version__


RUN_TIMESTAMP_PATTERN = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{6}Z$")
EVIDENCE_STAGES = (
    ("dna", "DNA"),
    ("short_read", "SHORT-READ RNA"),
    ("short_read_assembly", "SHORT-READ + ASSEMBLY"),
    ("long_read", "LONG-READ RNA"),
)
STATUS_STYLE = {
    "detected": ("DETECTED", "#176b4d", "#dff4e9"),
    "not_detected": ("NOT DETECTED", "#7a2938", "#f8e3e7"),
    "not_assessed": ("NOT ASSESSED", "#596776", "#eef2f5"),
    "unresolved": ("UNRESOLVED", "#755400", "#fff2c2"),
}


def _lines(text, width=42, limit=3):
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
        lines[-1] = lines[-1].rstrip(".") + "…"
    return lines


def _evidence_card(stage, label, evidence, x):
    status = evidence.get("status", "not_assessed")
    if status not in STATUS_STYLE:
        raise ValueError("Unknown evidence status %r for %s" % (status, stage))
    status_label, color, background = STATUS_STYLE[status]
    details = "".join(
        f'<text x="{x + 18}" y="{213 + index * 17}" class="card-body">'
        f'{escape(line)}</text>'
        for index, line in enumerate(_lines(evidence.get("detail", "")))
    )
    return f'''
<rect x="{x}" y="140" width="272" height="132" rx="8" fill="{background}"/>
<text x="{x + 18}" y="166" class="card-label">{escape(label)}</text>
<text x="{x + 18}" y="190" class="card-status" fill="{color}">{status_label}</text>
{details}'''


def _nucleotide_track(record):
    nucleotide = record["nucleotide"]
    left = nucleotide["left_sequence"][-42:]
    insert = nucleotide.get("insert_sequence", "")
    right = nucleotide["right_sequence"][:42]
    x = 185
    width = 9
    blocks = []
    for base in left:
        blocks.append(
            f'<rect x="{x}" y="363" width="8" height="24" rx="2" fill="#dceefa"/>'
            f'<text x="{x + 4}" y="380" class="base">{escape(base)}</text>')
        x += width
    junction_x = x
    for base in insert:
        blocks.append(
            f'<rect x="{x}" y="363" width="8" height="24" rx="2" fill="#eee4fa"/>'
            f'<text x="{x + 4}" y="380" class="base">{escape(base)}</text>')
        x += width
    for base in right:
        blocks.append(
            f'<rect x="{x}" y="363" width="8" height="24" rx="2" fill="#fde7df"/>'
            f'<text x="{x + 4}" y="380" class="base">{escape(base)}</text>')
        x += width
    insert_label = (
        f' · {len(insert)}-nt insert' if insert else ""
    )
    return "".join([
        '<text x="32" y="319" class="panel-title">1 · OBSERVED OR INFERRED JUNCTION SEQUENCE</text>',
        f'<text x="32" y="345" class="body">{escape(nucleotide["left_label"])} to '
        f'{escape(nucleotide["right_label"])}{escape(insert_label)}</text>',
        *blocks,
        f'<line x1="{junction_x - 1}" y1="355" x2="{junction_x - 1}" y2="397" stroke="#8b4a25" stroke-width="2"/>',
        f'<text x="{junction_x - 1}" y="416" text-anchor="middle" class="junction">junction</text>',
        f'<text x="185" y="439" class="legend">{escape(nucleotide["left_label"])} side</text>',
        f'<text x="{x}" y="439" text-anchor="end" class="legend">{escape(nucleotide["right_label"])} side</text>',
    ])


def render_evidence_svg(record):
    """Render one cross-platform sequence-evidence record as SVG."""
    evidence = record["evidence"]
    cards = "".join(
        _evidence_card(stage, label, evidence.get(stage, {}), 32 + index * 284)
        for index, (stage, label) in enumerate(EVIDENCE_STAGES)
    )
    protein = record.get("protein", {})
    protein_sequence = protein.get("sequence")
    if protein_sequence:
        protein_panel = (
            '<rect x="32" y="492" width="1136" height="86" rx="8" fill="#dff4e9"/>'
            '<text x="52" y="520" class="protein-status" fill="#176b4d">TRANSLATED</text>'
            f'<text x="52" y="552" class="protein-sequence">{escape(protein_sequence)}</text>')
    else:
        protein_panel = (
            '<rect x="32" y="492" width="1136" height="86" rx="8" fill="#fff2c2"/>'
            '<text x="52" y="520" class="protein-status" fill="#755400">PROTEIN WITHHELD</text>'
            f'<text x="52" y="552" class="body">{escape(protein.get("reason", "Coding frame not established."))}</text>')
    sources = " · ".join(record.get("source_urls", ()))
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="760" viewBox="0 0 1200 760">
<rect width="1200" height="760" fill="#ffffff"/>
<style>
text {{ font-family: Inter, Helvetica, Arial, sans-serif; }}
.title {{ font-size: 27px; font-weight: 700; fill: #17212b; }}
.subtitle {{ font-size: 14px; fill: #445364; }}
.conclusion {{ font-size: 13px; font-weight: 700; letter-spacing: 1px; fill: #185f8d; }}
.card-label, .panel-title {{ font-size: 11px; font-weight: 700; letter-spacing: 1px; fill: #526273; }}
.card-status {{ font-size: 14px; font-weight: 700; }}
.card-body, .body {{ font-size: 12px; fill: #445364; }}
.base {{ font-family: "SFMono-Regular", Consolas, monospace; font-size: 11px; font-weight: 600; text-anchor: middle; fill: #17212b; }}
.junction {{ font-size: 10px; font-weight: 700; fill: #8b4a25; }}
.legend {{ font-size: 10px; fill: #687786; }}
.protein-status {{ font-size: 12px; font-weight: 700; letter-spacing: 1px; }}
.protein-sequence {{ font-family: "SFMono-Regular", Consolas, monospace; font-size: 15px; font-weight: 600; fill: #17212b; }}
.interpretation {{ font-size: 13px; font-weight: 600; fill: #263544; }}
.footer {{ font-size: 9px; fill: #687786; }}
</style>
<text x="32" y="48" class="title">{escape(record["title"])}</text>
<text x="32" y="75" class="subtitle">{escape(record["kind"])} · {escape(record["locus"])}</text>
<rect x="32" y="94" width="1136" height="30" rx="7" fill="#dceefa"/>
<text x="50" y="114" class="conclusion">{escape(record["conclusion"])}</text>
{cards}
{_nucleotide_track(record)}
<text x="32" y="476" class="panel-title">2 · PROTEIN-LEVEL INTERPRETATION</text>
{protein_panel}
<rect x="32" y="600" width="1136" height="65" rx="8" fill="#eef3f7"/>
<text x="52" y="625" class="panel-title">SCIENTIFIC INTERPRETATION</text>
<text x="52" y="648" class="interpretation">{escape(record["interpretation"])}</text>
<line x1="32" y1="694" x2="1168" y2="694" stroke="#d8e0e7"/>
<text x="32" y="717" class="footer">{escape(sources)}</text>
<text x="1168" y="739" text-anchor="end" class="footer">Vaxrank {escape(__version__)} · sequence evidence, not a clinical recommendation</text>
</svg>'''


def generate_evidence_figures(
        input_json, output_root, formats=("svg", "pdf", "png"), timestamp=None,
        png_scale=3):
    """Create a timestamped run of cross-platform evidence figures."""
    input_json = Path(input_json)
    output_root = Path(output_root)
    timestamp = timestamp or datetime.now(timezone.utc).strftime("%Y-%m-%dT%H%M%SZ")
    if not RUN_TIMESTAMP_PATTERN.fullmatch(timestamp):
        raise ValueError("Timestamp must use UTC form YYYY-MM-DDTHHMMSSZ")
    run_directory = output_root / timestamp
    if run_directory.exists():
        raise FileExistsError("Figure run already exists: %s" % run_directory)
    records = json.loads(input_json.read_text())
    if not isinstance(records, list) or not records:
        raise ValueError("Evidence JSON must be a non-empty list")
    formats = tuple(dict.fromkeys(formats))
    unknown = set(formats).difference({"svg", "pdf", "png"})
    if unknown:
        raise ValueError("Unsupported figure formats: %s" % sorted(unknown))
    if not 0 < png_scale <= 10:
        raise ValueError("PNG scale must be greater than 0 and at most 10")

    output_root.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".%s-" % timestamp, dir=output_root))
    figures = []
    try:
        for record in records:
            slug = _slug(record["id"])
            directory = staging / slug
            directory.mkdir()
            svg = render_evidence_svg(record)
            files = []
            for output_format in formats:
                path = directory / ("evidence-context.%s" % output_format)
                if output_format == "svg":
                    path.write_text(svg)
                elif output_format == "pdf":
                    _write_pdf(svg, path)
                else:
                    _write_png(svg, path, png_scale)
                files.append(str(path.relative_to(staging)))
            record_path = directory / "record.json"
            record_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
            files.append(str(record_path.relative_to(staging)))
            figures.append({"id": record["id"], "directory": slug, "files": files})
        manifest = {
            "schema_version": 1,
            "created_utc": timestamp,
            "vaxrank_version": __version__,
            "source_json": str(input_json),
            "source_json_sha256": sha256(input_json.read_bytes()).hexdigest(),
            "formats": list(formats),
            "png_scale": png_scale if "png" in formats else None,
            "png_dimensions": ({"width": round(1200 * png_scale),
                                "height": round(760 * png_scale)}
                               if "png" in formats else None),
            "figures": figures,
        }
        (staging / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        staging.rename(run_directory)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return run_directory
