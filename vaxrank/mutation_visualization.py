"""Publication-ready views of annotation and RNA-assembled protein context."""

from __future__ import annotations

from datetime import datetime, timezone
from hashlib import sha256
from html import escape
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import tempfile

import pandas as pd

from .version import __version__


REQUIRED_COLUMNS = {
    "variant",
    "predicted_effect",
    "predicted_effect_gene_name",
    "predicted_effect_transcript_id",
    "predicted_effect_original_protein_sequence",
    "predicted_effect_mutant_protein_sequence",
    "predicted_effect_aa_mutation_start_offset",
    "protein_sequence",
    "protein_sequence_gene_names",
    "protein_sequence_gene_ids",
    "protein_sequence_transcript_names",
    "protein_sequence_transcript_ids",
    "num_alt_fragments",
    "num_ref_fragments",
    "num_other_fragments",
    "num_fragments_supporting_top_protein_sequence",
}

RUN_TIMESTAMP_PATTERN = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{6}Z$")


def _value(record, name, default=None):
    value = record.get(name, default)
    if value is None:
        return default
    try:
        if pd.isna(value):
            return default
    except (TypeError, ValueError):
        pass
    return value


def _text(record, name, default=""):
    value = _value(record, name, default)
    return str(value) if value is not None else default


def _integer(record, name, default=0):
    value = _value(record, name, default)
    return int(value) if value is not None else default


def _slug(text):
    slug = re.sub(r"[^A-Za-z0-9]+", "-", text).strip("-")
    return slug[:100] or "variant"


def _json_value(value):
    if value is None:
        return None
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if hasattr(value, "item"):
        value = value.item()
    return value


def _bounded_context(sequence, start, flank=24):
    """Return a bounded protein window centered on a mutation offset."""
    if not sequence:
        return ""
    start = min(max(0, start), len(sequence))
    left = max(0, start - flank)
    right = min(len(sequence), start + flank + 1)
    return sequence[left:right]


def _protein_context(record, sequence_column, flank=24):
    sequence = _text(record, sequence_column)
    start = _integer(record, "predicted_effect_aa_mutation_start_offset")
    return _bounded_context(sequence, start, flank=flank)


def _display_sequences(record):
    reference = _text(record, "trimmed_reference_protein_sequence")
    annotation = _text(record, "trimmed_predicted_mutant_protein_sequence")
    assembled = _text(record, "protein_sequence")
    if reference and annotation:
        start = _integer(record, "protein_sequence_mutation_start_idx")
        reference = _bounded_context(reference, start)
        annotation = _bounded_context(annotation, start)
        assembled = _bounded_context(assembled, start)
    else:
        reference = _protein_context(
            record, "predicted_effect_original_protein_sequence")
        annotation = _protein_context(
            record, "predicted_effect_mutant_protein_sequence")
        assembled = _bounded_context(
            assembled, _integer(record, "protein_sequence_mutation_start_idx"))
    return reference, annotation, assembled


def global_alignment(left, right):
    """Deterministic Needleman-Wunsch alignment for short protein contexts."""
    match, mismatch, gap = 2, -1, -2
    scores = [[0] * (len(right) + 1) for _ in range(len(left) + 1)]
    moves = [[""] * (len(right) + 1) for _ in range(len(left) + 1)]
    for i in range(1, len(left) + 1):
        scores[i][0] = i * gap
        moves[i][0] = "up"
    for j in range(1, len(right) + 1):
        scores[0][j] = j * gap
        moves[0][j] = "left"
    for i, left_residue in enumerate(left, 1):
        for j, right_residue in enumerate(right, 1):
            candidates = [
                (scores[i - 1][j - 1] + (
                    match if left_residue == right_residue else mismatch), "diag"),
                (scores[i - 1][j] + gap, "up"),
                (scores[i][j - 1] + gap, "left"),
            ]
            scores[i][j], moves[i][j] = max(candidates, key=lambda item: item[0])
    aligned_left = []
    aligned_right = []
    i, j = len(left), len(right)
    while i or j:
        move = moves[i][j]
        if move == "diag":
            aligned_left.append(left[i - 1])
            aligned_right.append(right[j - 1])
            i -= 1
            j -= 1
        elif move == "up":
            aligned_left.append(left[i - 1])
            aligned_right.append("–")
            i -= 1
        else:
            aligned_left.append("–")
            aligned_right.append(right[j - 1])
            j -= 1
    return "".join(reversed(aligned_left)), "".join(reversed(aligned_right))


def _sequence_row(label, sequence, comparison, y, x=252, residue_width=14):
    parts = [
        f'<text x="32" y="{y + 15}" class="track-label">{escape(label)}</text>'
    ]
    for index, residue in enumerate(sequence):
        different = residue != comparison[index]
        fill = "#fde7df" if different else "#f4f7fa"
        color = "#a23b21" if different else "#17212b"
        parts.append(
            f'<rect x="{x + index * residue_width}" y="{y}" width="13" height="22" '
            f'rx="2" fill="{fill}"/>'
        )
        parts.append(
            f'<text x="{x + index * residue_width + 6.5}" y="{y + 16}" '
            f'class="residue" fill="{color}">{escape(residue)}</text>'
        )
    return "".join(parts)


def _comparison_panel(title, left_label, left, right_label, right, y):
    aligned_left, aligned_right = global_alignment(left, right)
    return "".join([
        f'<text x="32" y="{y}" class="panel-title">{escape(title)}</text>',
        _sequence_row(left_label, aligned_left, aligned_right, y + 18),
        _sequence_row(right_label, aligned_right, aligned_left, y + 48),
    ])


def _classification(annotation, assembled):
    if not assembled:
        return (
            "WITHHELD",
            "No RNA-assembled protein sequence",
            "#755400",
            "#fff2c2",
        )
    if annotation == assembled:
        return (
            "CONFIRMED",
            "RNA assembly agrees with the displayed annotation-only sequence",
            "#176b4d",
            "#dff4e9",
        )
    return (
        "REFINED",
        "RNA assembly changes the local protein sequence",
        "#185f8d",
        "#dceefa",
    )


def _source_provenance(record, assembled):
    annotation = {
        "gene_names": _text(record, "predicted_effect_gene_name", "unavailable"),
        "gene_ids": _text(record, "predicted_effect_gene_id", "unavailable"),
        "transcript_names": _text(
            record, "predicted_effect_transcript_name", "unavailable"),
        "transcript_ids": _text(
            record, "predicted_effect_transcript_id", "unavailable"),
    }
    if assembled:
        assembly = {
            "gene_names": _text(
                record, "protein_sequence_gene_names", "unavailable"),
            "gene_ids": _text(record, "protein_sequence_gene_ids", "unavailable"),
            "transcript_names": _text(
                record, "protein_sequence_transcript_names", "unavailable"),
            "transcript_ids": _text(
                record, "protein_sequence_transcript_ids", "unavailable"),
        }
    else:
        assembly = {
            "gene_names": "no assembled protein",
            "gene_ids": "not applicable",
            "transcript_names": "not applicable",
            "transcript_ids": "not applicable",
        }
    return annotation, assembly


def _source_block(label, provenance, x):
    return "".join([
        f'<text x="{x}" y="88" class="source-label">{escape(label)}</text>',
        f'<text x="{x}" y="106" class="source-text">genes  '
        f'{escape(provenance["gene_names"])}  ·  transcripts  '
        f'{escape(provenance["transcript_ids"])}</text>',
        f'<text x="{x}" y="122" class="source-ids">gene IDs  '
        f'{escape(provenance["gene_ids"])}  ·  transcript names  '
        f'{escape(provenance["transcript_names"])}</text>',
    ])


def _no_assembly_message(record):
    num_alt_fragments = _integer(record, "num_alt_fragments")
    if num_alt_fragments:
        return (
            "%d alternate RNA fragments observed; no protein assembled."
            % num_alt_fragments,
            "Assembly or translation produced no top protein; the annotation-only "
            "sequence is not labeled RNA-supported.",
        )
    return (
        "No alternate RNA fragments; no protein could be assembled.",
        "The annotation-only mutant is not presented as RNA-supported sequence.",
    )


def render_mutation_svg(record):
    """Render one Isovar CSV record as a self-contained SVG string."""
    reference, annotation, assembled = _display_sequences(record)
    state, state_detail, state_color, state_background = _classification(
        annotation, assembled)
    gene = _text(record, "figure_label") or _text(
        record, "predicted_effect_gene_name", "Unknown gene")
    effect = _text(record, "predicted_effect", "Unknown effect")
    variant = _text(record, "variant", "Unknown variant")
    note = _text(record, "figure_note")
    source = _text(record, "figure_source_url")
    annotation_source, assembly_source = _source_provenance(record, assembled)
    no_assembly_title, no_assembly_detail = _no_assembly_message(record)
    second_panel = (
        _comparison_panel(
            "2 · WHAT RNA ASSEMBLY CHANGES",
            "Annotation only",
            annotation,
            "RNA assembled",
            assembled,
            470,
        )
        if assembled
        else "".join([
            '<text x="32" y="470" class="panel-title">2 · RNA ASSEMBLY DECISION</text>',
            '<rect x="32" y="491" width="1136" height="66" rx="7" fill="#fff7dc"/>',
            f'<text x="52" y="520" class="empty-title">{escape(no_assembly_title)}</text>',
            f'<text x="52" y="544" class="body">{escape(no_assembly_detail)}</text>',
        ])
    )
    evidence = (
        f"ALT fragments  {_integer(record, 'num_alt_fragments')}     "
        f"REF fragments  {_integer(record, 'num_ref_fragments')}     "
        f"OTHER  {_integer(record, 'num_other_fragments')}     "
        f"Top protein support  {_integer(record, 'num_fragments_supporting_top_protein_sequence')}"
    )
    footer = " · ".join(part for part in [note, source] if part)
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="700" viewBox="0 0 1200 700">
<rect width="1200" height="700" fill="#ffffff"/>
<style>
text {{ font-family: Inter, Helvetica, Arial, sans-serif; }}
.title {{ font-size: 26px; font-weight: 700; fill: #17212b; }}
.subtitle {{ font-size: 14px; fill: #445364; }}
.source-label {{ font-size: 10px; font-weight: 700; letter-spacing: 1px; fill: #526273; }}
.source-text {{ font-size: 11px; font-weight: 600; fill: #263544; }}
.source-ids {{ font-size: 9px; fill: #687786; }}
.status {{ font-size: 13px; font-weight: 700; letter-spacing: 1px; }}
.status-detail {{ font-size: 13px; font-weight: 600; }}
.panel-title {{ font-size: 12px; font-weight: 700; letter-spacing: 1.2px; fill: #526273; }}
.track-label {{ font-size: 13px; font-weight: 600; fill: #263544; }}
.residue {{ font-family: "SFMono-Regular", Consolas, monospace; font-size: 12px; font-weight: 600; text-anchor: middle; }}
.evidence {{ font-family: "SFMono-Regular", Consolas, monospace; font-size: 13px; fill: #253443; }}
.empty-title {{ font-size: 15px; font-weight: 700; fill: #755400; }}
.body {{ font-size: 13px; fill: #526273; }}
.footer {{ font-size: 10px; fill: #687786; }}
</style>
<text x="32" y="48" class="title">{escape(gene)} · {escape(effect)}</text>
<text x="32" y="75" class="subtitle">{escape(variant)}</text>
{_source_block("ANNOTATION SOURCE", annotation_source, 32)}
{_source_block("RNA-ASSEMBLY SOURCE", assembly_source, 620)}
<rect x="32" y="139" width="1136" height="52" rx="7" fill="{state_background}"/>
<text x="52" y="161" class="status" fill="{state_color}">{state}</text>
<text x="52" y="179" class="status-detail" fill="{state_color}">{escape(state_detail)}</text>
<rect x="32" y="211" width="1136" height="47" rx="7" fill="#eef3f7"/>
<text x="52" y="240" class="evidence">{escape(evidence)}</text>
{_comparison_panel("1 · WHAT ANNOTATION PREDICTS", "Reference", reference, "Annotation only", annotation, 294)}
{second_panel}
<line x1="32" y1="632" x2="1168" y2="632" stroke="#d8e0e7"/>
<rect x="32" y="653" width="12" height="12" rx="2" fill="#fde7df"/>
<text x="52" y="663" class="footer">mismatch or gap</text>
<text x="1168" y="663" text-anchor="end" class="footer">{escape(footer)}</text>
<text x="1168" y="682" text-anchor="end" class="footer">Vaxrank {escape(__version__)} · sequence evidence, not a clinical recommendation</text>
</svg>'''


def _configure_native_library_path():
    if platform.system() == "Darwin" and platform.machine() == "arm64":
        homebrew_lib = "/opt/homebrew/lib"
        existing = os.environ.get("DYLD_FALLBACK_LIBRARY_PATH", "")
        paths = [path for path in existing.split(":") if path]
        if homebrew_lib not in paths:
            os.environ["DYLD_FALLBACK_LIBRARY_PATH"] = ":".join(
                [homebrew_lib, *paths])


def _write_pdf(svg, path):
    _configure_native_library_path()
    from weasyprint import HTML

    html = (
        '<html><head><style>@page { size: 12.5in 7.3in; margin: 0; } '
        'body { margin: 0; background: white; } svg { width: 100%; height: auto; }'
        f'</style></head><body>{svg}</body></html>'
    )
    HTML(string=html).write_pdf(path)


def _write_png(svg, path, scale):
    _configure_native_library_path()
    from cairosvg import svg2png

    svg2png(
        bytestring=svg.encode("utf-8"),
        write_to=str(path),
        scale=scale,
        background_color="#ffffff",
    )


def _record_for_json(record):
    fields = [
        "variant",
        "predicted_effect",
        "predicted_effect_gene_name",
        "predicted_effect_gene_id",
        "predicted_effect_transcript_id",
        "predicted_effect_transcript_name",
        "protein_sequence",
        "protein_sequence_gene_names",
        "protein_sequence_gene_ids",
        "protein_sequence_transcript_names",
        "protein_sequence_transcript_ids",
        "protein_sequence_mutation_start_idx",
        "protein_sequence_mutation_end_idx",
        "trimmed_reference_protein_sequence",
        "trimmed_predicted_mutant_protein_sequence",
        "num_alt_fragments",
        "num_ref_fragments",
        "num_other_fragments",
        "num_fragments_supporting_top_protein_sequence",
        "figure_label",
        "figure_note",
        "figure_source_url",
        "figure_rna_sample",
    ]
    result = {name: _json_value(record.get(name)) for name in fields if name in record}
    for name in (
            "protein_sequence_mutation_start_idx",
            "protein_sequence_mutation_end_idx",
            "num_alt_fragments",
            "num_ref_fragments",
            "num_other_fragments",
            "num_fragments_supporting_top_protein_sequence"):
        if result.get(name) is not None:
            result[name] = int(result[name])
    reference, annotation, assembled = _display_sequences(record)
    result["displayed_sequences"] = {
        "reference": reference,
        "annotation_only": annotation,
        "rna_assembled": assembled or None,
    }
    annotation_source, assembly_source = _source_provenance(record, assembled)
    result["annotation_source"] = annotation_source
    result["rna_assembly_source"] = assembly_source if assembled else None
    result["assembly_outcome"] = _classification(annotation, assembled)[0].lower()
    return result


def generate_mutation_figures(
        input_csv, output_root, variants=(), formats=("svg", "pdf", "png"),
        timestamp=None, png_scale=3):
    """Create a timestamped figure run from a Vaxrank/Isovar CSV."""
    input_csv = Path(input_csv)
    output_root = Path(output_root)
    timestamp = timestamp or datetime.now(timezone.utc).strftime("%Y-%m-%dT%H%M%SZ")
    if not RUN_TIMESTAMP_PATTERN.fullmatch(timestamp):
        raise ValueError(
            "Timestamp must use UTC form YYYY-MM-DDTHHMMSSZ, got %r" % timestamp)
    run_directory = output_root / timestamp
    if run_directory.exists():
        raise FileExistsError("Figure run already exists: %s" % run_directory)
    dataframe = pd.read_csv(input_csv)
    missing = sorted(REQUIRED_COLUMNS.difference(dataframe.columns))
    if missing:
        raise ValueError("Isovar CSV is missing required columns: %s" % ", ".join(missing))
    if variants:
        search = dataframe.astype(str).apply(
            lambda row: " ".join(row.values).lower(), axis=1)
        keep = pd.Series(False, index=dataframe.index)
        for query in variants:
            keep |= search.str.contains(str(query).lower(), regex=False)
        dataframe = dataframe[keep]
    if dataframe.empty:
        raise ValueError("No variants matched the requested figure selection")
    formats = tuple(dict.fromkeys(formats))
    unknown_formats = set(formats).difference({"svg", "pdf", "png"})
    if unknown_formats:
        raise ValueError("Unsupported figure formats: %s" % sorted(unknown_formats))
    if not 0 < png_scale <= 10:
        raise ValueError("PNG scale must be greater than 0 and at most 10")

    output_root.mkdir(parents=True, exist_ok=True)
    staging_directory = Path(tempfile.mkdtemp(
        prefix=".%s-" % timestamp, dir=output_root))
    figures = []
    used_slugs = set()
    try:
        for _, series in dataframe.iterrows():
            record = series.to_dict()
            label = _text(record, "figure_label") or _text(
                record, "predicted_effect_gene_name", "variant")
            slug_base = _slug("%s-%s" % (label, _text(record, "variant")))
            slug = slug_base
            suffix = 2
            while slug in used_slugs:
                slug = "%s-%d" % (slug_base, suffix)
                suffix += 1
            used_slugs.add(slug)
            variant_directory = staging_directory / slug
            variant_directory.mkdir()
            svg = render_mutation_svg(record)
            files = []
            if "svg" in formats:
                svg_path = variant_directory / "protein-context.svg"
                svg_path.write_text(svg)
                files.append(str(svg_path.relative_to(staging_directory)))
            if "pdf" in formats:
                pdf_path = variant_directory / "protein-context.pdf"
                _write_pdf(svg, pdf_path)
                files.append(str(pdf_path.relative_to(staging_directory)))
            if "png" in formats:
                png_path = variant_directory / "protein-context.png"
                _write_png(svg, png_path, png_scale)
                files.append(str(png_path.relative_to(staging_directory)))
            record_path = variant_directory / "record.json"
            record_path.write_text(
                json.dumps(_record_for_json(record), indent=2, sort_keys=True) + "\n")
            files.append(str(record_path.relative_to(staging_directory)))
            figures.append({
                "label": label,
                "variant": _text(record, "variant"),
                "directory": slug,
                "files": files,
            })
        manifest = {
            "schema_version": 1,
            "created_utc": timestamp,
            "vaxrank_version": __version__,
            "source_csv": str(input_csv),
            "source_csv_sha256": sha256(input_csv.read_bytes()).hexdigest(),
            "formats": list(formats),
            "png_scale": png_scale if "png" in formats else None,
            "png_dimensions": ({
                "width": round(1200 * png_scale),
                "height": round(700 * png_scale),
            } if "png" in formats else None),
            "figures": figures,
        }
        (staging_directory / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        staging_directory.rename(run_directory)
    except BaseException:
        shutil.rmtree(staging_directory, ignore_errors=True)
        raise
    return run_directory
