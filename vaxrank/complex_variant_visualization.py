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
        while len(word) > width:
            if current:
                lines.append(" ".join(current))
                current = []
            lines.append(word[:width])
            word = word[width:]
        if not word:
            continue
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
{_text_lines('%s | %s' % (rna['platform'], rna['source']), 636, 248, css_class='strong', width=76, limit=1)}
{_text_lines(counts, 636, 272, width=82, limit=2, step=15)}
<text x="636" y="307" class="body">Vaxrank RNA gates: </text>
<text x="771" y="307" class="strong" fill="{filter_color}">{filters}</text>

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

{_text_lines(record.get('limitation', ''), 52, 646, css_class='muted', width=168, limit=2, step=15)}

<line x1="32" y1="690" x2="1168" y2="690" stroke="#d8e0e7"/>
<text x="32" y="714" class="footer">{escape(sources)}</text>
<text x="1168" y="737" text-anchor="end" class="footer">Vaxrank {escape(__version__)} | research use only</text>
</svg>'''


def render_platform_comparison_svg(record):
    """Render matched long- versus short-read evidence for one variant."""
    comparison = record["platform_comparison"]
    rows = comparison["rows"]
    if len(rows) != 2:
        raise ValueError("Platform comparison pages require exactly two rows")
    outcome_label, outcome_color, outcome_background = OUTCOME_STYLE[
        record["outcome"]]

    cards = []
    for index, row in enumerate(rows):
        x = 32 + index * 584
        cards.append(f'''
<rect x="{x}" y="318" width="552" height="236" rx="8" fill="#eef3f7"/>
<text x="{x + 20}" y="347" class="strong">{escape(row["label"])}</text>
{_text_lines(row['evidence'], x + 20, 372, width=73, limit=3, step=16)}
<text x="{x + 20}" y="425" class="panel-title">TRANSCRIPT NUCLEOTIDE / STRUCTURE</text>
{_text_lines(row['transcript_nt'], x + 20, 448, css_class='sequence', width=67, limit=3, step=17)}
<text x="{x + 20}" y="507" class="panel-title">PROTEIN AND PEPTIDE SPACE</text>
{_text_lines(row['protein'], x + 20, 530, css_class='sequence', width=67, limit=1)}
<text x="{x + 20}" y="550" class="footer">{escape(row['peptide_pool'])}</text>''')

    selected = record.get("selected_long_peptide")
    selection = (
        "Selected 25-aa construct: %s" % selected
        if selected else "No evidence-backed vaccine construct selected.")
    top_epitopes = record.get("top_epitopes", ())
    if top_epitopes:
        top = top_epitopes[0]
        selection += " Top binder: %s / %s / %.2f nM." % (
            top["sequence"], top["allele"], top["ic50_nm"])
    sources = " | ".join(record.get("source_labels", ()))
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="760" viewBox="0 0 1200 760">
<rect width="1200" height="760" fill="#ffffff"/>
<style>{_base_style()}</style>
<text x="32" y="45" class="kicker">SEQUENCING-PLATFORM COMPARISON</text>
<text x="32" y="78" class="title">{escape(record["title"])}</text>
<text x="32" y="103" class="subtitle">{escape(record["variant"])} | {escape(record["sample"])}</text>
<rect x="32" y="122" width="1136" height="55" rx="8" fill="{outcome_background}"/>
<text x="52" y="146" class="outcome" fill="{outcome_color}">{outcome_label}</text>
{_text_lines(comparison['conclusion'], 238, 145, width=105, limit=2, step=17)}

<text x="32" y="207" class="panel-title">DNA EVENT AND CODING MODEL</text>
<rect x="32" y="221" width="1136" height="68" rx="8" fill="#f8fafb"/>
{_text_lines(comparison['dna_model'], 52, 246, width=158, limit=2, step=17)}

<text x="32" y="309" class="panel-title">LONG READ</text>
<text x="616" y="309" class="panel-title">SHORT READ</text>
{''.join(cards)}

<rect x="32" y="578" width="1136" height="89" rx="8" fill="#e8f3fa"/>
<text x="52" y="603" class="panel-title">VACCINE CONSEQUENCE</text>
{_text_lines(comparison['vaccine_impact'], 52, 627, width=158, limit=2, step=17)}
{_text_lines(selection, 52, 658, css_class='strong', width=158, limit=1)}

<line x1="32" y1="690" x2="1168" y2="690" stroke="#d8e0e7"/>
<text x="32" y="714" class="footer">{escape(sources)}</text>
<text x="1168" y="737" text-anchor="end" class="footer">Vaxrank {escape(__version__)} | research use only</text>
</svg>'''


def render_platform_overview_svg(records):
    """Summarize which platform differences change antigen conclusions."""
    compared = [record for record in records if record.get("platform_comparison")]
    if len(compared) > 5:
        raise ValueError("Platform overview supports at most five rows")
    rows = []
    for index, record in enumerate(compared):
        comparison = record["platform_comparison"]
        y = 242 + index * 91
        fill = "#f5f8fa" if index % 2 else "#eef3f7"
        rows.append(
            f'<rect x="32" y="{y - 27}" width="1136" height="81" rx="7" fill="{fill}"/>'
            + _text_lines(
                record["title"], 52, y - 8,
                css_class="strong", width=30, limit=2, step=15)
            + _text_lines(
                comparison["overview"]["long"], 270, y - 8,
                width=42, limit=3, step=15)
            + _text_lines(
                comparison["overview"]["short"], 570, y - 8,
                width=42, limit=3, step=15)
            + _text_lines(
                comparison["overview"]["vaccine"], 870, y - 8,
                css_class="strong", width=38, limit=3, step=15))
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="760" viewBox="0 0 1200 760">
<rect width="1200" height="760" fill="#ffffff"/>
<style>{_base_style()}</style>
<text x="32" y="47" class="kicker">MATCHED PLATFORM AUDIT</text>
<text x="32" y="84" class="title">Sequencing platform changes evidence - not always the antigen</text>
<text x="32" y="112" class="subtitle">Counts are within-product evidence, not a quantitative comparison of library sensitivity.</text>
<rect x="32" y="139" width="1136" height="48" rx="8" fill="#dceefa"/>
<text x="52" y="168" class="strong">A platform difference matters to Vaxrank only when it changes a defensible mutant sequence or construct window.</text>
<text x="52" y="207" class="panel-title">VARIANT</text>
<text x="270" y="207" class="panel-title">LONG-READ RNA</text>
<text x="570" y="207" class="panel-title">SHORT-READ RNA</text>
<text x="870" y="207" class="panel-title">VACCINE CONSEQUENCE</text>
{''.join(rows)}
<line x1="32" y1="690" x2="1168" y2="690" stroke="#d8e0e7"/>
<text x="32" y="714" class="footer">Matched regional audit; exact source identifiers and transcript coordinates are retained in record.json.</text>
<text x="1168" y="737" text-anchor="end" class="footer">Vaxrank {escape(__version__)} | research use only</text>
</svg>'''


def render_orf_platform_summary_svg(audit):
    """Render full-matrix RNA platform coverage and protein conclusions."""
    rows = []
    for index, record in enumerate(
            audit["full_matrix"]["platform_summary"]):
        y = 342 + index * 63
        fill = "#eef3f7" if index % 2 == 0 else "#f7f9fb"
        rows.append(
            f'<rect x="32" y="{y - 31}" width="1136" height="53" '
            f'rx="7" fill="{fill}"/>'
            f'<text x="52" y="{y}" class="strong">'
            f'{escape(record["platform"])}</text>'
            f'<text x="230" y="{y}" class="body">'
            f'{record["catalogued_products"]}</text>'
            f'<text x="380" y="{y}" class="body">'
            f'{record["completed_source_variant_rows"]}</text>'
            f'<text x="555" y="{y}" class="body">'
            f'{record["variants_with_alt_evidence"]} / 44</text>'
            f'<text x="720" y="{y}" class="strong">'
            f'{record["variants_with_validated_protein_window"]} / 44</text>'
            f'<text x="900" y="{y}" class="body">'
            f'{record["variants_exact_in_at_least_one_product"]}</text>'
            f'<text x="1070" y="{y}" class="body">'
            f'{record["variants_nonexact_in_at_least_one_product"]}</text>')
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="760" viewBox="0 0 1200 760">
<rect width="1200" height="760" fill="#ffffff"/>
<style>{_base_style()}</style>
<text x="32" y="47" class="kicker">FULL RNA-PRODUCT AUDIT</text>
<text x="32" y="84" class="title">Which platforms establish a translated local coding window?</text>
<text x="32" y="112" class="subtitle">44 vaccine-included loci × 164 catalogued RNA products; Isovar assembly enabled</text>
<rect x="32" y="133" width="1136" height="55" rx="8" fill="#fff2c2"/>
{_text_lines('A protein result is an independently validated local RNA-derived window, not a full-length ORF. Products are technical/data products, not independent biological samples or a sensitivity benchmark.', 52, 157, css_class='strong', width=160, limit=2, step=17)}

<rect x="32" y="213" width="1136" height="65" rx="8" fill="#e8f3fa"/>
{_text_lines("Varcode supplies an exact reference-plus-isolated-edit protein for all 44 loci. RNA asks a different question: was a sample-specific coding frame observed, and did its local protein agree?", 52, 239, width=160, limit=2, step=18)}

<text x="52" y="300" class="table-head">PLATFORM</text>
<text x="230" y="300" class="table-head">PRODUCTS</text>
<text x="380" y="300" class="table-head">COMPLETED ROWS</text>
<text x="555" y="300" class="table-head">ALT LOCI</text>
<text x="720" y="300" class="table-head">PROTEIN LOCI</text>
<text x="900" y="300" class="table-head">EXACT ≥1</text>
<text x="1070" y="300" class="table-head">NON-EXACT ≥1</text>
{''.join(rows)}

<rect x="32" y="531" width="355" height="112" rx="8" fill="#eef3f7"/>
<text x="52" y="558" class="panel-title">ILMN</text>
{_text_lines('Broadest product inventory and 39/44 loci with a validated local protein window.', 52, 583, width=42, limit=3, step=18)}
<rect x="423" y="531" width="355" height="112" rx="8" fill="#eef3f7"/>
<text x="443" y="558" class="panel-title">ONT</text>
{_text_lines('30/44 loci yield a protein window; NTF3 is the sole non-exact isolated-effect result.', 443, 583, width=42, limit=3, step=18)}
<rect x="814" y="531" width="354" height="112" rx="8" fill="#f8e9e5"/>
<text x="834" y="558" class="panel-title">PACBIO</text>
{_text_lines('Only one of nine catalogued products has genomic coordinates: three protein loci. This is a data-availability limit.', 834, 583, width=42, limit=3, step=18)}

<line x1="32" y1="690" x2="1168" y2="690" stroke="#d8e0e7"/>
<text x="32" y="714" class="footer">Pinned Isovar full-product audit; unavailable rows remain unavailable rather than zero.</text>
<text x="1168" y="737" text-anchor="end" class="footer">Vaxrank {escape(__version__)} | research use only</text>
</svg>'''


def render_orf_assembly_summary_svg(audit):
    """Render the strict assembly-on/off comparison on identical read subsets."""
    summary = audit["paired_corpus"]["platform_mode_summary"]
    rows = []
    for index, record in enumerate(summary):
        y = 285 + index * 42
        fill = "#eef3f7" if index % 2 == 0 else "#f7f9fb"
        rows.append(
            f'<rect x="32" y="{y - 26}" width="1136" height="36" '
            f'rx="6" fill="{fill}"/>'
            f'<text x="52" y="{y}" class="strong">'
            f'{escape(record["platform"])}</text>'
            f'<text x="200" y="{y}" class="body">'
            f'{"ON" if record["mode"] == "assembly_on" else "OFF"}</text>'
            f'<text x="335" y="{y}" class="body">{record["cases"]}</text>'
            f'<text x="500" y="{y}" class="strong">'
            f'{record["validated_local_protein_windows"]}</text>'
            f'<text x="690" y="{y}" class="body">'
            f'{record["exact_isolated_edit_windows"]}</text>'
            f'<text x="865" y="{y}" class="body">'
            f'{record["nonexact_windows"]}</text>'
            f'<text x="1045" y="{y}" class="body">'
            f'{record["no_protein_window"]}</text>')
    changed = audit["paired_corpus"]["assembly_effect_summary"][0]
    changed_names = ", ".join(
        record["variant_id"].split("-chr", 1)[0]
        for record in changed["changed_cases"])
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="760" viewBox="0 0 1200 760">
<rect width="1200" height="760" fill="#ffffff"/>
<style>{_base_style()}</style>
<text x="32" y="47" class="kicker">PAIRED ASSEMBLY AUDIT</text>
<text x="32" y="84" class="title">Same original reads, Isovar assembly on versus off</text>
<text x="32" y="112" class="subtitle">49 checksummed source/variant BAM subsets; current Isovar runtime</text>
<rect x="32" y="133" width="1136" height="61" rx="8" fill="#fff2c2"/>
{_text_lines('This selected corpus is enriched for informative loci. It supports paired software behavior, not population or platform sensitivity estimates. PacBio contributes one case.', 52, 158, css_class='strong', width=160, limit=2, step=18)}

<text x="52" y="235" class="table-head">PLATFORM</text>
<text x="200" y="235" class="table-head">ASSEMBLY</text>
<text x="335" y="235" class="table-head">CASES</text>
<text x="500" y="235" class="table-head">PROTEIN WINDOWS</text>
<text x="690" y="235" class="table-head">EXACT</text>
<text x="865" y="235" class="table-head">NON-EXACT</text>
<text x="1045" y="235" class="table-head">NO WINDOW</text>
{''.join(rows)}

<rect x="32" y="522" width="552" height="126" rx="8" fill="#dff4e9"/>
<text x="52" y="550" class="panel-title">WHAT ASSEMBLY CHANGED</text>
{_text_lines('ILMN assembly changed local protein length in 5/32 cases while preserving whether a window was established. Cases: %s.' % changed_names, 52, 578, width=70, limit=4, step=18)}
<rect x="616" y="522" width="552" height="126" rx="8" fill="#eef3f7"/>
<text x="636" y="550" class="panel-title">WHAT IT DID NOT CHANGE</text>
{_text_lines('The assembly-on and assembly-off ORF counts are identical in this selected corpus. All 14 paired ONT proteins are sequence-identical because individual long reads already carry the local context.', 636, 578, width=70, limit=4, step=18)}

<line x1="32" y1="690" x2="1168" y2="690" stroke="#d8e0e7"/>
<text x="32" y="714" class="footer">Every case retains source URL, BAM identity, Varcode effects, sequences, and transcript edits in the checksummed source JSON.</text>
<text x="1168" y="737" text-anchor="end" class="footer">Vaxrank {escape(__version__)} | research use only</text>
</svg>'''


def render_orf_attribution_svg(audit):
    """Render non-exact protein results and current attribution limits."""
    rows = []
    for index, record in enumerate(audit["full_matrix"]["nonexact_records"]):
        y = 228 + index * 57
        fill = "#eef3f7" if index % 2 == 0 else "#f7f9fb"
        variant = record["variant_id"].split("-chr", 1)[0]
        rows.append(
            f'<rect x="32" y="{y - 27}" width="1136" height="49" '
            f'rx="7" fill="{fill}"/>'
            f'<text x="52" y="{y - 3}" class="strong">{escape(variant)}</text>'
            f'<text x="180" y="{y - 3}" class="body">'
            f'{escape(" + ".join(record["platforms"]))}</text>'
            f'<text x="300" y="{y - 3}" class="body">'
            f'{len(record["rows"])}</text>'
            + _text_lines(
                record["attribution_note"], 380, y - 8,
                width=103, limit=2, step=16))
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="760" viewBox="0 0 1200 760">
<rect width="1200" height="760" fill="#ffffff"/>
<style>{_base_style()}</style>
<text x="32" y="47" class="kicker">ISOLATED DNA EFFECT VERSUS RNA</text>
<text x="32" y="84" class="title">Six loci have at least one non-exact RNA protein window</text>
<text x="32" y="112" class="subtitle">Sixteen product rows; exact and non-exact products can coexist for one locus</text>
<rect x="32" y="133" width="1136" height="48" rx="8" fill="#e8f3fa"/>
<text x="52" y="162" class="strong">NTF3 is cross-platform and explained at sequence level; biological origin remains unresolved.</text>

<text x="52" y="193" class="table-head">LOCUS</text>
<text x="180" y="193" class="table-head">PLATFORM</text>
<text x="300" y="193" class="table-head">ROWS</text>
<text x="380" y="193" class="table-head">INTERPRETATION</text>
{''.join(rows)}

<rect x="32" y="559" width="552" height="104" rx="8" fill="#f8e9e5"/>
<text x="52" y="586" class="panel-title">GERMLINE ATTRIBUTION: NOT YET</text>
{_text_lines("The public run_isovar path has no matched-germline variant input. Additional transcript edits remain 'unexplained' even when they may be germline.", 52, 612, width=70, limit=3, step=17)}
<rect x="616" y="559" width="552" height="104" rx="8" fill="#fff2c2"/>
<text x="636" y="586" class="panel-title">CO-SOMATIC PHASING: CAPABLE, NOT TESTED HERE</text>
{_text_lines('Isovar can phase multiple supplied somatic variants by shared fragments, but the 44 nominated loci contain no nearby pair. MAP2 and CD109 use explicit external compound inputs.', 636, 612, width=70, limit=3, step=17)}

<line x1="32" y1="690" x2="1168" y2="690" stroke="#d8e0e7"/>
<text x="32" y="714" class="footer">Exact local sequences and transcript-relative edits are retained in the checksummed source JSON.</text>
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


def _balanced_chunks(records, maximum_size=7):
    """Split records into the fewest summary pages with balanced row counts."""
    page_count = (len(records) + maximum_size - 1) // maximum_size
    base_size, extra = divmod(len(records), page_count)
    sizes = [base_size + (index < extra) for index in range(page_count)]
    chunks = []
    start = 0
    for size in sizes:
        chunks.append(records[start:start + size])
        start += size
    return chunks


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
    orf_platform_audit = payload.get("orf_platform_audit")
    if not records:
        raise ValueError("Complex-variant JSON must contain records")

    output_root.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".%s-" % timestamp, dir=output_root))
    try:
        summary_chunks = _balanced_chunks(records)
        summary_svgs = [
            render_summary_svg(
                chunk, metadata, page_index=index,
                page_count=len(summary_chunks))
            for index, chunk in enumerate(summary_chunks, start=1)
        ]
        platform_overview_svg = (
            render_platform_overview_svg(records)
            if any(record.get("platform_comparison") for record in records)
            else None)
        provenance_svg = render_provenance_svg(records, metadata)
        page_entries = [
            ("summary" if index == 1 else "summary-%d" % index, svg, None)
            for index, svg in enumerate(summary_svgs, start=1)
        ]
        if platform_overview_svg:
            page_entries.append((
                "platform-overview", platform_overview_svg, None))
        if orf_platform_audit:
            page_entries.extend([
                ("orf-platform-summary",
                 render_orf_platform_summary_svg(orf_platform_audit),
                 None),
                ("orf-assembly-comparison",
                 render_orf_assembly_summary_svg(orf_platform_audit),
                 None),
                ("orf-attribution",
                 render_orf_attribution_svg(orf_platform_audit),
                 None),
            ])
        page_entries.extend(
            (
                _slug(record["id"]),
                render_platform_comparison_svg(record)
                if record.get("platform_comparison")
                else render_complex_result_svg(record),
                record,
            )
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
