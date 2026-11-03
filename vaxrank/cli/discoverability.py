# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Help a user find the right command (openvax/vaxrank#772).

Worked examples for ``--help``, a short screen for a bare ``vaxrank``, a
"did you mean" for mistyped flags, and the closing list of what a run wrote.
"""

import difflib
import os
import sys
from argparse import HelpFormatter

EXAMPLES = """\
Examples:
  Rank variants against tumor RNA and write peptide and mRNA constructs:
    vaxrank --vcf somatic.vcf --bam rnaseq.bam \\
        --mhc-predictor netmhcpan --mhc-alleles HLA-A*02:01,HLA-B*07:02 \\
        --output-dir vaccine_out/

  Import a pVACseq or LENS report (no BAM needed):
    vaxrank --input-pvacseq patient.all_epitopes.tsv --output-dir vaccine_out/
    vaxrank --input-lens patient.lens.tsv --output-dir vaccine_out/

  Replay a run saved earlier with --output-json-file:
    vaxrank --input-json-file ranked.json --output-dir replay_out/

Run `vaxrank --help` for every option, or `vaxrank --list-policies` for the
bundled selection policies."""

BARE_SUMMARY = """\
vaxrank ranks personalized cancer neoantigens from somatic variants, tumor RNA
and patient HLA type, and writes peptide or mRNA vaccine constructs plus reports.

%s
""" % EXAMPLES


class ExamplesHelpFormatter(HelpFormatter):
    """Wrap descriptions as usual but keep the examples epilog as written."""

    def _fill_text(self, text, width, indent):
        if text.startswith("Examples:"):
            return text
        return super()._fill_text(text, width, indent)


def suggest_flag(flag, known_flags, cutoff=0.6):
    """Return the closest known option to ``flag``, or None."""
    matches = difflib.get_close_matches(flag, known_flags, n=1, cutoff=cutoff)
    return matches[0] if matches else None


def check_unknown_flags(parser, args_list):
    """Fail with a suggestion when a long option is not one the parser knows.

    argparse only reports an unknown option after it has reported missing
    required ones, so ``--bamm x.bam`` surfaces as "required: --bam". Looking
    at the flags first names the mistake. Unique prefixes stay legal
    (``allow_abbrev``).
    """
    known = {name for action in parser._actions for name in action.option_strings}
    longs = sorted(name for name in known if name.startswith("--"))
    for token in args_list:
        if token == "--":
            return
        if not token.startswith("--"):
            continue
        flag = token.split("=", 1)[0]
        if flag in known:
            continue
        if parser.allow_abbrev and any(name.startswith(flag) for name in longs):
            continue
        suggestion = suggest_flag(flag, longs)
        hint = " (did you mean %s?)" % suggestion if suggestion else ""
        parser.error("unrecognized argument: %s%s" % (flag, hint))


def print_bare_summary(stream=None):
    (stream or sys.stdout).write(BARE_SUMMARY)


OUTPUT_LABELS = (
    ("ranked candidates", "output_csv"),
    ("ranked workbook", "output_xlsx_report"),
    ("prediction table", "output_neoepitope_report"),
    ("native evidence", "output_epitopes"),
    ("antigen export", "output_antigens"),
    ("policy exports", "output_policy_dir"),
    ("input predictions", "output_input_predictions"),
    ("JSON data", "output_json_file"),
    ("ASCII report", "output_ascii_report"),
    ("HTML report", "output_html_report"),
    ("PDF report", "output_pdf_report"),
)


def written_outputs(args, construct_dirs=()):
    """(label, path) pairs for the requested outputs that exist on disk."""
    found = []
    for label, attr in OUTPUT_LABELS:
        path = getattr(args, attr, "") or ""
        if path and os.path.exists(path):
            found.append((label, path))
    for label, path in construct_dirs:
        if path and os.path.exists(path):
            found.append((label, path))
    output_dir = getattr(args, "output_dir", "") or ""
    summary = os.path.join(output_dir, "run_summary.txt") if output_dir else ""
    if summary and os.path.exists(summary):
        found.append(("run summary", summary))
    return found


def format_outputs_block(entries):
    if not entries:
        return ""
    lines = ["Outputs written:"]
    lines += ["  %-18s %s" % (label + ":", path) for label, path in entries]
    return "\n".join(lines) + "\n"
