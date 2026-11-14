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
    vaxrank replay ranked.json --output-dir replay_out/

Running a design takes flags only; `vaxrank run FLAGS` means the same thing.
Other activities start with a word:
    vaxrank policies [show NAME]   list the bundled selection policies, or show one
    vaxrank config                 print the default YAML configuration
    vaxrank replay FILE            re-emit a saved run (ranked .json or native evidence)
    vaxrank help [TOPIC]           options for one topic (below)

Help topics: %s, all (every option)."""

BARE_SUMMARY = """\
vaxrank ranks personalized cancer neoantigens from somatic variants, tumor RNA
and patient HLA type, and writes peptide or mRNA vaccine constructs plus reports.

%s
"""


# topic -> (what it covers, help-group titles it prints)
HELP_TOPICS = {
    "variants": ("variant input: VCF, MAF, genome", ("Variants",)),
    "isovar": ("RNA reads, assembly, translation and filtering (Isovar)", (
        "RNA", "Determine coding sequence from RNA",
        "Translation from cDNA to protein sequence",
        "Protein sequences (grouping equivalent translations)", "Filtering",
        "Matched germline", "Read phasing")),
    "cta": ("patient cancer-testis antigen expression input",
            ("Patient CTA expression input",)),
    "predict": ("MHC / T-cell epitope prediction and scoring", (
        "Prediction Options", "T-cell epitope prediction options",
        "Attribution of peptide-level evidence to alleles",
        "Default predictor methods (for unqualified DSL references)")),
    "design": ("what goes in the vaccine: limits, antigen types, LENS labels", (
        "Vaccine design options", "Vaccine peptide options",
        "Antigen type selection", "LENS source labels")),
    "peptide": ("peptide construct options", ("Peptide construct options",)),
    "mrna": ("mRNA construct options", ("mRNA construct options",)),
    "annotation": ("processing-aware annotation",
                   ("Processing-aware annotation options",)),
    "output": ("output files and reports",
               ("Output options", "Supplemental report options")),
    "advanced": ("advanced options", ("Advanced options",)),
}

ACTIVITIES = ("run", "policies", "config", "replay", "help")

EXAMPLES = EXAMPLES % ", ".join(HELP_TOPICS)
BARE_SUMMARY = BARE_SUMMARY % EXAMPLES


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


SHORT_USAGE = (
    "usage: vaxrank [run] FLAGS   run a design from variants and RNA, "
    "or an imported report\n"
    "       vaxrank policies [show NAME] | config | replay FILE | help [TOPIC]\n\n")


def print_short_help(stream=None):
    """The one-screen ``--help``: modes, activities and help topics."""
    (stream or sys.stdout).write(SHORT_USAGE + BARE_SUMMARY)


def _command_error(message):
    """One line on stderr, exit 2: the manners of an argparse usage error."""
    print("vaxrank: error: %s" % message, file=sys.stderr)
    raise SystemExit(2)


def topic_help(topic):
    """Help text for one topic, from the full parser's own argument groups."""
    from .arg_parser import make_vaxrank_arg_parser
    if topic not in HELP_TOPICS and topic != "all":
        _command_error("unknown help topic %r; choose from %s, all"
                       % (topic, ", ".join(HELP_TOPICS)))
    parser = make_vaxrank_arg_parser()
    if topic == "all":
        return parser.format_help()
    description, titles = HELP_TOPICS[topic]
    groups = {group.title: group for group in parser._action_groups}
    formatter = parser._get_formatter()
    formatter.add_text("vaxrank help %s: %s" % (topic, description))
    for title in titles:
        group = groups[title]
        formatter.start_section(title)
        formatter.add_text(group.description)
        formatter.add_arguments(group._group_actions)
        formatter.end_section()
    return formatter.format_help()


def resolve_command_line(args_list):
    """Interpret a bare first word and the help flags.

    Returns ``(handled, args_list)``. Running a design takes flags only. A
    first word that is not a flag names another activity: it is rewritten to
    the equivalent flags so one code path serves both spellings, or answered
    here (``help``). ``run`` is an optional spelling of the default. Help
    flags print the short screen; ``vaxrank help all`` lists every option.
    """
    tokens = list(args_list)
    if tokens and tokens[0] == "run":
        # Optional spelling of the default: what follows is an ordinary run.
        tokens = tokens[1:]
        if not tokens:
            print_bare_summary()
            return True, tokens
        if not tokens[0].startswith("-"):
            _command_error("run takes flags, starting with --; the other "
                           "commands are %s" % ", ".join(ACTIVITIES[1:]))
    end = tokens.index("--") if "--" in tokens else len(tokens)
    if tokens and not tokens[0].startswith("-"):
        word, rest = tokens[0], tokens[1:]
        if word == "help":
            if not rest:
                print_short_help()
            elif len(rest) == 1:
                sys.stdout.write(topic_help(rest[0]))
            else:
                _command_error("help takes one topic; choose from %s, all"
                               % ", ".join(HELP_TOPICS))
            return True, tokens
        if word == "policies":
            if not rest:
                return False, ["--list-policies"]
            if rest[0] == "show" and len(rest) == 2:
                return False, ["--show-policy", rest[1]]
            _command_error("policies takes nothing, or `show NAME`")
        if word == "config":
            if rest:
                _command_error("vaxrank config takes no arguments; use "
                                 "--config FILE on a run to apply one")
            return False, ["--print-default-config"]
        if word == "replay":
            if not rest or rest[0].startswith("-"):
                _command_error("replay needs the saved run first: vaxrank replay FILE [flags]")
            flag = ("--input-json-file" if rest[0].lower().endswith(".json")
                    else "--input-epitopes")
            return False, [flag, rest[0]] + rest[1:]
        _command_error(
            "unknown command %r; the commands are %s (a run takes flags, "
            "starting with --)" % (word, ", ".join(ACTIVITIES[1:])))
    head = tokens[:end]
    if "-h" in head or "--help" in head:
        print_short_help()
        raise SystemExit(0)
    return False, tokens


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
