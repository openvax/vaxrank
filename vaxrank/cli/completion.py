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

"""Shell tab completion (argcomplete) helpers for the vaxrank command line."""

import argparse
from importlib.resources import files
import os

import argcomplete

from .discoverability import ACTIVITIES, HELP_TOPICS
from .help import HelpArgumentParser


def is_completing(environ=None):
    """Return whether this process is answering a shell completion request."""
    environ = os.environ if environ is None else environ
    return "_ARGCOMPLETE" in environ


def _split_completion_line(environ=None):
    """Return ``(words, current_prefix)`` from argcomplete environment variables."""
    environ = os.environ if environ is None else environ
    comp_line = environ.get("COMP_LINE", "")
    comp_point = int(environ.get("COMP_POINT", len(comp_line)))
    _, current_prefix, _, words, _ = argcomplete.split_line(comp_line, comp_point)
    start_index = max(0, int(environ.get("_ARGCOMPLETE", "1")) - 1)
    return words[start_index:], current_prefix


def _bundled_policy_names():
    """Return bundled policy names without importing YAML/msgspec."""
    names = []
    for path in files("vaxrank.config").iterdir():
        if path.name.endswith(".yaml") and path.name != "default.yaml":
            names.append(path.name[:-5])
    return sorted(names)


def _attach_subparser(subparsers, name, sub_parser, help_text=""):
    """Attach an existing parser (including ones with nested subparsers)."""
    has_nested = any(
        isinstance(action, argparse._SubParsersAction)
        for action in sub_parser._actions
    )
    if has_nested:
        subparsers._name_parser_map[name] = sub_parser
        subparsers._choices_actions.append(
            subparsers._ChoicesPseudoAction(name, (), help_text)
        )
    else:
        subparsers.add_parser(
            name,
            help=help_text,
            add_help=False,
            parents=[sub_parser],
        )


def build_completion_parser(environ=None):
    """Build an import-light parser that answers a vaxrank tab-completion request.

    Inspects ``COMP_LINE`` / ``COMP_POINT`` so that:
    - Flag-first invocations (``vaxrank --...``, ``vaxrank --vcf a.vcf --...``,
      ``vaxrank --input-lens x.tsv --...``) complete against the mode-appropriate
      vaxrank argument parser directly.
    - Activity invocations (``vaxrank <TAB>``, ``vaxrank data ...``,
      ``vaxrank compare ...``, ``vaxrank mhc-status ...``, ``vaxrank help ...``,
      ``vaxrank policies ...``, ``vaxrank replay ...``, ``vaxrank run ...``)
      complete activity names and only build the selected activity's parser.
    """
    words, current_prefix = _split_completion_line(environ)
    tail = words[1:]

    # Option-first mode: `vaxrank --...` or `vaxrank --vcf foo.vcf --...`
    if (not tail and current_prefix.startswith("-")) or (
        tail and tail[0].startswith("-")
    ):
        from .arg_parser import choose_arg_parser

        return choose_arg_parser(tail)

    parser = HelpArgumentParser(prog="vaxrank")
    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND")
    selected = tail[0] if tail and tail[0] in ACTIVITIES else None

    for name in ACTIVITIES:
        if name != selected:
            subparsers.add_parser(name, add_help=False)
            continue

        if name == "data":
            from .data_cli import _build_parser as build_data_parser

            _attach_subparser(subparsers, name, build_data_parser())
        elif name == "compare":
            from .compare_cli import _build_parser as build_compare_parser

            _attach_subparser(subparsers, name, build_compare_parser())
        elif name == "mhc-status":
            from .mhc_status_cli import _build_parser as build_mhc_status_parser

            _attach_subparser(subparsers, name, build_mhc_status_parser())
        elif name == "run":
            from .arg_parser import choose_arg_parser

            _attach_subparser(subparsers, name, choose_arg_parser(tail[1:]))
        elif name == "help":
            p_help = subparsers.add_parser(name, add_help=False)
            p_help.add_argument(
                "topic",
                nargs="?",
                choices=tuple(HELP_TOPICS) + ("all",),
            )
        elif name == "policies":
            p_pol = subparsers.add_parser(name, add_help=False)
            pol_sub = p_pol.add_subparsers(dest="policies_action")
            p_show = pol_sub.add_parser("show", add_help=False)
            p_show.add_argument("policy_name", choices=_bundled_policy_names())
        elif name == "replay":
            from .arg_parser import cached_run_arg_parser, external_input_arg_parser

            base = (
                external_input_arg_parser()
                if len(tail) > 1 and not tail[1].lower().endswith(".json")
                else cached_run_arg_parser()
            )
            p_replay = subparsers.add_parser(
                name, add_help=False, parents=[base]
            )
            p_replay.add_argument("saved_run")
        else:
            subparsers.add_parser(name, add_help=False)

    return parser


def autocomplete_vaxrank(environ=None):
    """Run ``argcomplete.autocomplete`` on the vaxrank completion parser."""
    argcomplete.autocomplete(build_completion_parser(environ))
