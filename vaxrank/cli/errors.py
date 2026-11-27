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

"""Which failures are user mistakes, and how the CLI reports them.

argparse already gets this right for the choices it validates: one line on
stderr, a non-zero exit, no traceback. Everything vaxrank validates itself
reached the user as a traceback instead — a missing ``--vcf``, an unparseable
allele, a typo'd ``--config`` key, and even vaxrank's own carefully worded
"No output path specified", which was delivered as the last frame of a stack
(openvax/vaxrank#457).

Nothing is hidden by reporting one line. ``logging.conf`` gives the ``vaxrank``
logger a DEBUG file handler and an INFO console handler, so the full traceback
still reaches ``--log-path`` without reaching the terminal, and ``--verbose``
prints it.
"""

import logging
import sys
import traceback


logger = logging.getLogger(__name__)

_OPTIONAL_USER_ERROR_SPECS = (
    ("mhcgnomes.errors", "ParseError"),
    ("msgspec", "DecodeError"),
    ("pyensembl.download_cache", "MissingLocalFile"),
)


def user_error_types():
    """Return user-error exception classes from already-loaded modules.

    Looking up optional exception classes from ``sys.modules`` keeps CLI startup
    and tab completion import-light while still catching every domain error that
    can actually be raised at runtime.
    """
    types = [OSError, ValueError]
    for module_name, attr_name in _OPTIONAL_USER_ERROR_SPECS:
        module = sys.modules.get(module_name)
        if module is not None:
            err_type = getattr(module, attr_name, None)
            if isinstance(err_type, type) and issubclass(err_type, Exception):
                types.append(err_type)
    return tuple(types)


def __getattr__(name):
    """Lazily resolve ``USER_ERRORS`` for callers inspecting the full tuple."""
    if name == "USER_ERRORS":
        from mhcgnomes.errors import ParseError as AlleleParseError
        from msgspec import DecodeError
        from pyensembl.download_cache import MissingLocalFile

        value = (
            OSError,
            ValueError,
            AlleleParseError,
            DecodeError,
            MissingLocalFile,
        )
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def _verbose(args_list):
    """Whether the run asked for the traceback.

    Read from the argument list rather than the parsed namespace: a failure
    during parsing or during input loading has no namespace to consult, and
    this has to behave the same either way.
    """
    argv = sys.argv[1:] if args_list is None else args_list
    return "--verbose" in argv or "-v" in argv


def exit_with_user_error(error, args_list=None, prog="vaxrank", traceback_hint=True):
    """Print one line on stderr, keep the traceback available, exit 1.

    ``prog`` names the command in the message. The companion figure commands
    have no ``--verbose`` or ``--log-path``, so they pass ``traceback_hint=False``
    rather than point at options that do not exist.
    """
    logger.debug("Exiting on user error: %s", error, exc_info=True)
    message = str(error).strip() or error.__class__.__name__
    print("%s: error: %s" % (prog, message), file=sys.stderr)
    if _verbose(args_list):
        traceback.print_exc()
    elif traceback_hint:
        print("%s: re-run with --verbose, or read the file named by "
              "--log-path, for the full traceback." % prog, file=sys.stderr)
    raise SystemExit(1)


def run_reporting_user_errors(prog, function, *args, **kwargs):
    """Run a companion command's work, reporting an input mistake as one line.

    Gives companion commands the manners ``vaxrank`` already has: a missing or
    malformed input file is a message and exit status 1, not a traceback.
    """
    try:
        return function(*args, **kwargs)
    except user_error_types() as error:
        # Companion commands expose no verbosity option. Do not inherit a
        # caller's argv (for example pytest -v or an embedding application).
        exit_with_user_error(error, args_list=(), prog=prog, traceback_hint=False)
