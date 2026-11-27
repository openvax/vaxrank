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

"""Tests for CLI argcomplete tab completion, lazy startup, and help styling (#794)."""

import importlib
import io
import os
from pathlib import Path
import re
import subprocess
import sys
import tomllib

import pytest

from vaxrank.cli.arg_parser import (
    cached_run_arg_parser,
    external_input_arg_parser,
    make_vaxrank_arg_parser,
)
from vaxrank.cli.compare_cli import _build_parser as build_compare_parser
from vaxrank.cli.data_cli import _build_parser as build_data_parser
from vaxrank.cli.discoverability import (
    ACTIVITIES,
    HELP_TOPICS,
    print_short_help,
    topic_help,
)
from vaxrank.cli.help import (
    HelpArgumentParser,
    color_enabled,
    style_help,
    style_table_header,
)
from vaxrank.cli.mhc_status_cli import _build_parser as build_mhc_status_parser
from vaxrank.testdata.download import make_arg_parser as build_testdata_parser

REPO_ROOT = Path(__file__).resolve().parents[1]

# Heavy runtime libraries that must never load during CLI import or parser construction.
HEAVY_RUNTIME_MODULES = (
    "pandas",
    "numpy",
    "torch",
    "pysam",
    "matplotlib",
    "pyarrow",
    "weasyprint",
    "topiary",
    "pyensembl",
    "oncoref",
    "orfdata",
    "osteosarc",
    "hitlist",
    "dnachisel",
)


def _subprocess_env(**extra):
    env = os.environ.copy()
    existing = env.get("PYTHONPATH", "")
    parts = [str(REPO_ROOT)]
    if existing:
        parts.append(existing)
    env["PYTHONPATH"] = os.pathsep.join(parts)
    env.update(extra)
    return env


def _run_argcomplete(tmp_path, entry_expr, comp_line, comp_point=None):
    """Execute a console script entry point under argcomplete's subprocess protocol."""
    out_file = tmp_path / "argcomplete_out.txt"
    point = len(comp_line) if comp_point is None else comp_point
    env = _subprocess_env(
        _ARGCOMPLETE="1",
        _ARGCOMPLETE_IFS="\n",
        _ARGCOMPLETE_STDOUT_FILENAME=str(out_file),
        COMP_LINE=comp_line,
        COMP_POINT=str(point),
    )
    proc = subprocess.run(
        [sys.executable, "-c", entry_expr],
        cwd=str(REPO_ROOT),
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, (
        f"argcomplete failed (rc={proc.returncode}) for {comp_line!r}:\n"
        f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
    )
    assert out_file.exists(), f"argcomplete produced no output for {comp_line!r}"
    return [
        line.strip()
        for line in out_file.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]


def test_every_console_script_entry_module_has_argcomplete_marker():
    pyproject = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    scripts = pyproject["project"]["scripts"]
    assert scripts, "Expected at least one console script in pyproject.toml"
    for script_name, target in scripts.items():
        module_name = target.split(":", 1)[0]
        rel_path = Path(*module_name.split(".")).with_suffix(".py")
        head = (REPO_ROOT / rel_path).read_bytes()[:1024]
        assert b"PYTHON_ARGCOMPLETE_OK" in head, (
            f"Console script {script_name!r} entry module {rel_path} is missing "
            "PYTHON_ARGCOMPLETE_OK within its first 1024 bytes"
        )


def test_requirements_declare_argcomplete_and_datacache_cli_help_floor():
    lines = [
        line.split("#", 1)[0].strip()
        for line in (REPO_ROOT / "requirements.txt").read_text(encoding="utf-8").splitlines()
        if line.split("#", 1)[0].strip()
    ]
    assert any(re.match(r"^argcomplete\s*>=\s*3\b", line) for line in lines), (
        f"Expected argcomplete>=3 in requirements.txt, got: {lines}"
    )
    assert any(re.match(r"^datacache\s*>=\s*1\.3[89]\b|^datacache\s*>=\s*1\.[4-9]\d\b", line) for line in lines), (
        f"Expected datacache>=1.38.0 in requirements.txt, got: {lines}"
    )


@pytest.mark.parametrize(
    "snippet",
    [
        "import vaxrank",
        "import argcomplete, vaxrank.cli.entry_point",
        (
            "import argcomplete, vaxrank.cli.arg_parser as ap; "
            "ap.make_vaxrank_arg_parser(); "
            "ap.external_input_arg_parser(); "
            "ap.cached_run_arg_parser()"
        ),
        "import argcomplete, vaxrank.cli.data_cli as dc; dc._build_parser()",
        "import argcomplete, vaxrank.cli.compare_cli as cc; cc._build_parser()",
        "import argcomplete, vaxrank.cli.mhc_status_cli as mc; mc._build_parser()",
        "import argcomplete, vaxrank.testdata.download as dl; dl.make_arg_parser()",
    ],
)
def test_cli_import_and_parser_construction_avoid_heavy_modules(snippet):
    code = (
        "import sys\n"
        f"{snippet}\n"
        f"loaded = [m for m in {HEAVY_RUNTIME_MODULES!r} if m in sys.modules]\n"
        "assert not loaded, f'Heavy modules loaded: {loaded}'\n"
    )
    proc = subprocess.run(
        [sys.executable, "-c", code],
        cwd=str(REPO_ROOT),
        env=_subprocess_env(),
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr


def test_top_level_package_lazy_exports_preserve_public_api():
    import vaxrank

    assert isinstance(vaxrank.__version__, str) and vaxrank.__version__
    assert "run_vaxrank" in vaxrank.__all__
    assert "VaccineConfig" in vaxrank.__all__
    assert "EpitopeConfig" in vaxrank.__all__
    assert set(vaxrank.__all__) <= set(dir(vaxrank))
    # Resolving a lazy export succeeds and caches it on the module.
    assert vaxrank.STANDARD_AMINO_ACIDS
    with pytest.raises(AttributeError):
        getattr(vaxrank, "definitely_not_a_real_vaxrank_export")


@pytest.mark.parametrize(
    "comp_line,expected_items",
    [
        ("vaxrank ", set(ACTIVITIES)),
        ("vaxrank da", {"data"}),
        ("vaxrank --", {"--vcf", "--bam", "--mhc-predictor", "--output-dir", "--vaccine-type", "--input-lens", "--version"}),
        ("vaxrank --vaccine-type ", {"peptide", "mrna"}),
        ("vaxrank --scoring-mode ", {"affinity", "percentile_rank"}),
        ("vaxrank run --", {"--vcf", "--bam", "--mhc-predictor", "--output-dir"}),
        ("vaxrank data ", {"dirs", "status", "path", "prune"}),
        ("vaxrank data prune --", {"--keep", "--kmer-indexes", "--yes", "--json"}),
        ("vaxrank compare --", {"--output-dir", "--label-a", "--label-b", "--top-k", "--pdf", "--json"}),
        ("vaxrank mhc-status --", {"--output", "--config", "--tumor-purity", "--from-t1k", "--allele-fasta", "--scrna-counts"}),
        ("vaxrank help ", set(HELP_TOPICS) | {"all"}),
        ("vaxrank policies ", {"show"}),
        ("vaxrank policies show ", {"openvax-v1", "lens-v1"}),
        ("vaxrank replay run.json --", {"--output-dir", "--output-csv", "--output-ascii-report"}),
    ],
)
def test_vaxrank_argcomplete_completions(tmp_path, comp_line, expected_items):
    completions = set(
        _run_argcomplete(
            tmp_path,
            "from vaxrank.cli.entry_point import main; main()",
            comp_line,
        )
    )
    assert expected_items <= completions, (
        f"Expected {expected_items} in completions {completions} for {comp_line!r}"
    )


def test_vaxrank_test_data_argcomplete_completions(tmp_path):
    completions = set(
        _run_argcomplete(
            tmp_path,
            "from vaxrank.testdata.download import main; main()",
            "vaxrank-test-data --",
        )
    )
    expected = {
        "--manifest",
        "--output",
        "--cache-root",
        "--offline",
        "--repair-cache",
        "--verify-only",
        "--progress",
        "--timeout",
        "--max-retries",
        "--version",
        "--help",
    }
    assert expected <= completions, (
        f"Expected {expected} in completions {completions}"
    )


def test_all_cli_parsers_use_help_argument_parser():
    parsers = [
        make_vaxrank_arg_parser(),
        external_input_arg_parser(),
        cached_run_arg_parser(),
        build_data_parser(),
        build_compare_parser(),
        build_mhc_status_parser(),
        build_testdata_parser(),
    ]
    for parser in parsers:
        assert isinstance(parser, HelpArgumentParser), (
            f"Parser {parser.prog!r} is {type(parser).__name__}, expected HelpArgumentParser"
        )
        plain = parser.format_help()
        assert "\033[" not in plain
        assert "``" not in plain


def test_help_styling_plain_in_format_help_and_respects_tty_and_no_color(monkeypatch):
    class _TTYBuffer(io.StringIO):
        def isatty(self):
            return True

    parser = build_testdata_parser()
    plain = parser.format_help()
    assert "\033[" not in plain
    assert "``" not in plain
    assert "\033[" not in topic_help("variants")

    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "xterm-256color")

    non_tty = io.StringIO()
    assert not color_enabled(non_tty)
    parser.print_help(file=non_tty)
    assert "\033[" not in non_tty.getvalue()

    short_non_tty = io.StringIO()
    print_short_help(stream=short_non_tty)
    assert "\033[" not in short_non_tty.getvalue()

    tty = _TTYBuffer()
    assert color_enabled(tty)
    parser.print_help(file=tty)
    styled = tty.getvalue()
    assert "\033[1musage:\033[0m" in styled
    assert "\033[36m--offline\033[0m" in styled
    assert "\033[1;36m" in styled
    assert style_table_header("name  size\nfoo   10\n", tty).startswith("\033[1;36m")

    short_tty = _TTYBuffer()
    print_short_help(stream=short_tty)
    short_styled = short_tty.getvalue()
    assert "\033[1musage:\033[0m" in short_styled
    assert "\033[1;36mExamples:\033[0m" in short_styled
    assert "\033[36m--vcf\033[0m" in short_styled

    monkeypatch.setenv("NO_COLOR", "1")
    assert not color_enabled(tty)
    assert "\033[" not in style_help(plain, tty)

    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "dumb")
    assert not color_enabled(tty)
    assert "\033[" not in style_help(plain, tty)


def test_console_script_entry_targets_are_importable():
    pyproject = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    for target in pyproject["project"]["scripts"].values():
        module_name, attr_name = target.split(":", 1)
        module = importlib.import_module(module_name)
        assert callable(getattr(module, attr_name))
