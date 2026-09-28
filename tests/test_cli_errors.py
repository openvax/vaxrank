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

"""Ordinary mistakes get one line, not a traceback (openvax/vaxrank#457).

argparse already behaved this way for the choices it validates. These pin the
same manners for everything vaxrank validates itself, and pin that the
traceback is still recoverable rather than discarded.
"""

import sys

import pytest

from vaxrank.cli import entry_point
from vaxrank.cli.errors import USER_ERRORS

from .testing_helpers import data_path


B16 = "b16.f10"


def _pipeline_args(tmp_path, *extra):
    return [
        "--vcf", data_path("%s/b16.vcf" % B16),
        "--bam", data_path("%s/b16.combined.sorted.bam" % B16),
        "--mhc-predictor", "random",
        "--mhc-alleles", "H2-Kb",
        "--mhc-peptide-lengths", "9",
        "--log-path", str(tmp_path / "python.log"),
    ] + list(extra)


# One representative per class of mistake in #457's table: a missing input
# file, an unparseable allele, a bad --config key, and vaxrank's own
# validation message.
MISTAKES = {
    "missing_vcf": lambda tmp: [
        "--vcf", str(tmp / "nope.vcf"),
        "--bam", data_path("%s/b16.combined.sorted.bam" % B16),
        "--mhc-predictor", "random", "--mhc-alleles", "H2-Kb",
        "--output-csv", str(tmp / "out.csv"),
        "--log-path", str(tmp / "python.log")],
    "unparseable_allele": lambda tmp: _pipeline_args(
        tmp, "--mhc-alleles", "NOTANALLELE", "--output-csv", str(tmp / "out.csv")),
    "unknown_config_key": lambda tmp: _pipeline_args(
        tmp, "--config-value", "nosuchkey=1", "--output-csv", str(tmp / "out.csv")),
    "no_output_flag": lambda tmp: _pipeline_args(tmp),
    "missing_lens_report": lambda tmp: [
        "--input-lens", str(tmp / "nope.tsv"),
        "--output-csv", str(tmp / "out.csv"),
        "--log-path", str(tmp / "python.log")],
}


@pytest.mark.parametrize("name", sorted(MISTAKES))
def test_user_mistakes_print_one_line_and_exit_nonzero(name, tmp_path, capsys):
    """No traceback, a non-zero exit, and a message an operator can act on."""
    with pytest.raises(SystemExit) as exit_info:
        entry_point.main(MISTAKES[name](tmp_path))
    assert exit_info.value.code == 1

    captured = capsys.readouterr()
    assert "vaxrank: error: " in captured.err, captured.err
    assert "Traceback (most recent call last)" not in captured.err
    assert 'File "' not in captured.err
    # The message itself has to carry content; an empty str(error) would
    # otherwise satisfy every assertion above.
    message = captured.err.split("vaxrank: error: ", 1)[1].splitlines()[0]
    assert len(message.strip()) > 10, message


def test_vaxrank_own_validation_message_survives_intact(tmp_path, capsys):
    """The wording check_args chose is the wording the user reads.

    This message was written for an operator and then delivered as the last
    frame of a stack.
    """
    with pytest.raises(SystemExit):
        entry_point.main(_pipeline_args(tmp_path))
    captured = capsys.readouterr()
    assert "vaxrank: error: No output path specified" in captured.err


def test_verbose_still_shows_the_traceback(tmp_path, capsys):
    """--verbose is the documented escape hatch, so it must actually work."""
    with pytest.raises(SystemExit):
        entry_point.main(MISTAKES["missing_lens_report"](tmp_path) + ["--verbose"])
    captured = capsys.readouterr()
    assert "vaxrank: error: " in captured.err
    assert "Traceback (most recent call last)" in captured.err


def test_log_path_records_the_traceback_when_stderr_does_not(tmp_path, capsys):
    """Nothing is lost by printing one line: the log file keeps the stack.

    logging.conf gives the vaxrank logger a DEBUG file handler and an INFO
    console handler, which is what lets the traceback reach the file without
    reaching the terminal.
    """
    log_path = tmp_path / "python.log"
    with pytest.raises(SystemExit):
        entry_point.main(MISTAKES["missing_lens_report"](tmp_path))
    assert "Traceback (most recent call last)" not in capsys.readouterr().err
    logged = log_path.read_text()
    assert "Exiting on user error" in logged
    assert "Traceback (most recent call last)" in logged


def test_programming_errors_still_propagate(tmp_path, monkeypatch):
    """A bug in vaxrank is ours to fix, so it keeps its traceback.

    The boundary catches OSError/ValueError, which a defect can also raise —
    that is the deliberate cost. What it must never do is swallow the error
    classes that only a defect produces.
    """
    def broken(args_list=None):
        raise KeyError("internal invariant")

    monkeypatch.setattr(entry_point, "run_cli", broken)
    with pytest.raises(KeyError, match="internal invariant"):
        entry_point.main([])

    for never_caught in (KeyError, AttributeError, TypeError, NameError,
                         IndexError, RuntimeError):
        assert not issubclass(never_caught, USER_ERRORS)


def test_argparse_validated_choices_are_unchanged(tmp_path, capsys):
    """argparse's own errors were already right; leave them exactly alone."""
    with pytest.raises(SystemExit) as exit_info:
        entry_point.main(_pipeline_args(tmp_path, "--mhc-predictor", "notapredictor"))
    assert exit_info.value.code == 2
    assert "Traceback (most recent call last)" not in capsys.readouterr().err


def test_weasyprint_without_pango_names_the_fix(monkeypatch):
    """A dlopen failure inside cffi is not an actionable message.

    test.sh and AGENTS.md both know the macOS incantation; the CLI should say
    it rather than print the loader's stack.
    """
    from vaxrank import report

    # None in sys.modules makes `import weasyprint` raise ImportError, which is
    # the same shape as the dlopen OSError this guards.
    monkeypatch.setitem(sys.modules, "weasyprint", None)
    with pytest.raises(ValueError) as error:
        report._pdf_via_weasyprint("ignored.html", "ignored.pdf")
    message = str(error.value)
    assert "DYLD_FALLBACK_LIBRARY_PATH" in message
    assert "--pdf-backend pdfkit" in message
    assert isinstance(error.value, USER_ERRORS)


def test_empty_report_warns_and_names_the_unwritten_path(tmp_path, caplog):
    """An empty input used to exit 0 having written nothing that was asked for."""
    from vaxrank import report

    csv_path = tmp_path / "vaccine-peptides.csv"
    with caplog.at_level("WARNING", logger="vaxrank.report"):
        report.make_csv_report([], csv_report_path=str(csv_path))
    assert not csv_path.exists()
    warnings = [r.getMessage() for r in caplog.records if r.levelname == "WARNING"]
    assert any(str(csv_path) in message for message in warnings), warnings
