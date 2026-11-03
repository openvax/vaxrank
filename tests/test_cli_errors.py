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


def test_empty_ranked_table_writes_a_schema(tmp_path):
    import pandas as pd
    from vaxrank import report
    path = tmp_path / 'candidates.csv'
    report.make_csv_report([], csv_report_path=path)
    frame = pd.read_csv(path)
    assert frame.empty and list(frame) == list(report.new_columns())


def test_argparse_errors_name_the_problem_without_the_usage_dump(capsys):
    """vaxrank's usage is over 150 lines; the one line that matters must not hide in it."""
    for argv in ([], ["--bam", "x.bam", "--mhc-predictor", "notapredictor"],
                 ["--input-pvacseq", "x.tsv", "--mhc-predictor", "notapredictor"]):
        with pytest.raises(SystemExit) as exit_info:
            entry_point.main(argv)
        assert exit_info.value.code == 2
        err = capsys.readouterr().err
        assert err.startswith("vaxrank: error: "), err
        assert "usage:" not in err
        assert "vaxrank --help" in err
        assert len(err.splitlines()) <= 3


def test_unparseable_allele_names_the_flag_and_the_entry_as_typed(tmp_path, capsys, monkeypatch):
    """mhctools drops the '*' from names it cannot parse, so a later check would
    quote 'HLA-Z99:99' for input the user typed as 'HLA-Z*99:99'."""
    def must_not_run(args):
        raise AssertionError("input loading started before the alleles were checked")

    monkeypatch.setattr(entry_point, "variant_collection_from_args", must_not_run)
    with pytest.raises(SystemExit) as exit_info:
        entry_point.main(_pipeline_args(
            tmp_path, "--mhc-alleles", "H2-Kb,HLA-Z*99:99", "--output-csv", str(tmp_path / "o.csv")))
    assert exit_info.value.code == 1
    assert "--mhc-alleles: cannot parse 'HLA-Z*99:99' as an MHC allele" in capsys.readouterr().err


def test_unwritable_output_destinations_fail_before_the_pipeline(tmp_path, capsys, monkeypatch):
    """A bad destination used to surface only after the whole analysis had run."""
    def must_not_run(args):
        raise AssertionError("pipeline started before the destination was checked")

    monkeypatch.setattr(entry_point, "variant_collection_from_args", must_not_run)
    blocker = tmp_path / "a_file"
    blocker.write_text("not a directory")
    for flag in ("--output-dir", "--output-csv"):
        target = blocker / "sub" if flag == "--output-dir" else blocker / "x.csv"
        with pytest.raises(SystemExit) as exit_info:
            entry_point.main(_pipeline_args(tmp_path, flag, str(target)))
        assert exit_info.value.code == 1
        assert "vaxrank: error: %s: " % flag in capsys.readouterr().err
    assert not (tmp_path / "newdir").exists()


def test_a_new_nested_output_dir_is_still_accepted_and_not_created_early(tmp_path):
    args = type("Args", (), {"output_dir": str(tmp_path / "new" / "nested")})()
    entry_point.populate_default_output_paths(args)
    entry_point.validate_output_destinations(args)
    assert not (tmp_path / "new").exists()
    assert args.output_csv.endswith("ranked_candidates.csv")


def test_figure_commands_report_input_mistakes_in_one_line(tmp_path, capsys):
    """The companion commands leaked a traceback for a missing or wrong-shaped file."""
    from vaxrank.cli import (
        complex_variant_visualization, evidence_visualization, mutation_visualization)

    wrong_shape = tmp_path / "wrong.json"
    wrong_shape.write_text('{"unrelated": 1}')
    cases = [
        ("vaxrank-mutation-figure", mutation_visualization.main, [str(tmp_path / "nope.csv")], "nope.csv"),
        ("vaxrank-evidence-figure", evidence_visualization.main, [str(tmp_path / "nope.json")], "nope.json"),
        ("vaxrank-complex-result-figure", complex_variant_visualization.main,
         [str(wrong_shape), "--output-root", str(tmp_path / "out")], "'records' and 'metadata'"),
    ]
    for prog, main, argv, expected in cases:
        with pytest.raises(SystemExit) as exit_info:
            main(argv)
        assert exit_info.value.code == 1
        err = capsys.readouterr().err
        assert err.startswith("%s: error: " % prog), err
        assert expected in err
        assert "Traceback" not in err and "--verbose" not in err
        assert len(err.splitlines()) == 1


def test_every_bundled_policy_has_its_own_description():
    from vaxrank.config.policies import list_policies

    # openvax-v1 is a byte-hash-frozen baseline, so its YAML cannot carry one.
    descriptions = [policy["description"] for policy in list_policies()
                    if policy["name"] != "openvax-v1"]
    assert not any("see --show-policy" in text for text in descriptions)
    assert len(set(descriptions)) == len(descriptions)


def test_config_errors_say_which_inputs_built_the_configuration(tmp_path):
    from vaxrank.config import load_vaxrank_config

    with pytest.raises(FileNotFoundError, match="config file not found: .*nope.yaml"):
        load_vaxrank_config(config_path=str(tmp_path / "nope.yaml"))
    # msgspec.ValidationError is not a ValueError on every supported msgspec release.
    with pytest.raises(USER_ERRORS, match=r"unknown field `nosuchkey` \(configuration built from: nosuchkey=1\)"):
        load_vaxrank_config(set_overrides=["nosuchkey=1"])
    # Overrides given as an ordered list are reported once, not once per list.
    with pytest.raises(USER_ERRORS) as error:
        load_vaxrank_config(ordered_overrides=[("set", "nosuchkey=1")], set_overrides=["nosuchkey=1"])
    assert str(error.value).count("nosuchkey=1") == 1


def test_fragments_without_epitopes_are_counted_not_logged_one_by_one():
    import inspect
    from vaxrank import core_logic

    source = inspect.getsource(core_logic)
    assert 'logger.info(\n                "No epitopes for mutant protein fragment' not in source
    assert "n_fragments_without_epitopes" in source


def _funnel(loaded, coding, rna, peptides):
    from types import SimpleNamespace
    return SimpleNamespace(
        num_somatic_variants=loaded, num_coding_effect_variants=coding,
        num_variants_with_rna_support=rna, num_variants_with_vaccine_peptides=peptides,
        antigen_summary=None, cta_expression_summary=None)


@pytest.mark.parametrize("funnel, expected", [
    ((0, 0, 0, 0), "no variants were loaded"),
    ((1, 0, 0, 0), "none has a coding effect"),
    ((1, 1, 0, 0), "none has RNA support"),
    ((1, 1, 1, 0), "none survived the RNA-evidence and peptide filters"),
])
def test_an_empty_result_is_warned_about_with_the_stage_that_emptied_it(funnel, expected, caplog):
    import logging
    with caplog.at_level(logging.WARNING, logger=entry_point.logger.name):
        entry_point.warn_if_no_vaccine_peptides(_funnel(*funnel))
    messages = [record.getMessage() for record in caplog.records if record.levelno == logging.WARNING]
    assert len(messages) == 1 and "No vaccine peptides were produced" in messages[0]
    assert expected in messages[0]


def test_a_run_with_peptides_or_without_funnel_counts_does_not_warn(caplog):
    import logging
    with caplog.at_level(logging.WARNING, logger=entry_point.logger.name):
        entry_point.warn_if_no_vaccine_peptides(_funnel(3, 3, 2, 1))
        entry_point.warn_if_no_vaccine_peptides(None)
    assert not caplog.records


def _without_pepsickle(monkeypatch):
    import importlib.util
    real = importlib.util.find_spec
    monkeypatch.setattr(importlib.util, "find_spec",
                        lambda name, *a, **k: None if name == "pepsickle" else real(name, *a, **k))


def test_missing_pepsickle_is_one_info_line_not_a_warning_by_default(monkeypatch, caplog):
    import logging
    _without_pepsickle(monkeypatch)
    entry_point.require_pepsickle_if_requested(["--mhc-predictor", "random"])  # not requested: no error
    with caplog.at_level(logging.INFO, logger=entry_point.logger.name):
        assert entry_point.processing_annotation_available() is False
    assert [r.levelno for r in caplog.records] == [logging.INFO]
    assert "mhctools[pepsickle]" in caplog.records[0].getMessage()


def test_asking_for_pepsickle_by_name_fails_when_it_is_missing(monkeypatch):
    _without_pepsickle(monkeypatch)
    with pytest.raises(ValueError, match=r"--processing-aware-annotation: pepsickle is not installed"):
        entry_point.require_pepsickle_if_requested(["--processing-aware-annotation"])


def test_installed_pepsickle_runs_the_annotation_silently(monkeypatch, caplog):
    import importlib.util
    monkeypatch.setattr(importlib.util, "find_spec", lambda name, *a, **k: object())
    entry_point.require_pepsickle_if_requested(["--processing-aware-annotation"])
    assert entry_point.processing_annotation_available() is True
    assert not caplog.records


def test_explicit_pepsickle_request_fails_before_the_pipeline(tmp_path, monkeypatch, capsys):
    _without_pepsickle(monkeypatch)

    def must_not_run(args):
        raise AssertionError("pipeline started before pepsickle availability was checked")

    monkeypatch.setattr(entry_point, "variant_collection_from_args", must_not_run)
    with pytest.raises(SystemExit) as exit_info:
        entry_point.main(_pipeline_args(
            tmp_path, "--processing-aware-annotation", "--output-csv", str(tmp_path / "o.csv")))
    assert exit_info.value.code == 1
    assert "--processing-aware-annotation: pepsickle is not installed" in capsys.readouterr().err


def test_pepsickle_in_a_configured_interpreter_counts_as_installed(monkeypatch):
    _without_pepsickle(monkeypatch)
    monkeypatch.delenv("PEPSICKLE_PYTHON", raising=False)
    monkeypatch.delenv("PEPSICKLE_GB_PYTHON", raising=False)
    assert entry_point.pepsickle_installed() is False
    monkeypatch.setenv("PEPSICKLE_PYTHON", "/opt/pepsickle-venv/bin/python")
    assert entry_point.pepsickle_installed() is True


def test_every_companion_command_answers_version_like_vaxrank(capsys):
    from vaxrank.testdata import download
    from vaxrank.cli import (
        complex_variant_visualization, evidence_visualization, mutation_visualization, window_exposure)
    from vaxrank.version import __version__

    for prog, main in [("vaxrank-mutation-figure", mutation_visualization.main),
                       ("vaxrank-evidence-figure", evidence_visualization.main),
                       ("vaxrank-complex-result-figure", complex_variant_visualization.main),
                       ("vaxrank-window-exposure", window_exposure.main),
                       ("vaxrank-test-data", download.main)]:
        with pytest.raises(SystemExit) as exit_info:
            main(["--version"])
        assert exit_info.value.code == 0
        assert capsys.readouterr().out.strip() == "%s %s" % (prog, __version__)


def test_companion_options_all_have_help():
    from vaxrank.cli import complex_variant_visualization, evidence_visualization, mutation_visualization

    for module in (mutation_visualization, evidence_visualization, complex_variant_visualization):
        for action in module.make_parser()._actions:
            assert action.help, "%s %s has no help" % (module.__name__, action.option_strings or action.dest)



def test_missing_input_paths_name_the_flag_before_any_work(tmp_path, monkeypatch, capsys):
    def must_not_run(args):
        raise AssertionError("input loading started before the paths were checked")

    monkeypatch.setattr(entry_point, "variant_collection_from_args", must_not_run)
    out = ["--output-csv", str(tmp_path / "o.csv")]
    cases = [
        ("--bam", str(tmp_path / "nope.bam"), "--bam: no such file"),
        ("--prediction-cache", str(tmp_path), "--prediction-cache: %r is a directory, expected a file" % str(tmp_path)),
        ("--mhc-alleles-file", str(tmp_path / "alleles.txt"), "--mhc-alleles-file: no such file"),
    ]
    for flag, value, expected in cases:
        with pytest.raises(SystemExit) as exit_info:
            entry_point.main(_pipeline_args(tmp_path, flag, value, *out))
        assert exit_info.value.code == 1
        assert expected in capsys.readouterr().err, flag


def test_remote_and_existing_inputs_pass_the_path_check(tmp_path):
    from types import SimpleNamespace
    present = tmp_path / "x.vcf"
    present.write_text("")
    entry_point.validate_input_paths(SimpleNamespace(
        vcf=[str(present)], bam="https://example.org/sample.bam", hitlist_evidence_bundle=str(tmp_path)))


def test_a_missing_directory_input_names_the_flag(tmp_path):
    from types import SimpleNamespace
    with pytest.raises(ValueError, match="--hitlist-evidence-bundle: no such directory"):
        entry_point.validate_input_paths(SimpleNamespace(hitlist_evidence_bundle=str(tmp_path / "nope")))
