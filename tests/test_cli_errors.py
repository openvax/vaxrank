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
from pathlib import Path

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
    for argv in (["--bam", "x.bam", "--mhc-predictor", "notapredictor"],
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

    monkeypatch.setattr("varcode.cli.variant_collection_from_args", must_not_run)
    with pytest.raises(SystemExit) as exit_info:
        entry_point.main(_pipeline_args(
            tmp_path, "--mhc-alleles", "H2-Kb,HLA-Z*99:99", "--output-csv", str(tmp_path / "o.csv")))
    assert exit_info.value.code == 1
    assert "--mhc-alleles: cannot parse 'HLA-Z*99:99' as an MHC allele" in capsys.readouterr().err


def test_unwritable_output_destinations_fail_before_the_pipeline(tmp_path, capsys, monkeypatch):
    """A bad destination used to surface only after the whole analysis had run."""
    def must_not_run(args):
        raise AssertionError("pipeline started before the destination was checked")

    monkeypatch.setattr("varcode.cli.variant_collection_from_args", must_not_run)
    blocker = tmp_path / "a_file"
    blocker.write_text("not a directory")
    for flag in ("--output-dir", "--output-policy-dir", "--output-csv"):
        target = blocker / "sub" if flag in ("--output-dir", "--output-policy-dir") else blocker / "x.csv"
        with pytest.raises(SystemExit) as exit_info:
            entry_point.main(_pipeline_args(tmp_path, flag, str(target)))
        assert exit_info.value.code == 1
        assert "vaxrank: error: %s: " % flag in capsys.readouterr().err
    assert not (tmp_path / "newdir").exists()


def test_every_primary_output_flag_is_validated_by_validate_output_destinations(tmp_path):
    """Every _PRIMARY_OUTPUT_FLAGS entry must be checked fail-fast (#781)."""
    from vaxrank.cli.arg_parser import _PRIMARY_OUTPUT_FLAGS

    blocker = tmp_path / "blocker_file"
    blocker.write_text("not a directory")
    for flag, attr, _ in _PRIMARY_OUTPUT_FLAGS:
        is_dir = attr in entry_point._OUTPUT_DIRECTORY_FLAGS
        target = str(blocker / "subdir") if is_dir else str(blocker / "out.dat")
        args = type("Args", (), {attr: target})()
        with pytest.raises(ValueError, match=flag):
            entry_point.validate_output_destinations(args)


def test_cta_admission_default_is_derived_before_destination_validation(tmp_path, monkeypatch):
    """Derived --output-cta-admission default is populated before validation (#781)."""
    seen = {}

    def capture_validation(args):
        seen["output_cta_admission"] = getattr(args, "output_cta_admission", None)
        raise ValueError("stop-after-validation")

    monkeypatch.setattr(entry_point, "validate_output_destinations", capture_validation)
    cta_csv = tmp_path / "cta.csv"
    cta_csv.write_text("gene,tpm\nMAGEA3,10.0\n")
    out_dir = tmp_path / "out"
    with pytest.raises(ValueError, match="stop-after-validation"):
        entry_point.run_cli([
            "--input-cta-expression", str(cta_csv),
            "--mhc-predictor", "random",
            "--mhc-alleles", "HLA-A*02:01",
            "--output-dir", str(out_dir),
            "--log-path", str(tmp_path / "cta.log"),
        ])
    assert seen["output_cta_admission"] == str(out_dir / "cta_admission.json")


def test_a_new_nested_output_dir_is_still_accepted_and_not_created_early(tmp_path):
    args = type("Args", (), {"output_dir": str(tmp_path / "new" / "nested")})()
    entry_point.populate_default_output_paths(args)
    entry_point.validate_output_destinations(args)
    assert not (tmp_path / "new").exists()
    assert args.output_csv.endswith("ranked_candidates.csv")


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


def test_processing_aware_prefix_and_override_respect_parsed_value(tmp_path, monkeypatch, capsys):
    """#783: --processing-aware prefix fails fast; trailing --no-processing-aware-annotation wins."""
    _without_pepsickle(monkeypatch)

    def must_not_run(args):
        raise AssertionError("pipeline started before pepsickle availability was checked")

    monkeypatch.setattr("varcode.cli.variant_collection_from_args", must_not_run)
    with pytest.raises(SystemExit) as exit_info:
        entry_point.main(_pipeline_args(
            tmp_path, "--processing-aware", "--output-csv", str(tmp_path / "o.csv")))
    assert exit_info.value.code == 1
    assert "--processing-aware-annotation: pepsickle is not installed" in capsys.readouterr().err

    # Passing --no-processing-aware-annotation afterward disables it; must not fail.
    from vaxrank.cli.arg_parser import parse_vaxrank_args
    parsed = parse_vaxrank_args(_pipeline_args(
        tmp_path, "--processing-aware-annotation", "--no-processing-aware-annotation",
        "--output-csv", str(tmp_path / "o.csv")))
    entry_point.require_pepsickle_if_requested(parsed)


def test_pepsickle_availability_single_source_of_truth_in_processing(monkeypatch):
    """#783: processing.pepsickle_installed is the shared check for CLI and processing layer."""
    from vaxrank import processing

    _without_pepsickle(monkeypatch)
    monkeypatch.delenv("PEPSICKLE_PYTHON", raising=False)
    monkeypatch.delenv("PEPSICKLE_GB_PYTHON", raising=False)
    assert processing.pepsickle_installed() is False
    assert entry_point.pepsickle_installed() is False
    assert processing.load_default_processing_predictor() is None


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

    monkeypatch.setattr("varcode.cli.variant_collection_from_args", must_not_run)
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
    from vaxrank.version import __version__

    for prog, main in [("vaxrank-test-data", download.main)]:
        with pytest.raises(SystemExit) as exit_info:
            main(["--version"])
        assert exit_info.value.code == 0
        assert capsys.readouterr().out.strip() == "%s %s" % (prog, __version__)


def test_vaxrank_test_data_version_does_not_import_cli_stack():
    """#783: vaxrank-test-data --version only needs vaxrank.version, not vaxrank.cli."""
    import subprocess
    import sys

    code = (
        "import sys\n"
        "from vaxrank.testdata import download\n"
        "try:\n"
        "    download.main(['--version'])\n"
        "except SystemExit as exc:\n"
        "    assert exc.code == 0\n"
        "assert 'vaxrank.cli' not in sys.modules, 'vaxrank.cli was imported'\n"
    )
    proc = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=False)
    assert proc.returncode == 0, proc.stderr


def test_missing_input_paths_name_the_flag_before_any_work(tmp_path, monkeypatch, capsys):
    def must_not_run(args):
        raise AssertionError("input loading started before the paths were checked")

    monkeypatch.setattr("varcode.cli.variant_collection_from_args", must_not_run)
    out = ["--output-csv", str(tmp_path / "o.csv")]
    cases = [
        ("--bam", str(tmp_path / "nope.bam"), "--bam: no such file"),
        ("--prediction-cache", str(tmp_path / "missing.parquet"),
         "--prediction-cache: no such file or directory"),
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
    parquet_dir = tmp_path / "partitioned_cache"
    parquet_dir.mkdir()
    entry_point.validate_input_paths(SimpleNamespace(
        vcf=[str(present)],
        bam="https://example.org/sample.bam",
        hitlist_evidence_bundle=str(tmp_path),
        prediction_cache=str(parquet_dir),
    ))


def test_a_missing_directory_input_names_the_flag(tmp_path):
    from types import SimpleNamespace
    with pytest.raises(ValueError, match="--hitlist-evidence-bundle: no such directory"):
        entry_point.validate_input_paths(SimpleNamespace(hitlist_evidence_bundle=str(tmp_path / "nope")))



def test_undeclared_antigen_scope_says_where_to_declare_it():
    from types import SimpleNamespace
    from vaxrank.antigen_input import AntigenDataset
    from vaxrank.input_scope import InputScope

    stand_in = SimpleNamespace(
        antigens={"a": object()}, source_evidence={},
        scope=InputScope(patient_id=None, reference_assembly="GRCh38", mhc_alleles=("HLA-A*02:01",)))
    with pytest.raises(ValueError) as error:
        AntigenDataset.validate_for_prediction(stand_in, genome=object())
    message = str(error.value)
    assert "requires declared patient_id" in message
    assert "--input-manifest" in message
    assert "--output-patient-id only labels the report" in message
    assert "--mhc-alleles" not in message

    # #783: when only mhc_alleles is missing, mention --mhc-alleles and omit --output-patient-id.
    missing_alleles = SimpleNamespace(
        antigens={"a": object()}, source_evidence={},
        scope=InputScope(patient_id="P1", reference_assembly="GRCh38", mhc_alleles=None))
    with pytest.raises(ValueError) as alleles_error:
        AntigenDataset.validate_for_prediction(missing_alleles, genome=object())
    alleles_msg = str(alleles_error.value)
    assert "requires declared mhc_alleles" in alleles_msg
    assert "--mhc-alleles" in alleles_msg
    assert "--output-patient-id" not in alleles_msg


def test_predictor_with_an_input_that_reuses_predictions_explains_the_choice(tmp_path, capsys):
    lens = Path(__file__).parent / "data" / "epitope_fixtures" / "lens_example.tsv"
    with pytest.raises(SystemExit) as exit_info:
        entry_point.main(["--input-lens", str(lens), "--mhc-predictor", "random",
                          "--mhc-alleles", "HLA-A*02:01", "--output-csv", str(tmp_path / "o.csv"),
                          "--log-path", str(tmp_path / "python.log")])
    assert exit_info.value.code == 1
    err = capsys.readouterr().err
    assert "reuses the predictions already in its tables" in err
    assert "--external-predictions fresh or additive" in err and "drop --mhc-predictor" in err


def test_bare_vaxrank_prints_examples_and_exits_cleanly(capsys):
    assert entry_point.main([]) is None
    out = capsys.readouterr().out
    assert "Examples:" in out and "--input-pvacseq" in out


def test_mistyped_flag_gets_a_suggestion(capsys):
    with pytest.raises(SystemExit) as exit_info:
        entry_point.main(["--bamm", "x.bam", "--mhc-predictor", "random"])
    assert exit_info.value.code == 2
    err = capsys.readouterr().err
    assert "did you mean --bam?" in err
    assert len(err.strip().splitlines()) == 2


def test_help_carries_worked_examples(capsys):
    for args in (["--vcf", "a.vcf", "--help"], ["--input-lens", "a.tsv", "--help"],
                 ["--input-json-file", "r.json", "--help"]):
        with pytest.raises(SystemExit) as exit_info:
            entry_point.main(args)
        assert exit_info.value.code == 0
        assert "Examples:" in capsys.readouterr().out


def test_unique_flag_prefixes_still_work():
    from vaxrank.cli.arg_parser import parse_vaxrank_args
    parsed = parse_vaxrank_args(["--vcf", "a.vcf", "--bam", "b.bam", "--mhc-predictor", "random",
                                 "--output-d", "out"])
    assert parsed.output_dir == "out"


def test_final_outputs_block_lists_what_was_written(tmp_path, capsys):
    out = tmp_path / "out"
    entry_point.main(_pipeline_args(tmp_path, "--output-dir", str(out),
                                    "--pdf-backend", "weasyprint"))
    err = capsys.readouterr().err
    assert "Outputs written:" in err
    assert str(out / "ranked_candidates.csv") in err
    assert str(out / "run_summary.txt") in err


def test_vaxrank_version_flag_outputs_vaxrank_version(capsys):
    """varcode>=12.0.1 omits --version in make_variants_parser (openvax/varcode#759)."""
    from vaxrank.cli import arg_parser
    from vaxrank.version import __version__

    parser = arg_parser.make_vaxrank_arg_parser()
    with pytest.raises(SystemExit) as exit_info:
        parser.parse_args(['--version'])
    assert exit_info.value.code == 0
    assert f"Vaxrank {__version__}" in capsys.readouterr().out
