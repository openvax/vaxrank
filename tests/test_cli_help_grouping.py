"""Mode inference in errors and design flags grouped by purpose (#766)."""

import sys

import pytest

from vaxrank.cli.arg_parser import infer_mode, parse_vaxrank_args


@pytest.mark.parametrize("argv, mode, trigger", [
    (["--input-json-file", "r.json"], "replay", "--input-json-file"),
    (["--input-json-file=r.json"], "replay", "--input-json-file"),
    (["--vcf", "a.vcf", "--bam", "x.bam"], "run", "--vcf"),
    (["--input-pvacseq", "p.tsv"], "import", "--input-pvacseq"),
    (["--input-lens=l.tsv"], "import", "--input-lens"),
    (["--input-epitopes", "e.tsv"], "import", "--input-epitopes"),
    ([], "run", None),
    (["--output-dir", "o"], "run", None),
])
def test_infer_mode_names_the_flag_that_decided_it(argv, mode, trigger):
    assert infer_mode(argv) == (mode, trigger)


def test_a_saved_run_wins_over_pipeline_flags():
    assert infer_mode(["--vcf", "a.vcf", "--input-json-file", "r.json"]) == ("replay", "--input-json-file")


def test_error_says_which_mode_a_flag_selected(capsys):
    """``--vcf`` with a LENS report fails on ``--bam``; say why the VCF mode applied."""
    with pytest.raises(SystemExit) as raised:
        parse_vaxrank_args(["--input-lens", "x.tsv", "--vcf", "y.vcf"])
    assert raised.value.code == 2
    err = capsys.readouterr().err
    assert "required: --bam" in err
    assert "Read as a variant and RNA run because --vcf was given." in err


def test_flag_from_another_mode_names_the_mode_that_rejected_it(capsys):
    with pytest.raises(SystemExit):
        parse_vaxrank_args(["--input-json-file", "r.json", "--vcf", "a.vcf"])
    err = capsys.readouterr().err
    assert "unrecognized argument: --vcf" in err
    assert "Read as the re-rendering of a saved run because --input-json-file was given." in err


def test_error_without_a_deciding_flag_stays_one_line_plus_pointer(capsys):
    with pytest.raises(SystemExit):
        parse_vaxrank_args(["--bamm", "x.bam", "--mhc-predictor", "random"])
    err = capsys.readouterr().err
    assert "Read as" not in err
    assert len(err.strip().splitlines()) == 2


def _option_sections(help_text):
    """Map each help section to the flags it defines (not flags its prose mentions)."""
    sections, current = {}, None
    for line in help_text.splitlines():
        if line and not line.startswith(" ") and line.endswith(":"):
            current = line[:-1]
            sections[current] = set()
        elif current and line.startswith("  -"):
            sections[current].update(word.rstrip(",") for word in line.split() if word.startswith("--"))
    return sections


def test_design_flags_are_not_listed_under_output_options(capsys):
    with pytest.raises(SystemExit):
        parse_vaxrank_args(["--vcf", "a.vcf", "--help"])
    sections = _option_sections(capsys.readouterr().out)
    for flag in ("--vaccine-type", "--antigen-content", "--max-sources", "--mrna-poly-a-length",
                 "--peptide-linker", "--pepsickle-threshold"):
        assert flag not in sections["Output options"], flag
    assert "--vaccine-type" in sections["Vaccine design options"]
    assert "--max-sources" in sections["Vaccine design options"]
    assert "--mrna-poly-a-length" in sections["mRNA construct options"]
    assert "--peptide-linker" in sections["Peptide construct options"]
    assert "--pepsickle-threshold" in sections["Processing-aware annotation options"]
    assert "--output-dir" in sections["Output options"]


# --- activities and help topics -------------------------------------------

def test_activity_words_rewrite_to_the_equivalent_flags():
    from vaxrank.cli.discoverability import resolve_command_line

    assert resolve_command_line(["policies"]) == (False, ["--list-policies"])
    assert resolve_command_line(["policies", "show", "lens-v1"]) == (
        False, ["--show-policy", "lens-v1"])
    assert resolve_command_line(["config"]) == (False, ["--print-default-config"])
    assert resolve_command_line(["replay", "ranked.json", "--output-dir", "o"]) == (
        False, ["--input-json-file", "ranked.json", "--output-dir", "o"])
    assert resolve_command_line(["replay", "design.tsv"]) == (
        False, ["--input-epitopes", "design.tsv"])
    # A run takes flags only: nothing is rewritten.
    run = ["--vcf", "a.vcf", "--bam", "b.bam"]
    assert resolve_command_line(run) == (False, run)
    # `run` is an optional spelling of the default.
    assert resolve_command_line(["run", "--vcf", "a.vcf"]) == (
        False, ["--vcf", "a.vcf"])


@pytest.mark.parametrize("argv", [
    ["frobnicate"], ["run", "policies"], ["policies", "bogus"], ["replay"], ["config", "x"],
    ["help", "nope"], ["help", "cta", "extra"]])
def test_activity_mistakes_are_one_line_usage_errors(argv, capsys):
    from vaxrank.cli.entry_point import main

    with pytest.raises(SystemExit) as exit_info:
        main(argv)
    assert exit_info.value.code == 2
    err = capsys.readouterr().err
    assert err.startswith("vaxrank: error: ")
    assert len(err.strip().splitlines()) == 1


def test_short_help_is_one_screen_and_help_all_topic_is_the_full_listing(capsys):
    from vaxrank.cli.entry_point import main

    with pytest.raises(SystemExit) as exit_info:
        main(["--help"])
    assert exit_info.value.code == 0
    short = capsys.readouterr().out
    assert short.startswith("usage: vaxrank")
    assert "Examples:" in short and "vaxrank help [TOPIC]" in short
    assert len(short.splitlines()) < 40
    assert main(["help", "all"]) is None
    full = capsys.readouterr().out
    assert len(full.splitlines()) > 400
    assert "--include-lens-source" in full


def test_policies_and_config_run_without_flags(capsys):
    from vaxrank.cli.entry_point import main

    with pytest.raises(SystemExit) as exit_info:
        main(["policies"])
    assert exit_info.value.code == 0
    assert "openvax-v1" in capsys.readouterr().out
    with pytest.raises(SystemExit):
        main(["config"])
    assert "vaccine_peptides" in capsys.readouterr().out


def test_every_help_group_belongs_to_a_topic():
    from vaxrank.cli.arg_parser import external_input_arg_parser, make_vaxrank_arg_parser
    from vaxrank.cli.discoverability import HELP_TOPICS, topic_help

    covered = [title for _, titles in HELP_TOPICS.values() for title in titles]
    assert len(covered) == len(set(covered))
    for parser in (make_vaxrank_arg_parser(), external_input_arg_parser()):
        for group in parser._action_groups:
            if group.title in ("positional arguments", "options"):
                continue
            assert group.title in covered, "%s has no help topic" % group.title
    for topic in HELP_TOPICS:
        text = topic_help(topic)
        assert text.startswith("vaxrank help %s:" % topic)
        for title in HELP_TOPICS[topic][1]:
            assert title + ":" in text


@pytest.mark.parametrize("argv", [
    ["config", "--help"], ["policies", "-h"], ["help", "--help"],
    ["replay", "x.json", "--help"]])
def test_help_flag_after_an_activity_word_shows_help(argv, capsys):
    from vaxrank.cli.entry_point import main

    with pytest.raises(SystemExit) as exit_info:
        main(argv) or sys.exit(0)
    assert exit_info.value.code == 0
    assert capsys.readouterr().out.startswith("usage: vaxrank")
