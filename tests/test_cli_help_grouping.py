"""Mode inference in errors and design flags grouped by purpose (#766)."""

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
