from tempfile import NamedTemporaryFile
from vaxrank.cli import main as run_shell_script

from .testing_helpers import data_path

cli_args_for_b16_seqdata = [
    "--vcf", data_path("b16.f10/b16.vcf"),
    "--bam", data_path("b16.f10/b16.combined.bam"),
    "--vaccine-peptide-length", "25",
    "--mhc-predictor", "random",
    "--mhc-alleles", "H2-Kb,H2-Db",
    "--padding-around-mutation", "5"
]

def test_csv_report():
    with NamedTemporaryFile(mode="r", delete=False) as f:
        csv_args = cli_args_for_b16_seqdata + ["--output-csv", f.name]
        run_shell_script(csv_args)
        contents = f.read()
        lines = contents.split("\n")
        assert len(lines) > 0

def test_html_report():
    with NamedTemporaryFile(mode="r") as f:
        html_args = cli_args_for_b16_seqdata + ["--output-html", f.name]
        run_shell_script(html_args)
        contents = f.read()
        lines = contents.split("\n")
        assert len(lines) > 0

if __name__ == "__main__":
    test_csv_report()
    test_html_report()
