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


from os.path import getsize
from mock import patch
from nose.plugins.attrib import attr
from tempfile import NamedTemporaryFile

import pandas as pd
from xlrd import open_workbook

from vaxrank.cli import main as run_shell_script

from .testing_helpers import data_path

cli_args_for_b16_seqdata = [
    "--vcf", data_path("b16.f10/b16.vcf"),
    "--bam", data_path("b16.f10/b16.combined.bam"),
    "--vaccine-peptide-length", "25",
    "--mhc-predictor", "random",
    "--mhc-alleles", "H2-Kb,H2-Db",
    "--padding-around-mutation", "5",
    "--count-mismatches-after-variant",
]

cli_args_for_b16_seqdata_real_predictor = [
    "--vcf", data_path("b16.f10/b16.vcf"),
    "--bam", data_path("b16.f10/b16.combined.bam"),
    "--vaccine-peptide-length", "25",
    "--mhc-predictor", "netmhcpan",
    "--mhc-alleles", "H2-Kb,H2-Db",
    "--mhc-epitope-lengths", "8",
    "--padding-around-mutation", "5",
    "--count-mismatches-after-variant"
]


def test_ascii_report():
    with NamedTemporaryFile(mode="r") as f:
        ascii_args = cli_args_for_b16_seqdata + ["--output-ascii-report", f.name]
        run_shell_script(ascii_args)
        contents = f.read()
        lines = contents.split("\n")
        assert len(lines) > 0


def test_ascii_report_real_netmhc_predictor():
    with NamedTemporaryFile(mode="r") as f:
        ascii_args = cli_args_for_b16_seqdata_real_predictor + [
            "--output-ascii-report", f.name]
        run_shell_script(ascii_args)
        contents = f.read()
        lines = contents.split("\n")
        assert len(lines) > 0
        no_variants_text = 'No variants'
        assert no_variants_text not in contents


def test_json_report():
    with NamedTemporaryFile(mode="r") as f:
        json_args = cli_args_for_b16_seqdata + ["--output-json-file", f.name]
        run_shell_script(json_args)
        contents = f.read()
        lines = contents.split("\n")
        assert len(lines) > 0


def test_csv_report():
    with NamedTemporaryFile(mode="r") as f:
        csv_args = cli_args_for_b16_seqdata + ["--output-csv", f.name]
        run_shell_script(csv_args)
        contents = f.read()
        lines = contents.split("\n")
        assert len(lines) > 1


def test_all_variant_csv_report():
    with NamedTemporaryFile(mode="r") as f:
        all_csv_args = cli_args_for_b16_seqdata + [
            "--output-passing-variants-csv", f.name,
            # TODO: make this flag not necessary
            "--output-csv", f.name + "ignored"]
        run_shell_script(all_csv_args)
        contents = f.read()
        lines = contents.split("\n")
        assert len(lines) > 1
        # make sure it can be a valid dataframe
        f.seek(0)
        df = pd.read_csv(f)
        assert len(df) > 1

def test_isovar_csv():
    with NamedTemporaryFile(mode="r") as f:
        isovar_csv_args = cli_args_for_b16_seqdata + [
            "--output-isovar-csv", f.name,
            # TODO: make this flag not necessary
            "--output-csv", f.name + "ignored"
        ]
        run_shell_script(isovar_csv_args)
        df = pd.read_csv(f)
        assert len(df) > 1

def test_xlsx_report():
    with NamedTemporaryFile(mode="r") as f:
        xlsx_args = cli_args_for_b16_seqdata + ["--output-xlsx-report", f.name]
        run_shell_script(xlsx_args)
        book = open_workbook(f.name)
        assert book.nsheets > 1




def test_html_report():
    with NamedTemporaryFile(mode="r") as f:
        html_args = cli_args_for_b16_seqdata + ["--output-html", f.name]
        run_shell_script(html_args)
        contents = f.read()
        lines = contents.split("\n")
        assert len(lines) > 1


@attr('skip')  # want the ability to skip this test on some machines
def test_pdf_report():
    with NamedTemporaryFile(mode="rb") as f:
        pdf_args = cli_args_for_b16_seqdata + ["--output-pdf-report", f.name]
        run_shell_script(pdf_args)
        assert getsize(f.name) > 1


@patch('vaxrank.core_logic.vaccine_peptides_for_variant')
def test_report_no_peptides(mock_vaccine_peptides_for_variant):
    # simulate case where we have no epitopes for any variant
    mock_vaccine_peptides_for_variant.return_value = []
    with NamedTemporaryFile(mode="r") as f:
        html_args = cli_args_for_b16_seqdata + ["--output-csv", f.name]
        # test that this doesn't crash and that the CSV output is empty
        run_shell_script(html_args)
        contents = f.read()
        assert contents == ''


if __name__ == "__main__":
    test_csv_report()
    test_html_report()
