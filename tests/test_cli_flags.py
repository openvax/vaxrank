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

"""
Tests for --ensembl-release and --verbose CLI flags.

See: https://github.com/openvax/vaxrank/issues/29
      https://github.com/openvax/vaxrank/issues/62
"""

import logging
import os
import subprocess
import sys
from unittest.mock import MagicMock, patch

from vaxrank.cli import make_vaxrank_arg_parser
from vaxrank.cli.entry_point import (
    resolve_ensembl_release,
    configure_logging,
)


def test_help_succeeds_under_ascii_locale():
    """CLI help must not depend on a UTF-8 terminal (issue #308)."""
    env = os.environ.copy()
    env.update({
        "LANG": "C",
        "LC_ALL": "C",
        "PYTHONCOERCECLOCALE": "0",
        "PYTHONUTF8": "0",
    })
    command = (
        "from vaxrank.cli.entry_point import main; "
        "main(['-h'])"
    )

    result = subprocess.run(
        [sys.executable, "-c", command],
        capture_output=True,
        text=True,
        check=False,
        env=env,
    )

    assert result.returncode == 0, result.stderr
    assert result.stdout.startswith("usage: vaxrank")
    assert result.stdout.isascii()


# ---- --ensembl-release (#29) ------------------------------------------------

def test_ensembl_release_arg_parsed():
    parser = make_vaxrank_arg_parser()
    args = parser.parse_args([
        "--vcf", "dummy.vcf",
        "--bam", "dummy.bam",
        "--mhc-predictor", "random",
        "--mhc-alleles", "HLA-A*02:01",
        "--ensembl-release", "75",
        "--output-csv", "out.csv",
    ])
    assert args.ensembl_release == 75


def test_ensembl_release_default_is_none():
    parser = make_vaxrank_arg_parser()
    args = parser.parse_args([
        "--vcf", "dummy.vcf",
        "--bam", "dummy.bam",
        "--mhc-predictor", "random",
        "--mhc-alleles", "HLA-A*02:01",
        "--output-csv", "out.csv",
    ])
    assert args.ensembl_release is None


def test_resolve_ensembl_release_sets_genome():
    args = MagicMock()
    args.ensembl_release = 75
    args.genome = None

    with patch(
        "pyensembl.EnsemblRelease",
    ) as mock_cls:
        mock_genome = MagicMock()
        mock_cls.return_value = mock_genome
        resolve_ensembl_release(args)

    mock_cls.assert_called_once_with(75)
    assert args.genome is mock_genome


def test_resolve_ensembl_release_noop_when_not_set():
    args = MagicMock()
    args.ensembl_release = None
    args.genome = "GRCh38"
    resolve_ensembl_release(args)
    assert args.genome == "GRCh38"


# ---- --verbose (#62) --------------------------------------------------------

def test_verbose_arg_parsed():
    parser = make_vaxrank_arg_parser()
    args = parser.parse_args([
        "--vcf", "dummy.vcf",
        "--bam", "dummy.bam",
        "--mhc-predictor", "random",
        "--mhc-alleles", "HLA-A*02:01",
        "--verbose",
        "--output-csv", "out.csv",
    ])
    assert args.verbose is True


def test_verbose_default_is_false():
    parser = make_vaxrank_arg_parser()
    args = parser.parse_args([
        "--vcf", "dummy.vcf",
        "--bam", "dummy.bam",
        "--mhc-predictor", "random",
        "--mhc-alleles", "HLA-A*02:01",
        "--output-csv", "out.csv",
    ])
    assert args.verbose is False


def test_verbose_sets_console_handlers_to_debug(tmp_path):
    """When --verbose is set, console handlers for vaxrank loggers
    should be set to DEBUG level."""
    args = MagicMock()
    args.verbose = True
    args.log_path = str(tmp_path / "test.log")

    configure_logging(args)

    vaxrank_logger = logging.getLogger("vaxrank")
    console_handlers = [
        h for h in vaxrank_logger.handlers
        if isinstance(h, logging.StreamHandler)
        and not isinstance(h, logging.FileHandler)
    ]
    assert len(console_handlers) > 0
    for h in console_handlers:
        assert h.level == logging.DEBUG


def test_no_verbose_keeps_console_handlers_at_info(tmp_path):
    """Without --verbose, console handlers should stay at INFO level."""
    args = MagicMock()
    args.verbose = False
    args.log_path = str(tmp_path / "test.log")

    configure_logging(args)

    vaxrank_logger = logging.getLogger("vaxrank")
    console_handlers = [
        h for h in vaxrank_logger.handlers
        if isinstance(h, logging.StreamHandler)
        and not isinstance(h, logging.FileHandler)
    ]
    assert len(console_handlers) > 0
    for h in console_handlers:
        assert h.level == logging.INFO


# ---- linker arg validation (#246) -------------------------------------------

def test_peptide_linker_rejects_unknown(capsys):
    """argparse should reject unknown --peptide-linker before assembly runs."""
    parser = make_vaxrank_arg_parser()
    try:
        parser.parse_args([
            "--vcf", "dummy.vcf", "--bam", "dummy.bam",
            "--mhc-predictor", "random", "--mhc-alleles", "HLA-A*02:01",
            "--output-csv", "out.csv",
            "--peptide-linker", "NOTREAL",
        ])
    except SystemExit:
        captured = capsys.readouterr()
        assert "NOTREAL" in captured.err
    else:
        raise AssertionError("argparse did not reject --peptide-linker NOTREAL")


def test_mrna_linker_rejects_unknown(capsys):
    parser = make_vaxrank_arg_parser()
    try:
        parser.parse_args([
            "--vcf", "dummy.vcf", "--bam", "dummy.bam",
            "--mhc-predictor", "random", "--mhc-alleles", "HLA-A*02:01",
            "--output-csv", "out.csv",
            "--mrna-linker", "NOTREAL",
        ])
    except SystemExit:
        captured = capsys.readouterr()
        assert "NOTREAL" in captured.err
    else:
        raise AssertionError("argparse did not reject --mrna-linker NOTREAL")


def test_peptide_linker_accepts_lowercase():
    """`--peptide-linker p2a` should work the same as `P2A` via str.upper."""
    parser = make_vaxrank_arg_parser()
    args = parser.parse_args([
        "--vcf", "dummy.vcf", "--bam", "dummy.bam",
        "--mhc-predictor", "random", "--mhc-alleles", "HLA-A*02:01",
        "--output-csv", "out.csv",
        "--peptide-linker", "p2a",
    ])
    assert args.peptide_linker == "P2A"


def test_mrna_linker_accepts_lowercase():
    parser = make_vaxrank_arg_parser()
    args = parser.parse_args([
        "--vcf", "dummy.vcf", "--bam", "dummy.bam",
        "--mhc-predictor", "random", "--mhc-alleles", "HLA-A*02:01",
        "--output-csv", "out.csv",
        "--mrna-linker", "gs3",
    ])
    assert args.mrna_linker == "GS3"


def test_peptide_mode_is_case_insensitive():
    """``--peptide-mode SLP`` / ``Slp`` / ``slp`` all parse to the
    same canonical lowercase ``'slp'`` so users don't get bitten by
    case mismatch in YAML configs or scripts."""
    parser = make_vaxrank_arg_parser()
    for raw in ['slp', 'SLP', 'Slp', 'sLP']:
        args = parser.parse_args([
            "--vcf", "dummy.vcf", "--bam", "dummy.bam",
            "--mhc-predictor", "random", "--mhc-alleles", "HLA-A*02:01",
            "--output-csv", "out.csv",
            "--peptide-mode", raw,
        ])
        assert args.peptide_mode == "slp"


def test_peptide_construct_config_normalizes_uppercase_mode():
    """YAML configs bypass argparse, so the normalization must also
    happen in ``PeptideConstructConfig.__post_init__``."""
    from vaxrank.peptide import PeptideConstructConfig
    cfg = PeptideConstructConfig(mode='SLP', antigens_per_construct=1)
    assert cfg.mode == 'slp'
    cfg2 = PeptideConstructConfig(
        mode='Minimal_Epitope', antigens_per_construct=1)
    assert cfg2.mode == 'minimal_epitope'
    assert cfg2.antigen_content == 'minimal_epitope'


def test_default_peptide_linker_is_unambiguous_and_canonical():
    """Default ``--peptide-linker G4Sx3`` is the self-documenting
    "(G4S) repeated 3x" form — resolves to GGGGSGGGGSGGGGS (15 aa,
    canonical FixVac flexible linker). The bare ``G4S3`` parses as
    a different 7aa literal (4 G + 3 S, GnSm grammar) and the
    legacy ``GS3`` reads ambiguously as "G + S + 3," so neither is
    used as the default."""
    from vaxrank.vaccine_library import get_linker
    parser = make_vaxrank_arg_parser()
    args = parser.parse_args([
        "--vcf", "dummy.vcf", "--bam", "dummy.bam",
        "--mhc-predictor", "random", "--mhc-alleles", "HLA-A*02:01",
        "--output-csv", "out.csv",
    ])
    assert args.peptide_linker == "G4SX3"   # _linker_arg uppercases
    linker = get_linker(args.peptide_linker)
    assert linker.amino_acids == "GGGGSGGGGSGGGGS"
    assert len(linker.amino_acids) == 15


def test_default_mrna_linker_is_unambiguous_and_canonical():
    """Default ``--mrna-linker G4Sx2`` resolves to GGGGSGGGGS (10 aa,
    canonical FixVac mRNA inter-antigen linker). Same self-
    documenting compositional form as the peptide-side default."""
    from vaxrank.vaccine_library import get_linker
    parser = make_vaxrank_arg_parser()
    args = parser.parse_args([
        "--vcf", "dummy.vcf", "--bam", "dummy.bam",
        "--mhc-predictor", "random", "--mhc-alleles", "HLA-A*02:01",
        "--output-csv", "out.csv",
    ])
    assert args.mrna_linker == "G4SX2"
    linker = get_linker(args.mrna_linker)
    assert linker.amino_acids == "GGGGSGGGGS"
    assert len(linker.amino_acids) == 10
