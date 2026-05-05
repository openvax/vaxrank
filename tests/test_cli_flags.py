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
from unittest.mock import MagicMock, patch

from vaxrank.cli import make_vaxrank_arg_parser
from vaxrank.cli.entry_point import (
    _resolve_ensembl_release,
    configure_logging,
)


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
        _resolve_ensembl_release(args)

    mock_cls.assert_called_once_with(75)
    assert args.genome is mock_genome


def test_resolve_ensembl_release_noop_when_not_set():
    args = MagicMock()
    args.ensembl_release = None
    args.genome = "GRCh38"
    _resolve_ensembl_release(args)
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


def test_default_peptide_linker_resolves_to_canonical_g4s3():
    """Default ``--peptide-linker GS3`` resolves to the canonical
    15aa flexible linker GGGGSGGGGSGGGGS — not the literal 7aa
    ``GGGGSSS`` that the bare ``G4S3`` form parses as via the
    GnSm grammar."""
    from vaxrank.vaccine_library import get_linker
    parser = make_vaxrank_arg_parser()
    args = parser.parse_args([
        "--vcf", "dummy.vcf", "--bam", "dummy.bam",
        "--mhc-predictor", "random", "--mhc-alleles", "HLA-A*02:01",
        "--output-csv", "out.csv",
    ])
    assert args.peptide_linker == "GS3"
    linker = get_linker(args.peptide_linker)
    assert linker.amino_acids == "GGGGSGGGGSGGGGS"
    assert len(linker.amino_acids) == 15
