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
