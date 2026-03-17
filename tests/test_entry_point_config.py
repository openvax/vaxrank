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
Tests for config override behavior in CLI entry points.
"""

from unittest.mock import MagicMock, patch

from vaxrank.cli import make_vaxrank_arg_parser, run_vaxrank_from_parsed_args


def test_config_overrides_defaults_for_protein_sequence_length(tmp_path):
    config_path = tmp_path / "vaxrank.yaml"
    config_path.write_text(
        "\n".join(
            [
                "vaccine_config:",
                "  vaccine_peptide_length: 31",
                "  padding_around_mutation: 9",
                "",
            ]
        )
    )

    parser = make_vaxrank_arg_parser()
    args = parser.parse_args(
        [
            "--config",
            str(config_path),
            "--vcf",
            "variants.vcf",
            "--bam",
            "tumor.bam",
            "--mhc-predictor",
            "random",
            "--mhc-alleles",
            "HLA-A*02:01",
        ]
    )

    captured = {}

    def fake_run_isovar(parsed_args):
        captured["protein_sequence_length"] = parsed_args.protein_sequence_length
        return []

    with patch(
        "vaxrank.cli.entry_point.run_isovar_from_parsed_args",
        side_effect=fake_run_isovar,
    ), patch(
        "vaxrank.cli.entry_point.mhc_binding_predictor_from_args",
        return_value=MagicMock(),
    ), patch(
        "vaxrank.cli.entry_point.run_vaxrank",
        return_value=MagicMock(),
    ):
        run_vaxrank_from_parsed_args(args)

    expected = 31 + 2 * 9
    assert captured["protein_sequence_length"] == expected
