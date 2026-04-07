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
Tests for filtering variants on non-canonical contigs.

See: https://github.com/openvax/vaxrank/issues/193
"""

from unittest.mock import MagicMock, patch

from varcode import VariantCollection

from vaxrank.cli.entry_point import _filter_variants_to_valid_contigs


CANONICAL_CONTIGS = [
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
    "11", "12", "13", "14", "15", "16", "17", "18", "19",
    "20", "21", "22", "X", "Y", "MT",
]


def _make_variant(contig):
    """Create a mock variant with a contig and a mock ensembl genome."""
    v = MagicMock()
    v.contig = contig
    mock_genome = MagicMock()
    mock_genome.contigs.return_value = CANONICAL_CONTIGS
    mock_genome.reference_name = "GRCh38"
    v.ensembl = mock_genome
    return v


def _make_collection(variants):
    mock = MagicMock()
    mock.__len__ = MagicMock(return_value=len(variants))
    mock.__iter__ = MagicMock(side_effect=lambda: iter(variants))
    mock.__getitem__ = MagicMock(side_effect=lambda i: variants[i])
    return mock


def test_empty_collection():
    variants = VariantCollection([])
    result = _filter_variants_to_valid_contigs(variants)
    assert len(result) == 0


def test_all_valid_contigs():
    v1 = _make_variant("1")
    v2 = _make_variant("X")
    collection = _make_collection([v1, v2])
    result = _filter_variants_to_valid_contigs(collection)
    # Should return the original collection unchanged
    assert result is collection


def test_filters_alt_contigs():
    v_valid = _make_variant("14")
    v_alt = _make_variant("chr14_GL000194v1_random")

    with patch(
        "vaxrank.cli.entry_point.VariantCollection",
        side_effect=lambda variants: variants,
    ):
        result = _filter_variants_to_valid_contigs(
            _make_collection([v_valid, v_alt])
        )
    assert len(result) == 1
    assert result[0].contig == "14"


def test_filters_multiple_alt_contigs():
    v1 = _make_variant("1")
    v2 = _make_variant("chr14_GL000194v1_random")
    v3 = _make_variant("chrUn_gl000220")
    v4 = _make_variant("22")

    with patch(
        "vaxrank.cli.entry_point.VariantCollection",
        side_effect=lambda variants: variants,
    ):
        result = _filter_variants_to_valid_contigs(
            _make_collection([v1, v2, v3, v4])
        )
    assert len(result) == 2
    assert result[0].contig == "1"
    assert result[1].contig == "22"
