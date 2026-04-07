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
Tests for filtering variants on unannotatable contigs.

See: https://github.com/openvax/vaxrank/issues/193
"""

from unittest.mock import MagicMock, patch, PropertyMock

from varcode import VariantCollection

from vaxrank.cli.entry_point import _filter_unannotatable_variants


def _make_variant(contig, annotatable=True):
    """
    Create a mock variant.

    If annotatable=False, accessing gene_names will raise ValueError,
    simulating a contig that varcode can't resolve.
    """
    v = MagicMock()
    v.contig = contig
    mock_genome = MagicMock()
    mock_genome.reference_name = "GRCh38"
    v.ensembl = mock_genome
    if annotatable:
        type(v).gene_names = PropertyMock(return_value=["BRCA1"])
    else:
        type(v).gene_names = PropertyMock(
            side_effect=ValueError(
                "Invalid contig name '%s' for reference 'GRCh38'" % contig
            )
        )
    return v


def _make_collection(variants):
    mock = MagicMock()
    mock.__len__ = MagicMock(return_value=len(variants))
    mock.__iter__ = MagicMock(side_effect=lambda: iter(variants))
    mock.__getitem__ = MagicMock(side_effect=lambda i: variants[i])
    return mock


def test_empty_collection():
    variants = VariantCollection([])
    result = _filter_unannotatable_variants(variants)
    assert len(result) == 0


def test_all_annotatable():
    v1 = _make_variant("1")
    v2 = _make_variant("X")
    collection = _make_collection([v1, v2])
    result = _filter_unannotatable_variants(collection)
    # Should return the original collection unchanged
    assert result is collection


def test_filters_unannotatable_contig():
    v_ok = _make_variant("14")
    v_bad = _make_variant("chr14_GL000194v1_random", annotatable=False)

    with patch(
        "vaxrank.cli.entry_point.VariantCollection",
        side_effect=lambda variants: variants,
    ):
        result = _filter_unannotatable_variants(
            _make_collection([v_ok, v_bad])
        )
    assert len(result) == 1
    assert result[0].contig == "14"


def test_filters_multiple_unannotatable_contigs():
    v1 = _make_variant("1")
    v2 = _make_variant("chr14_GL000194v1_random", annotatable=False)
    v3 = _make_variant("chrUn_gl000220", annotatable=False)
    v4 = _make_variant("22")

    with patch(
        "vaxrank.cli.entry_point.VariantCollection",
        side_effect=lambda variants: variants,
    ):
        result = _filter_unannotatable_variants(
            _make_collection([v1, v2, v3, v4])
        )
    assert len(result) == 2
    assert result[0].contig == "1"
    assert result[1].contig == "22"


def test_chr_prefixed_contigs_kept_if_annotatable():
    """Variants with chr-prefixed contigs should be kept if varcode can resolve them."""
    v_chr1 = _make_variant("chr1", annotatable=True)
    v_chrX = _make_variant("chrX", annotatable=True)
    collection = _make_collection([v_chr1, v_chrX])
    result = _filter_unannotatable_variants(collection)
    assert result is collection
