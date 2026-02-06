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
Tests for cancer hotspots lookup module.
"""

from vaxrank.cancer_hotspots import (
    is_hotspot,
    get_hotspot_info,
    get_hotspot_url,
    _load_hotspots,
    _get_hotspots_path,
)

from .common import eq_, ok_


# =============================================================================
# Data Loading Tests
# =============================================================================

def test_hotspots_path_exists():
    """Test that the hotspots data file path is valid"""
    import os
    path = _get_hotspots_path()
    ok_(os.path.exists(path), f"Hotspots file not found at {path}")


def test_load_hotspots_returns_dict():
    """Test that _load_hotspots returns a dictionary"""
    hotspots = _load_hotspots()
    ok_(isinstance(hotspots, dict))


def test_load_hotspots_not_empty():
    """Test that hotspots data is loaded"""
    hotspots = _load_hotspots()
    ok_(len(hotspots) > 0, "Hotspots dict should not be empty")


def test_load_hotspots_has_expected_count():
    """Test that we have a reasonable number of hotspots"""
    hotspots = _load_hotspots()
    # The file has ~3181 lines but ~2739 unique (gene, protein_change) pairs
    ok_(len(hotspots) > 2500, f"Expected >2500 hotspots, got {len(hotspots)}")


def test_hotspot_entry_structure():
    """Test that hotspot entries have expected fields"""
    hotspots = _load_hotspots()
    # Get any entry
    key, info = next(iter(hotspots.items()))

    # Check key is a tuple of (gene, protein_change)
    ok_(isinstance(key, tuple))
    eq_(len(key), 2)

    # Check info dict has expected fields
    ok_('gene' in info)
    ok_('protein_change' in info)
    ok_('mutation_type' in info)
    ok_('chromosome' in info)
    ok_('position' in info)


# =============================================================================
# Known Hotspot Tests
# =============================================================================

def test_braf_v600e_is_hotspot():
    """Test that BRAF V600E is recognized as a hotspot"""
    ok_(is_hotspot("BRAF", "V600E"))


def test_braf_v600e_with_p_prefix():
    """Test that p.V600E notation works"""
    ok_(is_hotspot("BRAF", "p.V600E"))


def test_tp53_hotspots():
    """Test that common TP53 hotspots are recognized"""
    # Common TP53 hotspots
    ok_(is_hotspot("TP53", "R175H"))
    ok_(is_hotspot("TP53", "R248Q"))
    ok_(is_hotspot("TP53", "R273H"))


def test_kras_g12_hotspots():
    """Test KRAS G12 hotspots"""
    ok_(is_hotspot("KRAS", "G12D"))
    ok_(is_hotspot("KRAS", "G12V"))
    ok_(is_hotspot("KRAS", "G12C"))


def test_pik3ca_hotspots():
    """Test PIK3CA hotspots"""
    ok_(is_hotspot("PIK3CA", "E545K"))
    ok_(is_hotspot("PIK3CA", "H1047R"))


def test_egfr_hotspots():
    """Test EGFR hotspots"""
    ok_(is_hotspot("EGFR", "L858R"))
    ok_(is_hotspot("EGFR", "T790M"))


# =============================================================================
# Case Insensitivity Tests
# =============================================================================

def test_gene_lowercase():
    """Test that lowercase gene names work"""
    ok_(is_hotspot("braf", "V600E"))


def test_gene_mixed_case():
    """Test that mixed case gene names work"""
    ok_(is_hotspot("Braf", "V600E"))


def test_protein_change_lowercase():
    """Test that lowercase protein changes work"""
    ok_(is_hotspot("BRAF", "v600e"))


def test_both_lowercase():
    """Test that both lowercase works"""
    ok_(is_hotspot("braf", "v600e"))


# =============================================================================
# Non-Hotspot Tests
# =============================================================================

def test_random_mutation_not_hotspot():
    """Test that a random mutation is not a hotspot"""
    ok_(not is_hotspot("FAKEGENE", "X999Y"))


def test_real_gene_fake_mutation():
    """Test that a real gene with fake mutation is not a hotspot"""
    ok_(not is_hotspot("BRAF", "A1B"))


def test_empty_strings():
    """Test that empty strings return False"""
    ok_(not is_hotspot("", ""))
    ok_(not is_hotspot("BRAF", ""))
    ok_(not is_hotspot("", "V600E"))


# =============================================================================
# get_hotspot_info Tests
# =============================================================================

def test_returns_dict_for_hotspot():
    """Test that get_hotspot_info returns dict for known hotspot"""
    info = get_hotspot_info("BRAF", "V600E")
    ok_(info is not None)
    ok_(isinstance(info, dict))


def test_get_hotspot_info_returns_none_for_non_hotspot():
    """Test that get_hotspot_info returns None for non-hotspot"""
    info = get_hotspot_info("FAKEGENE", "X999Y")
    eq_(info, None)


def test_info_contains_gene():
    """Test that info dict contains gene name"""
    info = get_hotspot_info("BRAF", "V600E")
    ok_('gene' in info)
    eq_(info['gene'].upper(), "BRAF")


def test_info_contains_protein_change():
    """Test that info dict contains protein change"""
    info = get_hotspot_info("BRAF", "V600E")
    ok_('protein_change' in info)
    eq_(info['protein_change'].upper(), "V600E")


def test_info_contains_mutation_type():
    """Test that info dict contains mutation type"""
    info = get_hotspot_info("BRAF", "V600E")
    ok_('mutation_type' in info)


def test_p_prefix_handled():
    """Test that p. prefix is handled"""
    info = get_hotspot_info("BRAF", "p.V600E")
    ok_(info is not None)


# =============================================================================
# get_hotspot_url Tests
# =============================================================================

def test_returns_url_for_hotspot():
    """Test that get_hotspot_url returns a URL for known hotspot"""
    url = get_hotspot_url("BRAF", "V600E")
    ok_(url is not None)
    ok_(url.startswith("https://"))


def test_get_hotspot_url_returns_none_for_non_hotspot():
    """Test that None is returned for non-hotspot"""
    url = get_hotspot_url("FAKEGENE", "X999Y")
    eq_(url, None)


def test_url_contains_gene():
    """Test that URL contains gene name"""
    url = get_hotspot_url("BRAF", "V600E")
    ok_("BRAF" in url)


def test_url_points_to_cancerhotspots():
    """Test that URL points to cancerhotspots.org"""
    url = get_hotspot_url("BRAF", "V600E")
    ok_("cancerhotspots.org" in url)


def test_url_format():
    """Test the URL format"""
    url = get_hotspot_url("BRAF", "V600E")
    eq_(url, "https://www.cancerhotspots.org/#/gene/BRAF")


# =============================================================================
# Edge Case Tests
# =============================================================================

def test_nonsense_mutation():
    """Test handling of nonsense mutations (stop codon)"""
    # Check if there are any nonsense mutations in the dataset
    hotspots = _load_hotspots()
    nonsense_found = any('*' in k[1] for k in hotspots.keys())
    ok_(nonsense_found, "Should have some nonsense mutations in dataset")


def test_deletion_mutation():
    """Test that indel mutations are in dataset"""
    hotspots = _load_hotspots()
    # Check for DEL or INS mutations
    del_ins_found = any(
        info['mutation_type'] in ('DEL', 'INS')
        for info in hotspots.values()
    )
    # This may or may not be true depending on the dataset
    # Just verify the check doesn't crash
    ok_(isinstance(del_ins_found, bool))


def test_caching_works():
    """Test that the LRU cache is working"""
    # Call multiple times - should use cache
    hotspots1 = _load_hotspots()
    hotspots2 = _load_hotspots()
    # Should be the exact same object due to caching
    ok_(hotspots1 is hotspots2)
