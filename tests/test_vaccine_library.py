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

"""Tests for the shared vaccine_library and the mrna_library compatibility shim."""

import pytest

from vaxrank import mrna_library, vaccine_library


def test_mrna_library_exports_old_string_constants():
    # 2.10.0 publicly exposed plain-string linker constants. Ensure the
    # 2.11+ refactor (Linker dataclass in vaccine_library) preserves
    # those imports.
    assert mrna_library.LINKER_GS == "GGGGS"
    assert mrna_library.LINKER_GS3 == "GGGGSGGGGSGGGGS"
    assert mrna_library.LINKER_AAY == "AAY"
    assert mrna_library.LINKER_GPGPG == "GPGPG"


def test_mrna_library_linkers_dict_returns_strings():
    # Pre-refactor consumers used LINKERS["GS3"] expecting a string,
    # not a Linker object; preserve that shape in mrna_library.
    for name, value in mrna_library.LINKERS.items():
        assert isinstance(value, str), \
            "mrna_library.LINKERS['%s'] must remain a string for compat" % name


def test_vaccine_library_exposes_linker_objects():
    # The richer dataclass form lives in vaccine_library.
    gs3 = vaccine_library.LINKERS["GS3"]
    assert isinstance(gs3, vaccine_library.Linker)
    assert gs3.amino_acids == "GGGGSGGGGSGGGGS"
    assert gs3.dna is None
    assert gs3.freeze_in_mrna is False


def test_vaccine_library_p2a_carries_blessed_dna():
    p2a = vaccine_library.LINKERS["P2A"]
    assert p2a.amino_acids == "GSGATNFSLLKQAGDVEENPGP"
    # AA / DNA length consistency: 22 aa * 3 = 66 nt
    assert len(p2a.dna) == len(p2a.amino_acids) * 3
    assert p2a.freeze_in_mrna is True
    assert p2a.inert_in_peptide_mode is True


def test_vaccine_library_get_linker_unknown_raises():
    with pytest.raises(ValueError):
        vaccine_library.get_linker("notalinker")


def test_signal_peptides_include_tcr_additions():
    # CD8A and CD28 added from openvax/tcrsift's vetted set; confirm
    # they're discoverable through the existing SIGNAL_PEPTIDES dict.
    assert mrna_library.SIGNAL_PEPTIDES["CD8A"] == "MALPVTALLLPLALLLHAARP"
    assert mrna_library.SIGNAL_PEPTIDES["CD28"] == "MLRLLLALNLFPSIQVTG"
