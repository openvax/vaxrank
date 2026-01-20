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
Tests for the ReferenceProteome set-based kmer index.
"""

import os
import pickle
import tempfile
from unittest.mock import MagicMock, patch

from vaxrank.reference_proteome import (
    ReferenceProteome,
    build_kmer_set_index,
    load_kmer_set_index,
    kmer_set_index_path,
    get_cache_dir,
    DEFAULT_MIN_KMER_LENGTH,
    DEFAULT_MAX_KMER_LENGTH,
    _kmer_set_cache,
)

from .common import eq_, ok_


# =============================================================================
# Mock Genome and Transcript for Testing
# =============================================================================

def create_mock_transcript(transcript_id, protein_sequence, is_protein_coding=True):
    """Create a mock transcript object"""
    transcript = MagicMock()
    transcript.id = transcript_id
    transcript.protein_sequence = protein_sequence
    transcript.is_protein_coding = is_protein_coding
    return transcript


def create_mock_genome(transcripts, species_name="test_species", release=100):
    """Create a mock genome object"""
    genome = MagicMock()
    genome.transcripts.return_value = transcripts
    genome.species = MagicMock()
    genome.species.latin_name = species_name
    genome.release = release
    return genome


# =============================================================================
# ReferenceProteome Basic Tests
# =============================================================================

class TestReferenceProteomeBasic:
    """Basic tests for ReferenceProteome class"""

    def test_none_genome(self):
        """Test ReferenceProteome with None genome"""
        ref = ReferenceProteome(None)
        ok_(not ref.contains("ANYSEQ"))
        ok_(not ref.contains("MADEUP"))
        eq_(len(ref._kmer_set), 0)

    def test_contains_basic(self):
        """Test basic contains functionality with mock genome"""
        # Create a simple protein: "ABCDEFGHIJKLMNOP"
        transcripts = [
            create_mock_transcript("T1", "ABCDEFGHIJKLMNOP"),
        ]
        genome = create_mock_genome(transcripts)

        with patch('vaxrank.reference_proteome.load_kmer_set_index') as mock_load:
            # Manually create the kmer set for this protein
            kmers = set()
            protein = "ABCDEFGHIJKLMNOP"
            for k in range(8, 16):
                for i in range(len(protein) - k + 1):
                    kmers.add(protein[i:i+k])
            mock_load.return_value = kmers

            ref = ReferenceProteome(genome)

            # 8-mers that should exist
            ok_(ref.contains("ABCDEFGH"))
            ok_(ref.contains("BCDEFGHI"))
            ok_(ref.contains("IJKLMNOP"))

            # 9-mers that should exist
            ok_(ref.contains("ABCDEFGHI"))
            ok_(ref.contains("HIJKLMNOP"))

            # Should not exist (not in protein)
            ok_(not ref.contains("ZZZZZZZZ"))
            ok_(not ref.contains("QRSTUVWX"))


class TestBuildKmerSetIndex:
    """Tests for build_kmer_set_index function"""

    def test_single_protein(self):
        """Test building kmer index from single protein"""
        transcripts = [
            create_mock_transcript("T1", "ABCDEFGHIJKLMNOP"),
        ]
        genome = create_mock_genome(transcripts)

        kmers = build_kmer_set_index(genome, min_len=8, max_len=10)

        # Check some expected 8-mers
        ok_("ABCDEFGH" in kmers)
        ok_("BCDEFGHI" in kmers)
        ok_("IJKLMNOP" in kmers)

        # Check some expected 9-mers
        ok_("ABCDEFGHI" in kmers)
        ok_("HIJKLMNOP" in kmers)

        # Check some expected 10-mers
        ok_("ABCDEFGHIJ" in kmers)
        ok_("GHIJKLMNOP" in kmers)

        # 11-mers should NOT be included (max_len=10)
        ok_("ABCDEFGHIJK" not in kmers)

    def test_multiple_proteins(self):
        """Test building kmer index from multiple proteins"""
        transcripts = [
            create_mock_transcript("T1", "AAAAAAAA"),  # Simple protein
            create_mock_transcript("T2", "BBBBBBBB"),  # Another simple protein
        ]
        genome = create_mock_genome(transcripts)

        kmers = build_kmer_set_index(genome, min_len=8, max_len=8)

        ok_("AAAAAAAA" in kmers)
        ok_("BBBBBBBB" in kmers)
        ok_("CCCCCCCC" not in kmers)

    def test_skip_non_coding_transcripts(self):
        """Test that non-coding transcripts are skipped"""
        transcripts = [
            create_mock_transcript("T1", "ABCDEFGH", is_protein_coding=True),
            create_mock_transcript("T2", "ZZZZZZZZ", is_protein_coding=False),
        ]
        genome = create_mock_genome(transcripts)

        kmers = build_kmer_set_index(genome, min_len=8, max_len=8)

        ok_("ABCDEFGH" in kmers)
        ok_("ZZZZZZZZ" not in kmers)  # Non-coding, should be skipped

    def test_skip_none_protein_sequence(self):
        """Test that transcripts with None protein sequence are skipped"""
        transcripts = [
            create_mock_transcript("T1", "ABCDEFGH"),
            create_mock_transcript("T2", None),  # No protein sequence
        ]
        genome = create_mock_genome(transcripts)

        kmers = build_kmer_set_index(genome, min_len=8, max_len=8)

        ok_("ABCDEFGH" in kmers)
        # Should not raise error, just skip the None protein

    def test_short_protein(self):
        """Test handling protein shorter than min kmer length"""
        transcripts = [
            create_mock_transcript("T1", "ABCDE"),  # Only 5 aa, shorter than min_len=8
        ]
        genome = create_mock_genome(transcripts)

        kmers = build_kmer_set_index(genome, min_len=8, max_len=10)

        # Should be empty since protein is too short
        eq_(len(kmers), 0)

    def test_custom_kmer_lengths(self):
        """Test building index with custom kmer length range"""
        transcripts = [
            create_mock_transcript("T1", "ABCDEFGHIJKLMNOPQRST"),  # 20 aa
        ]
        genome = create_mock_genome(transcripts)

        # Only 5-mers
        kmers = build_kmer_set_index(genome, min_len=5, max_len=5)

        ok_("ABCDE" in kmers)
        ok_("FGHIJ" in kmers)
        ok_("ABCDEF" not in kmers)  # 6-mer, not included

    def test_overlapping_kmers(self):
        """Test that all overlapping kmers are captured"""
        transcripts = [
            create_mock_transcript("T1", "ABCDEFGHIJ"),  # 10 aa
        ]
        genome = create_mock_genome(transcripts)

        kmers = build_kmer_set_index(genome, min_len=8, max_len=8)

        # There should be 3 8-mers: positions 0-7, 1-8, 2-9
        ok_("ABCDEFGH" in kmers)
        ok_("BCDEFGHI" in kmers)
        ok_("CDEFGHIJ" in kmers)
        eq_(len(kmers), 3)


class TestLoadKmerSetIndex:
    """Tests for load_kmer_set_index function with caching"""

    def test_builds_index_when_not_cached(self):
        """Test that index is built when cache doesn't exist"""
        transcripts = [
            create_mock_transcript("T1", "ABCDEFGHIJ"),
        ]
        genome = create_mock_genome(transcripts, species_name="test_species", release=100)

        with tempfile.TemporaryDirectory() as tmpdir:
            with patch('vaxrank.reference_proteome.get_cache_dir', return_value=tmpdir):
                kmers = load_kmer_set_index(genome, min_len=8, max_len=8)

                # Should contain expected kmers
                ok_("ABCDEFGH" in kmers)

                # Cache file should now exist
                cache_path = os.path.join(tmpdir, "test_species_100_kmer_set_8_8.pkl")
                ok_(os.path.exists(cache_path))

    def test_loads_from_cache_when_exists(self):
        """Test that index is loaded from disk cache when it exists"""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Pre-create a cache file
            cached_kmers = {"CACHED01", "CACHED02"}
            cache_path = os.path.join(tmpdir, "test_species_100_kmer_set_8_8.pkl")
            with open(cache_path, 'wb') as f:
                pickle.dump(cached_kmers, f)

            genome = create_mock_genome([], species_name="test_species", release=100)

            # Clear in-memory cache to test disk cache loading
            cache_key = ("test_species", 100, 8, 8)
            _kmer_set_cache.pop(cache_key, None)

            with patch('vaxrank.reference_proteome.get_cache_dir', return_value=tmpdir):
                # Should load from disk cache, not build new index
                kmers = load_kmer_set_index(genome, min_len=8, max_len=8)

                ok_("CACHED01" in kmers)
                ok_("CACHED02" in kmers)

    def test_force_reload_rebuilds_index(self):
        """Test that force_reload=True rebuilds index even if cached"""
        transcripts = [
            create_mock_transcript("T1", "NEWPROTEIN"),  # New protein
        ]
        genome = create_mock_genome(transcripts, species_name="test_species", release=100)

        with tempfile.TemporaryDirectory() as tmpdir:
            # Pre-create a cache file with different content
            cached_kmers = {"OLDKMERS"}
            cache_path = os.path.join(tmpdir, "test_species_100_kmer_set_8_8.pkl")
            with open(cache_path, 'wb') as f:
                pickle.dump(cached_kmers, f)

            with patch('vaxrank.reference_proteome.get_cache_dir', return_value=tmpdir):
                # With force_reload, should rebuild from genome
                kmers = load_kmer_set_index(genome, min_len=8, max_len=8, force_reload=True)

                # Should have new kmers, not old cached ones
                ok_("OLDKMERS" not in kmers)
                ok_("NEWPROTE" in kmers)  # 8-mer from "NEWPROTEIN"


class TestKmerSetIndexPath:
    """Tests for kmer_set_index_path function"""

    def test_path_format(self):
        """Test that path is formatted correctly"""
        genome = create_mock_genome([], species_name="homo_sapiens", release=104)

        with patch('vaxrank.reference_proteome.get_cache_dir', return_value="/cache"):
            path = kmer_set_index_path(genome, min_len=8, max_len=15)

            eq_(path, "/cache/homo_sapiens_104_kmer_set_8_15.pkl")

    def test_different_kmer_lengths(self):
        """Test that different kmer lengths produce different paths"""
        genome = create_mock_genome([], species_name="test", release=1)

        with patch('vaxrank.reference_proteome.get_cache_dir', return_value="/cache"):
            path1 = kmer_set_index_path(genome, min_len=8, max_len=15)
            path2 = kmer_set_index_path(genome, min_len=9, max_len=12)

            ok_(path1 != path2)
            ok_("8_15" in path1)
            ok_("9_12" in path2)


class TestReferenceProteomeContains:
    """Tests for ReferenceProteome.contains() method"""

    def test_contains_exact_kmer(self):
        """Test contains with exact kmer match"""
        with patch('vaxrank.reference_proteome.load_kmer_set_index') as mock_load:
            mock_load.return_value = {"PEPTIDE1", "PEPTIDE2", "PEPTIDE3"}

            genome = create_mock_genome([])
            ref = ReferenceProteome(genome)

            ok_(ref.contains("PEPTIDE1"))
            ok_(ref.contains("PEPTIDE2"))
            ok_(not ref.contains("NOTFOUND"))

    def test_contains_empty_string(self):
        """Test contains with empty string"""
        with patch('vaxrank.reference_proteome.load_kmer_set_index') as mock_load:
            mock_load.return_value = {"PEPTIDE1"}

            genome = create_mock_genome([])
            ref = ReferenceProteome(genome)

            ok_(not ref.contains(""))

    def test_contains_case_sensitive(self):
        """Test that contains is case-sensitive"""
        with patch('vaxrank.reference_proteome.load_kmer_set_index') as mock_load:
            mock_load.return_value = {"PEPTIDE1"}

            genome = create_mock_genome([])
            ref = ReferenceProteome(genome)

            ok_(ref.contains("PEPTIDE1"))
            ok_(not ref.contains("peptide1"))
            ok_(not ref.contains("Peptide1"))


class TestReferenceProteomeCustomKmerLengths:
    """Tests for ReferenceProteome with custom kmer lengths"""

    def test_custom_min_max_kmer_length(self):
        """Test ReferenceProteome with custom min/max kmer lengths"""
        transcripts = [
            create_mock_transcript("T1", "ABCDEFGHIJKLMNOP"),
        ]
        genome = create_mock_genome(transcripts)

        with tempfile.TemporaryDirectory() as tmpdir:
            with patch('vaxrank.reference_proteome.get_cache_dir', return_value=tmpdir):
                # Only index 10-12 mers
                ref = ReferenceProteome(genome, min_kmer_length=10, max_kmer_length=12)

                # 10-mer should be found
                ok_(ref.contains("ABCDEFGHIJ"))

                # 8-mer should NOT be found (below min length)
                # Note: contains() just checks set membership, it doesn't validate length
                ok_(not ref.contains("ABCDEFGH"))

    def test_kmer_length_attributes(self):
        """Test that kmer length attributes are set correctly"""
        with patch('vaxrank.reference_proteome.load_kmer_set_index') as mock_load:
            mock_load.return_value = set()

            genome = create_mock_genome([])
            ref = ReferenceProteome(genome, min_kmer_length=9, max_kmer_length=14)

            eq_(ref.min_kmer_length, 9)
            eq_(ref.max_kmer_length, 14)


class TestReferenceProteomeDefaults:
    """Tests for default values"""

    def test_default_kmer_lengths(self):
        """Test that default kmer lengths are correct"""
        eq_(DEFAULT_MIN_KMER_LENGTH, 8)
        eq_(DEFAULT_MAX_KMER_LENGTH, 15)

    def test_default_kmer_lengths_used(self):
        """Test that defaults are used when not specified"""
        with patch('vaxrank.reference_proteome.load_kmer_set_index') as mock_load:
            mock_load.return_value = set()

            genome = create_mock_genome([])
            ref = ReferenceProteome(genome)

            eq_(ref.min_kmer_length, DEFAULT_MIN_KMER_LENGTH)
            eq_(ref.max_kmer_length, DEFAULT_MAX_KMER_LENGTH)


class TestKmerSetIndexUniqueness:
    """Tests for kmer set uniqueness and deduplication"""

    def test_duplicate_kmers_deduplicated(self):
        """Test that duplicate kmers from multiple transcripts are deduplicated"""
        # Two transcripts with overlapping sequences
        transcripts = [
            create_mock_transcript("T1", "ABCDEFGH"),
            create_mock_transcript("T2", "ABCDEFGH"),  # Same protein
        ]
        genome = create_mock_genome(transcripts)

        kmers = build_kmer_set_index(genome, min_len=8, max_len=8)

        # Should only have one copy of the kmer
        eq_(len(kmers), 1)
        ok_("ABCDEFGH" in kmers)

    def test_partial_overlap_both_captured(self):
        """Test that partially overlapping proteins both contribute kmers"""
        transcripts = [
            create_mock_transcript("T1", "ABCDEFGHIJ"),  # Ends with CDEFGHIJ
            create_mock_transcript("T2", "CDEFGHIJKL"),  # Starts with CDEFGHIJ, adds KL
        ]
        genome = create_mock_genome(transcripts)

        kmers = build_kmer_set_index(genome, min_len=8, max_len=8)

        # From T1
        ok_("ABCDEFGH" in kmers)
        # Shared
        ok_("CDEFGHIJ" in kmers)
        # From T2
        ok_("DEFGHIJK" in kmers)
        ok_("EFGHIJKL" in kmers)


class TestGetCacheDir:
    """Tests for get_cache_dir function"""

    def test_cache_dir_created(self):
        """Test that cache directory is created if it doesn't exist"""
        with tempfile.TemporaryDirectory() as tmpdir:
            new_cache_dir = os.path.join(tmpdir, "new_cache")
            ok_(not os.path.exists(new_cache_dir))

            with patch('vaxrank.reference_proteome.get_data_dir', return_value=new_cache_dir):
                result = get_cache_dir()
                eq_(result, new_cache_dir)
                ok_(os.path.exists(new_cache_dir))
