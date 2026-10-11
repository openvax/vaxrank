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
Pytest configuration for vaxrank tests.
"""

from pathlib import Path

import pytest
import shutil


@pytest.fixture
def platform_cache_home(tmp_path, monkeypatch):
    """Move datacache's platform cache directory under tmp_path; return the openvax root.

    HOME decides the Linux (~/.cache) and macOS (~/Library/Caches) defaults
    whichever library datacache uses to find them (appdirs before 1.25,
    platformdirs after), so nothing inside datacache is patched.
    """
    from datacache import get_data_dir
    from vaxrank.download_test_data import CACHE_NAMESPACE

    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.delenv("XDG_CACHE_HOME", raising=False)
    root = Path(get_data_dir(CACHE_NAMESPACE))
    assert root.is_relative_to(tmp_path), root  # Never a real cache.
    return root


@pytest.fixture(autouse=True)
def isolate_shared_read_cache(request, monkeypatch):
    """Tests declare their cache environment unless they use the real bundle."""
    if request.node.get_closest_marker("shared_read_cache") is None:
        for name in ("OSTEOSARC_CACHE", "OPENVAX_DATA_CACHE"):
            monkeypatch.delenv(name, raising=False)


@pytest.fixture(scope="session")
def mouse_genome():
    """Lazy session-scoped pyensembl GRCm38 handle.

    Module-level ``genome_for_reference_name(...)`` runs at *collection*
    time — paid by every xdist worker before any test runs, even
    workers that won't touch the genome. Deferring to a fixture defers
    the load until something actually asks; scope='session' shares the
    handle across tests within one worker process.
    """
    from pyensembl import genome_for_reference_name
    return genome_for_reference_name("GRCm38")


@pytest.fixture(scope="session")
def human_genome_grch37():
    """Lazy session-scoped pyensembl GRCh37 (Ensembl release 75) handle."""
    from pyensembl import EnsemblRelease
    return EnsemblRelease(75)


def pytest_configure(config):
    """Register custom markers."""
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (use -m 'not slow' to skip)"
    )
    config.addinivalue_line(
        "markers", "requires_netmhcpan: marks tests that require NetMHCpan"
    )
    config.addinivalue_line(
        "markers", "shared_read_cache: use the caller's configured openvax-v1 cache"
    )


def netmhcpan_available():
    """Check if NetMHCpan is available on the system."""
    return shutil.which("netMHCpan") is not None


# Skip condition for tests requiring NetMHCpan
requires_netmhcpan = pytest.mark.skipif(
    not netmhcpan_available(),
    reason="NetMHCpan not found in PATH"
)
