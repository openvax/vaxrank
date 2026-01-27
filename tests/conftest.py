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

import pytest
import shutil


def pytest_configure(config):
    """Register custom markers."""
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (use -m 'not slow' to skip)"
    )
    config.addinivalue_line(
        "markers", "requires_netmhcpan: marks tests that require NetMHCpan"
    )


def netmhcpan_available():
    """Check if NetMHCpan is available on the system."""
    return shutil.which("netMHCpan") is not None


# Skip condition for tests requiring NetMHCpan
requires_netmhcpan = pytest.mark.skipif(
    not netmhcpan_available(),
    reason="NetMHCpan not found in PATH"
)
