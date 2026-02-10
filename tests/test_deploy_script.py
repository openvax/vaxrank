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

import os
import subprocess
from pathlib import Path


def test_deploy_rejects_non_release_branch():
    repo_root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env["VAXRANK_DEPLOY_BRANCH"] = "yaml-config"
    result = subprocess.run(
        ["bash", "deploy.sh", "--dry-run"],
        cwd=repo_root,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 1
    assert "only allowed from main or master" in result.stderr
