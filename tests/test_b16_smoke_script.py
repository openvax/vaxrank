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

from pathlib import Path
import shlex

from vaxrank.cli.arg_parser import parse_vaxrank_args


def test_b16_smoke_script_uses_current_self_contained_cli_contract():
    repo_root = Path(__file__).resolve().parents[1]
    script = (repo_root / "run-vaxrank-b16-test-data.sh").read_text()
    assert script.startswith("#!/usr/bin/env bash\n")
    command = script[script.index("vaxrank "):].replace("\\\n", " ")
    parsed = parse_vaxrank_args(shlex.split(command)[1:])

    assert all((repo_root / path).is_file() for path in parsed.vcf)
    assert (repo_root / parsed.bam).is_file()
    assert parsed.mhc_predictor == [["random"]]
    assert parsed.num_epitopes_per_vaccine_peptide == 5
