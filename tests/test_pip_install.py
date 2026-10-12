"""Execute CI retries with a controlled pip process and unchanged constraints."""

import json
import os
from pathlib import Path
import subprocess
import sys

import pytest
import yaml


ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def pip_process(tmp_path):
    python = tmp_path / "selected python"
    log = tmp_path / "calls.jsonl"
    python.write_text(
        "#!" + sys.executable + "\n"
        "import json, os, sys\n"
        "from pathlib import Path\n"
        "p = Path(os.environ['INSTALL_LOG'])\n"
        "calls = [json.loads(line) for line in p.read_text().splitlines()] if p.exists() else []\n"
        "attempt = len(calls)\n"
        "with p.open('a') as f: f.write(json.dumps(sys.argv[1:]) + '\\n')\n"
        "target = os.environ.get('FAIL_REQUIREMENT')\n"
        "if target:\n"
        "    sys.exit(17 if target in sys.argv[1:] and not any(target in call for call in calls) else 0)\n"
        "sys.exit(int(os.environ['FAIL_STATUS']) if attempt < int(os.environ['FAIL_COUNT']) else 0)\n")
    python.chmod(0o755)
    # Record backoff instead of spending 20 + 40 seconds in each offline regression.
    sleep = tmp_path / "sleep"
    sleep.write_text('#!/bin/sh\nprintf "%s\\n" "$1" >> "$SLEEP_LOG"\n')
    sleep.chmod(0o755)
    env = dict(os.environ, PYTHON=str(python), INSTALL_LOG=str(log),
               SLEEP_LOG=str(tmp_path / "sleeps"), FAIL_COUNT="0", FAIL_STATUS="17",
               PIP_INSTALL_ATTEMPTS="3", PIP_INSTALL_RETRY_DELAY="20",
               PATH=str(tmp_path) + os.pathsep + os.environ["PATH"])
    return env, log


@pytest.mark.parametrize("failures", [0, 1, 2, 3, 4])
def test_install_retries_with_identical_requirements(pip_process, failures):
    env, log = pip_process
    env["FAIL_COUNT"] = str(failures)
    args = ["-r", "requirements with spaces.txt", ".[hitlist]", "msgspec==0.18.6", "topiary>=5.94.4,<6"]
    result = subprocess.run(["bash", str(ROOT / "scripts/pip_install.sh"), *args],
                            env=env, capture_output=True, text=True)
    calls = [json.loads(line) for line in log.read_text().splitlines()]
    assert len(calls) == min(failures + 1, 3)
    assert calls[0] == ["-m", "pip", "install", *args]
    assert all(call == ["-m", "pip", "install", "--no-cache-dir", *args] for call in calls[1:])
    sleeps = Path(env["SLEEP_LOG"])
    assert (sleeps.read_text().splitlines() if sleeps.exists() else []) == ["20", "40"][:len(calls) - 1]
    assert result.returncode == (0 if failures < 3 else 17)


@pytest.mark.parametrize("status", [126, 127, 130, 143])
def test_install_does_not_retry_interrupts_or_bad_interpreters(pip_process, status):
    env, log = pip_process
    env.update(FAIL_COUNT="3", FAIL_STATUS=str(status))
    result = subprocess.run(["bash", str(ROOT / "scripts/pip_install.sh"), "."], env=env)
    assert result.returncode == status
    assert len(log.read_text().splitlines()) == 1
    assert not Path(env["SLEEP_LOG"]).exists()


@pytest.mark.parametrize("name,value", [("PIP_INSTALL_ATTEMPTS", "0"),
                                       ("PIP_INSTALL_RETRY_DELAY", "-1"),
                                       ("PIP_INSTALL_ATTEMPTS", "oops")])
def test_invalid_retry_configuration_fails_before_install(pip_process, name, value):
    env, log = pip_process
    env[name] = value
    result = subprocess.run(["bash", str(ROOT / "scripts/pip_install.sh"), "."], env=env)
    assert result.returncode == 2 and not log.exists()


@pytest.mark.parametrize("job", ["build", "supported-msgspec"])
def test_workflow_routes_dependency_installs_through_retry_script(tmp_path, pip_process, job):
    env, log = pip_process
    workflow = yaml.safe_load((ROOT / ".github/workflows/tests.yml").read_text())
    steps = workflow["jobs"][job]["steps"]
    name = "Install dependencies" if job == "build" else "Install supported msgspec and application dependencies"
    command = next(step["run"] for step in steps if step.get("name") == name)
    command = command.replace("${{ matrix.msgspec-version }}", "0.18.6")
    env["FAIL_REQUIREMENT"] = "-r" if job == "build" else "."
    # Tooling-only installs also call this interpreter, but must remain ordinary installs.
    (tmp_path / "python").symlink_to(env["PYTHON"])
    result = subprocess.run(["bash", "-eo", "pipefail", "-c", command],
                            cwd=ROOT, env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    calls = [json.loads(line) for line in log.read_text().splitlines()]
    expected = ([["-m", "pip", "install", "--upgrade", "pip"],
                 ["-m", "pip", "install", "pytest", "pytest-cov", "ruff"],
                 ["-m", "pip", "install", "-r", "requirements.txt"],
                 ["-m", "pip", "install", "--no-cache-dir", "-r", "requirements.txt"],
                 ["-m", "pip", "install", ".[hitlist]"]]
                if job == "build" else [["-m", "pip", "install", ".", "pytest", "msgspec==0.18.6"],
                                         ["-m", "pip", "install", "--no-cache-dir", ".", "pytest", "msgspec==0.18.6"]])
    assert calls == expected
