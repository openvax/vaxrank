"""Tests for bounded CI wkhtmltopdf installation (openvax/vaxrank#806)."""

import os
from pathlib import Path
import stat
import subprocess
import sys


REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT_PATH = REPO_ROOT / "scripts" / "install_wkhtmltopdf.sh"
WORKFLOW_PATH = REPO_ROOT / ".github" / "workflows" / "tests.yml"


def _write_executable(path: Path, body: str) -> Path:
    path.write_text(body, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return path


def _base_env(tmp_path: Path, *, apt_script: Path, wget_script: Path, dpkg_script: Path):
    env = dict(os.environ)
    env.update({
        "PYTHON": sys.executable,
        "WKHTMLTOPDF_APT_GET_BIN": str(apt_script),
        "WKHTMLTOPDF_WGET_BIN": str(wget_script),
        "WKHTMLTOPDF_DPKG_BIN": str(dpkg_script),
        "WKHTMLTOPDF_DEB_PATH": str(tmp_path / "wkhtmltox.deb"),
    })
    return env


def test_install_wkhtmltopdf_passes_bounded_flags_and_pinned_release(tmp_path):
    """Happy path preserves pinned 0.12.6-1 package and passes bounded timeouts/retries."""
    log_file = tmp_path / "invocations.log"
    stub_body = f"#!/usr/bin/env bash\necho \"$0 $*\" >> {str(log_file)!r}\n"
    apt_stub = _write_executable(tmp_path / "apt-get", stub_body)
    wget_stub = _write_executable(tmp_path / "wget", stub_body)
    dpkg_stub = _write_executable(tmp_path / "dpkg", stub_body)

    proc = subprocess.run(
        ["bash", str(SCRIPT_PATH)],
        env=_base_env(tmp_path, apt_script=apt_stub, wget_script=wget_stub, dpkg_script=dpkg_stub),
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr
    lines = log_file.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 4
    assert "Acquire::Retries=3" in lines[0] and "update" in lines[0]
    assert "xfonts-base xfonts-75dpi" in lines[1]
    assert "--tries=3" in lines[2]
    assert "--connect-timeout=20" in lines[2]
    assert "--read-timeout=30" in lines[2]
    assert "--waitretry=5" in lines[2]
    assert "0.12.6-1/wkhtmltox_0.12.6-1.bionic_amd64.deb" in lines[2]
    assert "-i " in lines[3]


def test_install_wkhtmltopdf_timeout_produces_clear_failure(tmp_path):
    """A stalled download exits 124 with a clear diagnostic instead of hanging."""
    ok_stub = _write_executable(tmp_path / "ok", "#!/usr/bin/env bash\nexit 0\n")
    slow_wget = _write_executable(
        tmp_path / "slow-wget",
        "#!/usr/bin/env bash\nsleep 10\n",
    )
    env = _base_env(tmp_path, apt_script=ok_stub, wget_script=slow_wget, dpkg_script=ok_stub)
    env["WKHTMLTOPDF_DOWNLOAD_TIMEOUT"] = "1"

    proc = subprocess.run(
        ["bash", str(SCRIPT_PATH)],
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 124
    assert "timed out after 1s" in proc.stderr
    assert "wget wkhtmltox package" in proc.stderr


def test_install_wkhtmltopdf_propagates_nonzero_error(tmp_path):
    """Non-zero failure from dpkg or wget is never swallowed."""
    ok_stub = _write_executable(tmp_path / "ok", "#!/usr/bin/env bash\nexit 0\n")
    fail_dpkg = _write_executable(tmp_path / "fail-dpkg", "#!/usr/bin/env bash\nexit 23\n")
    env = _base_env(tmp_path, apt_script=ok_stub, wget_script=ok_stub, dpkg_script=fail_dpkg)

    proc = subprocess.run(
        ["bash", str(SCRIPT_PATH)],
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 23
    assert "dpkg install wkhtmltox" in proc.stderr


def test_ci_workflow_bounds_wkhtmltopdf_step_deadline():
    """Workflow declares an explicit step timeout and delegates to install_wkhtmltopdf.sh."""
    workflow = WORKFLOW_PATH.read_text(encoding="utf-8")
    assert "- name: Install wkhtmltopdf\n        timeout-minutes: 5\n        run: bash scripts/install_wkhtmltopdf.sh" in workflow
