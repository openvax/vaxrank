"""Exercise Vaxrank's PDF path through the actual supported Xvfb wrapper."""

import subprocess
from unittest.mock import MagicMock, call, patch

import pytest
from xvfbwrapper import Xvfb

from vaxrank.report import _pdf_via_pdfkit


def owned_display(tmp_path, original_display):
    environ = {"TEST_DISPLAY_ENV": "isolated"}
    if original_display is not None:
        environ["DISPLAY"] = original_display
    with patch.object(Xvfb, "_xvfb_exists", return_value=True):
        display = Xvfb(environ=environ)
    environ["DISPLAY"] = ":owned-pdf-display"
    display._lock_display_file = (tmp_path / "owned-display.lock").open("w")
    display.proc = MagicMock(name="owned_xvfb_process")
    display.proc.args = ["Xvfb", ":owned-pdf-display"]
    return display, environ


@pytest.mark.parametrize("original_display", [None, ":existing-user-display"])
@pytest.mark.parametrize("stop_times_out", [False, True])
@pytest.mark.parametrize("render_fails", [False, True])
def test_pdf_cleanup_reaps_owned_display_and_preserves_render_error(
        tmp_path, original_display, stop_times_out, render_fails):
    display, environ = owned_display(tmp_path, original_display)
    process = display.proc
    if stop_times_out:
        process.wait.side_effect = [subprocess.TimeoutExpired(process.args, 10), 0]
    render_error = ValueError("render failed")
    with patch("vaxrank.report.sys.platform", "linux"), \
            patch("xvfbwrapper.Xvfb", return_value=display), \
            patch.object(display, "start") as start, \
            patch("pdfkit.from_file", side_effect=render_error if render_fails else None) as render:
        if render_fails:
            with pytest.raises(ValueError, match="render failed") as raised:
                _pdf_via_pdfkit("input.html", "output.pdf")
            assert raised.value is render_error
        else:
            _pdf_via_pdfkit("input.html", "output.pdf")
    start.assert_called_once_with()
    render.assert_called_once()
    process.terminate.assert_called_once_with()
    if stop_times_out:
        process.kill.assert_called_once_with()
        assert process.wait.call_args_list == [call(10), call(10)]
    else:
        process.kill.assert_not_called()
        process.wait.assert_called_once_with(10)
    assert display.proc is None
    assert environ.get("DISPLAY") == original_display
    assert display._lock_display_file.closed
    assert not (tmp_path / "owned-display.lock").exists()


@pytest.mark.parametrize("render_fails", [False, True])
def test_pdf_failed_forced_reap_is_not_reported_as_success(tmp_path, render_fails):
    display, _ = owned_display(tmp_path, None)
    process = display.proc
    graceful_error = subprocess.TimeoutExpired(process.args, 10)
    forced_error = subprocess.TimeoutExpired(process.args, 10)
    process.wait.side_effect = [graceful_error, forced_error]
    render_error = ValueError("render failed")
    with patch("vaxrank.report.sys.platform", "linux"), \
            patch("xvfbwrapper.Xvfb", return_value=display), \
            patch.object(display, "start"), \
            patch("pdfkit.from_file", side_effect=render_error if render_fails else None):
        with pytest.raises(subprocess.TimeoutExpired) as raised:
            _pdf_via_pdfkit("input.html", "output.pdf")
    assert raised.value is forced_error
    process.kill.assert_called_once_with()
    assert process.wait.call_args_list == [call(10), call(10)]
    if render_fails:
        assert forced_error.__context__ is graceful_error
        assert graceful_error.__context__ is render_error
