"""PDF cleanup must reap only its owned display and preserve render errors."""

import subprocess
from unittest.mock import MagicMock, call, patch

import pytest

from vaxrank.report import _pdf_via_pdfkit


@pytest.mark.parametrize("render_fails", [False, True])
@pytest.mark.parametrize("stop_times_out", [False, True])
def test_owned_display_cleanup_and_render_error_propagation(render_fails, stop_times_out):
    display = MagicMock()
    display.__enter__.side_effect = display.start
    display.__exit__.side_effect = lambda *args: display.stop()
    process = display.proc
    process.args = ["Xvfb", ":1234"]
    if stop_times_out:
        display.stop.side_effect = subprocess.TimeoutExpired(process.args, 10)
    render_error = ValueError("render failed")
    with patch("vaxrank.report.sys.platform", "linux"), \
            patch("xvfbwrapper.Xvfb", return_value=display), \
            patch("pdfkit.from_file", side_effect=render_error if render_fails else None) as render:
        if render_fails:
            with pytest.raises(ValueError, match="render failed"):
                _pdf_via_pdfkit("input.html", "output.pdf")
        else:
            _pdf_via_pdfkit("input.html", "output.pdf")
    display.start.assert_called_once_with()
    display.stop.assert_called_once_with()
    render.assert_called_once()
    if stop_times_out:
        process.kill.assert_called_once_with()
        process.wait.assert_called_once_with(timeout=10)
        assert display.proc is None
    else:
        process.kill.assert_not_called()


def test_failed_forced_reap_is_not_reported_as_success():
    display = MagicMock()
    display.proc.args = ["Xvfb", ":1234"]
    timeout = subprocess.TimeoutExpired(display.proc.args, 10)
    display.stop.side_effect = timeout
    display.proc.wait.side_effect = timeout
    with patch("vaxrank.report.sys.platform", "linux"), \
            patch("xvfbwrapper.Xvfb", return_value=display), patch("pdfkit.from_file"):
        with pytest.raises(subprocess.TimeoutExpired):
            _pdf_via_pdfkit("input.html", "output.pdf")
    display.proc.kill.assert_called_once_with()


def test_display_start_failure_still_cleans_up_owned_child():
    display = MagicMock()
    display.start.side_effect = RuntimeError("display unavailable")
    with patch("vaxrank.report.sys.platform", "linux"), \
            patch("xvfbwrapper.Xvfb", return_value=display), patch("pdfkit.from_file") as render:
        with pytest.raises(RuntimeError, match="display unavailable"):
            _pdf_via_pdfkit("input.html", "output.pdf")
    display.stop.assert_called_once_with()
    render.assert_not_called()


@pytest.mark.parametrize("original_display", [None, ":existing-user-display"])
def test_real_wrapper_restores_environment_and_reaps_timed_out_child(tmp_path, original_display):
    from xvfbwrapper import Xvfb

    environ = {"TEST_DISPLAY_ENV": "isolated"}
    if original_display is not None:
        environ["DISPLAY"] = original_display
    with patch.object(Xvfb, "_xvfb_exists", return_value=True):
        display = Xvfb(environ=environ)
    environ["DISPLAY"] = ":owned-pdf-display"
    display._lock_display_file = (tmp_path / "owned-display.lock").open("w")
    process = display.proc = MagicMock()
    process.args = ["Xvfb", ":owned-pdf-display"]
    process.wait.side_effect = [subprocess.TimeoutExpired(process.args, 10), 0]
    with patch("vaxrank.report.sys.platform", "linux"), \
            patch("xvfbwrapper.Xvfb", return_value=display), \
            patch.object(display, "start"), patch("pdfkit.from_file"):
        _pdf_via_pdfkit("input.html", "output.pdf")
    process.terminate.assert_called_once_with()
    process.kill.assert_called_once_with()
    assert process.wait.call_args_list == [call(10), call(timeout=10)]
    assert display.proc is None
    assert environ.get("DISPLAY") == original_display
    assert display._lock_display_file.closed
    assert not (tmp_path / "owned-display.lock").exists()
