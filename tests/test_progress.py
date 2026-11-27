"""Status lines and progress bars for the slow steps of a run."""

import logging
import sys

from vaxrank.progress import progress, timed_step


def test_bars_are_off_when_stderr_is_not_a_terminal(capsys):
    """Redirected and CI logs get the INFO status lines, not carriage-return spam (#460)."""
    assert list(progress(range(3), "demo", "item")) == [0, 1, 2]
    assert capsys.readouterr().err == ""


def test_bars_draw_on_a_terminal(monkeypatch):
    class Terminal:
        def __init__(self):
            self.written = []

        def isatty(self):
            return True

        def write(self, text):
            self.written.append(text)

        def flush(self):
            pass

    terminal = Terminal()
    monkeypatch.setattr(sys, "stderr", terminal)
    assert list(progress(range(3), "demo bar", "item", ncols=80)) == [0, 1, 2]
    assert any("demo bar" in text for text in terminal.written)


def test_timed_step_logs_start_and_elapsed_time(caplog):
    logger = logging.getLogger("vaxrank.test_progress")
    with caplog.at_level(logging.INFO, logger=logger.name):
        with timed_step(logger, "Doing %d things", "Did %d things in %.0fs", 4):
            pass
    messages = [record.getMessage() for record in caplog.records]
    assert messages[0] == "Doing 4 things"
    assert messages[1].startswith("Did 4 things in ") and messages[1].endswith("s")
    # Reported from the caller's line, not from inside the helper.
    assert {record.funcName for record in caplog.records} == {"test_timed_step_logs_start_and_elapsed_time"}


def test_bars_are_enabled_in_jupyter_even_when_stderr_is_not_a_tty(monkeypatch):
    """#783: Jupyter notebooks have non-TTY stderr but should still enable bars."""
    from types import SimpleNamespace
    import tqdm.auto

    class ZMQInteractiveShell:
        pass

    captured = {}

    def fake_auto_tqdm(iterable, **kwargs):
        captured.update(kwargs)
        return list(iterable)

    monkeypatch.setitem(sys.modules, "IPython", SimpleNamespace(get_ipython=lambda: ZMQInteractiveShell()))
    monkeypatch.setattr(tqdm.auto, "tqdm", fake_auto_tqdm)
    assert progress(range(3), "notebook bar", "item") == [0, 1, 2]
    assert captured["disable"] is False
    assert captured["desc"] == "notebook bar"


def test_make_pdf_report_validates_backend_before_logging_and_reports_elapsed_time(tmp_path, caplog, monkeypatch):
    """#783: invalid PDF backend fails before logging; valid backend logs elapsed time."""
    import pytest
    from vaxrank import report

    pdf_path = tmp_path / "out.pdf"
    with caplog.at_level(logging.INFO, logger=report.logger.name):
        with pytest.raises(ValueError, match="Unknown PDF backend"):
            report.make_pdf_report({}, str(pdf_path), backend="bogus")
    assert not caplog.records

    monkeypatch.setattr(report, "render_report", lambda data, f, tpl: f.write("<html></html>"))
    monkeypatch.setattr(report, "_pdf_via_weasyprint", lambda html_path, out_path: None)
    with caplog.at_level(logging.INFO, logger=report.logger.name):
        report.make_pdf_report({}, str(pdf_path), backend="weasyprint")
    messages = [r.getMessage() for r in caplog.records]
    assert any("Rendering PDF report with weasyprint" in m for m in messages)
    assert any("Wrote PDF report with weasyprint" in m and m.endswith("s") for m in messages)
