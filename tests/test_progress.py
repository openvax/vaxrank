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
