"""Regression tests for scoped CLI logging lifecycle (openvax/vaxrank#882)."""

import io
import logging
from contextlib import redirect_stderr, redirect_stdout

import pytest

from vaxrank.cli import entry_point


def test_consecutive_captured_cli_runs_do_not_leak_closed_streams(tmp_path):
    """Two captured CLI invocations whose stdout streams are closed between
    calls must not leave stale StreamHandlers on pyensembl/datacache/vaxrank
    loggers, while still writing real diagnostics to each run's log file.
    """
    for idx in (1, 2):
        log_path = tmp_path / f"run_{idx}.log"
        missing_lens = tmp_path / f"missing_{idx}.tsv"
        out_dir = tmp_path / f"out_{idx}"
        captured_out = io.StringIO()
        captured_err = io.StringIO()
        with redirect_stdout(captured_out), redirect_stderr(captured_err):
            with pytest.raises(SystemExit) as exc_info:
                entry_point.main([
                    "--input-lens", str(missing_lens),
                    "--output-dir", str(out_dir),
                    "--log-path", str(log_path),
                ])
        assert exc_info.value.code == 1
        assert "vaxrank: error:" in captured_err.getvalue()
        assert "Exiting on user error" in log_path.read_text()
        # Simulate pytest capsys tearing down and closing the captured stream.
        captured_out.close()
        captured_err.close()

        # Subsequent library logging (e.g., pyensembl/datacache copying fixtures)
        # must not write into the closed stream or emit logging tracebacks.
        post_err = io.StringIO()
        with redirect_stderr(post_err):
            logging.getLogger("pyensembl").info("copying fixture %d", idx)
            logging.getLogger("datacache").info("caching fixture %d", idx)
            logging.getLogger("vaxrank").info("post-run diagnostic %d", idx)
        assert "I/O operation on closed file" not in post_err.getvalue()
        assert "Logging error" not in post_err.getvalue()


def test_direct_run_cli_restores_logging_and_closes_file_handlers(tmp_path):
    """Two consecutive direct run_cli() invocations restore loggers and close log files."""
    for idx in (1, 2):
        log_path = tmp_path / f"direct_{idx}.log"
        missing_lens = tmp_path / f"missing_{idx}.tsv"
        out_dir = tmp_path / f"out_{idx}"
        captured_out = io.StringIO()
        with redirect_stdout(captured_out):
            with pytest.raises(ValueError, match="--input-lens"):
                entry_point.run_cli([
                    "--input-lens", str(missing_lens),
                    "--output-dir", str(out_dir),
                    "--log-path", str(log_path),
                ])
        captured_out.close()

        for name in ("vaxrank", "pyensembl", "datacache", "isovar", "varcode", "mhctools"):
            for handler in logging.getLogger(name).handlers:
                stream = getattr(handler, "stream", None)
                if stream is not None:
                    assert not getattr(stream, "closed", False), (
                        f"Logger {name!r} retained a closed handler stream: {handler!r}"
                    )

        post_err = io.StringIO()
        with redirect_stderr(post_err):
            logging.getLogger("datacache").info("after direct run_cli %d", idx)
        assert "I/O operation on closed file" not in post_err.getvalue()


def test_existing_logger_handlers_preserved_across_cli_invocation(tmp_path):
    """Caller-owned handlers (such as pytest caplog) remain attached and enabled."""
    custom_logger = logging.getLogger("vaxrank.test_caller_owned")
    records = []

    class _Recorder(logging.Handler):
        def emit(self, record):
            records.append(record.getMessage())

    handler = _Recorder()
    custom_logger.addHandler(handler)
    custom_logger.setLevel(logging.INFO)
    custom_logger.disabled = False
    try:
        log_path = tmp_path / "caller.log"
        with pytest.raises(SystemExit):
            entry_point.main([
                "--input-lens", str(tmp_path / "missing.tsv"),
                "--output-dir", str(tmp_path / "out"),
                "--log-path", str(log_path),
            ])
        assert handler in custom_logger.handlers
        assert custom_logger.disabled is False
        custom_logger.info("still-alive")
        assert records == ["still-alive"]
    finally:
        custom_logger.removeHandler(handler)
