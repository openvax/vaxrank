"""Progress reporting for the slow steps of a run.

Bars go to stderr and switch themselves off when stderr is not a terminal, so a
redirected or CI log gets the INFO status lines but no carriage-return spam
(openvax/vaxrank#460). In Jupyter notebooks, ``tqdm.auto`` renders notebook
progress bars even though ``sys.stderr.isatty()`` is ``False`` (#783).
"""

import logging
import sys
import time
from contextlib import contextmanager

def _in_jupyter_notebook() -> bool:
    """Return True when running inside an active Jupyter notebook kernel."""
    ipython = sys.modules.get("IPython")
    if ipython is None:
        return False
    get_ipython = getattr(ipython, "get_ipython", None)
    if not callable(get_ipython):
        return False
    shell = get_ipython()
    return shell is not None and type(shell).__name__ == "ZMQInteractiveShell"


def progress(iterable, desc, unit, **kwargs):
    """``tqdm`` with vaxrank's defaults: stderr only, auto-disabled off a terminal,
    and enabled in Jupyter notebooks."""
    if _in_jupyter_notebook():
        from tqdm.auto import tqdm as auto_tqdm
        kwargs.setdefault("disable", False)
        return auto_tqdm(iterable, desc=desc, unit=unit, **kwargs)
    from tqdm import tqdm
    kwargs.setdefault("file", sys.stderr)
    kwargs.setdefault("disable", None)
    return tqdm(iterable, desc=desc, unit=unit, **kwargs)


@contextmanager
def redirect_logging_for_progress(loggers=None):
    """Route console log lines through ``tqdm.write`` so progress bars stay intact."""
    from tqdm.contrib.logging import logging_redirect_tqdm
    with logging_redirect_tqdm(loggers=loggers):
        yield


@contextmanager
def timed_step(logger, start_message, done_message, *args):
    """Log ``start_message`` now and ``done_message`` with the elapsed time after.

    ``done_message`` is a %-format string that receives ``args`` followed by the
    elapsed seconds, so it should end with ``%.0fs``.
    """
    logger.info(start_message, *args, stacklevel=3)
    started = time.monotonic()
    yield
    logger.log(logging.INFO, done_message, *args, time.monotonic() - started, stacklevel=3)
