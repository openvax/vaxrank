"""Progress reporting for the slow steps of a run.

Bars go to stderr and switch themselves off when stderr is not a terminal, so a
redirected or CI log gets the INFO status lines but no carriage-return spam
(openvax/vaxrank#460).
"""

import logging
import sys
import time
from contextlib import contextmanager

from tqdm import tqdm


def progress(iterable, desc, unit, **kwargs):
    """``tqdm`` with vaxrank's defaults: stderr only, auto-disabled off a terminal."""
    return tqdm(iterable, desc=desc, unit=unit, file=sys.stderr, disable=None, **kwargs)


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
