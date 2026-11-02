"""Behavior every ``vaxrank-*`` companion command shares with ``vaxrank`` itself."""

from ..version import __version__


def add_version_argument(parser):
    """Give a companion parser the ``--version`` that ``vaxrank`` already has.

    Only ``vaxrank`` answered ``--version``; the other commands treated it as a
    usage error (openvax/vaxrank#768).
    """
    parser.add_argument("--version", action="version",
                        version="%(prog)s " + __version__)
    return parser
