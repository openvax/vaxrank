"""SHA-256 identity and validation helpers shared by every provenance record."""

import hashlib
import json
from pathlib import Path


_HEX64 = frozenset("0123456789abcdef")


def is_sha256_digest(value) -> bool:
    """Return True iff ``value`` is a 64-character lowercase hexadecimal string."""
    return isinstance(value, str) and len(value) == 64 and set(value) <= _HEX64


def text_sha256(text: str) -> str:
    """Return the hexadecimal SHA-256 digest of a UTF-8 encoded string."""
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def canonical_json_sha256(value) -> str:
    """Return the SHA-256 digest of compact, sorted, strict JSON."""
    return text_sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    )


def file_sha256(path, chunk_size=1024 * 1024):
    """Return the hexadecimal SHA-256 digest of a file's bytes.

    The file is read in chunks, so large BAM or reference files are never held
    in memory.

    Parameters
    ----------
    path : str or os.PathLike
        File to digest.
    chunk_size : int, optional
        Number of bytes read per iteration.

    Returns
    -------
    str
        64-character lowercase hexadecimal digest.

    Raises
    ------
    OSError
        If the file cannot be opened or read.
    ValueError
        If chunk_size is not a positive integer.
    """
    if type(chunk_size) is not int or chunk_size <= 0:
        raise ValueError('chunk_size must be a positive integer')
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        while chunk := stream.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()
