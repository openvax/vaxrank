"""SHA-256 identity of files on disk, shared by every provenance record."""

import hashlib
from pathlib import Path


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
    """
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        while chunk := stream.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()
