"""
catcif_tools.py — high-level API for building and working with .catcif files.
"""

import re

from .structure import rename_structure, compress_structure
from .scores import write_scores

# Matches zero or more leading blank lines and # comment lines.
_LEADING_COMMENTS_RE = re.compile(r'\A(?:#[^\n]*\n|\n)*')


def get_tags_from_index(index):
    """Return the ordered list of canonical tags from a catcif index."""
    return list(index["index"].keys())


def get_tags(path):
    """Return the ordered list of canonical tags in a .catcif file."""
    from .cache import get_catcif_index
    return get_tags_from_index(get_catcif_index(path))


def to_catcif_string(cif_str, tag, scores=None, add_header=False, compress=False):
    """
    Return cif_str prepared for inclusion in a .catcif file.

    Renames the data_ block header and _entry.id to tag.  If scores is not
    None, replaces any existing _catcif_scores block with the new scores.

    Parameters
    ----------
    cif_str : str
        Raw CIF text.  May or may not already have a data_ header.
    tag : str
        The canonical tag to assign to this structure.
    scores : dict or None
        If not None, scores are written into the block via write_scores().
    add_header : bool
        If the CIF text has no data_ header and add_header is True, one is
        prepended automatically.  If False (default), a ValueError is raised
        instead — this prevents accidentally passing non-CIF content such as
        PDB text.
    compress : bool
        If True, return the structure as gzip-compressed bytes instead of str.

    Returns
    -------
    str or bytes
        Plain CIF text when compress is False; gzip bytes when compress is True.
    """
    cif_str = cif_str[_LEADING_COMMENTS_RE.match(cif_str).end():]
    if not cif_str.startswith('data_'):
        if not add_header:
            raise ValueError(
                "CIF string has no data_ block header. "
                "Pass add_header=True to insert one automatically. "
                "Note: PDB format is not supported by this function."
            )
        cif_str = f'data_{tag}\n' + cif_str

    cif_str = rename_structure(cif_str, tag)

    if scores is not None:
        cif_str = write_scores(cif_str, scores)

    if compress:
        return compress_structure(cif_str)
    return cif_str


def append_to_catcif_file(catcif_file, cif_str, tag, scores=None, add_header=False,
                           compress=False):
    """
    Append a single CIF structure to a .catcif file on disk.

    Opens catcif_file in append mode and delegates to
    append_to_catcif_file_open().

    Parameters
    ----------
    catcif_file : str
        Path to the .catcif file (created if it does not exist).
    cif_str : str
        Raw CIF text for the structure to append.
    tag : str
        Canonical tag to assign to the structure.
    scores : dict or None
        Passed through to to_catcif_string().
    add_header : bool
        Passed through to to_catcif_string().
    compress : bool
        If True, the structure is gzip-compressed before being written.
    """
    mode = 'ab' if compress else 'a'
    with open(catcif_file, mode) as f:
        append_to_catcif_file_open(f, cif_str, tag, scores=scores,
                                   add_header=add_header, compress=compress)


def append_to_catcif_file_open(catcif_file_pointer, cif_str, tag,
                                scores=None, add_header=False, compress=False):
    """
    Append a single CIF structure to an already-open .catcif file.

    Parameters
    ----------
    catcif_file_pointer : file-like
        Open file object in text append (or write) mode when compress is False;
        binary append (or write) mode when compress is True.
    cif_str : str
        Raw CIF text for the structure to append.
    tag : str
        Canonical tag to assign to the structure.
    scores : dict or None
        Passed through to to_catcif_string().
    add_header : bool
        Passed through to to_catcif_string().
    compress : bool
        If True, the structure is gzip-compressed before being written.
        The file pointer must have been opened in binary mode.
    """
    catcif_file_pointer.write(
        to_catcif_string(cif_str, tag, scores=scores, add_header=add_header,
                         compress=compress)
    )
