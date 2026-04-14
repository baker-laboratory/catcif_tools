"""
structure.py — reading and renaming individual CIF structures from .catcif files.
"""

import os
import re
import zlib

from .path import split_catcif_tag
from .cache import get_catcif_index as _get_catcif_index

_CHUNK = 65536
_GZ_MAGIC = b'\x1f\x8b'

# Matches the _entry.id line; value may be plain or single-quoted.
# Group 1 captures the key + whitespace so spacing is preserved on replacement.
_ENTRY_ID_RE = re.compile(r'^(_entry\.id\s+)\S+', re.MULTILINE)


def compress_structure(structure):
    """
    Return a gzip-compressed copy of structure as bytes.

    The result can be written directly into a .catcif file alongside
    plain-text members; _read_gz_structure will decompress it on demand.

    Parameters
    ----------
    structure : str
        Full CIF text for a single structure (already tagged/renamed).

    Returns
    -------
    bytes
        gzip-compressed CIF bytes.
    """
    c = zlib.compressobj(wbits=31)  # wbits=31 → gzip format
    data = structure.encode('ascii')
    return c.compress(data) + c.flush()


def rename_structure(structure, new_tag):
    """
    Return a copy of structure with every occurrence of the old tag replaced
    by new_tag.

    Two locations are updated:
      1. The ``data_<tag>`` block header line (always the first line).
      2. The ``_entry.id <tag>`` key-value item inside the block, if present.

    Parameters
    ----------
    structure : str
        Full CIF text for a single structure, starting with its data_ line.
    new_tag : str
        The replacement tag (the text that follows data_).

    Returns
    -------
    str
        CIF text with the tag updated in both locations.
    """
    nl = structure.index('\n')
    structure = 'data_' + new_tag + structure[nl:]
    structure = _ENTRY_ID_RE.sub(r'\g<1>' + new_tag, structure, count=1)
    return structure


def _read_plain_structure(f_open):
    """
    Read one plain-text CIF structure from the current file position.

    Reads until the next ``\\ndata_``, gzip magic (0x1f 0x8b), or EOF,
    whichever comes first.  The trailing newline before the next block is
    included; the delimiter itself is not consumed.
    """
    buf = bytearray()

    while True:
        # Keep 5 bytes of overlap so boundary patterns aren't missed.
        search_from = max(0, len(buf) - 5)
        chunk = f_open.read(_CHUNK)
        if not chunk:
            break
        buf.extend(chunk)

        nl_data  = buf.find(b'\ndata_', search_from)
        gz_magic = buf.find(_GZ_MAGIC,  search_from)

        ends = []
        if nl_data  != -1: ends.append(nl_data + 1)  # include the \n
        if gz_magic != -1: ends.append(gz_magic)      # exclude the magic

        if ends:
            end = min(ends)
            # Rewind so the next read starts at the delimiter.
            f_open.seek(end - len(buf), 1)
            return buf[:end].decode('ascii')

    return buf.decode('ascii')


def _read_gz_structure(f_open):
    """
    Decompress one gzip member from the current file position.

    Rewinds the file pointer to the first byte after the gzip member so
    subsequent reads are correctly positioned.
    """
    d = zlib.decompressobj(wbits=47)  # 47 = auto-detect gzip/zlib header
    out = bytearray()

    while True:
        chunk = f_open.read(_CHUNK)
        if not chunk:
            raise ValueError("Unexpected EOF inside gzip member")
        out.extend(d.decompress(chunk))
        if d.unused_data:
            f_open.seek(-len(d.unused_data), 1)
            break

    return out.decode('ascii')


def _get_structure(tag, catcif_file, catcif_index, f_open, preserve_tags=False):
    """
    Load and return the raw CIF text for tag from an open catcif file.

    Seeks to the correct offset, decompresses if needed, and renames the
    data_ line (and _entry.id) if deduplication changed the on-disk name,
    unless preserve_tags is True.
    """
    entry = catcif_index['index'][tag]

    if 'gz' in entry:
        f_open.seek(entry['gz'])
        structure = _read_gz_structure(f_open)
    else:
        f_open.seek(entry['o'])
        structure = _read_plain_structure(f_open)

    if not preserve_tags:
        orig_tag = catcif_index['orig_tags'][entry['idx']]
        if orig_tag != tag:
            structure = rename_structure(structure, tag)

    return structure


def get_all_structures(catcif_file, preserve_tags=False, dont_rename_structure=False):
    """
    Yield (structure, tag) for every structure in catcif_file, in file order.

    Designed for maximum throughput: structures are read sequentially with no
    backward seeks, and the rename step can be skipped entirely when the caller
    only needs to inspect structure content (e.g. score parsing).

    Parameters
    ----------
    catcif_file : str
        Path to the .catcif file.
    preserve_tags : bool
        If True, the returned tag is the original on-disk name and the
        structure text is not renamed.
    dont_rename_structure : bool
        If True, the structure string is returned exactly as read from disk
        (no rename_structure call) even when deduplication changed the name.
        The returned tag is still the canonical deduplicated name unless
        preserve_tags is also True.
    """
    index, f_open, caller_must_close = _get_catcif_index(catcif_file)
    try:
        orig_tags = index['orig_tags']
        for tag, entry in index['index'].items():
            if 'gz' in entry:
                f_open.seek(entry['gz'])
                structure = _read_gz_structure(f_open)
            else:
                f_open.seek(entry['o'])
                structure = _read_plain_structure(f_open)

            orig_tag = orig_tags[entry['idx']]
            returned_tag = orig_tag if preserve_tags else tag

            if (not preserve_tags
                    and not dont_rename_structure
                    and orig_tag != tag):
                structure = rename_structure(structure, tag)

            yield structure, returned_tag
    finally:
        if caller_must_close:
            f_open.close()


def get_structure(tag, catcif_file=None, no_cache=False, instant_cache=False,
                  preserve_tags=False):
    """
    Return the CIF text for a single structure.

    tag may be a bare tag name (requires catcif_file) or a
    ``path/to/file.catcif:tag`` string.  If both are supplied they must
    resolve to the same path.

    If preserve_tags is True, the data_ line and _entry.id are left as they
    appear on disk and no renaming is performed.
    """
    tag_catcif_file, tag = split_catcif_tag(tag)

    assert tag_catcif_file is not None or catcif_file is not None, \
        f"No catcif file specified for tag {tag!r}"

    if tag_catcif_file is not None and catcif_file is not None:
        assert os.path.realpath(tag_catcif_file) == os.path.realpath(catcif_file), \
            f"Conflicting catcif files: {tag_catcif_file!r} vs {catcif_file!r}"

    catcif_file = tag_catcif_file or catcif_file

    index, f_open, caller_must_close = _get_catcif_index(
        catcif_file, no_cache=no_cache, instant_cache=instant_cache)
    try:
        return _get_structure(tag, catcif_file, index, f_open,
                              preserve_tags=preserve_tags)
    finally:
        if caller_must_close:
            f_open.close()


def get_structures(tags, catcif_file=None, no_cache=False, instant_cache=False,
                   preserve_tags=False):
    """
    Yield CIF text for each tag in the iterable tags.

    All tags must reference the same catcif file (either via the
    ``path.catcif:tag`` syntax or via the catcif_file argument).
    The assertion is checked up-front before any structures are yielded.
    f_open is closed at the end if the caller is responsible for it.

    If preserve_tags is True, no renaming is performed on any structure.
    """
    # Parse all tags first so we can assert before yielding anything.
    parsed = []       # list of (file_path, clean_tag)
    real_paths = {}   # realpath -> original path (for the assertion message)

    if catcif_file is not None:
        real_paths[os.path.realpath(catcif_file)] = catcif_file

    for t in tags:
        tag_catcif_file, clean_tag = split_catcif_tag(t)
        file_for_tag = tag_catcif_file if tag_catcif_file is not None else catcif_file

        assert file_for_tag is not None, \
            f"No catcif file for tag {t!r}"

        rp = os.path.realpath(file_for_tag)
        real_paths[rp] = file_for_tag
        parsed.append((file_for_tag, clean_tag))

    assert len(real_paths) <= 1, \
        ("get_structures: all tags must come from the same catcif file, "
         f"got: {list(real_paths.values())}")

    if not parsed:
        return

    actual_catcif = parsed[0][0]
    index, f_open, caller_must_close = _get_catcif_index(
        actual_catcif, no_cache=no_cache, instant_cache=instant_cache)
    try:
        for _, clean_tag in parsed:
            yield _get_structure(clean_tag, actual_catcif, index, f_open,
                                 preserve_tags=preserve_tags)
    finally:
        if caller_must_close:
            f_open.close()
