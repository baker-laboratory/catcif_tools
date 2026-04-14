"""
build_index.py — build and load .idx files for .catcif files.

A .catcif file may contain plain CIF text and individually-gzipped CIF structures
interleaved arbitrarily. Each gzip member encodes exactly one CIF structure.

Index file format (JSON, compact):
{
  "version": 2,
  "orig_tags": ["orig1", "orig2", ...],   // original names as found in file, in order
  "index": {
    "tag1": {"o": <int>},          // plain: absolute byte offset of data_ line
    "tag2": {"gz": <int>},         // gzip: byte offset of gzip member in file
    ...
  }
}

The ordered tag list is derivable from index.keys() (Python 3.7+ dict insertion order),
so it is not stored separately.

At runtime (not written to disk), each entry also carries:
  "idx": <int>    // 0-based position; used to look up orig_tags[entry["idx"]]

Duplicate tag handling: if a tag already exists in the index the incoming
structure is renamed tag_1, then tag_2, etc. (counter starts at 1 and
increments per collision). The index is always keyed by the deduplicated name;
orig_tags[entry["idx"]] gives the original name for any entry.

Keys are kept short because the index can contain hundreds of thousands of entries.
"""

import mmap
import os
import json
import zlib

CATCIF_INDEX_VERSION = 2
_GZ_MAGIC = b'\x1f\x8b'
_DATA_PREFIX = b'data_'
_CHUNK = 65536  # 64KB decompression chunks


def get_index_path(catcif_path):
    """Return the .idx path for a given .catcif path."""
    return catcif_path + ".idx"


def get_uncached_catcif_index(catcif_path):
    """
    Return the index for a catcif file, building it if necessary.

    Rebuilds if the data file is newer than the index or the index is corrupt/stale.
    Runtime-only fields (idx) are added after loading.
    """
    index_path = get_index_path(catcif_path)

    if os.path.exists(index_path):
        if os.path.getmtime(catcif_path) <= os.path.getmtime(index_path):
            try:
                with open(index_path) as f:
                    index = json.load(f)
                if index.get("version") == CATCIF_INDEX_VERSION:
                    _add_runtime_fields(index)
                    return index
            except Exception:
                pass

    return build_catcif_index(catcif_path)


def _add_runtime_fields(index):
    """Add runtime-only fields to each index entry (not written to disk)."""
    for i, entry in enumerate(index["index"].values()):
        entry["idx"] = i


def _dedup_tag(orig_tag, seen):
    """
    Return a collision-free tag name and update the seen counter.

    First occurrence of a name is returned unchanged.
    Each subsequent occurrence becomes name_1, name_2, etc.
    """
    if orig_tag not in seen:
        seen[orig_tag] = 0
        return orig_tag
    seen[orig_tag] += 1
    return f"{orig_tag}_{seen[orig_tag]}"


def build_catcif_index(catcif_path):
    """
    Scan a .catcif file and write a .idx file next to it. Returns the index dict.

    Handles plain CIF text and individually-gzipped CIF structures in any order.
    Plain sections are scanned with a fast C-level mmap byte search.
    Gzip members are decompressed with zlib to extract the tag and find the
    member boundary; output is buffered only until the tag (first line) is found.

    Duplicate tag names are resolved: the second occurrence of 'foo' becomes
    'foo_1', the third 'foo_2', etc. orig_tags preserves the original names.
    """
    orig_tags = []
    entries = {}
    seen = {}

    with open(catcif_path, 'rb') as f:
        mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)

    try:
        pos = 0
        size = len(mm)

        while pos < size:
            gz_pos = mm.find(_GZ_MAGIC, pos)

            # Process the plain-text segment before the next gzip member (or EOF)
            plain_end = gz_pos if gz_pos != -1 else size
            _scan_plain_segment(mm, pos, plain_end, orig_tags, entries, seen)

            if gz_pos == -1:
                break

            # Process the gzip member
            orig_tag, gz_end = _consume_gz_member(mm, gz_pos)
            tag = _dedup_tag(orig_tag, seen)
            orig_tags.append(orig_tag)
            entries[tag] = {"gz": gz_pos}
            pos = gz_end
    finally:
        mm.close()

    index = {"version": CATCIF_INDEX_VERSION, "orig_tags": orig_tags, "index": entries}

    index_path = get_index_path(catcif_path)
    with open(index_path, 'w') as f:
        json.dump(index, f, separators=(',', ':'))

    _add_runtime_fields(index)
    return index


def _scan_plain_segment(mm, start, end, orig_tags, entries, seen):
    """
    Scan a plain-text byte range of the mmap for data_ lines.

    Uses mm.find() for a fast C-level search rather than iterating lines.
    Duplicate tag names are resolved via _dedup_tag before being stored.
    """
    if start >= end:
        return

    # Check if the segment begins with data_ (start of file or after a gz member)
    if mm[start:start + 5] == _DATA_PREFIX:
        nl = mm.find(b'\n', start, end)
        line_end = nl if nl != -1 else end
        orig_tag = mm[start + 5:line_end].rstrip().decode('ascii')
        tag = _dedup_tag(orig_tag, seen)
        orig_tags.append(orig_tag)
        entries[tag] = {"o": start}

    # Find subsequent data_ lines (must be preceded by a newline)
    pos = start
    while True:
        pos = mm.find(b'\ndata_', pos, end)
        if pos == -1:
            break
        tag_start = pos + 1  # skip the leading \n
        nl = mm.find(b'\n', tag_start, end)
        line_end = nl if nl != -1 else end
        orig_tag = mm[tag_start + 5:line_end].rstrip().decode('ascii')
        tag = _dedup_tag(orig_tag, seen)
        orig_tags.append(orig_tag)
        entries[tag] = {"o": tag_start}
        pos = line_end


def _consume_gz_member(mm, gz_start):
    """
    Decompress a single gzip member to extract its tag and locate the member end.

    Returns (tag, gz_end) where gz_end is the byte offset of the first byte
    AFTER the gzip member in the file.

    Decompressed output is buffered only until the first newline (sufficient to
    extract the data_ tag). After that, chunks are fed to the decompressor and
    discarded — we must consume all compressed bytes to find the member boundary.
    """
    d = zlib.decompressobj(wbits=47)  # 47 = auto-detect gzip or zlib header

    buf = bytearray()
    tag = None
    pos = gz_start

    while True:
        chunk = bytes(mm[pos:pos + _CHUNK])
        if not chunk:
            # EOF — this is fine if the gzip member ended exactly at the file
            # boundary (d.unused_data stays b"" because there are no trailing
            # bytes). Verify the stream was complete by flushing; zlib raises
            # zlib.error if the stream was truncated.
            try:
                d.flush()
            except zlib.error:
                raise ValueError(
                    f"Unexpected end of file inside gzip member at offset {gz_start}"
                )
            if tag is None:
                raise ValueError(
                    f"Gzip member at offset {gz_start} does not start with data_"
                )
            return tag, pos

        out = d.decompress(chunk)

        if tag is None:
            buf.extend(out)
            nl = buf.find(b'\n')
            if nl != -1:
                first_line = buf[:nl]
                if not first_line.startswith(_DATA_PREFIX):
                    raise ValueError(
                        f"Gzip member at offset {gz_start} does not start with data_"
                    )
                tag = first_line[5:].rstrip().decode('ascii')
                buf = None  # release buffer; discard all further output

        if d.unused_data:
            # The gzip member ended somewhere within `chunk`.
            # d.unused_data holds bytes from `chunk` that belong to the next segment.
            gz_end = pos + len(chunk) - len(d.unused_data)
            return tag, gz_end

        pos += len(chunk)
