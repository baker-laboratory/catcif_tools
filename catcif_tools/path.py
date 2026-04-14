"""
path.py — helpers for parsing catcif path:tag strings.
"""


def is_catcif_path_tag(tag):
    """Return True if tag contains '.catcif:'. or '.cif:'"""
    return '.catcif:' in tag or '.cif:' in tag


def split_catcif_tag(tag):
    """
    Split a path:tag string into (path, tag).

    If tag contains '.catcif:', splits on the first ':' and returns
    (path, ':'.join(remaining_parts)).

    Otherwise returns (None, tag).
    """
    if not is_catcif_path_tag(tag):
        return None, tag
    parts = tag.split(':')
    return parts[0], ':'.join(parts[1:])
