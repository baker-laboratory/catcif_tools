"""
scores.py — reading and writing per-structure scores in .catcif files.

Scores are stored as a custom CIF key-value block within the structure's
data block, using the _catcif_scores category:

    #
    _catcif_scores.total_score   -123.45
    _catcif_scores.rmsd          1.234
    #

write_scores() places the block immediately after the data_ line (as high
as possible) so that score-only scans can stop reading early.
"""

import re

from .structure import get_all_structures

# Removes a catcif_scores block, including the bracketing # comment lines
# when they are immediately adjacent to the score lines.
# Pattern: optional leading #, one or more score lines, optional trailing #.
# re.MULTILINE makes ^ anchor to the start of each line.
_SCORES_BLOCK_RE = re.compile(
    r'(?:^#\n)?(?:^_catcif_scores\.[^\n]+\n?)+(?:^#\n)?',
    re.MULTILINE,
)

# Fast key-value scanner: captures (name, value) from _catcif_scores lines.
_SCORES_SCAN_RE = re.compile(
    r'^_catcif_scores\.(\S+)\s+(\S+)',
    re.MULTILINE,
)


def write_scores(structure, scores={}):
    """
    Return structure with the _catcif_scores block replaced by scores.

    The new block is inserted immediately after the data_ line.
    Any existing _catcif_scores block (including its bracketing # lines) is
    removed first.  If scores is empty the block is removed with nothing
    inserted.

    str() is called on every value.

    Parameters
    ----------
    structure : str
        Full CIF text for a single structure.
    scores : dict
        Mapping of score name to value.

    Returns
    -------
    str
    """
    # Only pay for the regex if scores already exist.
    if '_catcif_scores.' in structure:
        structure = _SCORES_BLOCK_RE.sub('', structure)

    if not scores:
        return structure

    # Build the new block.
    lines = ['#']
    for name, value in scores.items():
        lines.append(f'_catcif_scores.{name}   {value}')
    lines.append('#')
    block = '\n'.join(lines) + '\n'

    # Insert immediately after the data_ line.
    nl = structure.index('\n')
    return structure[:nl + 1] + block + structure[nl + 1:]


def get_scores(structure):
    """
    Parse and return the _catcif_scores block from structure.

    Values are returned as strings.  Designed to be as fast as possible
    since this function is called on every structure in score-heavy workloads.

    Parameters
    ----------
    structure : str
        Full CIF text for a single structure.

    Returns
    -------
    dict[str, str]
        Mapping of score name to value string.  Empty dict if no scores found.
    """
    return {m.group(1): m.group(2) for m in _SCORES_SCAN_RE.finditer(structure)}


def parse_score_file(catcif_file, preserve_tags=False):
    """
    Yield one scores dict per structure in catcif_file.

    Each dict contains the _catcif_scores values as strings, with the
    structure's tag stored under the key 'tag' as the last entry.

    The structure text is never renamed (dont_rename_structure=True) since
    only the score values and tag are needed.

    Parameters
    ----------
    catcif_file : str
        Path to the .catcif file.
    preserve_tags : bool
        If True, 'tag' reflects the original on-disk name rather than the
        canonical deduplicated name.

    Yields
    ------
    dict[str, str]
    """
    for structure, tag in get_all_structures(
            catcif_file, preserve_tags=preserve_tags, dont_rename_structure=True):
        scores = get_scores(structure)
        scores['tag'] = tag
        yield scores
