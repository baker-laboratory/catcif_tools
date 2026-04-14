"""catcifslice — extract a subset of structures from a .catcif file."""

import argparse
import os
import stat
import sys
import textwrap
from collections import defaultdict
from signal import signal, SIGPIPE, SIG_DFL

signal(SIGPIPE, SIG_DFL)

from .cache import get_catcif_index
from .structure import compress_structure, get_structure, get_structures
from .path import split_catcif_tag


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def _get_offset(entry):
    return entry.get('o', entry.get('gz', 0))


def _load_index(path, index_cache):
    """Return (original_path, index) for path, loading and caching it if needed."""
    rp = os.path.realpath(path)
    if rp not in index_cache:
        index, f_open, caller_must_close = get_catcif_index(path, instant_cache=True)
        if caller_must_close:
            f_open.close()
        index_cache[rp] = (path, index)
    return index_cache[rp]


def _find_in_index(clean_tag, index, g_flag):
    """
    Return list of (canonical_tag, offset) matches for clean_tag in index.
    With g_flag, clean_tag is compared against original names.
    """
    results = []
    if g_flag:
        for canonical_tag, entry in index['index'].items():
            if index['orig_tags'][entry['idx']] == clean_tag:
                results.append((canonical_tag, _get_offset(entry)))
    else:
        if clean_tag in index['index']:
            entry = index['index'][clean_tag]
            results.append((clean_tag, _get_offset(entry)))
    return results


def _build_parser(prog, description, epilog):
    """
    Return the shared argument parser with flags common to all catcif commands.

    Callers add their own command-specific flags to the returned parser before
    calling parse_args().
    """
    parser = argparse.ArgumentParser(
        prog=prog,
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(epilog),
    )
    parser.add_argument(
        '-g', action='store_true',
        help='slice by original name instead of deduplicated canonical names',
    )
    parser.add_argument(
        '-m', action='store_true',
        help='multi-catcif: all positional args are catcif files; tags must come via stdin',
    )
    parser.add_argument(
        '-z', action='store_true',
        help='write gzip-compressed output',
    )
    parser.add_argument(
        'positional', nargs='*', metavar='FILE_OR_TAG',
        help='catcif file followed by tags, or just tags (with embedded paths)',
    )
    return parser


def _collect_and_resolve(opts, prog):
    """
    Collect tags from stdin and CLI, resolve each to (orig_path, canonical_tag, offset).

    Returns the resolved list, or calls sys.exit(1) with usage if no tags were given.
    """
    # Collect tags from stdin
    stdin_tags = []
    if stat.S_ISFIFO(os.stat('/dev/stdin').st_mode):
        for line in sys.stdin:
            line = line.strip()
            if line:
                stdin_tags.extend(line.split())

    # Split positional args into catcif files and tag args
    if opts.m:
        catcif_files = opts.positional
        cli_tags = []
    elif opts.positional and (opts.positional[0].endswith('.catcif') or
                              opts.positional[0].endswith('.cif')):
        catcif_files = [opts.positional[0]]
        cli_tags = opts.positional[1:]
    else:
        catcif_files = []
        cli_tags = opts.positional

    all_tags = stdin_tags + cli_tags

    if not all_tags:
        # Re-build parser just to print help cleanly
        return None

    if len(all_tags) != len(set(all_tags)):
        eprint(f'{prog}: Warning: duplicate tags specified')

    # Load indexes for explicitly-passed catcif files
    index_cache = {}
    for cf in catcif_files:
        _load_index(cf, index_cache)

    # Resolve each tag to (orig_path, canonical_tag, offset)
    resolved = []

    for raw_tag in all_tags:
        tag_catcif_file, clean_tag = split_catcif_tag(raw_tag)

        if tag_catcif_file is not None:
            # Path tag — use the embedded file, ignore passed catcif files.
            orig_path, index = _load_index(tag_catcif_file, index_cache)
            matches = _find_in_index(clean_tag, index, opts.g)
            if not matches:
                eprint(f'{prog}: Unable to find tag: {raw_tag}')
                continue
            for canonical_tag, offset in matches:
                resolved.append((orig_path, canonical_tag, offset))

        else:
            # Bare tag — search the explicitly-passed catcif files in order.
            found = False
            for cf in catcif_files:
                orig_path, index = index_cache[os.path.realpath(cf)]
                matches = _find_in_index(clean_tag, index, opts.g)
                if matches:
                    if found and not opts.g:
                        eprint(f'{prog}: Tag found multiple times: {clean_tag}')
                    for canonical_tag, offset in matches:
                        resolved.append((orig_path, canonical_tag, offset))
                    found = True
            if not found:
                eprint(f'{prog}: Unable to find tag: {raw_tag}')

    return resolved


def _iter_structures(resolved, opts):
    """
    Yield (canonical_tag, structure) in the order determined by opts.e.

    With -e: input order, one get_structure() call per tag.
    Without -e: file order, get_structures() called once per catcif file.
    """
    if getattr(opts, 'e', False):
        for orig_path, canonical_tag, _ in resolved:
            yield canonical_tag, get_structure(canonical_tag, catcif_file=orig_path)
    else:
        file_order = []
        by_file = defaultdict(list)
        for orig_path, canonical_tag, offset in resolved:
            rp = os.path.realpath(orig_path)
            if rp not in by_file:
                file_order.append(orig_path)
            by_file[rp].append((offset, canonical_tag))

        for orig_path in file_order:
            rp = os.path.realpath(orig_path)
            tags_for_file = sorted(by_file[rp], key=lambda x: x[0])
            canonical_tags = [t for _, t in tags_for_file]
            for canonical_tag, structure in zip(
                    canonical_tags, get_structures(canonical_tags, catcif_file=orig_path)):
                yield canonical_tag, structure


def main():
    parser = _build_parser(
        prog='catcifslice',
        description='Extract a subset of structures from a .catcif file.',
        epilog="""\
            examples:
              cat tags.list | catcifslice my.catcif > out.catcif
              catcifslice my.catcif tag1 tag2 tag3 > out.catcif
              cat tags.list | catcifslice > out.catcif    # tags contain embedded paths
        """,
    )
    parser.add_argument(
        '-e', action='store_true',
        help='preserve input order of tags (default: file order, which is faster)',
    )
    opts = parser.parse_args()
    resolved = _collect_and_resolve(opts, prog='catcifslice')
    if not resolved:
        parser.print_help(sys.stderr)
        sys.exit(1)
    if opts.z:
        for _, structure in _iter_structures(resolved, opts):
            sys.stdout.buffer.write(compress_structure(structure))
            sys.stdout.buffer.flush()
    else:
        for _, structure in _iter_structures(resolved, opts):
            sys.stdout.write(structure)
            sys.stdout.flush()
