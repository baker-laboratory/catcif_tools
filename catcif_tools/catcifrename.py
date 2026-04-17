"""catcifrename — rename the tags of structures in a .catcif file."""

import os
import argparse
import stat
import sys
import textwrap
from signal import signal, SIGPIPE, SIG_DFL

signal(SIGPIPE, SIG_DFL)

from .cache import get_catcif_index
from .structure import compress_structure, get_all_structures, rename_structure


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def main():
    parser = argparse.ArgumentParser(
        prog='catcifrename',
        description='Rename the tags of structures in a .catcif file.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            examples:
              catcifls my.catcif > tags.list
              vi tags.list
              cat tags.list | catcifrename my.catcif > new.catcif

              catcifls my.catcif | sed 's/_0001$//g' | catcifrename my.catcif > new.catcif
        """),
    )
    parser.add_argument('catcif_file', metavar='FILE', help='.catcif file to rename')
    parser.add_argument(
        'new_tags', nargs='*', metavar='TAG',
        help='new tag names (may also be supplied via stdin)',
    )
    parser.add_argument(
        '-b', action='store_true',
        help='broadcast mode: duplicate a single-structure file with arbitrarily many names',
    )
    parser.add_argument(
        '-z', action='store_true',
        help='write gzip-compressed output',
    )
    opts = parser.parse_args()

    # Collect new tags from stdin and positional args
    tag_lines = []
    if stat.S_ISFIFO(os.stat('/dev/stdin').st_mode):
        tag_lines += sys.stdin.readlines()
    tag_lines += opts.new_tags

    new_tags = []
    any_spaces = False
    for line in tag_lines:
        line = line.strip()
        if not line:
            continue
        new_tags.append(line)
        if ' ' in line:
            any_spaces = True

    if any_spaces:
        eprint('catcifrename: Warning: spaces in tag names may break downstream tools')

    if not new_tags:
        parser.print_help(sys.stderr)
        sys.exit(1)

    # Check structure count up front
    index = get_catcif_index(opts.catcif_file)
    n_structures = len(index['index'])

    if opts.b:
        if n_structures != 1:
            eprint(f'catcifrename: -b requires exactly 1 structure in the file, '
                   f'found {n_structures}')
            sys.exit(1)
    else:
        if len(new_tags) != n_structures:
            eprint(f'catcifrename: number of new tags ({len(new_tags)}) does not match '
                   f'number of structures in file ({n_structures})')
            sys.exit(1)

    def _write(structure):
        if opts.z:
            sys.stdout.buffer.write(compress_structure(structure))
            sys.stdout.buffer.flush()
        else:
            sys.stdout.write(structure)
            sys.stdout.flush()

    if opts.b:
        # Read the single structure once, output it N times with different names
        structure, _ = next(iter(get_all_structures(
            opts.catcif_file, dont_rename_structure=True)))
        for new_tag in new_tags:
            _write(rename_structure(structure, new_tag))
    else:
        for (structure, _), new_tag in zip(
                get_all_structures(opts.catcif_file, dont_rename_structure=True),
                new_tags):
            _write(rename_structure(structure, new_tag))
