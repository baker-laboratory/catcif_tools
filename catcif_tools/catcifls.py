"""catcifls — list tags in a .catcif file."""

import argparse
import os
import sys
from signal import signal, SIGPIPE, SIG_DFL

signal(SIGPIPE, SIG_DFL)

from .cache import get_catcif_index


def main():
    parser = argparse.ArgumentParser(
        prog='catcifls',
        description='List structure tags in one or more .catcif files.',
    )
    parser.add_argument(
        'catcif_files', nargs='+', metavar='FILE',
        help='.catcif file(s) to list',
    )
    parser.add_argument(
        '-l', action='store_true',
        help='long format: prepend each tag with the real path of its file '
             '(output is suitable for use as catcif path:tag strings)',
    )
    parser.add_argument(
        '-g', action='store_true',
        help='show original names instead of deduplicated canonical names',
    )

    opts = parser.parse_args()

    for path in opts.catcif_files:
        index, f_open, caller_must_close = get_catcif_index(path)
        if caller_must_close:
            f_open.close()

        prefix = os.path.realpath(path) + ':' if opts.l else ''
        tags = index['orig_tags'] if opts.g else index['index']

        for tag in tags:
            print(prefix + tag)
