"""catcifsequence — extract per-chain sequences from a .catcif file."""

import argparse
import sys
import textwrap
from signal import signal, SIGPIPE, SIG_DFL

signal(SIGPIPE, SIG_DFL)

from .biology import get_sequence
from .structure import get_all_structures


def main():
    parser = argparse.ArgumentParser(
        prog='catcifsequence',
        description='Extract the amino-acid sequence of each structure in a .catcif file.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            examples:
              catcifsequence my.catcif
              catcifsequence -c my.catcif
        """),
    )
    parser.add_argument('catcif_file', metavar='FILE', help='.catcif file to read')
    parser.add_argument(
        '-c', action='store_true',
        help='prefix each chain sequence with its chain ID and a colon (e.g. A:ACDE)',
    )
    opts = parser.parse_args()

    for structure, tag in get_all_structures(opts.catcif_file):
        chains = get_sequence(structure)

        if not chains:
            continue

        if opts.c:
            chunks = [f"{ch['chain_id']}:{ch['sequence']}" for ch in chains]
        else:
            chunks = [ch['sequence'] for ch in chains]

        sys.stdout.write(' '.join(chunks) + ' ' + tag + '\n')
        sys.stdout.flush()
