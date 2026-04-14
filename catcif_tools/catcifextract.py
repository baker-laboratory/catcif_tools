"""catcifextract — write each structure from a .catcif file to its own .cif file."""

import os
import sys
from signal import signal, SIGPIPE, SIG_DFL

signal(SIGPIPE, SIG_DFL)

from .catcifslice import _build_parser, _collect_and_resolve, _iter_structures
from .structure import compress_structure, get_all_structures


def _write_structure(canonical_tag, structure, out_dir, compress):
    if compress:
        out_path = os.path.join(out_dir, f'{canonical_tag}.cif.gz')
        with open(out_path, 'wb') as f:
            f.write(compress_structure(structure))
    else:
        out_path = os.path.join(out_dir, f'{canonical_tag}.cif')
        with open(out_path, 'w') as f:
            f.write(structure)


def main():
    parser = _build_parser(
        prog='catcifextract',
        description='Write each structure from a .catcif file to its own .cif file.',
        epilog="""\
            examples:
              catcifextract my.catcif
              cat tags.list | catcifextract my.catcif
              catcifextract my.catcif tag1 tag2 tag3
              catcifextract -o out_dir my.catcif tag1 tag2 tag3
              cat tags.list | catcifextract    # tags contain embedded paths
        """,
    )
    parser.add_argument(
        '-o', metavar='DIR', default='.',
        help='output directory for .cif files (default: current directory)',
    )
    opts = parser.parse_args()
    os.makedirs(opts.o, exist_ok=True)

    resolved = _collect_and_resolve(opts, prog='catcifextract')
    if resolved is None:
        # No tags given — extract everything from the specified catcif file(s).
        if opts.m:
            catcif_files = opts.positional
        elif opts.positional and (opts.positional[0].endswith('.catcif') or
                                  opts.positional[0].endswith('.cif')):
            catcif_files = [opts.positional[0]]
        else:
            catcif_files = []

        if not catcif_files:
            parser.print_help(sys.stderr)
            sys.exit(1)

        for catcif_file in catcif_files:
            for structure, tag in get_all_structures(catcif_file):
                _write_structure(tag, structure, opts.o, opts.z)
    else:
        for canonical_tag, structure in _iter_structures(resolved, opts):
            _write_structure(canonical_tag, structure, opts.o, opts.z)
