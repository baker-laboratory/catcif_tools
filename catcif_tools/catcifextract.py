"""catcifextract — write each structure from a .catcif file to its own .cif file."""

import os
import sys
from signal import signal, SIGPIPE, SIG_DFL

signal(SIGPIPE, SIG_DFL)

from .catcifslice import _build_parser, _collect_and_resolve, _iter_structures
from .structure import compress_structure


def main():
    parser = _build_parser(
        prog='catcifextract',
        description='Write each structure from a .catcif file to its own .cif file.',
        epilog="""\
            examples:
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
    resolved = _collect_and_resolve(opts, prog='catcifextract')
    if not resolved:
        parser.print_help(sys.stderr)
        sys.exit(1)
    os.makedirs(opts.o, exist_ok=True)
    for canonical_tag, structure in _iter_structures(resolved, opts):
        if opts.z:
            out_path = os.path.join(opts.o, f'{canonical_tag}.cif.gz')
            with open(out_path, 'wb') as f:
                f.write(compress_structure(structure))
        else:
            out_path = os.path.join(opts.o, f'{canonical_tag}.cif')
            with open(out_path, 'w') as f:
                f.write(structure)
