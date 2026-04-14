"""catciffromfiles — build a .catcif file from individual structure files."""

import argparse
import gzip
import os
import stat
import sys
import textwrap
from signal import signal, SIGPIPE, SIG_DFL

signal(SIGPIPE, SIG_DFL)

from .catcif_tools import to_catcif_string
from .conversion import pdb_to_cif
from .structure import compress_structure


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


# Ordered longest-first so .cif.gz is checked before .cif
_SUPPORTED_EXTS = ('.cif.gz', '.cif', '.pdb.gz', '.pdb')
_PDB_EXTS       = ('.pdb.gz', '.pdb')


def _tag_from_path(path):
    """
    Derive the structure tag from a file path by stripping the file extension.

    Returns the tag string, or raises ValueError for unsupported extensions.
    """
    basename = os.path.basename(path)
    for ext in _SUPPORTED_EXTS:
        if basename.endswith(ext):
            return basename[:-len(ext)]
    raise ValueError(
        f'catciffromfiles: unsupported file type: {path!r}\n'
        f'  supported extensions: {", ".join(_SUPPORTED_EXTS)}'
    )


def _read_text(path):
    """Read a plain or gzip-compressed text file and return its contents as a string."""
    if path.endswith('.gz'):
        with gzip.open(path, 'rt', encoding='ascii') as f:
            return f.read()
    else:
        with open(path, encoding='ascii') as f:
            return f.read()


def _process_file(path):
    """
    Read path and return a CIF string suitable for appending to a .catcif file.

    The data_ block header is renamed to the basename-derived tag.
    Raises ImportError if a .pdb/.pdb.gz is given but biopython is not installed.
    Raises ValueError for unsupported extensions or empty files.
    """
    tag = _tag_from_path(path)
    content = _read_text(path)

    if not content.strip():
        raise ValueError(f'catciffromfiles: empty file: {path!r}')

    basename = os.path.basename(path)
    for ext in _PDB_EXTS:
        if basename.endswith(ext):
            content = pdb_to_cif(content, tag)
            return to_catcif_string(content, tag)

    return to_catcif_string(content, tag, add_header=True)


def main():
    parser = argparse.ArgumentParser(
        prog='catciffromfiles',
        description='Build a .catcif file from individual .cif or .cif.gz files.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            examples:
              catciffromfiles a.cif b.cif.gz > out.catcif
              ls *.cif | catciffromfiles > out.catcif
              ls *.cif | catciffromfiles extra.cif > out.catcif
        """),
    )
    parser.add_argument(
        '-z', action='store_true',
        help='write gzip-compressed output',
    )
    parser.add_argument(
        'files', nargs='*', metavar='FILE',
        help='.cif or .cif.gz files to include',
    )
    opts = parser.parse_args()

    # Collect file paths from stdin (if piped) and positional args
    stdin_files = []
    if stat.S_ISFIFO(os.stat('/dev/stdin').st_mode):
        for line in sys.stdin:
            line = line.strip()
            if line:
                stdin_files.extend(line.split())

    all_files = stdin_files + opts.files

    if not all_files:
        parser.print_help(sys.stderr)
        sys.exit(1)

    for path in all_files:
        try:
            structure = _process_file(path)
        except (ImportError, ValueError) as exc:
            eprint(exc)
            sys.exit(1)
        if opts.z:
            sys.stdout.buffer.write(compress_structure(structure))
            sys.stdout.buffer.flush()
        else:
            sys.stdout.write(structure)
            sys.stdout.flush()
