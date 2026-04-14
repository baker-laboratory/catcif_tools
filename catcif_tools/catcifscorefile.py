"""catcifscorefile — extract scores from a .catcif file to a score file."""

import argparse
import csv
import os
import sys
from signal import signal, SIGPIPE, SIG_DFL

signal(SIGPIPE, SIG_DFL)

from .scores import parse_score_file


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def _write_scorefile(path, cols, rows, sep, space_rep):
    with open(path, 'w', newline='') as f:
        if sep == ',':
            writer = csv.writer(f)
            writer.writerow(cols)
            for row_dict in rows:
                writer.writerow([row_dict.get(col, '') for col in cols])
        else:
            f.write(sep.join(cols) + '\n')
            for row_dict in rows:
                values = []
                for col in cols:
                    val = row_dict.get(col, '')
                    if space_rep and ' ' in val:
                        val = val.replace(' ', space_rep)
                    values.append(val)
                f.write(sep.join(values) + '\n')


def main():
    parser = argparse.ArgumentParser(
        prog='catcifscorefile',
        description='Extract scores from a .catcif file to a score file.',
        epilog='Output is written to <basename>.sc (or .csv) in the current directory.',
    )
    parser.add_argument(
        'catcif_file',
        help='.catcif file to read scores from',
    )
    parser.add_argument(
        '-o', metavar='DIR', default=None,
        help='write scorefile to this directory instead of the current directory',
    )
    parser.add_argument(
        '-g', action='store_true',
        help='write scorefile using original (pre-deduplication) tags',
    )
    parser.add_argument(
        '-c', action='store_true',
        help='write a CSV file instead of a space-separated scorefile',
    )
    parser.add_argument(
        '--space-rep', metavar='CHAR', default='_',
        help='replace spaces in values with this character in space-separated '
             'output (default: "_"; pass empty string to disable)',
    )

    opts = parser.parse_args()

    if not os.path.isfile(opts.catcif_file):
        eprint(f'catcifscorefile: File not found: {opts.catcif_file}')
        sys.exit(1)

    # Collect all score rows, grouping by the tuple of column names.
    # Structures with different score sets end up in separate groups.
    score_groups = {}  # tuple(col names) -> list of row dicts
    for scores in parse_score_file(opts.catcif_file, preserve_tags=opts.g):
        cols = tuple(scores.keys())
        if cols not in score_groups:
            score_groups[cols] = []
        score_groups[cols].append(scores)

    # Determine output directory and base filename.
    out_dir = opts.o if opts.o is not None else os.getcwd()
    os.makedirs(out_dir, exist_ok=True)

    basename = os.path.basename(opts.catcif_file)
    if basename.endswith('.catcif'):
        basename = basename[:-len('.catcif')]
    ext = '.csv' if opts.c else '.sc'
    sep = ',' if opts.c else ' '

    # Write output file(s).
    if len(score_groups) == 1:
        cols, rows = next(iter(score_groups.items()))
        out_path = os.path.join(out_dir, basename + ext)
        _write_scorefile(out_path, cols, rows, sep, opts.space_rep)
    else:
        eprint('catcifscorefile: Warning: multiple score column sets found; '
               'producing multiple files')
        for i, (cols, rows) in enumerate(score_groups.items(), 1):
            eprint(f'  set {i}: {" ".join(cols)}')
            out_path = os.path.join(out_dir, basename + str(i) + ext)
            _write_scorefile(out_path, cols, rows, sep, opts.space_rep)
