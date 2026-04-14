"""catcifsplit — split a .catcif file into multiple smaller files."""

import argparse
import itertools
import os
import string
import sys
import textwrap
from signal import signal, SIGPIPE, SIG_DFL

signal(SIGPIPE, SIG_DFL)

from .biology import chain_token_lengths
from .cache import get_catcif_index
from .structure import get_all_structures

_LETTERS = string.ascii_lowercase
_MAX_OPEN = 512 - 30  # 482


def _suffix_gen():
    """Yield xaa, xab, ..., xzz, xaaa, xaab, ..., xzzz, ... indefinitely."""
    length = 2
    while True:
        for combo in itertools.product(_LETTERS, repeat=length):
            yield 'x' + ''.join(combo)
        length += 1


def _outpath(out_dir, prefix, suffix):
    return os.path.join(out_dir, prefix + suffix)


# ---------------------------------------------------------------------------
# Plain mode
# ---------------------------------------------------------------------------

def _split_plain(catcif_file, n, out_dir, prefix=''):
    suffixes = _suffix_gen()
    f_out = None
    count = 0

    for structure, _ in get_all_structures(catcif_file):
        if count % n == 0:
            if f_out is not None:
                f_out.close()
            f_out = open(_outpath(out_dir, prefix, next(suffixes)), 'w')
        f_out.write(structure)
        count += 1

    if f_out is not None:
        f_out.close()


# ---------------------------------------------------------------------------
# Shuffle mode (-s)
# ---------------------------------------------------------------------------

def _split_shuffle(catcif_file, n, out_dir):
    index, f_open, caller_must_close = get_catcif_index(catcif_file,
                                                         instant_cache=True)
    if caller_must_close:
        f_open.close()

    total = len(index['index'])
    if total == 0:
        return

    k = (total + n - 1) // n  # ceil(total / n)
    suffixes = list(itertools.islice(_suffix_gen(), k))

    if k <= _MAX_OPEN:
        files = [open(_outpath(out_dir, '', suffixes[j]), 'w') for j in range(k)]
        try:
            for i, (structure, _) in enumerate(get_all_structures(catcif_file)):
                files[i % k].write(structure)
        finally:
            for f in files:
                f.close()
    else:
        for batch_start in range(0, k, _MAX_OPEN):
            batch_end = min(batch_start + _MAX_OPEN, k)
            files = [open(_outpath(out_dir, '', suffixes[j]), 'w')
                     for j in range(batch_start, batch_end)]
            try:
                for i, (structure, _) in enumerate(
                        get_all_structures(catcif_file)):
                    file_idx = i % k
                    if batch_start <= file_idx < batch_end:
                        files[file_idx - batch_start].write(structure)
            finally:
                for f in files:
                    f.close()


# ---------------------------------------------------------------------------
# Target-length mode (-t) and bucket mode (-b)
# ---------------------------------------------------------------------------

def _split_by_prefix(catcif_file, n, out_dir, prefix_fn):
    """
    Generic single-pass split where each structure is assigned a filename
    prefix by calling prefix_fn(structure).

    Maintains one open file handle per active prefix; opens new files on
    demand when a prefix is first seen or when its current file fills up.
    """
    # Per-prefix state: suffix generator and current (file_handle, count)
    suffix_gens = {}   # prefix -> iterator from _suffix_gen()
    open_files  = {}   # prefix -> open file handle

    try:
        for structure, _ in get_all_structures(catcif_file):
            pfx = prefix_fn(structure)

            if pfx not in suffix_gens:
                suffix_gens[pfx] = _suffix_gen()
                open_files[pfx] = (open(_outpath(out_dir, pfx,
                                                  next(suffix_gens[pfx])), 'w'),
                                   0)

            f, count = open_files[pfx]
            if count == n:
                f.close()
                f = open(_outpath(out_dir, pfx, next(suffix_gens[pfx])), 'w')
                count = 0

            f.write(structure)
            open_files[pfx] = (f, count + 1)
    finally:
        for f, _ in open_files.values():
            f.close()


def _target_prefix(structure):
    chains = chain_token_lengths(structure)
    return 'len' + '_'.join(str(ch['tokens']) for ch in chains)


def _bucket_prefix(structure, bucket_step):
    chains = chain_token_lengths(structure)
    total = sum(ch['tokens'] for ch in chains)
    bucketed = (total // bucket_step) * bucket_step
    return 'len' + str(bucketed)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        prog='catcifsplit',
        description='Split a .catcif file into multiple smaller files.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            examples:
              catcifsplit my.catcif 1000
              catcifsplit -s my.catcif 1000
              catcifsplit -t my.catcif 1000
              catcifsplit -b my.catcif 1000
              catcifsplit -b --bucket-step 32 -o out_dir/ my.catcif 1000
        """),
    )
    parser.add_argument('catcif_file', metavar='FILE', help='.catcif file to split')
    parser.add_argument('n', metavar='N', type=int,
                        help='maximum number of structures per output file')
    parser.add_argument('-o', metavar='DIR', default='.',
                        help='output directory (default: current directory)')

    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('-s', action='store_true',
                      help='round-robin shuffle across output files')
    mode.add_argument('-t', action='store_true',
                      help='prefix output files with per-chain token lengths '
                           '(requires biopython)')
    mode.add_argument('-b', action='store_true',
                      help='prefix output files with total bucketed token count '
                           '(requires biopython)')

    parser.add_argument('--bucket-step', type=int, default=16, metavar='STEP',
                        help='rounding step for -b mode (default: 16)')

    opts = parser.parse_args()

    os.makedirs(opts.o, exist_ok=True)

    if opts.s:
        _split_shuffle(opts.catcif_file, opts.n, opts.o)
    elif opts.t:
        _split_by_prefix(opts.catcif_file, opts.n, opts.o, _target_prefix)
    elif opts.b:
        _split_by_prefix(opts.catcif_file, opts.n, opts.o,
                         lambda s: _bucket_prefix(s, opts.bucket_step))
    else:
        _split_plain(opts.catcif_file, opts.n, opts.o)
