"""
cache.py — in-memory cache for catcif indexes and open file pointers.

The module-level `catcif_cache` dict is a singleton: Python executes each
module exactly once per interpreter session, so all importers share the same
dict regardless of how catcif_tools was imported.

Cache entry layout
------------------
catcif_cache[real_path] = {
    'first_open': float,   # time of first cache entry creation (ms)
    'last_open':  float,   # time of most recent access (ms)
    'num_open':   int,     # number of times this path has been accessed
    'index':      dict,    # catcif index from build_index
    'f_open':     file | None,  # cached open binary file object, if any
}
"""

import os
import time

from .build_index import get_uncached_catcif_index
from .settings import catcif_settings

# Singleton cache dict shared by all importers.
catcif_cache = {}


def get_catcif_index(catcif_file, no_cache=False, instant_cache=False):
    """
    Return the index, a file object, and a close-responsibility flag.

    Parameters
    ----------
    catcif_file : str
        Path to the .catcif file.
    no_cache : bool
        If True, skip all caching entirely. Always builds index from disk.
        Caller is always responsible for closing the returned file object.
    instant_cache : bool
        If True, cache the file pointer immediately without waiting for the
        fast/slow access-count thresholds to be reached.

    Returns
    -------
    index : dict
        The catcif index.
    f_open : file object
        Open binary file object positioned at the start of catcif_file.
    caller_must_close : bool
        True when the caller is responsible for closing f_open.
        False when the file pointer is owned by the cache.
    """
    real_path = os.path.realpath(catcif_file)

    # ── no_cache: build fresh, skip all caching ───────────────────────────────
    if no_cache:
        index = get_uncached_catcif_index(catcif_file)
        return index, open(catcif_file, 'rb'), True

    # ── get or build index ─────────────────────────────────────────────────────
    if real_path in catcif_cache:
        index = catcif_cache[real_path]['index']
    else:
        index = get_uncached_catcif_index(catcif_file)

        if catcif_settings.cache_indexes:
            now = _now_ms()
            catcif_cache[real_path] = {
                'first_open': now,
                'last_open': now,
                'num_open': 0,
                'index': index,
                'f_open': None,
            }

    # ── get file pointer and note access ──────────────────────────────────────
    if real_path in catcif_cache:
        entry = catcif_cache[real_path]
        f_open = entry['f_open'] if entry['f_open'] is not None else open(catcif_file, 'rb')

        # instant_cache bypasses access-count thresholds.
        if instant_cache and entry['f_open'] is None:
            entry['f_open'] = f_open

        cached = _note_access(real_path, f_open)
    else:
        # cache_indexes was False; open without caching.
        f_open = open(catcif_file, 'rb')
        cached = False

    # ── enforce max_caches ─────────────────────────────────────────────────────
    if len(catcif_cache) > catcif_settings.max_caches:
        _trim_cache()

    return index, f_open, not cached


def _now_ms():
    return time.time() * 1000


def _note_access(real_path, f_open):
    """
    Update access statistics for a cached entry; decide whether to own f_open.

    Returns True if f_open is now owned by the cache (caller must NOT close it).
    Returns False otherwise (caller must close it).
    """
    if real_path not in catcif_cache:
        return False

    entry = catcif_cache[real_path]
    entry['num_open'] += 1
    entry['last_open'] = _now_ms()

    if not catcif_settings.cache_indexes or not catcif_settings.cache_file_pointers:
        return False

    # File pointer already owned by cache — caller must not close.
    if entry['f_open'] is not None:
        return True

    elapsed = entry['last_open'] - entry['first_open']

    if (entry['num_open'] == catcif_settings.fast_cache_num
            and elapsed < catcif_settings.fast_cache_time):
        entry['f_open'] = f_open
        return True

    if (entry['num_open'] == catcif_settings.slow_cache_num
            and elapsed < catcif_settings.slow_cache_time):
        entry['f_open'] = f_open
        return True

    return False


def _trim_cache():
    """
    Remove the entry with the oldest last_open until len <= max_caches.
    Closes any cached file pointer for the evicted entry.
    """
    while len(catcif_cache) > catcif_settings.max_caches:
        oldest = min(catcif_cache, key=lambda k: catcif_cache[k]['last_open'])
        entry = catcif_cache.pop(oldest)
        if entry['f_open'] is not None:
            entry['f_open'].close()


def clear_cache():
    """
    Empty the entire cache, closing any owned file pointers.
    """
    for entry in catcif_cache.values():
        if entry['f_open'] is not None:
            entry['f_open'].close()
    catcif_cache.clear()
