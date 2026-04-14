"""
settings.py — global configuration for catcif_tools.

Usage
-----
    import catcif_tools
    catcif_tools.catcif_settings.cache_indexes = False

or

    from catcif_tools.settings import catcif_settings
    catcif_settings.cache_indexes = False

The module-level `catcif_settings` instance is a singleton: Python executes
each module exactly once per interpreter session, so every importer shares the
same object regardless of how catcif_tools is imported.
"""


class CatCifSettings:
    """Global settings for catcif_tools behaviour."""

    def __init__(self):
        # Whether to cache loaded indexes in memory.
        self.cache_indexes = True

        # Whether to cache open file pointers for random-access reads.
        self.cache_file_pointers = True

        # Fast cache: small number of entries, short TTL (milliseconds).
        # Intended for files accessed repeatedly in a tight loop.
        self.fast_cache_num = 3
        self.fast_cache_time = 5000

        # Slow cache: larger number of entries, longer TTL (milliseconds).
        # Intended for files that are opened occasionally but repeatedly.
        self.slow_cache_num = 10
        self.slow_cache_time = 100000

        # Maximum number of files to keep in the cache at once.
        self.max_caches = 10


# Singleton instance — shared by all importers.
catcif_settings = CatCifSettings()
