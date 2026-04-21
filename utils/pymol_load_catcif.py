"""
pymol_load_catcif.py — PyMOL plugin for loading structures from .catcif archives.

Install by adding to ~/.pymolrc:
    run /path/to/catcif_tools/utils/pymol_load_catcif.py

Or as a PyMOL plugin via the Plugin Manager.

Usage:
    cmd.load('/path/to/archive.catcif:tag')
    load /path/to/archive.catcif:tag          # PyMOL command line
    load /path/to/archive.catcif              # loads all structures (drag-and-drop)
"""

import os
import sys
import inspect
_here = os.path.dirname(os.path.abspath(inspect.getframeinfo(inspect.currentframe()).filename))
sys.path.append(os.path.abspath(os.path.join(_here, '..')))
del _here

max_structures = 100


def catcif_max_load(n, *, _self=None):
    '''
DESCRIPTION

    Set the maximum number of structures loaded when opening a .catcif
    archive without a specific tag (e.g. drag-and-drop).

USAGE

    catcif_max_load n

ARGUMENTS

    n = integer: maximum number of structures to load, or 0 for no limit

EXAMPLES

    catcif_max_load 500
    catcif_max_load 0
    '''
    global max_structures
    max_structures = int(n)


def _load_catcif(filename, object='', state=0, finish=1, discrete=-1,
                 quiet=1, multiplex=None, zoom=-1, **kw):
    """Load a single structure from a .catcif archive into PyMOL."""
    import catcif_tools
    from pymol import cmd

    cif_text = catcif_tools.get_structure(filename)

    if not object:
        _, tag = catcif_tools.split_catcif_tag(filename)
        object = tag.replace('/', '_').replace(' ', '_')

    return cmd.load_raw(cif_text, 'cif', object, state,
                        finish=finish, discrete=discrete,
                        quiet=quiet, multiplex=multiplex, zoom=zoom)


def _load_all_catcif(filename, state=0, finish=1, discrete=-1,
                     quiet=1, multiplex=None, zoom=-1, **kw):
    """Load every structure from a .catcif archive into PyMOL."""
    import catcif_tools
    from pymol import cmd

    for i, (cif_text, tag) in enumerate(catcif_tools.get_all_structures(filename)):
        if max_structures and i >= max_structures:
            print(f' catcif: limit of {max_structures} structures reached; use catcif_max_load to increase')
            break
        object = tag.replace('/', '_').replace(' ', '_')
        cmd.load_raw(cif_text, 'cif', object, state,
                     finish=finish, discrete=discrete,
                     quiet=quiet, multiplex=multiplex, zoom=zoom)


def _make_load_wrapper(original_load):
    def load(filename, object='', state=0, format='', finish=1,
             discrete=-1, quiet=1, multiplex=None, zoom=-1, partial=0,
             mimic=1, object_props=None, atom_props=None, **kw):
        if isinstance(filename, str) and '.catcif:' in filename:
            return _load_catcif(filename, object=object, state=state,
                                finish=finish, discrete=discrete, quiet=quiet,
                                multiplex=multiplex, zoom=zoom)
        if isinstance(filename, str) and filename.endswith('.catcif'):
            return _load_all_catcif(filename, state=state, finish=finish,
                                    discrete=discrete, quiet=quiet,
                                    multiplex=multiplex, zoom=zoom)
        return original_load(filename, object=object, state=state,
                             format=format, finish=finish, discrete=discrete,
                             quiet=quiet, multiplex=multiplex, zoom=zoom,
                             partial=partial, mimic=mimic,
                             object_props=object_props, atom_props=atom_props,
                             **kw)
    load.__doc__ = original_load.__doc__
    return load


def __init_plugin__(app=None):
    _install()


def _install():
    from pymol import cmd

    if getattr(cmd.load, '_catcif_patched', False):
        return

    wrapped = _make_load_wrapper(cmd.load)
    wrapped._catcif_patched = True
    cmd.load = wrapped
    cmd.extend('load', wrapped)
    cmd.extend('catcif_max_load', catcif_max_load)


_install()
