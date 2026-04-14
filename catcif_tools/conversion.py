"""
conversion.py — format conversion utilities for catcif_tools.

Heavy dependencies (e.g. biopython) are imported lazily inside each function
so they are never required at import time.
"""


def pdb_to_cif(pdb_string, tag):
    """
    Convert a PDB-format string to a CIF string using biopython.

    Parameters
    ----------
    pdb_string : str
        PDB-format structure text.
    tag : str
        Tag used as the structure identifier; becomes the data_ block name.

    Returns
    -------
    str
        CIF text with a data_ block header.

    Raises
    ------
    ImportError
        If biopython is not installed.
    """
    try:
        from Bio.PDB import PDBParser, MMCIFIO
    except ImportError:
        raise ImportError(
            "biopython is required for PDB-to-CIF conversion: pip install biopython"
        )
    import io

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(tag, io.StringIO(pdb_string))

    mmcif_io = MMCIFIO()
    mmcif_io.set_structure(structure)
    buf = io.StringIO()
    mmcif_io.save(buf)

    return buf.getvalue()
