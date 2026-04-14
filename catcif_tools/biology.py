"""
biology.py — biological sequence utilities for catcif_tools.

Heavy dependencies (e.g. biopython) are imported lazily inside each function
so they are never required at import time.
"""

_CANONICAL_3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}


def get_sequence(structure):
    """
    Return the 1-letter amino-acid sequence for each chain in a CIF structure.

    A new chain is started whenever ``label_asym_id`` (chain ID) or
    ``label_entity_id`` changes relative to the previous atom row.  A gap
    in residue numbering within the same chain ID and entity ID is treated as
    a chain break in the structure but NOT as a new chain here — the sequence
    is kept contiguous.

    Residues outside the canonical 20 amino acids are represented as ``X``.
    Atoms with ``label_seq_id == '.'`` (typically solvent) are skipped.

    Parameters
    ----------
    structure : str
        CIF text for a single structure.

    Returns
    -------
    list of dict
        Each dict has keys:
          ``chain_id``  — ``label_asym_id`` value
          ``entity_id`` — ``label_entity_id`` value
          ``sequence``  — 1-letter sequence string

    Raises
    ------
    ImportError
        If biopython is not installed.
    """
    try:
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    except ImportError:
        raise ImportError(
            "biopython is required for get_sequence: pip install biopython"
        )
    import io

    mmcif_dict = MMCIF2Dict(io.StringIO(structure))

    entity_ids = mmcif_dict.get('_atom_site.label_entity_id', [])
    chain_ids  = mmcif_dict.get('_atom_site.label_asym_id',   [])
    comp_ids   = mmcif_dict.get('_atom_site.label_comp_id',   [])
    seq_ids    = mmcif_dict.get('_atom_site.label_seq_id',    [])

    # Each chain is a consecutive run of identical (entity_id, chain_id).
    # Within a chain, residues are deduplicated by label_seq_id (first
    # occurrence wins); atoms with seq_id '.' (solvent) are skipped.
    chains = []   # list of {'entity_id', 'chain_id', 'residues': OrderedDict}
    prev_key = None

    for entity_id, chain_id, comp_id, seq_id in zip(
            entity_ids, chain_ids, comp_ids, seq_ids):

        if seq_id == '.':
            continue

        key = (entity_id, chain_id)
        if key != prev_key:
            chains.append({
                'entity_id': entity_id,
                'chain_id':  chain_id,
                'residues':  {},   # seq_id -> comp_id, insertion-ordered
            })
            prev_key = key

        residues = chains[-1]['residues']
        if seq_id not in residues:
            residues[seq_id] = comp_id

    result = []
    for chain in chains:
        seq = ''.join(
            _CANONICAL_3TO1.get(comp_id, 'X')
            for comp_id in chain['residues'].values()
        )
        result.append({
            'chain_id':  chain['chain_id'],
            'entity_id': chain['entity_id'],
            'sequence':  seq,
        })

    return result


def chain_token_lengths(structure):
    """
    Return the number of tokens in each chain of a CIF structure.

    Token definition:
      - Canonical amino acid residue: 1 token.
      - Anything else (non-standard residue, ligand, solvent): the number of
        non-hydrogen atoms in that residue/molecule.

    Chain splitting follows the same rule as :func:`get_sequence`: a new chain
    starts whenever ``label_asym_id`` or ``label_entity_id`` changes.

    For polymer residues (``label_seq_id != '.'``) atoms are grouped by
    ``label_seq_id``.  For non-polymer entries (``label_seq_id == '.'``) atoms
    are grouped by ``auth_seq_id``; if that is also ``'.'`` each atom is
    treated as its own residue.

    Parameters
    ----------
    structure : str
        CIF text for a single structure.

    Returns
    -------
    list of dict
        Each dict has keys:
          ``chain_id``  — ``label_asym_id`` value
          ``entity_id`` — ``label_entity_id`` value
          ``tokens``    — total token count for the chain (int)

    Raises
    ------
    ImportError
        If biopython is not installed.
    """
    try:
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    except ImportError:
        raise ImportError(
            "biopython is required for chain_token_lengths: pip install biopython"
        )
    import io

    mmcif_dict = MMCIF2Dict(io.StringIO(structure))

    entity_ids   = mmcif_dict.get('_atom_site.label_entity_id', [])
    chain_ids    = mmcif_dict.get('_atom_site.label_asym_id',   [])
    comp_ids     = mmcif_dict.get('_atom_site.label_comp_id',   [])
    seq_ids      = mmcif_dict.get('_atom_site.label_seq_id',    [])
    auth_seq_ids = mmcif_dict.get('_atom_site.auth_seq_id',     [])
    elements     = mmcif_dict.get('_atom_site.type_symbol',     [])

    # Pad auth_seq_ids and elements to the atom count so zip is safe when absent.
    n = len(entity_ids)
    if len(auth_seq_ids) < n:
        auth_seq_ids = auth_seq_ids + ['.'] * (n - len(auth_seq_ids))
    if len(elements) < n:
        elements = elements + [''] * (n - len(elements))

    chains = []
    prev_chain_key = None
    atom_index = 0  # fallback unique key when both seq_ids are '.'

    for entity_id, chain_id, comp_id, seq_id, auth_seq_id, element in zip(
            entity_ids, chain_ids, comp_ids, seq_ids, auth_seq_ids, elements):

        chain_key = (entity_id, chain_id)
        if chain_key != prev_chain_key:
            chains.append({
                'entity_id': entity_id,
                'chain_id':  chain_id,
                'residues':  {},  # res_key -> {'comp_id': str, 'non_h_atoms': int}
            })
            prev_chain_key = chain_key

        # Determine the residue grouping key for this atom.
        if seq_id != '.':
            res_key = seq_id
        elif auth_seq_id != '.':
            res_key = ('dot', auth_seq_id)
        else:
            res_key = ('atom', atom_index)

        residues = chains[-1]['residues']
        if res_key not in residues:
            residues[res_key] = {'comp_id': comp_id, 'non_h_atoms': 0}

        if element.upper() not in ('H', 'D'):
            residues[res_key]['non_h_atoms'] += 1

        atom_index += 1

    result = []
    for chain in chains:
        tokens = 0
        for res_info in chain['residues'].values():
            if res_info['comp_id'] in _CANONICAL_3TO1:
                tokens += 1
            else:
                tokens += res_info['non_h_atoms']
        result.append({
            'chain_id':  chain['chain_id'],
            'entity_id': chain['entity_id'],
            'tokens':    tokens,
        })

    return result
