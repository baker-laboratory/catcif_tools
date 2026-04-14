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

    Fast pure-Python implementation that parses the CIF text directly without
    biopython. Only the ``_atom_site`` loop is scanned; all other categories
    are skipped entirely. Each data row is expected to occupy a single line
    (true for all standard mmCIF files produced by structure-prediction tools
    and the PDB). Rows with whitespace-containing quoted values in earlier
    columns are skipped gracefully rather than crashing.

    Signifcantly faster than biopython

    Returns the same structure as :func:`get_sequence`.
    """
    lines = structure.splitlines()
    n = len(lines)
    i = 0

    # ── 1. Find the loop_ that opens the _atom_site category ─────────────────
    while i < n:
        if lines[i].strip() == 'loop_':
            j = i + 1
            while j < n and not lines[j].strip():
                j += 1
            if j < n and lines[j].strip().startswith('_atom_site.'):
                i = j   # start collecting headers from here
                break
        i += 1
    else:
        return []

    # ── 2. Collect column headers ─────────────────────────────────────────────
    col_names = []
    while i < n:
        stripped = lines[i].strip()
        if stripped.startswith('_atom_site.'):
            col_names.append(stripped)
            i += 1
        elif not stripped or stripped.startswith('#'):
            i += 1
        else:
            break   # first data row

    if not col_names:
        return []

    # ── 3. Locate the four target column indices ──────────────────────────────
    col_map = {name: idx for idx, name in enumerate(col_names)}
    try:
        idx_entity = col_map['_atom_site.label_entity_id']
        idx_chain  = col_map['_atom_site.label_asym_id']
        idx_comp   = col_map['_atom_site.label_comp_id']
        idx_seq    = col_map['_atom_site.label_seq_id']
    except KeyError:
        return []

    ncols = len(col_names)

    # ── 4. Iterate data rows ──────────────────────────────────────────────────
    # split() is a fast C-level call; all four target columns (entity_id,
    # asym_id, comp_id, seq_id) are always unquoted tokens in real mmCIF files.
    # If a row has a different token count (quoted value with spaces in an
    # earlier column), we skip it — the residue will appear in another atom row.
    chains   = []
    prev_key = None

    while i < n:
        line = lines[i]
        i += 1

        if not line:
            continue

        c0 = line[0]
        if c0 in (' ', '\t'):
            line = line.lstrip()
            if not line:
                continue
            c0 = line[0]

        if c0 == '#':
            continue

        # End-of-loop markers: new key, new loop_, new data_ block, stop_
        if c0 == '_' or line.startswith('loop_') or line.startswith('data_') or line.startswith('stop_'):
            break

        # Semicolon-delimited multi-line value (essentially never in _atom_site)
        if c0 == ';':
            while i < n and not lines[i].startswith(';'):
                i += 1
            i += 1
            continue

        parts = line.split()
        if len(parts) != ncols:
            continue

        seq_id = parts[idx_seq]
        if seq_id == '.':
            continue

        entity_id = parts[idx_entity]
        chain_id  = parts[idx_chain]
        comp_id   = parts[idx_comp]

        key = (entity_id, chain_id)
        if key != prev_key:
            chains.append({
                'entity_id': entity_id,
                'chain_id':  chain_id,
                'residues':  {},
            })
            prev_key = key

        residues = chains[-1]['residues']
        if seq_id not in residues:
            residues[seq_id] = comp_id

    # ── 5. Convert residues to 1-letter sequences ─────────────────────────────
    return [
        {
            'chain_id':  chain['chain_id'],
            'entity_id': chain['entity_id'],
            'sequence':  ''.join(
                _CANONICAL_3TO1.get(comp, 'X') for comp in chain['residues'].values()
            ),
        }
        for chain in chains
    ]


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

    Fast pure-Python implementation; does not require biopython.
    """
    lines = structure.splitlines()
    n = len(lines)
    i = 0

    # ── 1. Find the loop_ that opens the _atom_site category ─────────────────
    while i < n:
        if lines[i].strip() == 'loop_':
            j = i + 1
            while j < n and not lines[j].strip():
                j += 1
            if j < n and lines[j].strip().startswith('_atom_site.'):
                i = j
                break
        i += 1
    else:
        return []

    # ── 2. Collect column headers ─────────────────────────────────────────────
    col_names = []
    while i < n:
        stripped = lines[i].strip()
        if stripped.startswith('_atom_site.'):
            col_names.append(stripped)
            i += 1
        elif not stripped or stripped.startswith('#'):
            i += 1
        else:
            break

    if not col_names:
        return []

    # ── 3. Locate target column indices (auth_seq_id and type_symbol optional) ─
    col_map = {name: idx for idx, name in enumerate(col_names)}
    try:
        idx_entity = col_map['_atom_site.label_entity_id']
        idx_chain  = col_map['_atom_site.label_asym_id']
        idx_comp   = col_map['_atom_site.label_comp_id']
        idx_seq    = col_map['_atom_site.label_seq_id']
    except KeyError:
        return []

    idx_auth_seq = col_map.get('_atom_site.auth_seq_id')
    idx_element  = col_map.get('_atom_site.type_symbol')
    ncols = len(col_names)

    # ── 4. Iterate data rows ──────────────────────────────────────────────────
    chains     = []
    prev_key   = None
    atom_index = 0

    while i < n:
        line = lines[i]
        i += 1

        if not line:
            continue

        c0 = line[0]
        if c0 in (' ', '\t'):
            line = line.lstrip()
            if not line:
                continue
            c0 = line[0]

        if c0 == '#':
            continue

        if c0 == '_' or line.startswith('loop_') or line.startswith('data_') or line.startswith('stop_'):
            break

        if c0 == ';':
            while i < n and not lines[i].startswith(';'):
                i += 1
            i += 1
            continue

        parts = line.split()
        if len(parts) != ncols:
            atom_index += 1
            continue

        entity_id = parts[idx_entity]
        chain_id  = parts[idx_chain]
        comp_id   = parts[idx_comp]
        seq_id    = parts[idx_seq]

        auth_seq_id = parts[idx_auth_seq] if idx_auth_seq is not None else '.'
        element     = parts[idx_element]  if idx_element  is not None else ''

        key = (entity_id, chain_id)
        if key != prev_key:
            chains.append({
                'entity_id': entity_id,
                'chain_id':  chain_id,
                'residues':  {},
            })
            prev_key = key

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

    # ── 5. Sum tokens per chain ───────────────────────────────────────────────
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
