"""
conversion.py — format conversion utilities for catcif_tools.

Heavy dependencies (e.g. biopython) are imported lazily inside each function
so they are never required at import time.
"""


def pdb_to_cif(pdb_string, tag):
    """
    Convert a PDB-format string to a CIF string without biopython.

    Parses ATOM/HETATM records using fixed-width column slicing (no tokenizer,
    no object hierarchy) and writes only the _atom_site loop and _entry.id —
    the categories catcif_tools actually uses.

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
    """
    # ── 1. Parse ATOM / HETATM records ───────────────────────────────────────
    groups    = []
    serials   = []
    elems     = []
    atom_nms  = []
    alt_locs  = []
    res_names = []
    chain_ids = []
    res_seqs  = []
    ins_codes = []
    xs = []; ys = []; zs = []
    occs      = []
    bfacs     = []
    model_nums = []

    model_num = 1

    for line in pdb_string.splitlines():
        if len(line) < 6:
            continue
        rec = line[:6]
        if rec == 'MODEL ':
            try:
                model_num = int(line[6:].strip())
            except ValueError:
                model_num += 1
            continue
        if rec not in ('ATOM  ', 'HETATM'):
            continue
        if len(line) < 54:          # need at least through Z coordinate
            continue

        groups.append(rec.strip())
        serials.append(line[6:11].strip() or '0')
        atom_nms.append(line[12:16].strip())
        alt_locs.append(line[16].strip())
        res_names.append(line[17:20].strip())
        chain_ids.append(line[21].strip() or 'A')
        res_seqs.append(line[22:26].strip() or '0')
        ins_codes.append(line[26].strip() if len(line) > 26 else '')
        xs.append(line[30:38].strip())
        ys.append(line[38:46].strip())
        zs.append(line[46:54].strip())
        occs.append( line[54:60].strip() if len(line) > 54 else '')
        bfacs.append(line[60:66].strip() if len(line) > 60 else '')

        if len(line) > 76:
            elem = line[76:78].strip()
        else:
            # Derive from atom name: drop leading digits, take first letter.
            nm = line[12:16].strip().lstrip('0123456789')
            elem = nm[:1] if nm else ''
        elems.append(elem or '?')

        model_nums.append(str(model_num))

    if not groups:
        return f'data_{tag}\n#\n_entry.id {tag}\n#\n'

    # ── 2. Assign entity IDs — one per unique chain, in appearance order ──────
    entity_map: dict[str, str] = {}
    for cid in chain_ids:
        if cid not in entity_map:
            entity_map[cid] = str(len(entity_map) + 1)

    # ── 3. Assign label_seq_id — sequential within each chain, increments
    #        whenever (res_seq, ins_code) changes ──────────────────────────────
    label_seq_ids: list[str] = []
    seq_state: dict[str, tuple] = {}   # chain_id -> (res_seq, ins_code, seq_id)

    for cid, rseq, ins in zip(chain_ids, res_seqs, ins_codes):
        if cid not in seq_state:
            seq_state[cid] = (rseq, ins, 1)
            label_seq_ids.append('1')
        else:
            last_rseq, last_ins, last_id = seq_state[cid]
            if (rseq, ins) == (last_rseq, last_ins):
                label_seq_ids.append(str(last_id))
            else:
                new_id = last_id + 1
                seq_state[cid] = (rseq, ins, new_id)
                label_seq_ids.append(str(new_id))

    # ── 4. Normalise optional fields ─────────────────────────────────────────
    alt_display = [v or '.' for v in alt_locs]
    ins_display = [v or '?' for v in ins_codes]
    occs        = [v or '1.00' for v in occs]
    bfacs       = [v or '0.00' for v in bfacs]
    entity_ids  = [entity_map[cid] for cid in chain_ids]

    # ── 5. Write output ───────────────────────────────────────────────────────
    out = [
        f'data_{tag}',
        '#',
        f'_entry.id {tag}',
        '#',
        'loop_',
        '_atom_site.group_PDB',
        '_atom_site.id',
        '_atom_site.type_symbol',
        '_atom_site.label_atom_id',
        '_atom_site.label_alt_id',
        '_atom_site.label_comp_id',
        '_atom_site.label_asym_id',
        '_atom_site.label_entity_id',
        '_atom_site.label_seq_id',
        '_atom_site.pdbx_PDB_ins_code',
        '_atom_site.Cartn_x',
        '_atom_site.Cartn_y',
        '_atom_site.Cartn_z',
        '_atom_site.occupancy',
        '_atom_site.B_iso_or_equiv',
        '_atom_site.auth_seq_id',
        '_atom_site.auth_asym_id',
        '_atom_site.auth_atom_id',
        '_atom_site.auth_comp_id',
        '_atom_site.pdbx_PDB_model_num',
    ]

    for i in range(len(groups)):
        out.append(
            f'{groups[i]} {serials[i]} {elems[i]} {atom_nms[i]} '
            f'{alt_display[i]} {res_names[i]} {chain_ids[i]} '
            f'{entity_ids[i]} {label_seq_ids[i]} {ins_display[i]} '
            f'{xs[i]} {ys[i]} {zs[i]} '
            f'{occs[i]} {bfacs[i]} '
            f'{res_seqs[i]} {chain_ids[i]} {atom_nms[i]} {res_names[i]} '
            f'{model_nums[i]}'
        )

    out.append('#')
    return '\n'.join(out) + '\n'


def pdb_to_cif_bio(pdb_string, tag):
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
