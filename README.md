# catcif_tools

Tools for working with concatenated CIF files (`.catcif`).

## Motivation

Computational biology often relies on storing thousands to millions of structures. Written as individual
files, these structures can clog up file systems. Easy-to-use archive formats are needed.

Here we present the catcif file format. It's effectively just catted cif files. If all structures are
uncompressed, a catcif file is a valid cif file (thanks RCSB for making that one possible).

Now, instead of storing thousands of files from a run, a single catcif file can be stored.

## The catcif File Format

A `.catcif` file is a plain concatenation of standard CIF structures. This works
without any special framing because every CIF block begins with a `data_<tag>` line —
that line acts as a natural delimiter. A reader can split the file on `\ndata_` and
recover each structure exactly.

**Plain-text structures** are written as-is. **Compressed structures** are written as
raw gzip members (identifiable by the two-byte magic `\x1f\x8b`). Plain and compressed
members can be interleaved freely in the same file; catcif_tools detects the format of
each member automatically and decompresses on demand.

**Combining files** is as simple as `cat`:
```bash
cat runs/*/out.catcif > combined.catcif
```

**Note!** Many tools do not output cif files with the `data_` header set to the name
of the structure. Solve this with:
```bash
catciffromfiles *.cif > combined.catcif
```

### Index files

The first time a `.catcif` file is opened, catcif_tools builds an in-memory index that
maps each structure's tag to its byte offset in the file. This makes random access
(fetching a single structure out of millions) fast: one seek, one read.

If two structures share the same tag, the second occurrence is renamed `foo_1`, the
third `foo_2`, and so on. The on-disk names are preserved in `orig_tags` alongside the
deduplicated canonical names in the index.

### Score files

Scores are stored directly inside each CIF structure as a `_catcif_scores` key-value
block:

```
data_my_design
#
_catcif_scores.total_score   -123.45
_catcif_scores.rmsd          1.234
#
_atom_site.group_PDB ...
```

This means a single `.catcif` file holds both structures and scores; there is no
separate `.sc` file to keep synchronized. The score block is written at the very top of
each structure so that score-only scans (e.g., `catcifscorefile`) can stop reading each
member as soon as they find it.


## Quick API Demo

```python
import catcif_tools

# Read one structure by name (returns a cif string)
structure = catcif_tools.get_structure("my_design", catcif_file="designs.catcif")

# The same call using the embedded path:tag syntax (no catcif_file argument needed)
structure = catcif_tools.get_structure("designs.catcif:my_design")

# Append a structure to a .catcif file
catcif_tools.append_to_catcif_file("designs.catcif", cif_str, tag="my_design")

# Append with scores
catcif_tools.append_to_catcif_file(
    "designs.catcif", cif_str, tag="my_design",
    scores={"total_score": -123.45, "rmsd": 1.23},
)

# Append a gzip-compressed structure
catcif_tools.append_to_catcif_file(
    "designs.catcif", cif_str, tag="my_design", compress=True
)

# Read multiple structures; tags carry the file path so no catcif_file arg is needed
tags = catcif_tools.get_tags("designs.catcif")                      # see full API below
path_tags = [f"designs.catcif:{t}" for t in tags[:10]]
for structure in catcif_tools.get_structures(path_tags):
    scores = catcif_tools.get_scores(structure)
    print(scores["total_score"])
```


## Command-line Demo

```bash
# Make a catcif file from .cif files (also accepts .pdb with biopython installed)
catciffromfiles *.cif > my.catcif

# Ask what's in a catcif file
catcifls my.catcif

# Ask how many structures are in a catcif file
catcifls my.catcif | wc -l

# Extract all structures from a catcif file to individual .cif files
catcifextract my.catcif

# Extract the first 10 structures from a catcif file
catcifls my.catcif | head -n 10 | catcifextract my.catcif

# Extract a random 10 structures from a catcif file
catcifls my.catcif | shuf | head -n 10 | catcifextract my.catcif

# Extract a specific structure from a catcif file
catcifextract my.catcif name_of_structure_0001

# Produce a scorefile from a catcif file
catcifscorefile my.catcif

# Combine catcif files
cat 1.catcif 2.catcif 3.catcif > my.catcif

# Ensure all structures in a catcif file have unique names
catcifls my.catcif | catcifrename my.catcif > uniq.catcif

# Remove _0001 from all the names in a catcif file
catcifls my.catcif | sed 's/_0001$//g' | catcifrename my.catcif > renamed.catcif

# Make a new catcif file with the first 10 structures
catcifls my.catcif | head -n 10 | catcifslice my.catcif > new.catcif

# Get all sequences from a catcif file
catcifsequence my.catcif > my.seq

# Get sequences with chain IDs labeled
catcifsequence -c my.catcif > my.seq

# Split a catcif file into groups of 100
catcifsplit my.catcif 100

# Split a catcif file, shuffling structures round-robin across output files
catcifsplit -s my.catcif 100

# Split by per-chain token length (requires biopython)
catcifsplit -t my.catcif 1000

# Split into token-length buckets of size 32 (requires biopython)
catcifsplit -b --bucket-step 32 my.catcif 1000
```

### A realistic filter-and-extract workflow

```bash
# Collect results from many parallel runs
cat runs/*/out.catcif > combined.catcif

# Extract a scorefile
catcifscorefile combined.catcif   # writes combined.sc

# Filter by score and extract matching structures
cat combined.sc | awk '{if ($1 < -500) {print $NF}}' \
    | catcifextract combined.catcif

# Or slice into a new catcif file for downstream processing
cat combined.sc | awk '{if ($1 < -500) {print $NF}}' \
    | catcifslice combined.catcif > filtered.catcif
```

### Using the long path:tag format

`catcifls -l` emits tags in `path:tag` form. These can be passed directly to any tool
that reads tags — no need to specify the catcif file separately:

```bash
# List tags with embedded paths from multiple catcif files
catcifls -l run1.catcif run2.catcif > all_tags.list

# Slice across multiple files in a single command (tags carry the file path)
cat all_tags.list | shuf | head -n 50 | catcifslice > subset.catcif

# Extract specific structures from wherever they live
cat all_tags.list | grep "my_design" | catcifextract
```


## API Reference

All public symbols are importable directly from `catcif_tools`:

```python
import catcif_tools
catcif_tools.get_structure(...)
```

---

### Reading structures

#### `get_structure(tag, catcif_file=None, no_cache=False, instant_cache=False, preserve_tags=False) → str`

Return the CIF text for a single structure.

`tag` may be a bare tag name (requires `catcif_file`) or a `path/to/file.catcif:tag`
string. If both are supplied they must resolve to the same file.

```python
structure = catcif_tools.get_structure("my_design", catcif_file="designs.catcif")
structure = catcif_tools.get_structure("designs.catcif:my_design")
```

---

#### `get_structures(tags, catcif_file=None, no_cache=False, instant_cache=False, preserve_tags=False) → Iterator[str]`

Yield CIF text for each tag in the iterable `tags`. All tags must reference the same
catcif file (either via the `path.catcif:tag` syntax or the `catcif_file` argument).
Opens the file once and reads structures in a single forward pass.

```python
tags = ["alpha", "beta", "gamma"]
for structure in catcif_tools.get_structures(tags, catcif_file="designs.catcif"):
    print(catcif_tools.get_scores(structure))
```

---

#### `get_all_structures(catcif_file, preserve_tags=False, dont_rename_structure=False) → Iterator[tuple[str, str]]`

Yield `(structure, tag)` for every structure in `catcif_file`, in file order.

Reads sequentially with no backward seeks — use this when you want to process every
structure in a file rather than a specific subset.

```python
for structure, tag in catcif_tools.get_all_structures("designs.catcif"):
    print(tag, catcif_tools.get_scores(structure))
```

---

### Writing structures

#### `append_to_catcif_file(catcif_file, cif_str, tag, scores=None, add_header=False, compress=False)`

Append a single CIF structure to a `.catcif` file on disk (created if it does not
exist).

```python
catcif_tools.append_to_catcif_file(
    "out.catcif", cif_str, tag="design_0001",
    scores={"total_score": -512.3},
)
```

---

#### `append_to_catcif_file_open(catcif_file_pointer, cif_str, tag, scores=None, add_header=False, compress=False)`

Same as `append_to_catcif_file` but writes to an already-open file object. Use this
when you want to manage the file handle yourself (e.g., to keep it open across many
appends without repeatedly opening and closing it).

```python
with open("out.catcif", "a") as f:
    for cif_str, tag in my_structures:
        catcif_tools.append_to_catcif_file_open(f, cif_str, tag)
```

---

#### `to_catcif_string(cif_str, tag, scores=None, add_header=False, compress=False) → str | bytes`

Prepare a CIF string for inclusion in a `.catcif` file without writing to disk. Returns
plain text by default; returns gzip bytes when `compress=True`.

Renames the `data_` block header and `_entry.id` to `tag`, and optionally inserts a
`_catcif_scores` block.

If the CIF text has no `data_` header, pass `add_header=True` to insert one
automatically; otherwise a `ValueError` is raised.

---

#### `rename_structure(structure, new_tag) → str`

Return a copy of `structure` with the tag replaced by `new_tag` in both the `data_`
line and the `_entry.id` field.

---

### Scores

#### `get_scores(structure) → dict[str, str]`

Parse and return the `_catcif_scores` block from a CIF structure string. Values are
returned as strings. Returns an empty dict if no scores are present.

```python
structure = catcif_tools.get_structure("my_design", catcif_file="designs.catcif")
scores = catcif_tools.get_scores(structure)
print(scores["total_score"])
```

---

#### `parse_score_file(catcif_file, preserve_tags=False) → Iterator[dict[str, str]]`

Yield one scores dict per structure in `catcif_file`. Each dict contains the
`_catcif_scores` values plus a `'tag'` key as the last entry.

Reads only as much of each structure as needed to find the score block.

```python
for row in catcif_tools.parse_score_file("designs.catcif"):
    print(row["tag"], row.get("total_score"))
```

---

### Tags and paths

#### `get_tags(path) → list[str]`

Return the ordered list of canonical tags in a `.catcif` file.

```python
tags = catcif_tools.get_tags("designs.catcif")
```

---

#### `get_tags_from_index(index) → list[str]`

Return the ordered list of canonical tags from an already-built catcif index dict.

```python
index = catcif_tools.get_catcif_index("designs.catcif")
tags = catcif_tools.get_tags_from_index(index)
```

---

#### `get_catcif_index(catcif_file, no_cache=False, instant_cache=False, return_f=False) → dict | tuple`

Return the index for a catcif file. The index maps each canonical tag to its byte
offset. By default (`return_f=False`) only the index dict is returned and any file
opened for index-building is closed internally. Pass `return_f=True` to get the full
`(index, f_open, caller_must_close)` 3-tuple when you need to seek and read structures
directly.

Normally you do not need this directly — use `get_structure` or `get_all_structures`
instead.

---

#### `split_catcif_tag(tag) → tuple[str | None, str]`

Split a `path/to/file.catcif:tag` string into `(path, tag)`. Returns `(None, tag)` if
the string does not contain `.catcif:` or `.cif:`.

```python
path, tag = catcif_tools.split_catcif_tag("designs.catcif:my_design")
# path = "designs.catcif", tag = "my_design"
```

---

#### `is_catcif_path_tag(tag) → bool`

Return `True` if `tag` is a `path:tag` compound string (contains `.catcif:` or
`.cif:`).

---

### Biology (requires `pip install catcif_tools[pdb]`)

#### `get_sequence(structure) → list[dict]`

Return the 1-letter amino-acid sequence for each chain in a CIF structure.

Each dict in the returned list has keys `chain_id`, `entity_id`, and `sequence`.
Residues outside the canonical 20 amino acids are represented as `X`. Solvent atoms
(`label_seq_id == '.'`) are skipped.

Requires `pip install catcif_tools[pdb]`.

```python
for chain in catcif_tools.get_sequence(structure):
    print(chain["chain_id"], chain["sequence"])
```

---

#### `chain_token_lengths(structure) → list[dict]`

Return the token count for each chain in a CIF structure.

Token counts follow the AlphaFold/RoseTTAFold convention:
- Canonical amino acid residue → 1 token.
- Non-standard residue or ligand → number of non-hydrogen atoms in that residue.

Each dict has keys `chain_id`, `entity_id`, and `tokens`.

Requires `pip install catcif_tools[pdb]`.

```python
for chain in catcif_tools.chain_token_lengths(structure):
    print(chain["chain_id"], chain["tokens"])
```

---

### Cache control

#### `clear_cache()`

Empty the entire in-memory index cache and close any owned file pointers. Call this if
catcif files on disk have been modified since they were last read.

---

### Settings

`catcif_tools.catcif_settings` is a module-level singleton with the following
attributes:

| Attribute | Default | Description |
|---|---|---|
| `cache_indexes` | `True` | Cache loaded indexes in memory |
| `cache_file_pointers` | `True` | Cache open file pointers for random access |
| `fast_cache_num` | `3` | Access count threshold for fast caching |
| `fast_cache_time` | `5000` | Time window (ms) for fast caching |
| `slow_cache_num` | `10` | Access count threshold for slow caching |
| `slow_cache_time` | `100000` | Time window (ms) for slow caching |
| `max_caches` | `10` | Maximum number of files to keep in cache |

```python
import catcif_tools
catcif_tools.catcif_settings.cache_indexes = False   # disable all caching
```


## Installation

```bash
pip install catcif-tools
```

Or install from source:

```bash
git clone https://github.com/bcov77/catcif_tools
cd catcif_tools
pip install -e .
```
