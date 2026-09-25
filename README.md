# peff_digest

A PEFF-aware protein digest tool. Given a PEFF file and enzymatic digestion parameters, produces a CSV of peptides with their sequence (ProForma notation), variant annotation, length, and monoisotopic mass.

Each PEFF `VariantSimple` and `VariantComplex` annotation is applied independently (not combined). PEFF PTMs (`ModResPsi` / `ModResUnimod`) are applied combinatorially up to a configurable limit per peptide. Fixed and variable user-defined modifications — including terminal modifications — are also supported.

## Installation

```bash
uv sync
# or
just install
```

Requires Python 3.12+.

## Usage

### Via config file

```bash
peff-digest --config config.toml
```

### Via flags

```bash
peff-digest --peff-file human.peff --output-file peptides.csv
```

Flags override any values set in `--config`. Run `peff-digest --help` for the full flag reference.

### Output

The output CSV has five columns:

| Column | Description |
|---|---|
| `protein_id` | `db_unique_id` from the PEFF header |
| `sequence` | ProForma-annotated peptide sequence (includes mods) |
| `variant` | PEFF variant notation, e.g. `(42\|R)`, or empty for canonical |
| `length` | Peptide length in residues |
| `mass` | Monoisotopic mass in Da, or empty if not computable |

## Config reference

All options can be set in a TOML or JSON config file. TOML example:

```toml
peff_file = "human.peff"
output_file = "peptides.csv"

cleave_on = "KR"
missed_cleavages = 2
semi_enzymatic = false
max_ptm_per_peptide = 2
min_length = 7
max_length = 40
restrict_after = "P"
restrict_before = ""
cterminal = true
min_mass = 400.0
max_mass = 10000.0
drop_invalid_mass = false
annotate_variants = true

# Internal modifications — one [[internal_mods]] block per mod:
# [[internal_mods]]
# modification = "Carbamidomethyl"
# residue = "C"
# mod_type = "fixed"
#
# [[internal_mods]]
# modification = "Oxidation"
# residue = "M"
# mod_type = "variable"

# Terminal modifications — one [[terminal_mods]] block per mod:
# [[terminal_mods]]
# modification = "Acetyl"
# position = "nterm"
# mod_type = "variable"
# protein_terminus = true   # only the first peptide of each protein
#
# [[terminal_mods]]
# modification = "UNIMOD:737"
# position = "nterm"
# mod_type = "fixed"
# residue = "M"             # only if the terminal residue is M
#
# [[terminal_mods]]
# modification = "Amidated"
# position = "cterm"
# mod_type = "variable"
```

### `DigestConfig` fields

| Field | Type | Default | Description |
|---|---|---|---|
| `peff_file` | `str` | required | Path to the input PEFF file. Must exist. |
| `output_file` | `str` | `"peptides.csv"` | Path for the output CSV. |
| `cleave_on` | `str` | `"KR"` | Amino acids at which to cleave (e.g. `"KR"` for trypsin). |
| `missed_cleavages` | `int` | `2` | Maximum number of missed cleavage sites per peptide. Min 0. |
| `semi_enzymatic` | `bool` | `false` | Include semi-enzymatic peptides (one non-enzymatic terminus). |
| `max_ptm_per_peptide` | `int` | `2` | Maximum number of variable mods (PEFF + user) applied simultaneously per peptide. `0` disables all variable mods. Min 0. |
| `min_length` | `int` | `7` | Minimum peptide length in residues (inclusive). Min 1. |
| `max_length` | `int` | `40` | Maximum peptide length in residues (inclusive). Min 1. |
| `restrict_after` | `str` | `"P"` | Skip cleavage when the following residue is in this set (e.g. `"P"` for trypsin/Pro rule). |
| `restrict_before` | `str` | `""` | Skip cleavage when the preceding residue is in this set. |
| `cterminal` | `bool` | `true` | `true` = C-terminal cleavage (standard); `false` = N-terminal. |
| `internal_mods` | `list[InternalMod]` | `[]` | Per-residue modifications. See `InternalMod` fields below. |
| `terminal_mods` | `list[TerminalMod]` | `[]` | Terminal modifications. See `TerminalMod` fields below. |
| `min_mass` | `float \| None` | `None` | Minimum peptide mass in Da. Ignored if `None`. |
| `max_mass` | `float \| None` | `None` | Maximum peptide mass in Da. Ignored if `None`. |
| `drop_invalid_mass` | `bool` | `false` | If `true`, exclude peptides whose mass cannot be computed. |
| `annotate_variants` | `bool` | `true` | If `false`, do not set `peptide_name` on variant peptides. |
| `workers` | `int \| None` | `None` | Number of worker processes. Defaults to all available CPUs. Min 1. |

### `InternalMod` fields

| Field | Type | Default | Description |
|---|---|---|---|
| `modification` | `str` | required | Modification name (e.g. `"Carbamidomethyl"`, `"UNIMOD:21"`). |
| `residue` | `str` | required | One or more amino acids the mod applies to (e.g. `"C"` or `"KR"`). |
| `mod_type` | `"fixed" \| "variable"` | required | `"fixed"` = always applied; `"variable"` = enumerated combinatorially (counts against `max_ptm_per_peptide`). |

### `TerminalMod` fields

| Field | Type | Default | Description |
|---|---|---|---|
| `modification` | `str` | required | Modification name (e.g. `"Acetyl"`, `"UNIMOD:737"`). |
| `position` | `"nterm" \| "cterm"` | required | Which terminus to apply the mod to. |
| `mod_type` | `"fixed" \| "variable"` | required | `"fixed"` = always applied; `"variable"` = enumerated combinatorially (counts against `max_ptm_per_peptide`). |
| `residue` | `str \| None` | `None` | If set, the mod is only applied when the terminal residue is in this string (e.g. `"M"` or `"KR"`). |
| `protein_terminus` | `bool` | `false` | If `true`, only apply to the protein-level terminus (first peptide for N-term, last for C-term). |

## Python API

### Full digest → Polars DataFrame

```python
from peff_digest import DigestConfig, InternalMod, TerminalMod, digest

config = DigestConfig(
    input_file="human.peff",
    missed_cleavages=2,
    min_length=7,
    max_length=40,
    min_mass=400.0,
    max_mass=10000.0,
    drop_invalid_mass=True,
    internal_mods=[
        InternalMod(modification="Carbamidomethyl", residue="C", mod_type="fixed"),
        InternalMod(modification="Oxidation", residue="M", mod_type="variable"),
    ],
    terminal_mods=[
        TerminalMod(modification="Acetyl", position="nterm", mod_type="variable", protein_terminus=True),
    ],
)

df = digest(config)
print(df)
```

Returns a `polars.DataFrame` with columns `protein_id`, `sequence`, `variant`, `length`, `mass`. All filtering from the config (mass bounds, `drop_invalid_mass`) is applied before returning.

### Single-entry digest

```python
import pefftacular as pf
from peff_digest import DigestConfig, InternalMod, TerminalMod, digest_peff_sequence

with pf.PeffReader("human.peff") as reader:
    entry = next(iter(reader))

config = DigestConfig(
    cleave_on="KR",
    missed_cleavages=2,
    min_length=7,
    max_length=40,
    restrict_after="P",
    internal_mods=[
        InternalMod(modification="Carbamidomethyl", residue="C", mod_type="fixed"),
        InternalMod(modification="Oxidation", residue="M", mod_type="variable"),
    ],
    max_ptm_per_peptide=2,
    terminal_mods=[
        TerminalMod(modification="Acetyl", position="nterm", mod_type="variable", protein_terminus=True),
        TerminalMod(modification="Amidated", position="cterm", mod_type="variable"),
    ],
)

for peptide in digest_peff_sequence(entry, config):
    print(peptide.proforma.serialize(), peptide.sequence, peptide.mass)
```

Yields `Peptide` objects (a generator, each peptidoform once). Each has `.proforma` (a
`peptacular.ProFormaAnnotation`), `.sequence`, `.mass`, `.mod_map`, `.missed_cleavages`,
`.semi_enzymatic`, `.is_protein_nterm`/`.is_protein_cterm` and `.variant` (the applied PEFF
variant, or `None` for canonical). `ProFormaAnnotation` is not hashable in peptacular 5, so
key on `peptide.proforma.serialize()` to put peptides in a set or dict.

## Development

```bash
just lint      # ruff check
just format    # ruff format + import sort
just test      # pytest
just check     # lint + type check (ty) + test
```
