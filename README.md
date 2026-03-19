# peff_digest

A PEFF-aware protein digest tool. Given a PEFF file and enzymatic digestion parameters, produces a CSV of peptides with their sequence, variant annotation, length, and monoisotopic mass.

Each PEFF `VariantSimple` and `VariantComplex` annotation is applied independently (not combined). PEFF PTMs (`ModResPsi` / `ModResUnimod`) are applied combinatorially up to a configurable limit per peptide. Fixed and variable user-defined modifications are also supported.

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

[fixed_mods]
C = "Carbamidomethyl"

[variable_mods]
M = ["Oxidation"]
```

### `DigestConfig` fields

| Field | Type | Default | Description |
|---|---|---|---|
| `peff_file` | `str` | required | Path to the input PEFF file. Must exist. |
| `output_file` | `str` | `"peptides.csv"` | Path for the output CSV. |
| `cleave_on` | `str` | `"KR"` | Amino acids at which to cleave (e.g. `"KR"` for trypsin). |
| `missed_cleavages` | `int` | `2` | Maximum number of missed cleavage sites per peptide. Min 0. |
| `semi_enzymatic` | `bool` | `false` | Include semi-enzymatic peptides (one non-enzymatic terminus). |
| `max_ptm_per_peptide` | `int` | `2` | Maximum number of PEFF PTM annotations applied simultaneously. Set to `0` to skip PEFF PTMs. Min 0. |
| `min_length` | `int` | `7` | Minimum peptide length in residues (inclusive). Min 1. |
| `max_length` | `int` | `40` | Maximum peptide length in residues (inclusive). Min 1. |
| `restrict_after` | `str` | `"P"` | Skip cleavage when the following residue is in this set (e.g. `"P"` for trypsin/Pro rule). |
| `restrict_before` | `str` | `""` | Skip cleavage when the preceding residue is in this set. |
| `cterminal` | `bool` | `true` | `true` = C-terminal cleavage (standard); `false` = N-terminal. |
| `fixed_mods` | `dict[str, str]` | `{"C": "Carbamidomethyl"}` | Static modifications applied to all matching residues. Keys are single-letter amino acids, values are modification names. |
| `variable_mods` | `dict[str, list[str]]` | `{"M": ["Oxidation"]}` | Variable modifications enumerated combinatorially. Keys are single-letter amino acids, values are lists of modification names. |
| `min_mass` | `float \| None` | `None` | Minimum peptide mass in Da (exclusive lower bound). Ignored if `None`. |
| `max_mass` | `float \| None` | `None` | Maximum peptide mass in Da (exclusive upper bound). Ignored if `None`. |
| `drop_invalid_mass` | `bool` | `false` | If `true`, exclude peptides whose mass cannot be computed. |
| `workers` | `int \| None` | `None` | Number of worker processes. Defaults to all available CPUs. Min 1. |

## Python API

### Full digest → Polars DataFrame

```python
from peff_digest import DigestConfig, digest

config = DigestConfig(
    peff_file="human.peff",
    missed_cleavages=2,
    min_length=7,
    max_length=40,
    min_mass=400.0,
    max_mass=10000.0,
    drop_invalid_mass=True,
)

df = digest(config)
print(df)
```

Returns a `polars.DataFrame` with columns `protein_id`, `sequence`, `variant`, `length`, `mass`. All filtering from the config (mass bounds, `drop_invalid_mass`) is applied before returning.

### Single-entry digest

```python
import pefftacular as pf
from peff_digest import digest_peff_sequence

entry = next(iter(pf.PeffReader("human.peff")))

peptides = digest_peff_sequence(
    entry,
    cleave_on="KR",
    missed_cleavages=2,
    min_length=7,
    max_length=40,
    restrict_after="P",
    fixed_mods={"C": "Carbamidomethyl"},
    variable_mods={"M": ["Oxidation"]},
    max_ptm_per_peptide=2,
)

for peptide in peptides:
    print(str(peptide), len(peptide), peptide.mass())
```

Returns a `set[peptacular.ProFormaAnnotation]`. Each element supports `len()`, `.mass()`, `str()`, and `.peptide_name` (PEFF variant notation, or `None` for canonical).

## Development

```bash
just lint      # ruff check
just format    # ruff format + import sort
just test      # pytest
just check     # lint + type check (ty) + test
```
