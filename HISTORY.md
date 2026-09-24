# History

## Unreleased

* Requires `peptacular>=5.0,<6` and now lists `tacular>=2.0,<3` directly. peptacular 5
  no longer re-exports `AA_LOOKUP`; the residue masses come from `tacular.AA_LOOKUP`.
  Masses are unchanged.
* README: the single-entry example matches the current API (`DigestConfig`, a generator
  of `Peptide`), opens `PeffReader` as a context manager, and the full-digest example
  uses `input_file=`.

## 0.2.0 (2026-09-23)

* Requires the 1.0 releases of the PEFF and modification libraries:
  `pefftacular>=1.0,<2`, `psimodpy>=1.0,<2`, `unimodpy>=1.0,<2`, and `peptacular>=4.2,<5`.
  The unused `uniprotptmpy` dependency is dropped.
* New `DigestConfig` options: `consolidate_proteins` merges rows that share a peptidoform
  (keeping the row with the fewest missed cleavages, fully enzymatic first), and
  `fixed_mod_overrides_peff` (default on: a fixed mod replaces PEFF-annotated mods at
  the same residue).
* CSV and DataFrame output gain `peptidoform_id`, `neutral_mass`, `missed_cleavages`,
  `semi_enzymatic`, `mass_array` and `n_peptidoforms` columns.
* PEFF input: a malformed entry is now skipped and reading continues with the next entry.
  Before, the first malformed entry silently ended the read. Text before the first `>`
  line is skipped with a warning. An invalid file header now raises `PeffParseError`
  instead of giving an empty result.
* Fixed: with `use_mod_names`, PSI-MOD modifications whose Unimod xref has a `#site`
  suffix were silently dropped. Both UniMod conversion paths now parse xrefs the same way.

## 0.1.0 (2026-03-19)

* First release on PyPI.
