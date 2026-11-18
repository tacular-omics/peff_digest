# History

## Unreleased

* The variant arithmetic now comes from pefftacular rather than from a private copy here.
  `SequenceEntry.apply_variants` (new in pefftacular 1.4.0, corrected in 1.4.1) builds each
  sequence form, so the floor is raised to `pefftacular>=1.4.1,<2`. `_Variant`,
  `_canonical`, `_apply_simple`, `_apply_complex` and the position-mapping half of
  `_mods_in_span` are gone; cut sites, missed cleavages, semi-enzymatic spans, PTM
  combinatorics, PSI-MOD origin validation, terminal promotion and ProForma output are
  unchanged. The public API, the CLI and every output column stay as they were.
* **Changed, and it is an interpretation rather than a conformance fix — Ty to decide.** A
  `\VariantSimple` whose new residue is `*` now gives a form truncated before the stop,
  because that is what pefftacular 1.4.1's `apply_variants` does and this package now calls
  it. PEFF 1.0 §3.3.8 requires only that "The asterisk is to be interpreted as a nonsense
  mutation (stop codon) over which a peptide sequence MUST NOT span", which 0.3.0 already
  satisfied; §3.2's amino acid table glosses `*` as "Sequence interruption (stop codon,
  unknown linkage)", which reads against truncation. The specification does not settle it.
  Who sees a difference: anyone digesting a file with a `*` variant. Measured on
  `tests/data/UniProtExport_3prot.peff` (a UniProt export with ten such annotations),
  17 peptides are added under the default settings and 38 with semi-enzymatic spans, none
  removed; each addition is the peptide ending at the new C-terminus, carrying the variant.
  No fallback or switch for the old behaviour has been added.
* A modification on a residue a variant substituted is now dropped by the library (spec
  §3.3.10) as well as by this package's PSI-MOD origin check. On the repository's test
  files the two rules agree and no output changes; they can differ only for a modification
  whose PSI-MOD origin is `ANY` or whose accession is not in the PSI-MOD database.
* A variant that cannot be applied — a non-numeric position such as `?`, or one outside the
  sequence — now contributes no form at all, instead of a second copy of the reference
  sequence that the deduplication then absorbed. The output is unchanged.

## 0.3.0 (2026-11-14)

* Requires `peptacular>=5.0,<6` and now lists `tacular>=2.0,<3` directly. peptacular 5
  no longer re-exports `AA_LOOKUP`; the residue masses come from `tacular.AA_LOOKUP`.
  Masses are unchanged.
* The `pydantic` floor is raised to `>=2.9`: pydantic 2.0 cannot be installed on
  Python 3.13, which this package declares as supported, because there is no wheel
  for it and the pydantic-core 2.0.1 sdist fails to build there.
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
