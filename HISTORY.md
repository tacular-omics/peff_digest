# History

## Unreleased

* Works with pefftacular 1.0, psimodpy 1.0 and unimodpy 1.0 as well as the 0.x releases.
  Dependency pins are now ranges: `pefftacular>=0.4,<2`, `peptacular>=4.0,<5`,
  `psimodpy>=0.2,<2`, `unimodpy>=0.2,<2`. The unused `uniprotptmpy` dependency is dropped.
* PEFF input: a malformed entry is now skipped and reading continues with the next entry.
  Before, the first malformed entry silently ended the read. Text before the first `>`
  line is skipped with a warning. An invalid file header now raises `PeffParseError`
  instead of giving an empty result.

## 0.1.0 (2026-03-19)

* First release on PyPI.
