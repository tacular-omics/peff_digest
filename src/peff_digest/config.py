from __future__ import annotations

import json
import tomllib
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, Field, field_validator


class InternalMod(BaseModel, frozen=True):
    modification: str
    residue: str  # one or more AAs (e.g. "C" or "KR") — mod applies to each
    mod_type: Literal["fixed", "variable"]


class TerminalMod(BaseModel, frozen=True):
    modification: str
    position: Literal["nterm", "cterm"]
    mod_type: Literal["fixed", "variable"]
    residue: str | None = None  # AA that must be at the terminus; None = any
    protein_terminus: bool = False  # True = only the protein N- or C-terminal peptide


class DigestConfig(BaseModel):
    input_file: str = ""
    output_file: str = "peptides.csv"
    cleave_on: str = "KR"
    missed_cleavages: int = Field(default=2, ge=0)
    semi_enzymatic: bool = False
    max_ptm_per_peptide: int = Field(default=2, ge=0)
    min_length: int = Field(default=7, ge=1)
    max_length: int = Field(default=40, ge=1)
    restrict_after: str = "P"
    restrict_before: str = ""
    cterminal: bool = True
    internal_mods: list[InternalMod] = Field(default_factory=list)
    terminal_mods: list[TerminalMod] = Field(default_factory=list)
    min_mass: float | None = Field(default=None, gt=0)
    max_mass: float | None = Field(default=None, gt=0)
    drop_invalid_mass: bool = False
    annotate_variants: bool = True
    use_mod_names: bool = False
    use_psi_mods: bool = True
    use_unimod_output: bool = False
    include_simple_variants: bool = True
    include_complex_variants: bool = True
    workers: int | None = Field(default=None, ge=1)
    batch_size: int = Field(default=1, ge=1)

    @field_validator("input_file")
    @classmethod
    def input_file_must_exist(cls, v: str) -> str:
        if not v:
            return v
        if not Path(v).exists():
            raise ValueError(f"input_file does not exist: {v}")
        return v

    @classmethod
    def _parse_file(cls, path: Path) -> dict[str, Any]:
        if path.suffix == ".toml":
            with open(path, "rb") as f:
                return tomllib.load(f)
        with open(path) as f:
            return json.load(f)

    @classmethod
    def from_file(cls, path: str | Path, **overrides: Any) -> DigestConfig:
        data = cls._parse_file(Path(path))
        return cls(**{**data, **overrides})

    def to_file(self, path: str | Path) -> None:
        p = Path(path)
        data = self.model_dump()
        if p.suffix == ".toml":
            import tomli_w

            data = {k: v for k, v in data.items() if v is not None}
            with open(p, "wb") as f:
                tomli_w.dump(data, f)
        else:
            with open(p, "w") as f:
                json.dump(data, f, indent=2)
