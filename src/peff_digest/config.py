from __future__ import annotations

from pathlib import Path
from typing import Literal

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
    input_file: str
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
    workers: int | None = Field(default=None, ge=1)

    @field_validator("input_file")
    @classmethod
    def input_file_must_exist(cls, v: str) -> str:
        if not Path(v).exists():
            raise ValueError(f"input_file does not exist: {v}")
        return v
