from __future__ import annotations

from pathlib import Path

from pydantic import BaseModel, Field, field_validator


class DigestConfig(BaseModel):
    peff_file: str
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
    fixed_mods: dict[str, str] = {"C": "Carbamidomethyl"}
    variable_mods: dict[str, list[str]] = {"M": ["Oxidation"]}
    min_mass: float | None = Field(default=None, gt=0)
    max_mass: float | None = Field(default=None, gt=0)
    drop_invalid_mass: bool = False
    workers: int | None = Field(default=None, ge=1)

    @field_validator("peff_file")
    @classmethod
    def peff_file_must_exist(cls, v: str) -> str:
        if not Path(v).exists():
            raise ValueError(f"peff_file does not exist: {v}")
        return v
