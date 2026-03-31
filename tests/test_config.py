from __future__ import annotations

from pathlib import Path

import pytest
from pydantic import ValidationError

from peff_digest import DigestConfig, InternalMod

from conftest import _cfg


def test_config_round_trip_toml(tmp_path: Path) -> None:
    """DigestConfig should round-trip through TOML serialization."""
    config = DigestConfig(
        input_file="", cleave_on="KR", missed_cleavages=2,
        min_length=5, max_length=30, max_ptm_per_peptide=3,
    )
    path = tmp_path / "test.toml"
    config.to_file(path)
    loaded = DigestConfig.from_file(path)
    assert loaded.cleave_on == "KR"
    assert loaded.missed_cleavages == 2
    assert loaded.min_length == 5
    assert loaded.max_length == 30
    assert loaded.max_ptm_per_peptide == 3


def test_config_round_trip_json(tmp_path: Path) -> None:
    """DigestConfig should round-trip through JSON serialization."""
    config = DigestConfig(
        input_file="", cleave_on="R", missed_cleavages=1,
        internal_mods=[InternalMod(modification="Oxidation", residue="M", mod_type="variable")],
    )
    path = tmp_path / "test.json"
    config.to_file(path)
    loaded = DigestConfig.from_file(path)
    assert loaded.cleave_on == "R"
    assert loaded.missed_cleavages == 1
    assert len(loaded.internal_mods) == 1
    assert loaded.internal_mods[0].modification == "Oxidation"


def test_config_from_file_with_overrides(tmp_path: Path) -> None:
    """CLI overrides should take precedence over file values."""
    config = DigestConfig(input_file="", cleave_on="KR", missed_cleavages=2)
    path = tmp_path / "test.json"
    config.to_file(path)
    loaded = DigestConfig.from_file(path, missed_cleavages=5, min_length=10)
    assert loaded.missed_cleavages == 5
    assert loaded.min_length == 10
    assert loaded.cleave_on == "KR"  # from file


def test_digest_config_nonexistent_input_file_raises() -> None:
    with pytest.raises(ValidationError, match="input_file does not exist"):
        _cfg(input_file="/nonexistent/path/to/file.fasta")
