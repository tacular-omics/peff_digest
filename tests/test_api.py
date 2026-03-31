from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from peff_digest import DigestConfig
from peff_digest.api import _to_dataframe, _write_csv, digest


# ---------------------------------------------------------------------------
# _write_csv
# ---------------------------------------------------------------------------


def test_write_csv_creates_file_with_correct_headers_and_rows(tmp_path: Path) -> None:
    results: list[list[tuple[str, str, str | None, int, float | None]]] = [
        [("P001", "ACDK", None, 4, 432.17), ("P001", "LMNR", "(5|W)", 4, 530.24)],
        [("P002", "STVWYK", None, 6, 782.39)],
    ]
    output_file = tmp_path / "peptides.csv"
    n = _write_csv(results, str(output_file))

    assert output_file.exists()
    assert n == 3
    lines = output_file.read_text().splitlines()
    assert lines[0] == "protein_id,sequence,variant,length,mass"
    assert len(lines) == 4  # header + 3 data rows


def test_write_csv_row_values_match_input(tmp_path: Path) -> None:
    results: list[list[tuple[str, str, str | None, int, float | None]]] = [
        [("P001", "ACDK", "(2|W)", 4, 432.17)],
    ]
    output_file = tmp_path / "peptides.csv"
    _write_csv(results, str(output_file))

    lines = output_file.read_text().splitlines()
    assert lines[1] == "P001,ACDK,(2|W),4,432.17"


# ---------------------------------------------------------------------------
# _to_dataframe
# ---------------------------------------------------------------------------


def test_to_dataframe_returns_correct_schema() -> None:
    results: list[list[tuple[str, str, str | None, int, float | None]]] = [
        [("P001", "ACDK", None, 4, 432.17)],
    ]
    df = _to_dataframe(results)

    assert isinstance(df, pl.DataFrame)
    assert df.schema["protein_id"] == pl.String
    assert df.schema["sequence"] == pl.String
    assert df.schema["variant"] == pl.String
    assert df.schema["length"] == pl.Int64
    assert df.schema["mass"] == pl.Float64


def test_to_dataframe_row_values_match_input() -> None:
    results: list[list[tuple[str, str, str | None, int, float | None]]] = [
        [("P001", "ACDK", "(2|W)", 4, 432.17), ("P002", "LMR", None, 3, 390.20)],
    ]
    df = _to_dataframe(results)

    assert df.shape == (2, 5)
    assert df["protein_id"].to_list() == ["P001", "P002"]
    assert df["variant"].to_list() == ["(2|W)", None]


def test_to_dataframe_empty_results_returns_empty_frame() -> None:
    df = _to_dataframe([])
    assert df.shape == (0, 5)


# ---------------------------------------------------------------------------
# digest()
# ---------------------------------------------------------------------------


def test_digest_returns_dataframe_by_default(tiny_fasta: Path) -> None:
    cfg = DigestConfig(
        input_file=str(tiny_fasta), min_length=1, max_length=40, missed_cleavages=0, max_ptm_per_peptide=0
    )
    result = digest(cfg, write_file=False)

    assert isinstance(result, pl.DataFrame)
    assert set(result.columns) == {"protein_id", "sequence", "variant", "length", "mass"}


def test_digest_write_file_creates_csv_and_returns_none(tiny_fasta: Path, tmp_path: Path) -> None:
    output = tmp_path / "out.csv"
    cfg = DigestConfig(
        input_file=str(tiny_fasta),
        output_file=str(output),
        min_length=1,
        max_length=40,
        missed_cleavages=0,
        max_ptm_per_peptide=0,
    )
    result = digest(cfg, write_file=True)

    assert result is None
    assert output.exists()
    header = output.read_text().splitlines()[0]
    assert header == "protein_id,sequence,variant,length,mass"
