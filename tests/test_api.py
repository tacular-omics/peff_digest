from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from peff_digest import DigestConfig
from peff_digest.api import _consolidate_proteins, _to_dataframe, _write_csv, digest

# ---------------------------------------------------------------------------
# _write_csv
# ---------------------------------------------------------------------------


def _row(
    protein_id: str = "P001",
    sequence: str = "ACDK",
    peptidoform_id: str = "ACDK|",
    variant: str | None = None,
    length: int = 4,
    mass: float | None = 432.17,
    missed_cleavages: int = 0,
    semi_enzymatic: bool = False,
    mass_array: list[float] | None = None,
    n_peptidoforms: int = 1,
) -> tuple:
    return (
        protein_id,
        sequence,
        peptidoform_id,
        variant,
        length,
        mass,
        missed_cleavages,
        semi_enzymatic,
        mass_array,
        n_peptidoforms,
    )


def test_write_csv_creates_file_with_correct_headers_and_rows(tmp_path: Path) -> None:
    results = [
        [_row("P001", "ACDK"), _row("P001", "LMNR", "LMNR|", "(5|W)", 4, 530.24)],
        [_row("P002", "STVWYK", "STVWYK|", None, 6, 782.39)],
    ]
    output_file = tmp_path / "peptides.csv"
    n = _write_csv(results, str(output_file))

    assert output_file.exists()
    assert n == 3
    lines = output_file.read_text().splitlines()
    assert (
        lines[0]
        == "protein_id,sequence,peptidoform_id,variant,length,neutral_mass,missed_cleavages,semi_enzymatic,mass_array,n_peptidoforms"
    )
    assert len(lines) == 4  # header + 3 data rows


def test_write_csv_row_values_match_input(tmp_path: Path) -> None:
    results = [
        [_row("P001", "ACDK", "ACDK|", "(2|W)", 4, 432.17, 1, False, n_peptidoforms=2)],
    ]
    output_file = tmp_path / "peptides.csv"
    _write_csv(results, str(output_file))

    lines = output_file.read_text().splitlines()
    assert lines[1] == "P001,ACDK,ACDK|,(2|W),4,432.17,1,False,,2"  # mass_array=None → empty string


# ---------------------------------------------------------------------------
# _to_dataframe
# ---------------------------------------------------------------------------


def test_to_dataframe_returns_correct_schema() -> None:
    results = [[_row("P001", "ACDK")]]
    df = _to_dataframe(results)

    assert isinstance(df, pl.DataFrame)
    assert df.schema["protein_id"] == pl.String
    assert df.schema["sequence"] == pl.String
    assert df.schema["peptidoform_id"] == pl.String
    assert df.schema["variant"] == pl.String
    assert df.schema["length"] == pl.Int64
    assert df.schema["neutral_mass"] == pl.Float64
    assert df.schema["missed_cleavages"] == pl.Int32
    assert df.schema["semi_enzymatic"] == pl.Boolean
    assert df.schema["n_peptidoforms"] == pl.Int32


def test_to_dataframe_row_values_match_input() -> None:
    results = [
        [_row("P001", "ACDK", "ACDK|", "(2|W)", 4, 432.17), _row("P002", "LMR", "LMR|", None, 3, 390.20)],
    ]
    df = _to_dataframe(results)

    assert df.shape == (2, 10)
    assert df["protein_id"].to_list() == ["P001", "P002"]
    assert df["variant"].to_list() == ["(2|W)", None]


def test_to_dataframe_empty_results_returns_empty_frame() -> None:
    df = _to_dataframe([])
    assert df.shape == (0, 10)


# ---------------------------------------------------------------------------
# digest()
# ---------------------------------------------------------------------------


def test_digest_returns_dataframe_by_default(tiny_fasta: Path) -> None:
    cfg = DigestConfig(
        input_file=str(tiny_fasta), min_length=1, max_length=40, missed_cleavages=0, max_ptm_per_peptide=0
    )
    result = digest(cfg, write_file=False)

    assert isinstance(result, pl.DataFrame)
    assert set(result.columns) == {
        "protein_id",
        "sequence",
        "peptidoform_id",
        "variant",
        "length",
        "neutral_mass",
        "missed_cleavages",
        "semi_enzymatic",
        "mass_array",
        "n_peptidoforms",
    }


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
    assert (
        header
        == "protein_id,sequence,peptidoform_id,variant,length,neutral_mass,missed_cleavages,semi_enzymatic,mass_array,n_peptidoforms"
    )


# ---------------------------------------------------------------------------
# _consolidate_proteins
# ---------------------------------------------------------------------------


def test_consolidate_merges_shared_peptidoform() -> None:
    results = [[_row("P1", "ACDK", "ACDK|")], [_row("P2", "ACDK", "ACDK|")]]
    out = _consolidate_proteins(results)
    assert len(out[0]) == 1
    assert out[0][0][0] == "P1;P2"


def test_consolidate_best_row_by_mc() -> None:
    results = [[_row("P1", "ACDK", "ACDK|", missed_cleavages=2), _row("P2", "ACDK", "ACDK|", missed_cleavages=0)]]
    out = _consolidate_proteins(results)
    assert out[0][0][6] == 0


def test_consolidate_same_sequence_merged_different_kept() -> None:
    # Same ProForma sequence from two proteins → merged into one row, n_peptidoforms=1
    results = [
        [_row("P1", "ACDK", "ACDK|")],
        [_row("P2", "ACDK", "ACDK|")],
    ]
    out = _consolidate_proteins(results)
    assert len(out[0]) == 1
    assert out[0][0][0] == "P1;P2"
    assert out[0][0][9] == 1

    # Different ProForma sequences (same peptidoform_id) → kept as separate rows
    results2 = [
        [_row("P1", "AC[Phospho]DK", "ACDK|")],
        [_row("P2", "ACD[Phospho]K", "ACDK|")],
    ]
    out2 = _consolidate_proteins(results2)
    assert len(out2[0]) == 2  # two distinct sequences kept
    assert out2[0][0][9] == 2  # each row knows there are 2 sequences for "ACDK|"
    assert out2[0][1][9] == 2


def test_consolidate_no_duplicate_protein_ids() -> None:
    results = [[_row("P1", "ACDK", "ACDK|")], [_row("P1", "ACDK", "ACDK|")]]
    out = _consolidate_proteins(results)
    assert out[0][0][0] == "P1"


def test_consolidate_returns_nested_list() -> None:
    results = [[_row("P1", "ACDK", "ACDK|")]]
    out = _consolidate_proteins(results)
    assert isinstance(out, list)
    assert isinstance(out[0], list)


def test_consolidate_distinct_peptides_not_merged() -> None:
    results = [[_row("P1", "ACDK", "ACDK|"), _row("P1", "LMNR", "LMNR|")]]
    out = _consolidate_proteins(results)
    assert len(out[0]) == 2


def test_digest_consolidate_proteins_merges_across_proteins(tmp_path: Path) -> None:
    fasta = tmp_path / "shared.fasta"
    fasta.write_text(">sp|P1|PROT1\nACDEFGRACDEFGR\n>sp|P2|PROT2\nACDEFGRKKK\n")
    cfg = DigestConfig(
        input_file=str(fasta),
        min_length=1,
        max_length=40,
        missed_cleavages=0,
        max_ptm_per_peptide=0,
        consolidate_proteins=True,
    )
    df = digest(cfg)
    acdefgr = df.filter(pl.col("sequence") == "ACDEFGR")
    assert acdefgr.shape[0] == 1
    assert ";" in acdefgr["protein_id"][0]
