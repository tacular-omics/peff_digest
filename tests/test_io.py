from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pefftacular as pf
import peptacular as pt
import psimodpy
import pytest
from conftest import _cfg, _make_entry

from peff_digest import InternalMod
from peff_digest.io import (
    _dedup_rows,
    _digest_batch_worker,
    _digest_worker,
    _format_variant,
    _get_uni_db,
    read_sequences,
)

# ---------------------------------------------------------------------------
# _format_variant
# ---------------------------------------------------------------------------


def test_format_variant_simple() -> None:
    assert _format_variant(pf.VariantSimple(position=5, new_amino_acid="W")) == "(5|W)"


def test_format_variant_complex() -> None:
    v = pf.VariantComplex(start_pos=3, end_pos=7, new_sequence="DDK")
    assert _format_variant(v) == "(3|7|DDK)"


def test_format_variant_deletion() -> None:
    v = pf.VariantComplex(start_pos=2, end_pos=5, new_sequence="")
    assert _format_variant(v) == "(2|5|)"


def test_format_variant_none() -> None:
    assert _format_variant(None) is None


# ---------------------------------------------------------------------------
# _digest_worker mass filtering
# ---------------------------------------------------------------------------


def test_digest_worker_min_mass_filters() -> None:
    """Peptides below min_mass should be excluded."""
    entry = _make_entry("GR")  # Gly-Arg: very small peptide (~231 Da)
    config = _cfg(cleave_on="R", min_mass=5000.0)
    rows = _digest_worker(entry, config)
    assert len(rows) == 0


def test_digest_worker_max_mass_filters() -> None:
    """Peptides above max_mass should be excluded."""
    entry = _make_entry("ACDEFGHIKLMNPQRSTVWY")  # long, heavy peptide
    config = _cfg(cleave_on="X", max_mass=100.0)  # X won't cleave, one huge peptide
    rows = _digest_worker(entry, config)
    assert len(rows) == 0


def test_digest_worker_mass_within_range_kept() -> None:
    """Peptides within mass range should be kept."""
    entry = _make_entry("ACDEFGR")
    config = _cfg(cleave_on="R", min_mass=100.0, max_mass=2000.0)
    rows = _digest_worker(entry, config)
    assert len(rows) > 0
    for row in rows:
        assert row[5] is not None
        assert 100.0 <= row[5] <= 2000.0


def test_digest_worker_variant_column_from_peptide_variant() -> None:
    """The variant column should come from Peptide.variant, not peptide_name."""
    entry = _make_entry(
        "AAAKBBBR",
        variant_simple=(pf.VariantSimple(position=2, new_amino_acid="C"),),
    )
    config = _cfg(cleave_on="KR")
    rows = _digest_worker(entry, config)
    # Find the variant row for "ACAK"
    acak_rows = [r for r in rows if "ACAK" in r[1]]
    assert len(acak_rows) == 1
    assert acak_rows[0][3] == "(2|C)"
    # Canonical rows should have None variant
    canonical_rows = [r for r in rows if r[3] is None]
    assert len(canonical_rows) >= 2  # AAAK and BBBR from canonical


def test_digest_worker_min_mass_filters_light_peptides() -> None:
    """Peptides below min_mass must be excluded."""
    # AAK ~288 Da, R ~174 Da — set min_mass=200 to exclude R
    entry = pf.SequenceEntry(prefix="sp", db_unique_id="TEST", sequence="AAKR")
    cfg = _cfg(min_mass=200.0)

    rows = _digest_worker(entry, cfg)
    sequences = [r[1] for r in rows]

    assert "AAK" in sequences
    assert "R" not in sequences


def test_digest_worker_max_mass_filters_heavy_peptides() -> None:
    """Peptides above max_mass must be excluded."""
    # AAK ~288 Da, R ~174 Da — set max_mass=200 to exclude AAK
    entry = pf.SequenceEntry(prefix="sp", db_unique_id="TEST", sequence="AAKR")
    cfg = _cfg(max_mass=200.0)

    rows = _digest_worker(entry, cfg)
    sequences = [r[1] for r in rows]

    assert "R" in sequences
    assert "AAK" not in sequences


def test_digest_worker_use_unimod_output_converts_mod_tags() -> None:
    """With use_unimod_output=True, MOD: tags must be replaced by UNIMOD: tags."""
    entry = pf.SequenceEntry(
        prefix="sp",
        db_unique_id="TEST",
        sequence="ACDEFGR",
        mod_res_psi=(pf.ModResPsi(positions=(3,), accession="MOD:00696", name="phosphorylated residue"),),
    )
    cfg = _cfg(use_unimod_output=True, use_psi_mods=True, max_ptm_per_peptide=1)

    rows = _digest_worker(entry, cfg)
    sequences = [r[1] for r in rows]

    assert any("UNIMOD:21" in s for s in sequences)
    assert not any("MOD:00696" in s for s in sequences)


def test_digest_worker_use_unimod_output_drops_unconvertible_mods() -> None:
    """Peptides carrying a PSI-MOD with no UniMod xref must be dropped."""
    entry = pf.SequenceEntry(
        prefix="sp",
        db_unique_id="TEST",
        sequence="AAKR",
        mod_res_psi=(pf.ModResPsi(positions=(1,), accession="MOD:00050", name="N-acetyl-L-alanine"),),
    )
    cfg = _cfg(use_unimod_output=True, use_psi_mods=True, max_ptm_per_peptide=1)

    rows = _digest_worker(entry, cfg)
    # Modified peptides (with MOD:00050 that has no UniMod xref) must be absent
    assert not any("MOD:00050" in r[1] for r in rows)
    assert not any("UNIMOD:" in r[1] for r in rows)


def test_digest_worker_drop_invalid_mass_excludes_none_mass_rows() -> None:
    """Peptides whose mass cannot be computed must be dropped when drop_invalid_mass=True."""
    entry = pf.SequenceEntry(prefix="sp", db_unique_id="TEST", sequence="AKR")
    cfg = _cfg(drop_invalid_mass=True)

    with patch.object(pt.ProFormaAnnotation, "mass", side_effect=ValueError("no mass")):
        rows = _digest_worker(entry, cfg)

    assert rows == []


def test_digest_worker_keep_invalid_mass_rows_when_flag_off() -> None:
    """Peptides with None mass must be kept when drop_invalid_mass=False."""
    entry = pf.SequenceEntry(prefix="sp", db_unique_id="TEST", sequence="AKR")
    cfg = _cfg(drop_invalid_mass=False)

    with patch.object(pt.ProFormaAnnotation, "mass", side_effect=ValueError("no mass")):
        rows = _digest_worker(entry, cfg)

    assert len(rows) == 2
    assert all(r[5] is None for r in rows)


# ---------------------------------------------------------------------------
# _digest_batch_worker
# ---------------------------------------------------------------------------


def test_digest_batch_worker_multiple_sequences() -> None:
    """Batch worker should process multiple sequences and return combined rows."""
    entries = [_make_entry("AKB"), _make_entry("CKD")]
    config = _cfg(cleave_on="K")
    rows = _digest_batch_worker(entries, config)
    protein_ids = {r[0] for r in rows}
    assert protein_ids == {"sp|TEST_ID"}  # both have same ID in test helper
    assert len(rows) >= 4  # at least "AK", "B", "CK", "D"


def test_digest_batch_worker_empty_batch() -> None:
    """Empty batch should return no rows."""
    rows = _digest_batch_worker([], _cfg(cleave_on="K"))
    assert rows == []


# ---------------------------------------------------------------------------
# _get_uni_db caching
# ---------------------------------------------------------------------------


def test_get_uni_db_returns_database_and_caches() -> None:
    import peff_digest.io as io_mod

    # Reset cache to force a fresh load
    io_mod._UNI_DB = None
    db1 = _get_uni_db()
    db2 = _get_uni_db()

    assert db1 is not None
    assert db1 is db2  # same object — cached


# ---------------------------------------------------------------------------
# _dedup_rows
# ---------------------------------------------------------------------------


def _r(pid: str, pform_id: str, mc: int, semi: bool) -> tuple:
    # Use pform_id as part of the sequence so rows with distinct pform_ids are distinct sequences
    return (pid, f"SEQ_{pform_id}", pform_id, None, 3, 300.0, mc, semi)


def test_dedup_rows_no_duplicates_returns_all() -> None:
    rows = [_r("P1", "A|", 0, False), _r("P1", "B|", 0, False)]
    result = _dedup_rows(rows)
    assert len(result) == 2
    assert {r[2] for r in result} == {"A|", "B|"}
    assert all(r[8] == 1 for r in result)


def test_dedup_rows_lower_mc_wins() -> None:
    # Same ProForma sequence, two rows → 1 unique sequence
    rows = [_r("P1", "A|", 2, False), _r("P1", "A|", 1, False)]
    result = _dedup_rows(rows)
    assert len(result) == 1
    assert result[0][6] == 1
    assert result[0][8] == 1  # both rows had the same sequence "SEQ"


def test_dedup_rows_keeps_distinct_sequences() -> None:
    # Different ProForma sequences → both rows are kept, each gets n_peptidoforms=2
    row1 = ("P1", "SEQ[Phospho]", "A|", None, 3, 300.0, 1, False)
    row2 = ("P1", "S[Phospho]EQ", "A|", None, 3, 300.0, 0, False)
    result = _dedup_rows([row1, row2])
    assert len(result) == 2
    assert result[0][8] == 2
    assert result[1][8] == 2


def test_dedup_rows_same_mc_prefers_non_semi() -> None:
    rows = [_r("P1", "A|", 1, True), _r("P1", "A|", 1, False)]
    result = _dedup_rows(rows)
    assert len(result) == 1
    assert result[0][7] is False
    assert result[0][8] == 1  # same sequence "SEQ" in both rows


def test_dedup_rows_preserves_first_appearance_order() -> None:
    # Third row is an exact duplicate of first (same sequence) and should be dropped
    row_a1 = _r("P1", "A|", 0, False)
    row_b = _r("P1", "B|", 0, False)
    row_a2 = ("P1", row_a1[1], "A|", None, 3, 300.0, 1, False)  # same seq as row_a1, worse MC
    result = _dedup_rows([row_a1, row_b, row_a2])
    assert [r[2] for r in result] == ["A|", "B|"]


def test_dedup_rows_empty_returns_empty() -> None:
    assert _dedup_rows([]) == []


# ---------------------------------------------------------------------------
# read_sequences with malformed FASTA entry
# ---------------------------------------------------------------------------


def test_read_sequences_fasta_counts_malformed_entries(tmp_path: Path) -> None:
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">sp|P12345|GOOD\nACDEFGR\n>malformed_entry\nACDEFGR\n")

    with patch(
        "peff_digest.io._fasta_to_entry",
        side_effect=[
            pf.SequenceEntry(prefix="sp", db_unique_id="P12345", sequence="ACDEFGR"),
            ValueError("simulated malformed"),
        ],
    ):
        seqs, n_malformed = read_sequences(str(fasta))

    assert n_malformed == 1
    assert len(seqs) == 1


# ---------------------------------------------------------------------------
# read_sequences with malformed PEFF entries
# ---------------------------------------------------------------------------

_PEFF_HEADER = "# PEFF 1.0\n# //\n# DbName=testdb\n# Prefix=sp\n# DbVersion=1\n# DbSource=http://example.com\n# NumberOfEntries=4\n# SequenceType=AA\n# //\n"


def _ids(seqs: list[pf.SequenceEntry]) -> list[str]:
    return [s.db_unique_id for s in seqs]


def test_read_sequences_peff_skips_malformed_entry_and_keeps_reading(tmp_path: Path) -> None:
    peff = tmp_path / "test.peff"
    peff.write_text(
        _PEFF_HEADER + ">sp:P1 \\Length=5\nACDEF\n"
        ">sp:P2 \\Length=notanumber\nACDEF\n"  # bad integer: malformed on every pefftacular version
        ">sp:P3 \\Length=5\nGHIKL\n"
        ">sp:P4 \\Length=5\nMNPQR\n"
    )
    seqs, n_malformed = read_sequences(str(peff))
    assert _ids(seqs) == ["P1", "P3", "P4"]
    assert n_malformed == 1


def test_read_sequences_peff_malformed_last_entry(tmp_path: Path) -> None:
    peff = tmp_path / "test.peff"
    peff.write_text(_PEFF_HEADER + ">sp:P1 \\Length=5\nACDEF\n>sp:P2 \\Length=x\nACDEF\n")
    seqs, n_malformed = read_sequences(str(peff))
    assert _ids(seqs) == ["P1"]
    assert n_malformed == 1


def test_read_sequences_peff_consecutive_malformed_entries(tmp_path: Path) -> None:
    peff = tmp_path / "test.peff"
    peff.write_text(_PEFF_HEADER + ">sp:P1 \\Length=x\nACDEF\n>sp:P2 \\Length=y\nACDEF\n>sp:P3 \\Length=5\nGHIKL\n")
    seqs, n_malformed = read_sequences(str(peff))
    assert _ids(seqs) == ["P3"]
    assert n_malformed == 2


def test_read_sequences_peff_empty_sequence_entry(tmp_path: Path) -> None:
    # pefftacular >= 1.0 rejects an entry with no sequence; older versions accept it.
    peff = tmp_path / "test.peff"
    peff.write_text(_PEFF_HEADER + ">sp:P1\nACDEF\n>sp:P2\n>sp:P3\nGHIKL\n")
    seqs, n_malformed = read_sequences(str(peff))
    ids = _ids(seqs)
    assert ids[0] == "P1" and ids[-1] == "P3"
    assert len(ids) + n_malformed == 3


def test_read_sequences_peff_text_before_first_entry_is_skipped(tmp_path: Path) -> None:
    peff = tmp_path / "test.peff"
    peff.write_text(_PEFF_HEADER + "STRAYSEQ\n>sp:P1\nACDEF\n")
    seqs, n_malformed = read_sequences(str(peff))
    assert _ids(seqs) == ["P1"]
    assert n_malformed == 1


def test_read_sequences_peff_invalid_header_raises(tmp_path: Path) -> None:
    peff = tmp_path / "test.peff"
    peff.write_text("# PEFF\n>sp:P1\nACDEF\n")
    with pytest.raises(pf.PeffParseError):
        read_sequences(str(peff))
