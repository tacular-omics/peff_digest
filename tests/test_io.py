from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pefftacular as pf
import peptacular as pt
import psimodpy
import pytest

from peff_digest import InternalMod
from peff_digest.io import (
    _digest_batch_worker,
    _digest_worker,
    _format_variant,
    _get_uni_db,
    _try_convert_psimod_to_unimod,
    read_sequences,
)

from conftest import _cfg, _make_entry


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
        assert row[4] is not None
        assert 100.0 <= row[4] <= 2000.0


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
    assert acak_rows[0][2] == "(2|C)"
    # Canonical rows should have None variant
    canonical_rows = [r for r in rows if r[2] is None]
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
    assert all(r[4] is None for r in rows)


# ---------------------------------------------------------------------------
# _digest_batch_worker
# ---------------------------------------------------------------------------


def test_digest_batch_worker_multiple_sequences() -> None:
    """Batch worker should process multiple sequences and return combined rows."""
    entries = [_make_entry("AKB"), _make_entry("CKD")]
    config = _cfg(cleave_on="K")
    rows = _digest_batch_worker(entries, config)
    protein_ids = {r[0] for r in rows}
    assert protein_ids == {"TEST_ID"}  # both have same ID in test helper
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
# _try_convert_psimod_to_unimod
# ---------------------------------------------------------------------------


def test_try_convert_psimod_to_unimod_converts_mod_tag(psi_db: psimodpy.PsiModDatabase) -> None:
    """MOD:00696 (phospho) has a UniMod xref and must be converted to UNIMOD:21."""
    uni_db = _get_uni_db()
    ann = pt.parse("ACD[MOD:00696]EFGR")
    result = _try_convert_psimod_to_unimod(ann, psi_db, uni_db)

    assert result is not None
    assert "UNIMOD:21" in str(result)
    assert "MOD:00696" not in str(result)


def test_try_convert_psimod_to_unimod_returns_none_when_no_unimod_xref(
    psi_db: psimodpy.PsiModDatabase,
) -> None:
    """MOD:00050 has no UniMod xref — conversion must return None."""
    uni_db = _get_uni_db()
    ann = pt.parse("[MOD:00050]-AAK")
    result = _try_convert_psimod_to_unimod(ann, psi_db, uni_db)

    assert result is None


def test_try_convert_psimod_to_unimod_passes_through_non_mod_tags(
    psi_db: psimodpy.PsiModDatabase,
) -> None:
    """Tags that do not start with 'MOD:' must pass through unchanged."""
    uni_db = _get_uni_db()
    ann = pt.parse("AC[Carbamidomethyl]K")
    result = _try_convert_psimod_to_unimod(ann, psi_db, uni_db)

    assert result is not None
    assert "Carbamidomethyl" in str(result)


# ---------------------------------------------------------------------------
# read_sequences with malformed FASTA entry
# ---------------------------------------------------------------------------


def test_read_sequences_fasta_counts_malformed_entries(tmp_path: Path) -> None:
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">sp|P12345|GOOD\nACDEFGR\n>malformed_entry\nACDEFGR\n")

    with patch("peff_digest.io._fasta_to_entry", side_effect=[
        pf.SequenceEntry(prefix="sp", db_unique_id="P12345", sequence="ACDEFGR"),
        ValueError("simulated malformed"),
    ]):
        seqs, n_malformed = read_sequences(str(fasta))

    assert n_malformed == 1
    assert len(seqs) == 1
