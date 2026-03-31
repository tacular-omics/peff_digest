from __future__ import annotations

from peff_digest import digest_peff_sequence

from conftest import _cfg, _make_entry


def test_semi_enzymatic_right_open() -> None:
    """Enzymatic N-term with non-enzymatic C-term should be yielded."""
    # "AAKBBB": cut sites [0, 3, 6]; enzymatic: "AAK", "BBB"
    # right-open from pos 0: ends 1, 2 (no internal cuts crossed) → "A" (filtered), "AA"
    # right-open from pos 3: ends 4, 5 → "B" (filtered), "BB"
    entry = _make_entry("AAKBBB")
    result = list(digest_peff_sequence(
        entry, _cfg(cleave_on="K", min_length=2, semi_enzymatic=True),
    ))
    seqs = {str(p.proforma) for p in result}
    # Right-open: starts at enzymatic cut (pos 0 or pos 3), ends at non-cut position
    assert "AA" in seqs   # right-open: pos 0 → pos 2, no internal cuts
    assert "BB" in seqs   # right-open: pos 3 → pos 5, no internal cuts
    semi_peptides = [p for p in result if p.semi_enzymatic]
    assert len(semi_peptides) > 0
    assert any(p.sequence == "AA" and p.semi_enzymatic for p in semi_peptides)


def test_semi_enzymatic_left_open() -> None:
    """Non-enzymatic N-term with enzymatic C-term should be yielded."""
    # "AAKBBB": left-open ending at position 3: starts 1, 2 → "AK", "AK" wait...
    # Actually cut_set={0,3,6}; left-open: start not in cut_set, end in cut_set
    # end=3: start=2→"AK", start=1→"AAK" (start=0 is in cut_set, skip)
    # end=6: start=5→"B", start=4→"BB", start=3→"BBB" (enzymatic), start=2→"AKBBB"
    entry = _make_entry("AAKBBB")
    result = list(digest_peff_sequence(
        entry, _cfg(cleave_on="K", min_length=2, semi_enzymatic=True),
    ))
    seqs = {str(p.proforma) for p in result}
    assert "AK" in seqs  # left-open, ends at enzymatic cut after K
    semi_peptides = [p for p in result if p.semi_enzymatic]
    assert any(p.sequence == "AK" and p.semi_enzymatic for p in semi_peptides)


def test_semi_enzymatic_includes_fully_enzymatic() -> None:
    """Semi-enzymatic mode must still yield fully-enzymatic peptides."""
    entry = _make_entry("AAKBBB")
    result = list(digest_peff_sequence(
        entry, _cfg(cleave_on="K", semi_enzymatic=True),
    ))
    seqs = {str(p.proforma) for p in result}
    assert "AAK" in seqs
    assert "BBB" in seqs
    fully_enzymatic = [p for p in result if not p.semi_enzymatic]
    assert any(p.sequence == "AAK" for p in fully_enzymatic)
    assert any(p.sequence == "BBB" for p in fully_enzymatic)


def test_semi_enzymatic_with_missed_cleavages() -> None:
    """Semi-enzymatic should work correctly with missed cleavages > 0."""
    # "AAKBBKCC": cuts at [0, 3, 6, 8]
    # With missed_cleavages=1 and semi_enzymatic=True:
    # Fully enzymatic (mc=0): "AAK", "BBK", "CC"
    # Fully enzymatic (mc=1): "AAKBBK", "BBKCC"
    # Plus semi-enzymatic spans
    entry = _make_entry("AAKBBKCC")
    result = list(digest_peff_sequence(
        entry, _cfg(cleave_on="K", missed_cleavages=1, semi_enzymatic=True, min_length=2),
    ))
    seqs = {p.sequence for p in result}
    # Fully enzymatic should be present
    assert "AAK" in seqs
    assert "BBK" in seqs
    assert "CC" in seqs
    assert "AAKBBK" in seqs  # 1 missed
    assert "BBKCC" in seqs   # 1 missed
    # Semi-enzymatic examples
    assert "AA" in seqs   # right-open from pos 0
    assert "BB" in seqs   # right-open from pos 3
    semi = [p for p in result if p.semi_enzymatic]
    assert len(semi) > 0
