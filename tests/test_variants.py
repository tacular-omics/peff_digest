from __future__ import annotations

import pefftacular as pf
import pytest
from conftest import _cfg, _make_entry

from peff_digest import digest_peff_sequence
from peff_digest.digest import (
    _mods_in_span,
    _variant_in_span,
    _yield_mod_variants,
)


def _form(sequence: str, variant: pf.VariantSimple | pf.VariantComplex) -> pf.VariantForm:
    """The form pefftacular builds for *variant*, which is what the digest now walks."""
    key = "variant_simple" if isinstance(variant, pf.VariantSimple) else "variant_complex"
    entry = _make_entry(sequence, **{key: (variant,)})
    return entry.apply_variants([variant])


def test_variant_simple_produces_substituted_peptide() -> None:
    # Position 2 (1-based) is the second AA in "AAAKBBBR" → 'A' replaced by 'C'.
    # Canonical digest on K,R → "AAAK", "BBBR".
    # Variant sequence "ACAKBBBR" → "ACAK", "BBBR".
    entry = _make_entry(
        "AAAKBBBR",
        variant_simple=(pf.VariantSimple(position=2, new_amino_acid="C"),),
    )
    result = list(digest_peff_sequence(entry, _cfg(cleave_on="KR")))
    seqs = {str(p.proforma) for p in result}
    # Canonical peptide has no variant notation
    assert "AAAK" in seqs
    # Variant peptide carries PEFF notation in str() — check the AA sequence is present
    assert any("ACAK" in s for s in seqs)

    variant_peptides = [p for p in result if "ACAK" in str(p.proforma) and str(p.proforma) != "AAAK"]
    assert len(variant_peptides) == 1
    assert variant_peptides[0].variant is not None


def test_variant_complex_substitution() -> None:
    # "AABBBKCCR": positions 3-5 (1-based, "BBB") replaced with "DD"
    # → variant sequence "AADDKCCR".
    # Canonical digest on K,R: "AABBBK", "CCR".
    # Variant digest: "AADDK", "CCR".
    entry = _make_entry(
        "AABBBKCCR",
        variant_complex=(pf.VariantComplex(start_pos=3, end_pos=5, new_sequence="DD"),),
    )
    result = list(digest_peff_sequence(entry, _cfg(cleave_on="KR")))
    seqs = {str(p.proforma) for p in result}
    # Canonical peptide has no variant notation
    assert "AABBBK" in seqs
    # Variant peptide carries PEFF notation in str() — check the AA sequence is present
    assert any("AADDK" in s for s in seqs)

    variant_peptides = [p for p in result if "AADDK" in str(p.proforma) and str(p.proforma) != "AABBBK"]
    assert len(variant_peptides) == 1
    assert variant_peptides[0].variant is not None


def test_variant_complex_deletion() -> None:
    """VariantComplex with new_sequence='' deletes residues, producing a shorter peptide."""
    # "AABBBKR": canonical digest → "AABBBK", "R"
    # VariantComplex deletes positions 3-5 ("BBB") → variant sequence "AAKR"
    # Variant digest → "AAK", "R"
    entry = _make_entry(
        "AABBBKR",
        variant_complex=(pf.VariantComplex(start_pos=3, end_pos=5, new_sequence=""),),
    )
    result = list(digest_peff_sequence(entry, _cfg(cleave_on="KR")))
    seqs = {str(p.proforma) for p in result}
    assert "AABBBK" in seqs, "Canonical peptide must still be present"
    assert any("AAK" in s for s in seqs), "Deletion variant must produce shorter peptide"


def test_variant_complex_insertion_longer_sequence() -> None:
    """VariantComplex with new_sequence longer than replaced region."""
    # "AAKR": replace position 2 (1-based, "A") with "MMMM"
    # → variant "AMMMMKR" (length 7 vs original 4)
    entry = _make_entry(
        "AAKR",
        variant_complex=(pf.VariantComplex(start_pos=2, end_pos=2, new_sequence="MMMM"),),
    )
    result = list(digest_peff_sequence(entry, _cfg(cleave_on="KR")))
    seqs = {p.sequence for p in result}
    assert "AAK" in seqs  # canonical
    assert "AMMMMK" in seqs  # variant (longer)
    assert "R" in seqs  # both canonical and variant C-term


def test_no_simple_variants_skips_substituted_peptides() -> None:
    entry = _make_entry(
        "AAAKBBBR",
        variant_simple=(pf.VariantSimple(position=2, new_amino_acid="C"),),
    )
    result = list(
        digest_peff_sequence(
            entry,
            _cfg(cleave_on="KR", include_simple_variants=False),
        )
    )
    seqs = {str(p.proforma) for p in result}
    assert not any("ACAK" in s for s in seqs)
    assert "AAAK" in seqs


def test_no_complex_variants_skips_complex_peptides() -> None:
    entry = _make_entry(
        "AABBBKCCR",
        variant_complex=(pf.VariantComplex(start_pos=3, end_pos=5, new_sequence="DD"),),
    )
    result = list(
        digest_peff_sequence(
            entry,
            _cfg(cleave_on="KR", include_complex_variants=False),
        )
    )
    seqs = {str(p.proforma) for p in result}
    assert not any("AADDK" in s for s in seqs)
    assert "AABBBK" in seqs


def test_variant_field_set_when_peptide_contains_simple_variant() -> None:
    """Peptide.variant must be set when the span contains the VariantSimple site."""
    # "AAAKBBBR": variant at position 2 (A→C). Digest → "AAAK" contains pos 2, "BBBR" does not.
    entry = _make_entry(
        "AAAKBBBR",
        variant_simple=(pf.VariantSimple(position=2, new_amino_acid="C"),),
    )
    result = list(digest_peff_sequence(entry, _cfg(cleave_on="KR")))
    # Variant peptide "ACAK" should have .variant set
    acak = [p for p in result if p.sequence == "ACAK"]
    assert len(acak) == 1
    assert acak[0].variant is not None
    assert isinstance(acak[0].variant, pf.VariantSimple)
    assert acak[0].variant.position == 2


def test_variant_field_none_when_peptide_does_not_contain_simple_variant() -> None:
    """Peptide.variant must be None for peptides from a variant sequence that don't overlap the site."""
    # "AAAKBBBR": variant at position 2 (A→C). "BBBR" does not contain pos 2.
    entry = _make_entry(
        "AAAKBBBR",
        variant_simple=(pf.VariantSimple(position=2, new_amino_acid="C"),),
    )
    result = list(digest_peff_sequence(entry, _cfg(cleave_on="KR")))
    # "BBBR" peptides from variant sequence should have variant=None
    bbbr_from_variant = [p for p in result if p.sequence == "BBBR"]
    assert all(p.variant is None for p in bbbr_from_variant)


def test_variant_field_none_for_canonical_peptides() -> None:
    """Peptide.variant must be None for all canonical peptides."""
    entry = _make_entry(
        "AAAKBBBR",
        variant_simple=(pf.VariantSimple(position=2, new_amino_acid="C"),),
    )
    result = list(digest_peff_sequence(entry, _cfg(cleave_on="KR")))
    canonical = [p for p in result if p.sequence in ("AAAK", "BBBR") and p.variant is None]
    # There should be canonical "AAAK" and "BBBR" with no variant
    assert any(p.sequence == "AAAK" for p in canonical)
    assert any(p.sequence == "BBBR" for p in canonical)


def test_variant_field_set_for_complex_variant_in_span() -> None:
    """Peptide.variant must be set for complex variants when the span overlaps."""
    # "AABBBKCCR": positions 3-5 replaced with "DD" → variant "AADDKCCR"
    # "AADDK" overlaps the variant site, "CCR" does not.
    entry = _make_entry(
        "AABBBKCCR",
        variant_complex=(pf.VariantComplex(start_pos=3, end_pos=5, new_sequence="DD"),),
    )
    result = list(digest_peff_sequence(entry, _cfg(cleave_on="KR")))
    aaddk = [p for p in result if p.sequence == "AADDK"]
    assert len(aaddk) == 1
    assert aaddk[0].variant is not None
    assert isinstance(aaddk[0].variant, pf.VariantComplex)
    # "CCR" from variant sequence should have variant=None (doesn't overlap)
    ccr_all = [p for p in result if p.sequence == "CCR"]
    assert all(p.variant is None for p in ccr_all)


def test_variant_field_set_for_complex_deletion_in_span() -> None:
    """Peptide.variant should be set for a deletion variant when the span contains the site."""
    # "AABBBKR": delete positions 3-5 → variant "AAKR"
    # Variant peptide "AAK" (span [0,3)) contains the deletion site (start0=2).
    entry = _make_entry(
        "AABBBKR",
        variant_complex=(pf.VariantComplex(start_pos=3, end_pos=5, new_sequence=""),),
    )
    result = list(digest_peff_sequence(entry, _cfg(cleave_on="KR")))
    aak_variants = [p for p in result if p.sequence == "AAK" and p.variant is not None]
    assert len(aak_variants) == 1
    assert isinstance(aak_variants[0].variant, pf.VariantComplex)
    assert aak_variants[0].variant.new_sequence == ""


def test_variant_in_span_deletion_overlaps() -> None:
    """Deletion at positions 3-5 (1-based) should be detected in span that contains the site."""
    f = _form("AABBBKR", pf.VariantComplex(start_pos=3, end_pos=5, new_sequence=""))
    assert f.sequence == "AAKR"
    assert _variant_in_span(f, 0, 3)


def test_variant_in_span_deletion_outside() -> None:
    """Deletion should NOT be detected in a span that doesn't contain the site."""
    f = _form("AABBBKR", pf.VariantComplex(start_pos=3, end_pos=5, new_sequence=""))
    assert not _variant_in_span(f, 3, 5)


def test_variant_in_span_complex_insertion() -> None:
    """Long insertion should properly detect span overlap."""
    # "AAKR": replace position 2 (1-based, 'A') with "MMMMM"
    # → form "AMMMMMKR". Inserted region: [1, 1+5) = [1, 6) in form coordinates.
    f = _form("AAKR", pf.VariantComplex(start_pos=2, end_pos=2, new_sequence="MMMMM"))
    assert f.sequence == "AMMMMMKR"
    assert _variant_in_span(f, 0, 4)  # [0,4) overlaps [1,6)
    assert not _variant_in_span(f, 0, 1)  # [0,1) does NOT overlap [1,6)


def test_variant_in_span_boundary_exact() -> None:
    """Variant at exact span boundary (end-exclusive) should not overlap."""
    f = _form("AAABKR", pf.VariantComplex(start_pos=4, end_pos=4, new_sequence="X"))
    # Replaced slot sits at start0=3, new_len=1, so the site is [3, 4).
    # Span [0,3): 0 < 4 AND 3 < 3 → False
    assert not _variant_in_span(f, 0, 3)
    # Span [3, 6) should overlap: 3 < 4 AND 3 < 6 → True
    assert _variant_in_span(f, 3, 6)


def test_variant_in_span_simple_mapped_none_returns_false() -> None:
    """A position_map entry of None at the variant's position must give no overlap.

    Built by hand: pefftacular never deletes a residue for a VariantSimple except through
    a ``*`` stop, which ``truncated_at`` covers, so this guard has no natural input.
    """
    form = pf.VariantForm(
        entry=_make_entry("AAKR"),
        applied=(pf.VariantSimple(position=2, new_amino_acid="C"),),
        position_map=(None, 1, None, 3, 4),
    )
    assert not _variant_in_span(form, 0, 4)


def test_simple_variant_out_of_range_is_skipped() -> None:
    """A position beyond the sequence yields no extra form, so no extra peptide."""
    v = pf.VariantSimple(position=99, new_amino_acid="C")
    entry = _make_entry("AAAKBBBR", variant_simple=(v,))
    result = list(digest_peff_sequence(entry, _cfg(cleave_on="KR")))

    assert {p.sequence for p in result} == {"AAAK", "BBBR"}
    assert all(p.variant is None for p in result)


def test_complex_variant_out_of_range_is_skipped() -> None:
    """start/end beyond the sequence yields no extra form, so no extra peptide."""
    v = pf.VariantComplex(start_pos=10, end_pos=20, new_sequence="DD")
    entry = _make_entry("AAAKBBBR", variant_complex=(v,))
    result = list(digest_peff_sequence(entry, _cfg(cleave_on="KR")))

    assert {p.sequence for p in result} == {"AAAK", "BBBR"}
    assert all(p.variant is None for p in result)


def test_non_numeric_variant_position_is_skipped() -> None:
    """A '?' position cannot locate a variant, so it contributes no form."""
    v = pf.VariantSimple(position="?", new_amino_acid="C")
    entry = _make_entry("AAAKBBBR", variant_simple=(v,))
    result = list(digest_peff_sequence(entry, _cfg(cleave_on="KR")))

    assert {p.sequence for p in result} == {"AAAK", "BBBR"}


def test_nonsense_variant_truncates_and_drops_downstream_peptides() -> None:
    """A '*' substitution ends the protein (spec 3.3.8): nothing downstream survives.

    The reference form still yields the whole protein's peptides. The variant form is the
    sequence truncated before the stop, so its only new peptide is the one ending at the
    new C-terminus, and it carries the variant.
    """
    # Stop at position 7, inside the peptide "BBBBR", so the truncated form is "AAAKBB".
    v = pf.VariantSimple(position=7, new_amino_acid="*")
    entry = _make_entry("AAAKBBBBRCCCK", variant_simple=(v,))
    result = list(digest_peff_sequence(entry, _cfg(cleave_on="KR")))

    by_seq = {p.sequence: p for p in result}
    assert set(by_seq) == {"AAAK", "BBBBR", "CCCK", "BB"}
    # "BB" is the truncated form's new C-terminal peptide and carries the variant.
    assert by_seq["BB"].variant is v
    # The reference form's peptides are unaffected and unattributed.
    assert by_seq["BBBBR"].variant is None
    assert by_seq["CCCK"].variant is None
    # No peptide anywhere contains the stop character.
    assert all("*" not in p.sequence for p in result)


def test_mod_on_substituted_residue_is_dropped() -> None:
    """Spec 3.3.10: a modification on a residue a variant substituted does not carry over."""
    v = pf.VariantSimple(position=2, new_amino_acid="G")
    mod = pf.ModResPsi(positions=(2,), accession="MOD:00696", name="phosphorylated residue")
    entry = _make_entry("ASAKBBBR", variant_simple=(v,), mod_res_psi=(mod,))
    result = list(digest_peff_sequence(entry, _cfg(cleave_on="KR", max_ptm_per_peptide=1)))

    variant_peptides = [p for p in result if p.variant is not None]
    assert variant_peptides, "the variant form must still produce a peptide"
    assert all("MOD:00696" not in str(p.proforma) for p in variant_peptides)


def test_mods_in_span_skips_unknown_position() -> None:
    """A ModResPsi with position '?' must be silently skipped."""
    mod = pf.ModResPsi(positions=("?",), accession="MOD:00696", name="phosphorylated residue")
    result = _mods_in_span("ACDEFGR", 0, 7, [mod])

    assert result == []


def test_mods_in_span_skips_position_outside_sequence() -> None:
    """An annotation position past the end of this form's sequence is dropped, not indexed."""
    mod = pf.ModResPsi(positions=(681,), accession="MOD:00696", name="phosphorylated residue")
    result = _mods_in_span("ACDEFGR", 0, 7, [mod])

    assert result == []


def test_yield_mod_variants_skips_duplicate_position_combos() -> None:
    """Combinations where two mods share the same position must be skipped."""
    # Both mods target position 1 — the n=2 combo (pos1, pos1) must not appear.
    applicable = [(1, "Carbamidomethyl"), (1, "Oxidation")]
    result = _yield_mod_variants("ACK", applicable, max_ptm=2)

    sequences = [str(r) for r in result]
    # Unmodified + each single mod = 3 variants; the duplicate-position combo is skipped.
    assert len(sequences) == 3
    assert "ACK" in sequences
    assert "AC[Carbamidomethyl]K" in sequences
    assert "AC[Oxidation]K" in sequences
