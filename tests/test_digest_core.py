from __future__ import annotations

from types import GeneratorType

import pytest
from conftest import _cfg, _make_entry

from peff_digest import DigestConfig, InternalMod, digest_peff_sequence


def test_digest_returns_generator() -> None:
    entry = _make_entry("ACDEFGHIKLM")
    result = digest_peff_sequence(entry, _cfg(cleave_on="K"))
    assert isinstance(result, GeneratorType)


def test_digest_min_length_filters_short_peptides() -> None:
    # "AAKBBBBBBBB" → peptides: "AAK" (len 3) and "BBBBBBBB" (len 8)
    entry = _make_entry("AAKBBBBBBBB")
    result = digest_peff_sequence(entry, _cfg(cleave_on="K", min_length=5))
    sequences = {str(p.proforma) for p in result}
    assert "AAK" not in sequences
    assert "BBBBBBBB" in sequences


def test_digest_max_length_filters_long_peptides() -> None:
    entry = _make_entry("AAAAAAAAAAAK")  # 11 As + K = 12 chars, one peptide "AAAAAAAAAAAK"
    result = digest_peff_sequence(entry, _cfg(cleave_on="K", max_length=5))
    # The only enzymatic peptide is 12 chars long — should be filtered out
    assert len(list(result)) == 0


def test_digest_missed_cleavages() -> None:
    # "AAKBBR" → 0 missed: ["AAK", "BBR"]; 1 missed: also "AAKBBR"
    entry = _make_entry("AAKBBR")
    result_0 = digest_peff_sequence(entry, _cfg(cleave_on="KR"))
    result_1 = digest_peff_sequence(entry, _cfg(cleave_on="KR", missed_cleavages=1))
    seqs_0 = {str(p.proforma) for p in result_0}
    seqs_1 = {str(p.proforma) for p in result_1}
    assert "AAK" in seqs_0
    assert "BBR" in seqs_0
    assert "AAKBBR" not in seqs_0
    assert "AAKBBR" in seqs_1


def test_digest_no_ptms_by_default_returns_unmodified() -> None:
    entry = _make_entry("ACDEFGR")
    result = digest_peff_sequence(entry, _cfg(cleave_on="R"))
    seqs = {str(p.proforma) for p in result}
    assert "ACDEFGR" in seqs


def test_digest_fixed_mod_applied() -> None:
    entry = _make_entry("ACAR")
    cam = InternalMod(modification="Carbamidomethyl", residue="C", mod_type="fixed")
    result = digest_peff_sequence(entry, _cfg(cleave_on="R", internal_mods=[cam]))
    seqs = {str(p.proforma) for p in result}
    # The peptide with the fixed mod should be present, unmodified should not
    assert not any(s == "ACAR" for s in seqs)
    assert any("Carbamidomethyl" in s for s in seqs)


def test_digest_variable_mods_applied() -> None:
    # "AMKR": digest on K,R → "AMK" (M at peptide index 1), "R".
    # variable_mods={"M": ["Oxidation"]}, max_ptm_per_peptide=1.
    # "AMK" should appear both unmodified and with Oxidation on M.
    entry = _make_entry("AMKR")
    result = digest_peff_sequence(
        entry,
        _cfg(
            cleave_on="KR",
            max_ptm_per_peptide=1,
            internal_mods=[InternalMod(modification="Oxidation", residue="M", mod_type="variable")],
        ),
    )
    seqs = {str(p.proforma) for p in result}
    assert "AMK" in seqs
    assert any("Oxidation" in s for s in seqs)


def test_digest_fixed_mods_applied() -> None:
    # "ACKR": digest on K,R → "ACK" (C at peptide index 1), "R".
    # fixed_mods={"C": "Carbamidomethyl"} → every "ACK" must carry the mod.
    # Unmodified "ACK" must not appear.
    entry = _make_entry("ACKR")
    result = digest_peff_sequence(
        entry,
        _cfg(
            cleave_on="KR",
            internal_mods=[InternalMod(modification="Carbamidomethyl", residue="C", mod_type="fixed")],
        ),
    )
    seqs = {str(p.proforma) for p in result}
    assert "ACK" not in seqs
    assert any("Carbamidomethyl" in s for s in seqs)


def test_digest_peptide_protein_terminus_flags() -> None:
    """is_protein_nterm / is_protein_cterm must reflect protein position, not just peptide content."""
    entry = _make_entry("AAKBBR")
    result = list(digest_peff_sequence(entry, _cfg(cleave_on="KR")))
    nterm_peptides = [p for p in result if p.is_protein_nterm]
    cterm_peptides = [p for p in result if p.is_protein_cterm]
    assert all(p.sequence == "AAK" for p in nterm_peptides)
    assert all(p.sequence == "BBR" for p in cterm_peptides)
    assert all(not p.is_protein_cterm for p in nterm_peptides)
    assert all(not p.is_protein_nterm for p in cterm_peptides)


def test_digest_consecutive_cleavage_sites_yields_single_aa_peptides() -> None:
    # "AKKR" with min_length=1 → "AK", "K", "R"
    entry = _make_entry("AKKR")
    result = list(digest_peff_sequence(entry, _cfg(cleave_on="KR")))
    seqs = {str(p.proforma) for p in result}
    assert "AK" in seqs
    assert "K" in seqs
    assert "R" in seqs


def test_digest_missed_cleavages_two() -> None:
    # "AAKBBRCCQ" → 0 missed: ["AAK","BBR","CCQ"]; 2 missed also includes "AAKBBRCCQ"
    entry = _make_entry("AAKBBRCCQ")
    result = list(digest_peff_sequence(entry, _cfg(cleave_on="KR", missed_cleavages=2)))
    seqs = {str(p.proforma) for p in result}
    assert "AAK" in seqs
    assert "BBR" in seqs
    assert "CCQ" in seqs
    assert "AAKBBR" in seqs
    assert "BBRCCQ" in seqs
    assert "AAKBBRCCQ" in seqs


def test_digest_peptide_mass_returns_positive_float() -> None:
    """Peptide.mass must return a positive float for all digested peptides."""
    entry = _make_entry("ACDEFGR")
    result = list(
        digest_peff_sequence(
            entry,
            _cfg(
                cleave_on="R",
                internal_mods=[InternalMod(modification="Carbamidomethyl", residue="C", mod_type="fixed")],
            ),
        )
    )
    assert len(result) > 0
    for p in result:
        m = p.mass
        assert isinstance(m, float), f"mass should be float, got {type(m)}"
        assert m > 0, f"mass should be positive, got {m}"
