from __future__ import annotations

from peff_digest import digest_peff_sequence
from peff_digest.config import TerminalMod

from conftest import _cfg, _make_entry


def test_terminal_mod_fixed_nterm_all_peptides() -> None:
    """Fixed nterm mod with protein_terminus=False must appear on every peptide."""
    entry = _make_entry("AAKBBR")
    result = list(digest_peff_sequence(
        entry, _cfg(cleave_on="KR",
                     terminal_mods=[TerminalMod(modification="Acetyl", position="nterm", mod_type="fixed")]),
    ))
    seqs = [str(p.proforma) for p in result]
    assert len(seqs) == 2
    assert all(s.startswith("[Acetyl]") for s in seqs), f"Not all peptides have nterm Acetyl: {seqs}"


def test_terminal_mod_fixed_nterm_protein_terminus_only() -> None:
    """Fixed nterm mod with protein_terminus=True must appear only on the first peptide."""
    entry = _make_entry("AAKBBR")
    result = list(digest_peff_sequence(
        entry, _cfg(cleave_on="KR",
                     terminal_mods=[TerminalMod(
                         modification="Acetyl", position="nterm", mod_type="fixed", protein_terminus=True)]),
    ))
    seqs_by_seq = {p.sequence: str(p.proforma) for p in result}
    assert seqs_by_seq["AAK"].startswith("[Acetyl]"), "Protein N-term peptide must have Acetyl"
    assert not seqs_by_seq["BBR"].startswith("[Acetyl]"), "Internal peptide must not have Acetyl"


def test_terminal_mod_variable_nterm() -> None:
    """Variable nterm mod must produce both unmodified and modified forms."""
    entry = _make_entry("AAKBBR")
    result = list(digest_peff_sequence(
        entry, _cfg(cleave_on="KR", max_ptm_per_peptide=1,
                     terminal_mods=[TerminalMod(modification="Acetyl", position="nterm", mod_type="variable")]),
    ))
    aak_peptides = [str(p.proforma) for p in result if p.sequence == "AAK"]
    assert "AAK" in aak_peptides, "Unmodified form must be present"
    assert any(s.startswith("[Acetyl]") for s in aak_peptides), "Acetylated form must be present"


def test_terminal_mod_fixed_nterm_residue_filter() -> None:
    """Fixed nterm mod with residue='A' must only fire when the terminal AA is A."""
    # "AAKBBR": "AAK" starts with A → gets mod; "BBR" starts with B → no mod
    entry = _make_entry("AAKBBR")
    result = list(digest_peff_sequence(
        entry, _cfg(cleave_on="KR",
                     terminal_mods=[TerminalMod(
                         modification="Acetyl", position="nterm", mod_type="fixed", residue="A")]),
    ))
    seqs_by_seq = {p.sequence: str(p.proforma) for p in result}
    assert seqs_by_seq["AAK"].startswith("[Acetyl]"), "AAK starts with A, must have Acetyl"
    assert not seqs_by_seq["BBR"].startswith("[Acetyl]"), "BBR starts with B, must not have Acetyl"
