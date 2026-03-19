from __future__ import annotations

import pefftacular as pf

from peff_digest import digest_peff_sequence, get_cut_sites


def _make_entry(sequence: str, **kwargs) -> pf.SequenceEntry:
    return pf.SequenceEntry(prefix="sp", db_unique_id="TEST_ID", sequence=sequence, **kwargs)


# ---------------------------------------------------------------------------
# get_cut_sites
# ---------------------------------------------------------------------------


def test_get_cut_sites_basic_trypsin():
    # K and R are cleavage sites; P after cut suppresses it
    cuts = get_cut_sites("PEPTKPEPTIDE", cleave_on={"K", "R"}, restrict_after={"P"}, restrict_before=set(), cterminal=True)
    # K is at index 4, next char is 'P' → restricted, no cut there
    assert cuts == [0, 12]


def test_get_cut_sites_cleaves_on_kr():
    cuts = get_cut_sites("ACKR", cleave_on={"K", "R"}, restrict_after=set(), restrict_before=set(), cterminal=True)
    # K at index 2 → cut after → position 3; R at index 3 → terminal, skipped
    assert cuts == [0, 3, 4]


def test_get_cut_sites_nterminal():
    # N-terminal cleavage: cut is placed at the K/R position itself
    cuts = get_cut_sites("AAKB", cleave_on={"K"}, restrict_after=set(), restrict_before=set(), cterminal=False)
    # K at index 2 → cut_pos = 2 (not 0, not 4)
    assert cuts == [0, 2, 4]


def test_get_cut_sites_no_cleavage_sites():
    cuts = get_cut_sites("AAAA", cleave_on={"K", "R"}, restrict_after=set(), restrict_before=set(), cterminal=True)
    assert cuts == [0, 4]


def test_get_cut_sites_restrict_before():
    # restrict_before checks protein_sequence[cut_pos - 1], which for cterminal=True
    # is the cleavage residue itself (cut_pos = i+1, so cut_pos-1 = i = K's position).
    # restrict_before={"K"} therefore suppresses all K cleavages.
    cuts = get_cut_sites("AAKR", cleave_on={"K"}, restrict_after=set(), restrict_before={"K"}, cterminal=True)
    assert cuts == [0, 4]


# ---------------------------------------------------------------------------
# digest_peff_sequence
# ---------------------------------------------------------------------------


def test_digest_returns_set():
    entry = _make_entry("ACDEFGHIKLM")
    result = digest_peff_sequence(entry, cleave_on="K", missed_cleavages=0, max_ptm_per_peptide=0)
    assert isinstance(result, set)


def test_digest_min_length_filters_short_peptides():
    # "AAKBBBBBBBB" → peptides: "AAK" (len 3) and "BBBBBBBB" (len 8)
    entry = _make_entry("AAKBBBBBBBB")
    result = digest_peff_sequence(
        entry,
        cleave_on="K",
        missed_cleavages=0,
        max_ptm_per_peptide=0,
        min_length=5,
        max_length=40,
        fixed_mods=None,
        variable_mods=None,
    )
    sequences = {str(ann) for ann in result}
    assert "AAK" not in sequences
    assert "BBBBBBBB" in sequences


def test_digest_max_length_filters_long_peptides():
    entry = _make_entry("AAAAAAAAAAAK")  # 11 As + K = 12 chars, one peptide "AAAAAAAAAAAK"
    result = digest_peff_sequence(
        entry,
        cleave_on="K",
        missed_cleavages=0,
        max_ptm_per_peptide=0,
        min_length=1,
        max_length=5,
        fixed_mods=None,
        variable_mods=None,
    )
    # The only enzymatic peptide is 12 chars long — should be filtered out
    assert len(result) == 0


def test_digest_missed_cleavages():
    # "AAKBBR" → 0 missed: ["AAK", "BBR"]; 1 missed: also "AAKBBR"
    entry = _make_entry("AAKBBR")
    result_0 = digest_peff_sequence(
        entry,
        cleave_on="KR",
        missed_cleavages=0,
        max_ptm_per_peptide=0,
        min_length=1,
        max_length=40,
        fixed_mods=None,
        variable_mods=None,
    )
    result_1 = digest_peff_sequence(
        entry,
        cleave_on="KR",
        missed_cleavages=1,
        max_ptm_per_peptide=0,
        min_length=1,
        max_length=40,
        fixed_mods=None,
        variable_mods=None,
    )
    seqs_0 = {str(a) for a in result_0}
    seqs_1 = {str(a) for a in result_1}
    assert "AAK" in seqs_0
    assert "BBR" in seqs_0
    assert "AAKBBR" not in seqs_0
    assert "AAKBBR" in seqs_1


def test_digest_no_ptms_by_default_returns_unmodified():
    entry = _make_entry("ACDEFGR")
    result = digest_peff_sequence(
        entry,
        cleave_on="R",
        missed_cleavages=0,
        max_ptm_per_peptide=0,
        min_length=1,
        max_length=40,
        fixed_mods=None,
        variable_mods=None,
    )
    seqs = {str(a) for a in result}
    assert "ACDEFGR" in seqs


def test_digest_fixed_mod_applied():
    entry = _make_entry("ACAR")
    result = digest_peff_sequence(
        entry,
        cleave_on="R",
        missed_cleavages=0,
        max_ptm_per_peptide=0,
        min_length=1,
        max_length=40,
        fixed_mods={"C": "Carbamidomethyl"},
        variable_mods=None,
    )
    seqs = {str(a) for a in result}
    # The peptide with the fixed mod should be present, unmodified should not
    assert not any(s == "ACAR" for s in seqs)
    assert any("Carbamidomethyl" in s for s in seqs)
