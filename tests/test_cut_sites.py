from __future__ import annotations

from peff_digest import get_cut_sites


def test_get_cut_sites_basic_trypsin() -> None:
    # K and R are cleavage sites; P after cut suppresses it
    cuts = get_cut_sites(
        "PEPTKPEPTIDE", cleave_on={"K", "R"}, restrict_after={"P"}, restrict_before=set(), cterminal=True
    )
    # K is at index 4, next char is 'P' → restricted, no cut there
    assert cuts == [0, 12]


def test_get_cut_sites_cleaves_on_kr() -> None:
    cuts = get_cut_sites("ACKR", cleave_on={"K", "R"}, restrict_after=set(), restrict_before=set(), cterminal=True)
    # K at index 2 → cut after → position 3; R at index 3 → terminal, skipped
    assert cuts == [0, 3, 4]


def test_get_cut_sites_nterminal() -> None:
    # N-terminal cleavage: cut is placed at the K/R position itself
    cuts = get_cut_sites("AAKB", cleave_on={"K"}, restrict_after=set(), restrict_before=set(), cterminal=False)
    # K at index 2 → cut_pos = 2 (not 0, not 4)
    assert cuts == [0, 2, 4]


def test_get_cut_sites_no_cleavage_sites() -> None:
    cuts = get_cut_sites("AAAA", cleave_on={"K", "R"}, restrict_after=set(), restrict_before=set(), cterminal=True)
    assert cuts == [0, 4]


def test_get_cut_sites_restrict_before() -> None:
    # restrict_before checks protein_sequence[cut_pos - 1], which for cterminal=True
    # is the cleavage residue itself (cut_pos = i+1, so cut_pos-1 = i = K's position).
    # restrict_before={"K"} therefore suppresses all K cleavages.
    cuts = get_cut_sites("AAKR", cleave_on={"K"}, restrict_after=set(), restrict_before={"K"}, cterminal=True)
    assert cuts == [0, 4]


def test_get_cut_sites_consecutive_cleavage_sites() -> None:
    # "AKKR": K at index 1 → cut 2, K at index 2 → cut 3, R at index 3 is terminal → skipped
    cuts = get_cut_sites("AKKR", cleave_on={"K", "R"}, restrict_after=set(), restrict_before=set(), cterminal=True)
    assert cuts == [0, 2, 3, 4]
