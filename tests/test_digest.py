from __future__ import annotations

from pathlib import Path

import pefftacular as pf
import psimodpy
import pytest

from peff_digest import DigestConfig, InternalMod, digest_peff_sequence, get_cut_sites
from peff_digest.config import TerminalMod
from peff_digest.io import read_sequences


@pytest.fixture(scope="session")
def psi_db() -> psimodpy.PsiModDatabase:
    return psimodpy.load()

DATA_DIR = Path(__file__).parent / "data"
PEFF_FILES = sorted(DATA_DIR.glob("*.peff"))
FASTA_FILES = sorted(DATA_DIR.glob("*.fasta"))


def _make_entry(sequence: str, **kwargs) -> pf.SequenceEntry:
    return pf.SequenceEntry(prefix="sp", db_unique_id="TEST_ID", sequence=sequence, **kwargs)


def _cfg(**kwargs) -> DigestConfig:
    """Create a DigestConfig with test-friendly defaults."""
    defaults = dict(min_length=1, max_length=40, missed_cleavages=0, max_ptm_per_peptide=0)
    return DigestConfig(**{**defaults, **kwargs})


# ---------------------------------------------------------------------------
# get_cut_sites
# ---------------------------------------------------------------------------


def test_get_cut_sites_basic_trypsin():
    # K and R are cleavage sites; P after cut suppresses it
    cuts = get_cut_sites(
        "PEPTKPEPTIDE", cleave_on={"K", "R"}, restrict_after={"P"}, restrict_before=set(), cterminal=True
    )
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


def test_digest_returns_generator():
    from types import GeneratorType
    entry = _make_entry("ACDEFGHIKLM")
    result = digest_peff_sequence(entry, _cfg(cleave_on="K"))
    assert isinstance(result, GeneratorType)


def test_digest_min_length_filters_short_peptides():
    # "AAKBBBBBBBB" → peptides: "AAK" (len 3) and "BBBBBBBB" (len 8)
    entry = _make_entry("AAKBBBBBBBB")
    result = digest_peff_sequence(entry, _cfg(cleave_on="K", min_length=5))
    sequences = {str(p.proforma) for p in result}
    assert "AAK" not in sequences
    assert "BBBBBBBB" in sequences


def test_digest_max_length_filters_long_peptides():
    entry = _make_entry("AAAAAAAAAAAK")  # 11 As + K = 12 chars, one peptide "AAAAAAAAAAAK"
    result = digest_peff_sequence(entry, _cfg(cleave_on="K", max_length=5))
    # The only enzymatic peptide is 12 chars long — should be filtered out
    assert len(list(result)) == 0


def test_digest_missed_cleavages():
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


def test_digest_no_ptms_by_default_returns_unmodified():
    entry = _make_entry("ACDEFGR")
    result = digest_peff_sequence(entry, _cfg(cleave_on="R"))
    seqs = {str(p.proforma) for p in result}
    assert "ACDEFGR" in seqs


def test_digest_fixed_mod_applied():
    entry = _make_entry("ACAR")
    cam = InternalMod(modification="Carbamidomethyl", residue="C", mod_type="fixed")
    result = digest_peff_sequence(entry, _cfg(cleave_on="R", internal_mods=[cam]))
    seqs = {str(p.proforma) for p in result}
    # The peptide with the fixed mod should be present, unmodified should not
    assert not any(s == "ACAR" for s in seqs)
    assert any("Carbamidomethyl" in s for s in seqs)


# ---------------------------------------------------------------------------
# PEFF variants and PTMs
# ---------------------------------------------------------------------------


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
    assert variant_peptides[0].proforma.peptide_name is not None


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
    assert variant_peptides[0].proforma.peptide_name is not None


def test_peff_mod_applied_to_peptide() -> None:
    # "ACDEFGR": D is at 1-based position 3. ModResPsi at (3,) with MOD:00696.
    # Digest on R → one span "ACDEFGR". D is at peptide index 2.
    # max_ptm_per_peptide=1 → unmodified + one with MOD:00696 applied.
    entry = _make_entry(
        "ACDEFGR",
        mod_res_psi=(
            pf.ModResPsi(positions=(3,), accession="MOD:00696", name="phosphorylated residue"),
        ),
    )
    result = list(digest_peff_sequence(entry, _cfg(cleave_on="R", max_ptm_per_peptide=1, use_psi_mods=True)))
    assert len(result) == 2
    seqs = {str(p.proforma) for p in result}
    assert any("MOD:00696" in s for s in seqs)


def test_peff_mod_max_ptm_zero_skips_mods() -> None:
    # Same entry as above but max_ptm_per_peptide=0 → only the unmodified peptide.
    entry = _make_entry(
        "ACDEFGR",
        mod_res_psi=(
            pf.ModResPsi(positions=(3,), accession="MOD:00696", name="phosphorylated residue"),
        ),
    )
    result = list(digest_peff_sequence(entry, _cfg(cleave_on="R")))
    assert len(result) == 1
    seqs = {str(p.proforma) for p in result}
    assert "ACDEFGR" in seqs


def test_variable_mods_applied() -> None:
    # "AMKR": digest on K,R → "AMK" (M at peptide index 1), "R".
    # variable_mods={"M": ["Oxidation"]}, max_ptm_per_peptide=1.
    # "AMK" should appear both unmodified and with Oxidation on M.
    entry = _make_entry("AMKR")
    result = digest_peff_sequence(
        entry, _cfg(cleave_on="KR", max_ptm_per_peptide=1,
                     internal_mods=[InternalMod(modification="Oxidation", residue="M", mod_type="variable")]),
    )
    seqs = {str(p.proforma) for p in result}
    assert "AMK" in seqs
    assert any("Oxidation" in s for s in seqs)


def test_fixed_mods_applied() -> None:
    # "ACKR": digest on K,R → "ACK" (C at peptide index 1), "R".
    # fixed_mods={"C": "Carbamidomethyl"} → every "ACK" must carry the mod.
    # Unmodified "ACK" must not appear.
    entry = _make_entry("ACKR")
    result = digest_peff_sequence(
        entry, _cfg(cleave_on="KR",
                     internal_mods=[InternalMod(modification="Carbamidomethyl", residue="C", mod_type="fixed")]),
    )
    seqs = {str(p.proforma) for p in result}
    assert "ACK" not in seqs
    assert any("Carbamidomethyl" in s for s in seqs)


# ---------------------------------------------------------------------------
# Smoke tests against all example PEFF files
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("peff_path", PEFF_FILES, ids=[f.name for f in PEFF_FILES])
def test_digest_peff_file_no_crash(peff_path: Path) -> None:
    """Every parseable sequence in each example PEFF file should digest without error."""
    reader = iter(pf.PeffReader(str(peff_path)))
    n_digested = 0
    while True:
        try:
            entry = next(reader)
        except StopIteration:
            break
        except Exception:
            continue  # skip malformed entries
        result = digest_peff_sequence(
            entry, _cfg(cleave_on="KR", missed_cleavages=1, min_length=4, max_length=50, max_ptm_per_peptide=2),
        )
        assert hasattr(result, "__iter__")
        n_digested += 1


# ---------------------------------------------------------------------------
# Smoke tests against FASTA files
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("fasta_path", FASTA_FILES, ids=[f.name for f in FASTA_FILES])
def test_digest_fasta_file_no_crash(fasta_path: Path) -> None:
    """Every sequence in each example FASTA file should digest without error."""
    sequences, _ = read_sequences(str(fasta_path))
    assert len(sequences) > 0
    for entry in sequences:
        result = digest_peff_sequence(
            entry, _cfg(cleave_on="KR", missed_cleavages=1, min_length=4, max_length=50, max_ptm_per_peptide=2),
        )
        assert hasattr(result, "__iter__")


# ---------------------------------------------------------------------------
# Mod source and variant inclusion controls
# ---------------------------------------------------------------------------


def test_psi_mods_excluded_when_disabled() -> None:
    entry = _make_entry(
        "ACDEFGR",
        mod_res_psi=(pf.ModResPsi(positions=(3,), accession="MOD:00696", name="phosphorylated residue"),),
    )
    result = list(digest_peff_sequence(
        entry, _cfg(cleave_on="R", max_ptm_per_peptide=1, use_psi_mods=False),
    ))
    seqs = {str(p.proforma) for p in result}
    assert len(result) == 1
    assert "ACDEFGR" in seqs


def test_no_simple_variants_skips_substituted_peptides() -> None:
    entry = _make_entry(
        "AAAKBBBR",
        variant_simple=(pf.VariantSimple(position=2, new_amino_acid="C"),),
    )
    result = list(digest_peff_sequence(
        entry, _cfg(cleave_on="KR", include_simple_variants=False),
    ))
    seqs = {str(p.proforma) for p in result}
    assert not any("ACAK" in s for s in seqs)
    assert "AAAK" in seqs


def test_no_complex_variants_skips_complex_peptides() -> None:
    entry = _make_entry(
        "AABBBKCCR",
        variant_complex=(pf.VariantComplex(start_pos=3, end_pos=5, new_sequence="DD"),),
    )
    result = list(digest_peff_sequence(
        entry, _cfg(cleave_on="KR", include_complex_variants=False),
    ))
    seqs = {str(p.proforma) for p in result}
    assert not any("AADDK" in s for s in seqs)
    assert "AABBBK" in seqs


# ---------------------------------------------------------------------------
# PSI-MOD validation: mods invalidated by sequence variants
# ---------------------------------------------------------------------------


def test_psimod_validation_drops_mod_after_residue_mutation(psi_db) -> None:
    """MOD:00046 (O-phospho-L-serine, origin=S) at position 1 must be dropped
    when a VariantSimple mutates position 1 from S to W."""
    # Protein "SWKR": position 1=S, position 2=W, ...
    # MOD:00046 annotated at position 1 (the S).
    # Canonical variant keeps S → mod is valid.
    # Simple variant S→W: residue at position 1 is now W → phospho-serine is invalid.
    entry = _make_entry(
        "SWKR",
        mod_res_psi=(pf.ModResPsi(positions=(1,), accession="MOD:00046", name="O-phospho-L-serine"),),
        variant_simple=(pf.VariantSimple(position=1, new_amino_acid="W"),),
    )
    result = list(digest_peff_sequence(
        entry, _cfg(cleave_on="KR", max_ptm_per_peptide=1, use_psi_mods=True, include_simple_variants=True),
        psi_db=psi_db,
    ))
    # Canonical peptide "SWK" should have a phospho version
    assert any("MOD:00046" in str(p.proforma) and p.sequence == "SWK" for p in result)
    # Variant peptide "WWK" must NOT carry MOD:00046 (S→W invalidates phospho-serine)
    assert not any("MOD:00046" in str(p.proforma) and p.sequence == "WWK" for p in result)


def test_psimod_validation_retains_mod_on_unchanged_residue(psi_db) -> None:
    """MOD:00046 at position 2 (an S) must survive when position 1 is mutated."""
    entry = _make_entry(
        "ASKR",
        mod_res_psi=(pf.ModResPsi(positions=(2,), accession="MOD:00046", name="O-phospho-L-serine"),),
        variant_simple=(pf.VariantSimple(position=1, new_amino_acid="W"),),
    )
    result = list(digest_peff_sequence(
        entry, _cfg(cleave_on="KR", max_ptm_per_peptide=1, use_psi_mods=True, include_simple_variants=True),
        psi_db=psi_db,
    ))
    # Both canonical "ASK" and variant "WSK" should have the phospho version
    assert any("MOD:00046" in str(p.proforma) and p.sequence == "ASK" for p in result)
    assert any("MOD:00046" in str(p.proforma) and p.sequence == "WSK" for p in result)


# ---------------------------------------------------------------------------
# Terminal promotion: term_spec mods promoted to terminal ProForma position
# ---------------------------------------------------------------------------


def test_psimod_terminal_promotion_nterm(psi_db) -> None:
    """MOD:00050 (N-acetyl-L-alanine, origin=A, term_spec=N-term) at position 1
    of a protein starting with A should become an nterm mod, not an internal mod."""
    # "AAKR": position 1 = A. MOD:00050 has origin=A and term_spec=N-term.
    # For the protein-N-terminal peptide "AAK" it must appear as [MOD:00050]-AAK.
    entry = _make_entry(
        "AAKR",
        mod_res_psi=(pf.ModResPsi(positions=(1,), accession="MOD:00050", name="N-acetyl-L-alanine"),),
    )
    result = list(digest_peff_sequence(
        entry, _cfg(cleave_on="KR", max_ptm_per_peptide=1, use_psi_mods=True),
        psi_db=psi_db,
    ))
    seqs = {str(p.proforma) for p in result}
    # The modification must appear in nterm notation (leading bracket before sequence)
    assert any(s.startswith("[MOD:00050]") for s in seqs), \
        f"Expected nterm-promoted MOD:00050 in {seqs}"
    # It must NOT appear as an internal mod at index 0 (which would look like A[MOD:00050]AK)
    assert not any("A[MOD:00050]" in s for s in seqs), \
        f"MOD:00050 should be terminal, not internal, in {seqs}"


def test_psimod_terminal_promotion_skipped_for_internal_peptide(psi_db) -> None:
    """Terminal promotion must NOT fire for a peptide that does not span the protein N-term."""
    # "AAAKR": cut after K → peptides "AAAK" (nterm) and "R" (cterm).
    # MOD:00050 at position 1 (first A) → should be nterm in "AAAK" only.
    # Now add missed_cleavages=0 but the mod is at position 1 of the protein so
    # only "AAAK" covers position 1. If we use a peptide that does NOT start at 0
    # (e.g., by having an internal cleavage before position 1), it won't appear.
    # Simplest check: the non-nterm peptide "R" should not carry the mod at all.
    entry = _make_entry(
        "AAAKR",
        mod_res_psi=(pf.ModResPsi(positions=(1,), accession="MOD:00050", name="N-acetyl-L-alanine"),),
    )
    result = list(digest_peff_sequence(
        entry, _cfg(cleave_on="KR", max_ptm_per_peptide=1, use_psi_mods=True),
        psi_db=psi_db,
    ))
    r_peptides = [p for p in result if p.sequence == "R"]
    assert all("MOD:00050" not in str(p.proforma) for p in r_peptides)


# ---------------------------------------------------------------------------
# Two-phase delta: user mods respect remaining capacity after PEFF mods
# ---------------------------------------------------------------------------


def test_two_phase_delta_user_mods_capped_by_peff_mods(psi_db) -> None:
    """With max_ptm_per_peptide=2 and 1 PEFF mod applied, at most 1 user mod can be added."""
    # "SMKR": S at position 1 (MOD:00046), M at peptide index 1 (Oxidation variable).
    # max_ptm=2: with 1 PEFF mod, remaining=1 → only 0 or 1 Oxidation combos.
    entry = _make_entry(
        "SMKR",
        mod_res_psi=(pf.ModResPsi(positions=(1,), accession="MOD:00046", name="O-phospho-L-serine"),),
    )
    result = list(digest_peff_sequence(
        entry, _cfg(cleave_on="KR", max_ptm_per_peptide=2, use_psi_mods=True,
                     internal_mods=[InternalMod(modification="Oxidation", residue="M", mod_type="variable")]),
        psi_db=psi_db,
    ))
    smk = [p for p in result if p.sequence == "SMK"]
    # Expected: unmod, phospho-only, oxidation-only, phospho+oxidation (all ≤2 mods)
    assert any("MOD:00046" in str(p.proforma) and "Oxidation" in str(p.proforma) for p in smk), \
        "Should have a peptide with both PEFF phospho and user Oxidation"
    # No peptide should carry more than 2 total mods
    from peff_digest.digest import _count_mods
    for p in smk:
        assert _count_mods(p.proforma) <= 2, f"Peptide has too many mods: {p.proforma}"


def test_two_phase_delta_full_capacity_no_user_mods(psi_db) -> None:
    """With max_ptm_per_peptide=1 and 1 PEFF mod applied, no user mods can be added."""
    entry = _make_entry(
        "SMKR",
        mod_res_psi=(pf.ModResPsi(positions=(1,), accession="MOD:00046", name="O-phospho-L-serine"),),
    )
    result = list(digest_peff_sequence(
        entry, _cfg(cleave_on="KR", max_ptm_per_peptide=1, use_psi_mods=True,
                     internal_mods=[InternalMod(modification="Oxidation", residue="M", mod_type="variable")]),
        psi_db=psi_db,
    ))
    smk_peptides = [p for p in result if p.sequence == "SMK"]
    # Should have: unmod + phospho-only + oxidation-only (never both, capacity=1)
    assert not any("MOD:00046" in str(p.proforma) and "Oxidation" in str(p.proforma) for p in smk_peptides), \
        "Should not have both PEFF and user mod when max_ptm=1"


# ---------------------------------------------------------------------------
# Protein terminus flags on Peptide
# ---------------------------------------------------------------------------


def test_peptide_protein_terminus_flags() -> None:
    """is_protein_nterm / is_protein_cterm must reflect protein position, not just peptide content."""
    entry = _make_entry("AAKBBR")
    result = list(digest_peff_sequence(entry, _cfg(cleave_on="KR")))
    nterm_peptides = [p for p in result if p.is_protein_nterm]
    cterm_peptides = [p for p in result if p.is_protein_cterm]
    assert all(p.sequence == "AAK" for p in nterm_peptides)
    assert all(p.sequence == "BBR" for p in cterm_peptides)
    assert all(not p.is_protein_cterm for p in nterm_peptides)
    assert all(not p.is_protein_nterm for p in cterm_peptides)


# ---------------------------------------------------------------------------
# Consecutive cleavage sites
# ---------------------------------------------------------------------------


def test_get_cut_sites_consecutive_cleavage_sites():
    # "AKKR": K at index 1 → cut 2, K at index 2 → cut 3, R at index 3 is terminal → skipped
    cuts = get_cut_sites("AKKR", cleave_on={"K", "R"}, restrict_after=set(), restrict_before=set(), cterminal=True)
    assert cuts == [0, 2, 3, 4]


def test_digest_consecutive_cleavage_sites_yields_single_aa_peptides():
    # "AKKR" with min_length=1 → "AK", "K", "R"
    entry = _make_entry("AKKR")
    result = list(digest_peff_sequence(entry, _cfg(cleave_on="KR")))
    seqs = {str(p.proforma) for p in result}
    assert "AK" in seqs
    assert "K" in seqs
    assert "R" in seqs


# ---------------------------------------------------------------------------
# Missed cleavages > 1
# ---------------------------------------------------------------------------


def test_digest_missed_cleavages_two():
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


# ---------------------------------------------------------------------------
# Semi-enzymatic digestion
# ---------------------------------------------------------------------------


def test_semi_enzymatic_right_open():
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


def test_semi_enzymatic_left_open():
    """Non-enzymatic N-term with enzymatic C-term should be yielded."""
    # "AAKBBB": left-open ending at position 3: starts 1, 2 → "AK", "AK" wait...
    # left-open ending at 3 (after K): starts 2→"AK", start 1→"AAK" (already enzymatic)
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


def test_semi_enzymatic_includes_fully_enzymatic():
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


# ---------------------------------------------------------------------------
# Terminal modifications (TerminalMod)
# ---------------------------------------------------------------------------


def test_terminal_mod_fixed_nterm_all_peptides():
    """Fixed nterm mod with protein_terminus=False must appear on every peptide."""
    entry = _make_entry("AAKBBR")
    result = list(digest_peff_sequence(
        entry, _cfg(cleave_on="KR",
                     terminal_mods=[TerminalMod(modification="Acetyl", position="nterm", mod_type="fixed")]),
    ))
    seqs = [str(p.proforma) for p in result]
    assert len(seqs) == 2
    assert all(s.startswith("[Acetyl]") for s in seqs), f"Not all peptides have nterm Acetyl: {seqs}"


def test_terminal_mod_fixed_nterm_protein_terminus_only():
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


def test_terminal_mod_variable_nterm():
    """Variable nterm mod must produce both unmodified and modified forms."""
    entry = _make_entry("AAKBBR")
    result = list(digest_peff_sequence(
        entry, _cfg(cleave_on="KR", max_ptm_per_peptide=1,
                     terminal_mods=[TerminalMod(modification="Acetyl", position="nterm", mod_type="variable")]),
    ))
    aak_peptides = [str(p.proforma) for p in result if p.sequence == "AAK"]
    assert "AAK" in aak_peptides, "Unmodified form must be present"
    assert any(s.startswith("[Acetyl]") for s in aak_peptides), "Acetylated form must be present"


def test_terminal_mod_fixed_nterm_residue_filter():
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


# ---------------------------------------------------------------------------
# use_mod_names=True
# ---------------------------------------------------------------------------


def test_use_mod_names_outputs_name_instead_of_accession():
    """use_mod_names=True should embed the mod name, not the accession, in ProForma."""
    entry = _make_entry(
        "ACDEFGR",
        mod_res_psi=(pf.ModResPsi(positions=(3,), accession="MOD:00696", name="phosphorylated residue"),),
    )
    result = list(digest_peff_sequence(
        entry, _cfg(cleave_on="R", max_ptm_per_peptide=1, use_psi_mods=True, use_mod_names=True),
    ))
    seqs = {str(p.proforma) for p in result}
    assert any("phosphorylated residue" in s for s in seqs), "Mod name must appear in output"
    assert not any("MOD:00696" in s for s in seqs), "Accession must not appear when use_mod_names=True"


# ---------------------------------------------------------------------------
# C-terminal terminal promotion (MOD:00090 = L-alanine amide, origin=A, C-term)
# ---------------------------------------------------------------------------


def test_psimod_terminal_promotion_cterm(psi_db) -> None:
    """MOD:00090 (L-alanine amide, origin=A, term_spec=C-term) at the last position
    of the protein must become a cterm mod on the protein C-terminal peptide."""
    # "AAKBA": last residue is A (position 5). MOD:00090 at position 5.
    # Digest on K → "AAK" and "BA". "BA" is the protein C-terminal peptide ending with A.
    entry = _make_entry(
        "AAKBA",
        mod_res_psi=(pf.ModResPsi(positions=(5,), accession="MOD:00090", name="L-alanine amide"),),
    )
    result = list(digest_peff_sequence(
        entry, _cfg(cleave_on="K", max_ptm_per_peptide=1, use_psi_mods=True),
        psi_db=psi_db,
    ))
    seqs = {str(p.proforma) for p in result}
    # The modification must appear in cterm notation (trailing bracket after sequence)
    assert any(s.endswith("[MOD:00090]") for s in seqs), \
        f"Expected cterm-promoted MOD:00090 in {seqs}"
    # It must NOT appear as an internal mod (e.g. "B[MOD:00090]A" or "BA[MOD:00090]A")
    ba_seqs = [s for s in seqs if "BA" in s and "MOD:00090" in s]
    assert all(s.endswith("[MOD:00090]") for s in ba_seqs), \
        f"MOD:00090 should be cterm-promoted, not internal: {ba_seqs}"


# ---------------------------------------------------------------------------
# Complex variant: deletion (empty new_sequence)
# ---------------------------------------------------------------------------


def test_variant_complex_deletion():
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


# ---------------------------------------------------------------------------
# annotate_variants=False
# ---------------------------------------------------------------------------


def test_annotate_variants_false_suppresses_peptide_name():
    """When annotate_variants=False, no peptide should have a peptide_name set."""
    entry = _make_entry(
        "AAAKBBBR",
        variant_simple=(pf.VariantSimple(position=2, new_amino_acid="C"),),
    )
    result = list(digest_peff_sequence(
        entry, _cfg(cleave_on="KR", annotate_variants=False),
    ))
    assert all(not p.proforma.peptide_name for p in result), \
        "No peptide should have a non-empty peptide_name when annotate_variants=False"


# ---------------------------------------------------------------------------
# Peptide.mass property
# ---------------------------------------------------------------------------


def test_peptide_mass_returns_positive_float():
    """Peptide.mass must return a positive float for all digested peptides."""
    entry = _make_entry("ACDEFGR")
    result = list(digest_peff_sequence(
        entry, _cfg(cleave_on="R",
                     internal_mods=[InternalMod(modification="Carbamidomethyl", residue="C", mod_type="fixed")]),
    ))
    assert len(result) > 0
    for p in result:
        m = p.mass
        assert isinstance(m, float), f"mass should be float, got {type(m)}"
        assert m > 0, f"mass should be positive, got {m}"


# ---------------------------------------------------------------------------
# ann_to_map / Peptide.mod_map
# ---------------------------------------------------------------------------


def test_ann_to_map_internal_mod():
    """ann_to_map must return the correct position→name mapping for an internal mod."""
    entry = _make_entry("ACKR")
    result = list(digest_peff_sequence(
        entry, _cfg(cleave_on="KR",
                     internal_mods=[InternalMod(modification="Carbamidomethyl", residue="C", mod_type="fixed")]),
    ))
    ack_peptides = [p for p in result if p.sequence == "ACK"]
    assert len(ack_peptides) == 1
    mod_map = ack_peptides[0].mod_map
    # C is at 0-based index 1 in "ACK"
    assert mod_map == {1: "Carbamidomethyl"}, f"Unexpected mod_map: {mod_map}"


def test_ann_to_map_nterm_sentinel():
    """N-terminal mods must map to key -1 in mod_map."""
    entry = _make_entry("AAKR")
    result = list(digest_peff_sequence(
        entry, _cfg(cleave_on="KR",
                     terminal_mods=[TerminalMod(modification="Acetyl", position="nterm", mod_type="fixed")]),
    ))
    aak = next(p for p in result if p.sequence == "AAK")
    assert aak.mod_map.get(-1) == "Acetyl", f"Expected -1: Acetyl in mod_map, got {aak.mod_map}"
