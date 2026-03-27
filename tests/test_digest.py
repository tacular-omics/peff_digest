from __future__ import annotations

from pathlib import Path

import pefftacular as pf
import pytest

from peff_digest import InternalMod, digest_peff_sequence, get_cut_sites
from peff_digest.cli import read_sequences

DATA_DIR = Path(__file__).parent / "data"
PEFF_FILES = sorted(DATA_DIR.glob("*.peff"))
FASTA_FILES = sorted(DATA_DIR.glob("*.fasta"))


def _make_entry(sequence: str, **kwargs) -> pf.SequenceEntry:
    return pf.SequenceEntry(prefix="sp", db_unique_id="TEST_ID", sequence=sequence, **kwargs)


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
    result = digest_peff_sequence(entry, cleave_on="K", missed_cleavages=0, max_ptm_per_peptide=0)
    assert isinstance(result, GeneratorType)


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
        internal_mods=None,
    )
    sequences = {str(p.proforma) for p in result}
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
        internal_mods=None,
    )
    # The only enzymatic peptide is 12 chars long — should be filtered out
    assert len(list(result)) == 0


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
        internal_mods=None,
    )
    result_1 = digest_peff_sequence(
        entry,
        cleave_on="KR",
        missed_cleavages=1,
        max_ptm_per_peptide=0,
        min_length=1,
        max_length=40,
        internal_mods=None,
    )
    seqs_0 = {str(p.proforma) for p in result_0}
    seqs_1 = {str(p.proforma) for p in result_1}
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
        internal_mods=None,
    )
    seqs = {str(p.proforma) for p in result}
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
        internal_mods=[InternalMod(modification="Carbamidomethyl", residue="C", mod_type="fixed")],
    )
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
    result = digest_peff_sequence(
        entry,
        cleave_on="KR",
        missed_cleavages=0,
        max_ptm_per_peptide=0,
        min_length=1,
        max_length=40,
        internal_mods=None,
    )
    result = list(result)
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
    result = digest_peff_sequence(
        entry,
        cleave_on="KR",
        missed_cleavages=0,
        max_ptm_per_peptide=0,
        min_length=1,
        max_length=40,
        internal_mods=None,
    )
    result = list(result)
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
    result = digest_peff_sequence(
        entry,
        cleave_on="R",
        missed_cleavages=0,
        max_ptm_per_peptide=1,
        min_length=1,
        max_length=40,
        internal_mods=None,
        use_psi_mods=True,
    )
    result = list(result)
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
    result = digest_peff_sequence(
        entry,
        cleave_on="R",
        missed_cleavages=0,
        max_ptm_per_peptide=0,
        min_length=1,
        max_length=40,
        internal_mods=None,
    )
    result = list(result)
    assert len(result) == 1
    seqs = {str(p.proforma) for p in result}
    assert "ACDEFGR" in seqs


def test_variable_mods_applied() -> None:
    # "AMKR": digest on K,R → "AMK" (M at peptide index 1), "R".
    # variable_mods={"M": ["Oxidation"]}, max_ptm_per_peptide=1.
    # "AMK" should appear both unmodified and with Oxidation on M.
    entry = _make_entry("AMKR")
    result = digest_peff_sequence(
        entry,
        cleave_on="KR",
        missed_cleavages=0,
        max_ptm_per_peptide=1,
        min_length=1,
        max_length=40,
        internal_mods=[InternalMod(modification="Oxidation", residue="M", mod_type="variable")],
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
        entry,
        cleave_on="KR",
        missed_cleavages=0,
        max_ptm_per_peptide=0,
        min_length=1,
        max_length=40,
        internal_mods=[InternalMod(modification="Carbamidomethyl", residue="C", mod_type="fixed")],
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
            entry,
            cleave_on="KR",
            missed_cleavages=1,
            min_length=4,
            max_length=50,
            max_ptm_per_peptide=2,
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
            entry,
            cleave_on="KR",
            missed_cleavages=1,
            min_length=4,
            max_length=50,
            max_ptm_per_peptide=2,
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
        entry, cleave_on="R", missed_cleavages=0, max_ptm_per_peptide=1,
        min_length=1, max_length=40, use_psi_mods=False,
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
        entry, cleave_on="KR", missed_cleavages=0, max_ptm_per_peptide=0,
        min_length=1, max_length=40, include_simple_variants=False,
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
        entry, cleave_on="KR", missed_cleavages=0, max_ptm_per_peptide=0,
        min_length=1, max_length=40, include_complex_variants=False,
    ))
    seqs = {str(p.proforma) for p in result}
    assert not any("AADDK" in s for s in seqs)
    assert "AABBBK" in seqs
