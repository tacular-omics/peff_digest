from __future__ import annotations

import pefftacular as pf
import psimodpy
import pytest
from conftest import _cfg, _make_entry

from peff_digest import InternalMod, digest_peff_sequence


def test_peff_mod_applied_to_peptide() -> None:
    # "ACDEFGR": D is at 1-based position 3. ModResPsi at (3,) with MOD:00696.
    # Digest on R → one span "ACDEFGR". D is at peptide index 2.
    # max_ptm_per_peptide=1 → unmodified + one with MOD:00696 applied.
    entry = _make_entry(
        "ACDEFGR",
        mod_res_psi=(pf.ModResPsi(positions=(3,), accession="MOD:00696", name="phosphorylated residue"),),
    )
    result = list(digest_peff_sequence(entry, _cfg(cleave_on="R", max_ptm_per_peptide=1, use_psi_mods=True)))
    assert len(result) == 2
    seqs = {str(p.proforma) for p in result}
    assert any("MOD:00696" in s for s in seqs)


def test_peff_mod_max_ptm_zero_skips_mods() -> None:
    # Same entry as above but max_ptm_per_peptide=0 → only the unmodified peptide.
    entry = _make_entry(
        "ACDEFGR",
        mod_res_psi=(pf.ModResPsi(positions=(3,), accession="MOD:00696", name="phosphorylated residue"),),
    )
    result = list(digest_peff_sequence(entry, _cfg(cleave_on="R")))
    assert len(result) == 1
    seqs = {str(p.proforma) for p in result}
    assert "ACDEFGR" in seqs


def test_psi_mods_excluded_when_disabled() -> None:
    entry = _make_entry(
        "ACDEFGR",
        mod_res_psi=(pf.ModResPsi(positions=(3,), accession="MOD:00696", name="phosphorylated residue"),),
    )
    result = list(
        digest_peff_sequence(
            entry,
            _cfg(cleave_on="R", max_ptm_per_peptide=1, use_psi_mods=False),
        )
    )
    seqs = {str(p.proforma) for p in result}
    assert len(result) == 1
    assert "ACDEFGR" in seqs


def test_psimod_validation_drops_mod_after_residue_mutation(psi_db: psimodpy.PsiModDatabase) -> None:
    """MOD:00046 (O-phospho-L-serine, origin=S) at position 1 must be dropped
    when a VariantSimple mutates position 1 from S to W."""
    entry = _make_entry(
        "SWKR",
        mod_res_psi=(pf.ModResPsi(positions=(1,), accession="MOD:00046", name="O-phospho-L-serine"),),
        variant_simple=(pf.VariantSimple(position=1, new_amino_acid="W"),),
    )
    result = list(
        digest_peff_sequence(
            entry,
            _cfg(cleave_on="KR", max_ptm_per_peptide=1, use_psi_mods=True, include_simple_variants=True),
            psi_db=psi_db,
        )
    )
    # Canonical peptide "SWK" should have a phospho version
    assert any("MOD:00046" in str(p.proforma) and p.sequence == "SWK" for p in result)
    # Variant peptide "WWK" must NOT carry MOD:00046 (S→W invalidates phospho-serine)
    assert not any("MOD:00046" in str(p.proforma) and p.sequence == "WWK" for p in result)


def test_psimod_validation_retains_mod_on_unchanged_residue(psi_db: psimodpy.PsiModDatabase) -> None:
    """MOD:00046 at position 2 (an S) must survive when position 1 is mutated."""
    entry = _make_entry(
        "ASKR",
        mod_res_psi=(pf.ModResPsi(positions=(2,), accession="MOD:00046", name="O-phospho-L-serine"),),
        variant_simple=(pf.VariantSimple(position=1, new_amino_acid="W"),),
    )
    result = list(
        digest_peff_sequence(
            entry,
            _cfg(cleave_on="KR", max_ptm_per_peptide=1, use_psi_mods=True, include_simple_variants=True),
            psi_db=psi_db,
        )
    )
    # Both canonical "ASK" and variant "WSK" should have the phospho version
    assert any("MOD:00046" in str(p.proforma) and p.sequence == "ASK" for p in result)
    assert any("MOD:00046" in str(p.proforma) and p.sequence == "WSK" for p in result)


def test_psimod_terminal_promotion_nterm(psi_db: psimodpy.PsiModDatabase) -> None:
    """MOD:00050 (N-acetyl-L-alanine, origin=A, term_spec=N-term) at position 1
    of a protein starting with A should become an nterm mod, not an internal mod."""
    entry = _make_entry(
        "AAKR",
        mod_res_psi=(pf.ModResPsi(positions=(1,), accession="MOD:00050", name="N-acetyl-L-alanine"),),
    )
    result = list(
        digest_peff_sequence(
            entry,
            _cfg(cleave_on="KR", max_ptm_per_peptide=1, use_psi_mods=True),
            psi_db=psi_db,
        )
    )
    seqs = {str(p.proforma) for p in result}
    # The modification must appear in nterm notation (leading bracket before sequence)
    assert any(s.startswith("[MOD:00050]") for s in seqs), f"Expected nterm-promoted MOD:00050 in {seqs}"
    # It must NOT appear as an internal mod at index 0 (which would look like A[MOD:00050]AK)
    assert not any("A[MOD:00050]" in s for s in seqs), f"MOD:00050 should be terminal, not internal, in {seqs}"


def test_psimod_terminal_promotion_skipped_for_internal_peptide(psi_db: psimodpy.PsiModDatabase) -> None:
    """Terminal promotion must NOT fire for a peptide that does not span the protein N-term."""
    entry = _make_entry(
        "AAAKR",
        mod_res_psi=(pf.ModResPsi(positions=(1,), accession="MOD:00050", name="N-acetyl-L-alanine"),),
    )
    result = list(
        digest_peff_sequence(
            entry,
            _cfg(cleave_on="KR", max_ptm_per_peptide=1, use_psi_mods=True),
            psi_db=psi_db,
        )
    )
    r_peptides = [p for p in result if p.sequence == "R"]
    assert all("MOD:00050" not in str(p.proforma) for p in r_peptides)


def test_psimod_terminal_promotion_cterm(psi_db: psimodpy.PsiModDatabase) -> None:
    """MOD:00090 (L-alanine amide, origin=A, term_spec=C-term) at the last position
    of the protein must become a cterm mod on the protein C-terminal peptide."""
    entry = _make_entry(
        "AAKBA",
        mod_res_psi=(pf.ModResPsi(positions=(5,), accession="MOD:00090", name="L-alanine amide"),),
    )
    result = list(
        digest_peff_sequence(
            entry,
            _cfg(cleave_on="K", max_ptm_per_peptide=1, use_psi_mods=True),
            psi_db=psi_db,
        )
    )
    seqs = {str(p.proforma) for p in result}
    # The modification must appear in cterm notation (trailing bracket after sequence)
    assert any(s.endswith("[MOD:00090]") for s in seqs), f"Expected cterm-promoted MOD:00090 in {seqs}"
    # It must NOT appear as an internal mod
    ba_seqs = [s for s in seqs if "BA" in s and "MOD:00090" in s]
    assert all(s.endswith("[MOD:00090]") for s in ba_seqs), (
        f"MOD:00090 should be cterm-promoted, not internal: {ba_seqs}"
    )


def test_two_phase_delta_user_mods_capped_by_peff_mods(psi_db: psimodpy.PsiModDatabase) -> None:
    """With max_ptm_per_peptide=2 and 1 PEFF mod applied, at most 1 user mod can be added."""
    entry = _make_entry(
        "SMKR",
        mod_res_psi=(pf.ModResPsi(positions=(1,), accession="MOD:00046", name="O-phospho-L-serine"),),
    )
    result = list(
        digest_peff_sequence(
            entry,
            _cfg(
                cleave_on="KR",
                max_ptm_per_peptide=2,
                use_psi_mods=True,
                internal_mods=[InternalMod(modification="Oxidation", residue="M", mod_type="variable")],
            ),
            psi_db=psi_db,
        )
    )
    smk = [p for p in result if p.sequence == "SMK"]
    # Expected: unmod, phospho-only, oxidation-only, phospho+oxidation (all ≤2 mods)
    assert any("MOD:00046" in str(p.proforma) and "Oxidation" in str(p.proforma) for p in smk), (
        "Should have a peptide with both PEFF phospho and user Oxidation"
    )
    # No peptide should carry more than 2 total mods
    from peff_digest.digest import _count_mods

    for p in smk:
        assert _count_mods(p.proforma) <= 2, f"Peptide has too many mods: {p.proforma}"


def test_two_phase_delta_full_capacity_no_user_mods(psi_db: psimodpy.PsiModDatabase) -> None:
    """With max_ptm_per_peptide=1 and 1 PEFF mod applied, no user mods can be added."""
    entry = _make_entry(
        "SMKR",
        mod_res_psi=(pf.ModResPsi(positions=(1,), accession="MOD:00046", name="O-phospho-L-serine"),),
    )
    result = list(
        digest_peff_sequence(
            entry,
            _cfg(
                cleave_on="KR",
                max_ptm_per_peptide=1,
                use_psi_mods=True,
                internal_mods=[InternalMod(modification="Oxidation", residue="M", mod_type="variable")],
            ),
            psi_db=psi_db,
        )
    )
    smk_peptides = [p for p in result if p.sequence == "SMK"]
    # Should have: unmod + phospho-only + oxidation-only (never both, capacity=1)
    assert not any("MOD:00046" in str(p.proforma) and "Oxidation" in str(p.proforma) for p in smk_peptides), (
        "Should not have both PEFF and user mod when max_ptm=1"
    )


def test_use_mod_names_outputs_name_instead_of_accession() -> None:
    """use_mod_names=True should embed the mod name, not the accession, in ProForma."""
    entry = _make_entry(
        "ACDEFGR",
        mod_res_psi=(pf.ModResPsi(positions=(3,), accession="MOD:00696", name="phosphorylated residue"),),
    )
    result = list(
        digest_peff_sequence(
            entry,
            _cfg(cleave_on="R", max_ptm_per_peptide=1, use_psi_mods=True, use_mod_names=True),
        )
    )
    seqs = {str(p.proforma) for p in result}
    assert any("phosphorylated residue" in s for s in seqs), "Mod name must appear in output"
    assert not any("MOD:00696" in s for s in seqs), "Accession must not appear when use_mod_names=True"


@pytest.mark.parametrize("use_mod_names", [False, True])
def test_unimod_output_converts_psimod_tags(use_mod_names: bool) -> None:
    """MOD:00696 (xref "Unimod:21#S") converts to UniMod with and without use_mod_names.

    Regression: the M:<name> path passed the raw xref string to get_by_id, got None, and
    silently dropped every modified peptide.
    """
    entry = _make_entry(
        "ACDEFGR",
        mod_res_psi=(pf.ModResPsi(positions=(3,), accession="MOD:00696", name="phosphorylated residue"),),
    )
    cfg = _cfg(
        cleave_on="R",
        max_ptm_per_peptide=1,
        use_psi_mods=True,
        use_mod_names=use_mod_names,
        use_unimod_output=True,
    )
    seqs = {str(p.proforma) for p in digest_peff_sequence(entry, cfg)}
    assert len(seqs) == 2, seqs
    expected = "Phospho" if use_mod_names else "UNIMOD:21"
    assert any(expected in s for s in seqs), seqs
    assert not any("[MOD:" in s or "[M:" in s for s in seqs), seqs


@pytest.mark.parametrize(
    ("xref", "expected"),
    [("Unimod:21", 21), ("Unimod:21#S", 21), ("Unimod:", None), ("Unimod:x", None)],
)
def test_unimod_id_from_xref(xref: str, expected: int | None) -> None:
    from peff_digest.digest import _unimod_id_from_xref

    assert _unimod_id_from_xref(xref) == expected
