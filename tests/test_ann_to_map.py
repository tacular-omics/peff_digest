from __future__ import annotations

import peptacular as pt
import pytest
from conftest import _cfg, _make_entry

from peff_digest import InternalMod, digest_peff_sequence
from peff_digest.config import TerminalMod
from peff_digest.digest import ann_to_map


def test_ann_to_map_internal_mod() -> None:
    """ann_to_map must return the correct position→name mapping for an internal mod."""
    entry = _make_entry("ACKR")
    result = list(
        digest_peff_sequence(
            entry,
            _cfg(
                cleave_on="KR",
                internal_mods=[InternalMod(modification="Carbamidomethyl", residue="C", mod_type="fixed")],
            ),
        )
    )
    ack_peptides = [p for p in result if p.sequence == "ACK"]
    assert len(ack_peptides) == 1
    mod_map = ack_peptides[0].mod_map
    # C is at 0-based index 1 in "ACK"
    assert mod_map == {1: "Carbamidomethyl"}, f"Unexpected mod_map: {mod_map}"


def test_ann_to_map_nterm_sentinel() -> None:
    """N-terminal mods must map to key -1 in mod_map."""
    entry = _make_entry("AAKR")
    result = list(
        digest_peff_sequence(
            entry,
            _cfg(
                cleave_on="KR", terminal_mods=[TerminalMod(modification="Acetyl", position="nterm", mod_type="fixed")]
            ),
        )
    )
    aak = next(p for p in result if p.sequence == "AAK")
    assert aak.mod_map.get(-1) == "Acetyl", f"Expected -1: Acetyl in mod_map, got {aak.mod_map}"


def test_ann_to_map_cterm_sentinel() -> None:
    """C-terminal mods must map to key -2 in mod_map."""
    ann = pt.parse("AAKR")
    ann.append_cterm_mod("Amidated")
    _, mod_map = ann_to_map(ann)
    assert mod_map.get(-2) == "Amidated"


def test_ann_to_map_raises_on_cterm_multiplier() -> None:
    """A c-terminal mod with count > 1 must raise ValueError."""
    ann = pt.parse("A")
    ann.append_cterm_mod("Amide")
    ann.append_cterm_mod("Amide")  # multiplier = 2

    with pytest.raises(ValueError, match="modification multipliers"):
        ann_to_map(ann)


def test_ann_to_map_raises_on_nterm_multiplier() -> None:
    """An n-terminal mod with count > 1 must raise ValueError."""
    ann = pt.parse("A")
    ann.append_nterm_mod("Acetyl")
    ann.append_nterm_mod("Acetyl")  # multiplier = 2

    with pytest.raises(ValueError, match="modification multipliers"):
        ann_to_map(ann)


def test_ann_to_map_raises_on_multiple_mods_at_same_internal_site() -> None:
    """Two distinct mods at the same internal position must raise ValueError."""
    ann = pt.parse("ACK")
    ann.append_internal_mod_at_index(1, "Carbamidomethyl")
    ann.append_internal_mod_at_index(1, "Oxidation")  # two mods at position 1

    with pytest.raises(ValueError, match="multiple modifications at the same site"):
        ann_to_map(ann)


def test_ann_to_map_raises_on_internal_mod_multiplier() -> None:
    """An internal mod applied twice at the same position must raise ValueError."""
    ann = pt.parse("ACK")
    ann.append_internal_mod_at_index(1, "Carbamidomethyl")
    ann.append_internal_mod_at_index(1, "Carbamidomethyl")  # multiplier = 2

    with pytest.raises(ValueError, match="modification multipliers"):
        ann_to_map(ann)
