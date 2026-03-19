"""
PEFF-aware protein digest.

Produces all peptide ProFormaAnnotation variants from a PEFF SequenceEntry by:
  1. Generating one sequence per variant (canonical + each VariantSimple/VariantComplex
     applied independently — no combinations).
  2. Digesting each variant sequence with the standard enzymatic rules.
  3. For each peptide span, enumerating PEFF-annotated PTM combinations up to
     max_ptm_per_peptide (0 = no PEFF PTMs applied, only unmodified peptides returned).

Skipped: Processed entries (signal peptide / mature-chain trimming).
"""

from __future__ import annotations

import itertools
from bisect import bisect_left, bisect_right
from dataclasses import dataclass

import pefftacular as pf
import peptacular as pt


def get_cut_sites(
    protein_sequence: str,
    cleave_on: set[str],
    restrict_after: set[str],
    restrict_before: set[str],
    cterminal: bool,
) -> list[int]:
    seq_len = len(protein_sequence)
    cut_sites: list[int] = [0]
    for i, aa in enumerate(protein_sequence):
        if aa not in cleave_on:
            continue
        cut_pos = (i + 1) if cterminal else i
        if cut_pos == 0 or cut_pos == seq_len:
            continue
        if restrict_after and protein_sequence[cut_pos] in restrict_after:
            continue
        if restrict_before and protein_sequence[cut_pos - 1] in restrict_before:
            continue
        cut_sites.append(cut_pos)
    cut_sites.append(seq_len)
    return cut_sites


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


@dataclass
class _Variant:
    """A concrete protein sequence derived from one PEFF annotation event."""

    sequence: str
    # Maps 1-based PEFF position → 0-based position in this variant's sequence.
    # None means the position was deleted by a VariantComplex.
    pos_map: dict[int, int | None]
    source: pf.VariantSimple | pf.VariantComplex | None = None
    source_orig: str = ""  # original AA(s) at the variant site


def _canonical(sequence: str) -> _Variant:
    return _Variant(sequence, {i + 1: i for i in range(len(sequence))})


def _apply_simple(sequence: str, v: pf.VariantSimple) -> _Variant:
    """Single AA substitution — sequence length is unchanged."""
    pos0 = int(v.position) - 1  # 0-based
    if pos0 < 0 or pos0 >= len(sequence):
        return _canonical(sequence)  # out-of-range annotation → skip
    new_seq = sequence[:pos0] + v.new_amino_acid + sequence[pos0 + 1 :]
    pos_map = {i + 1: i for i in range(len(sequence))}  # identity
    return _Variant(new_seq, pos_map, source=v, source_orig=sequence[pos0])


def _apply_complex(sequence: str, v: pf.VariantComplex) -> _Variant:
    """Multi-residue variant — sequence length may change."""
    start0 = int(v.start_pos) - 1  # 0-based inclusive
    end0 = int(v.end_pos)  # 0-based exclusive  (PEFF end is 1-based inclusive)
    if start0 < 0 or end0 > len(sequence) or start0 > end0:
        return _canonical(sequence)  # out-of-range annotation → skip

    orig_seq = sequence[start0:end0]
    new_seq = sequence[:start0] + v.new_sequence + sequence[end0:]
    delta = len(v.new_sequence) - (end0 - start0)

    pos_map: dict[int, int | None] = {}
    for peff_pos in range(1, len(sequence) + 1):
        orig0 = peff_pos - 1
        if orig0 < start0:
            pos_map[peff_pos] = orig0
        elif orig0 < end0:
            pos_map[peff_pos] = None  # consumed by the variant
        else:
            pos_map[peff_pos] = orig0 + delta

    return _Variant(new_seq, pos_map, source=v, source_orig=orig_seq)


def _variant_name(variant: _Variant, pep_start: int, pep_end: int) -> str | None:
    """Return a PEFF-notation variant description, or None for canonical."""
    if variant.source is None:
        return None
    if isinstance(variant.source, pf.VariantSimple):
        v = variant.source
        return f"({v.position}|{v.new_amino_acid})"
    v = variant.source
    return f"({v.start_pos}|{v.end_pos}|{v.new_sequence})"


def _mods_in_span(
    variant: _Variant,
    span_start: int,  # 0-based in variant.sequence
    span_end: int,  # 0-based exclusive
    mod_entries: list[pf.ModResPsi | pf.ModResUnimod],
) -> list[tuple[int, str]]:
    """Return (0-based peptide-local position, ProForma tag) for each PEFF mod in the span."""
    result: list[tuple[int, str]] = []
    for mod in mod_entries:
        for peff_pos in mod.positions:
            local0 = variant.pos_map.get(int(peff_pos))
            if local0 is None:
                continue
            pep_local = local0 - span_start
            if 0 <= pep_local < (span_end - span_start):
                tag = mod.accession
                result.append((pep_local, tag))
    return result


def _yield_mod_variants(
    base_sequence: str,
    applicable: list[tuple[int, str]],
    max_ptm: int,
) -> list[pt.ProFormaAnnotation]:
    """
    Yield ProFormaAnnotation objects for the base peptide plus every combination
    of 1..max_ptm PEFF mods applied to it.
    """
    base = pt.parse(base_sequence)
    variants: list[pt.ProFormaAnnotation] = [base]

    if not applicable or max_ptm <= 0:
        return variants

    for n in range(1, min(max_ptm, len(applicable)) + 1):
        for combo in itertools.combinations(applicable, n):
            positions = [pos for pos, _ in combo]
            if len(positions) != len(set(positions)):
                continue  # multiple mods at the same site — skip this combo
            ann = pt.parse(base_sequence)
            for pos, tag in combo:
                ann.append_internal_mod_at_index(pos, tag)
            variants.append(ann)

    return variants


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def digest_peff_sequence(
    peff_entry: pf.SequenceEntry,
    cleave_on: str,
    min_length: int | None = None,
    max_length: int | None = None,
    restrict_after: str = "",
    restrict_before: str = "",
    cterminal: bool = True,
    missed_cleavages: int = 0,
    semi_enzymatic: bool = False,
    max_ptm_per_peptide: int = 2,
    fixed_mods: dict[str, str] | None = None,
    variable_mods: dict[str, list[str]] | None = None,
) -> set[pt.ProFormaAnnotation]:
    """
    Digest a PEFF SequenceEntry and return all peptide variants as ProFormaAnnotations.

    Each PEFF VariantSimple / VariantComplex is applied independently (not combined).
    PEFF PTMs (ModResPsi / ModResUnimod) are applied in combinations of up to
    max_ptm_per_peptide per peptide.  Pass 0 to skip PEFF PTMs entirely.

    Args:
        peff_entry:          Parsed PEFF SequenceEntry.
        cleave_on:           Set of amino acids to cleave on (e.g. {"K", "R"}).
        min_length:          Minimum peptide length (inclusive), or None.
        max_length:          Maximum peptide length (inclusive), or None.
        restrict_after:      Do not cleave when the next AA is in this set.
        restrict_before:     Do not cleave when the preceding AA is in this set.
        cterminal:           True = C-terminal cleavage (standard); False = N-terminal.
        missed_cleavages:    Maximum number of missed cleavage sites allowed.
        semi_enzymatic:      Include semi-enzymatic peptides (one non-enzymatic end).
        max_ptm_per_peptide: Max number of PEFF PTM annotations to apply simultaneously.
        fixed_mods:          Static modifications applied to all matching residues,
                             e.g. {"C": "UNIMOD:4"}.
        variable_mods:       Variable modifications enumerated combinatorially,
                             e.g. {"S": ["UNIMOD:21"], "T": ["UNIMOD:21"]}.

    Returns:
        Set of unique ProFormaAnnotation objects.
    """
    sequence = peff_entry.sequence
    _min = min_length if min_length is not None else 0
    _max = max_length if max_length is not None else len(sequence)

    restrict_after: set[str] = set(restrict_after)
    restrict_before: set[str] = set(restrict_before)
    cleave_on: set[str] = set(cleave_on)

    # Collect all PEFF mod annotations (PSI-MOD + Unimod)
    all_mods: list[pf.ModResPsi | pf.ModResUnimod] = [
        *peff_entry.mod_res_psi,
        *peff_entry.mod_res_unimod,
    ]

    # One variant per PEFF event (canonical + each simple/complex independently)
    variants: list[_Variant] = [_canonical(sequence)]
    for v in peff_entry.variant_simple:
        variants.append(_apply_simple(sequence, v))
    for v in peff_entry.variant_complex:
        variants.append(_apply_complex(sequence, v))

    seen: set[pt.ProFormaAnnotation] = set()

    for variant in variants:
        vseq = variant.sequence
        seq_len = len(vseq)

        cut_sites = get_cut_sites(vseq, cleave_on, restrict_after, restrict_before, cterminal)
        n_cuts = len(cut_sites)

        def _process_span(
            start: int,
            end: int,
            _vseq: str = vseq,
            _variant: _Variant = variant,
        ) -> None:
            length = end - start
            if length < _min or length > _max:
                return
            pep_seq = _vseq[start:end]

            # PEFF mods + user variable mods as (0-based peptide position, tag) pairs
            applicable = _mods_in_span(_variant, start, end, all_mods)
            if variable_mods:
                for aa, tags in variable_mods.items():
                    for i, res in enumerate(pep_seq):
                        if res == aa:
                            for tag in tags:
                                applicable.append((i, tag))

            # Drop any PEFF/variable mods at positions overridden by a fixed mod
            if fixed_mods:
                fixed_positions = {i for i, res in enumerate(pep_seq) if res in fixed_mods}
                applicable = [(pos, tag) for pos, tag in applicable if pos not in fixed_positions]

            name = _variant_name(_variant, start, end)
            for ann in _yield_mod_variants(pep_seq, applicable, max_ptm_per_peptide):
                if fixed_mods:
                    for aa, mod_str in fixed_mods.items():
                        for i, res in enumerate(pep_seq):
                            if res == aa:
                                ann.append_internal_mod_at_index(i, mod_str)
                if name is not None:
                    ann.peptide_name = name
                seen.add(ann)

        # Fully enzymatic peptides
        for i in range(n_cuts - 1):
            for j in range(i + 1, min(i + 2 + missed_cleavages, n_cuts)):
                _process_span(cut_sites[i], cut_sites[j])

        # Semi-enzymatic peptides (one free end)
        if semi_enzymatic:
            cut_set = set(cut_sites)
            for enz_idx, enz_pos in enumerate(cut_sites):
                # Right-open: enzymatic N-term, non-enzymatic C-term
                for end in range(enz_pos + _min, seq_len + 1):
                    if _max and (end - enz_pos) > _max:
                        break
                    if end in cut_set:
                        continue
                    mc = bisect_right(cut_sites, end - 1) - enz_idx - 1
                    if mc <= missed_cleavages:
                        _process_span(enz_pos, end)

                # Left-open: non-enzymatic N-term, enzymatic C-term
                for start in range(enz_pos - _min, -1, -1):
                    if _max and (enz_pos - start) > _max:
                        break
                    if start in cut_set:
                        continue
                    mc = enz_idx - bisect_left(cut_sites, start + 1)
                    if mc <= missed_cleavages:
                        _process_span(start, enz_pos)

    return seen
