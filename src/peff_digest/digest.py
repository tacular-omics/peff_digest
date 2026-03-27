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

from collections.abc import Generator

import pefftacular as pf
import peptacular as pt

from peff_digest.config import InternalMod, TerminalMod

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
    pos_map: dict[int, int | None] = {i + 1: i for i in range(len(sequence))}
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
    use_mod_names: bool = False,
) -> list[tuple[int, str]]:
    """Return (0-based peptide-local position, ProForma tag) for each PEFF mod in the span."""
    result: list[tuple[int, str]] = []
    for mod in mod_entries:
        for peff_pos in mod.positions:
            try:
                peff_pos_int = int(peff_pos)
            except (ValueError, TypeError):
                continue  # skip unknown positions (e.g. '?')
            local0 = variant.pos_map.get(peff_pos_int)
            if local0 is None:
                continue
            pep_local = local0 - span_start
            if 0 <= pep_local < (span_end - span_start):
                tag = mod.name if use_mod_names else mod.accession
                result.append((pep_local, tag))
    return result


_NTERM_POS = -1  # sentinel for N-terminal mods in the applicable list
_CTERM_POS = -2  # sentinel for C-terminal mods in the applicable list


def _apply_mod(ann: pt.ProFormaAnnotation, pos: int, tag: str, seq_len: int) -> None:
    if pos == _NTERM_POS:
        ann.append_nterm_mod(tag)
    elif pos == _CTERM_POS:
        ann.append_cterm_mod(tag)
    else:
        ann.append_internal_mod_at_index(pos, tag)


def _yield_mod_variants(
    base_sequence: str,
    applicable: list[tuple[int, str]],
    max_ptm: int,
) -> list[pt.ProFormaAnnotation]:
    """
    Yield ProFormaAnnotation objects for the base peptide plus every combination
    of 1..max_ptm mods applied to it.

    Sentinel positions: _NTERM_POS (-1) = N-terminal, _CTERM_POS (-2) = C-terminal.
    """
    base = pt.parse(base_sequence)
    variants: list[pt.ProFormaAnnotation] = [base]

    if not applicable or max_ptm <= 0:
        return variants

    seq_len = len(base_sequence)
    for n in range(1, min(max_ptm, len(applicable)) + 1):
        for combo in itertools.combinations(applicable, n):
            positions = [pos for pos, _ in combo]
            if len(positions) != len(set(positions)):
                continue  # multiple mods at the same site — skip this combo
            ann = pt.parse(base_sequence)
            for pos, tag in combo:
                _apply_mod(ann, pos, tag, seq_len)
            variants.append(ann)

    return variants


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


@dataclass
class Peptide:
    proforma: pt.ProFormaAnnotation
    missed_cleavages: int
    semi_enzymatic: bool

    @property
    def sequence(self) -> str:
        return self.proforma.stripped_sequence

    @property
    def mass(self) -> float:
        return self.proforma.neutral_mass()

    @property
    def mod_map(self) -> dict[int, str]:
        _, mod_map = ann_to_map(self.proforma)
        return mod_map


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
    internal_mods: list[InternalMod] | None = None,
    terminal_mods: list[TerminalMod] | None = None,
    annotate_variants: bool = True,
    use_mod_names: bool = False,
    use_psi_mods: bool = True,
    include_simple_variants: bool = True,
    include_complex_variants: bool = True,
) -> Generator[Peptide, None, None]:
    """
    Digest a PEFF SequenceEntry and return all peptide variants as ProFormaAnnotations.

    Each PEFF VariantSimple / VariantComplex is applied independently (not combined).
    PEFF PTMs (ModResPsi) are applied in combinations of up to
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
        internal_mods:       Per-residue modifications. Each InternalMod specifies a
                             residue string, modification name, and mod_type ("fixed"
                             or "variable").

    Returns:
        Set of unique ProFormaAnnotation objects.
    """
    sequence = peff_entry.sequence
    _min = min_length if min_length is not None else 0
    _max = max_length if max_length is not None else len(sequence)

    restrict_after: set[str] = set(restrict_after)
    restrict_before: set[str] = set(restrict_before)
    cleave_on: set[str] = set(cleave_on)

    fixed_mods: dict[str, str] = {}
    variable_mods: dict[str, list[str]] = {}
    for m in internal_mods or []:
        for aa in m.residue:
            if m.mod_type == "fixed":
                fixed_mods[aa] = m.modification
            else:
                variable_mods.setdefault(aa, []).append(m.modification)

    # Collect all PEFF mod annotations (PSI-MOD)
    all_mods: list[pf.ModResPsi] = []
    if use_psi_mods:
        all_mods.extend(peff_entry.mod_res_psi)

    # One variant per PEFF event (canonical + each simple/complex independently)
    variants: list[_Variant] = [_canonical(sequence)]
    if include_simple_variants:
        for v in peff_entry.variant_simple:
            variants.append(_apply_simple(sequence, v))
    if include_complex_variants:
        for v in peff_entry.variant_complex:
            variants.append(_apply_complex(sequence, v))

    for variant in variants:
        vseq = variant.sequence
        seq_len = len(vseq)

        cut_sites = get_cut_sites(
            vseq, cleave_on, restrict_after, restrict_before, cterminal
        )
        n_cuts = len(cut_sites)

        def _process_span(
            start: int,
            end: int,
            mc: int,
            is_semi: bool,
            _vseq: str = vseq,
            _variant: _Variant = variant,
            _seq_len: int = seq_len,
        ) -> list[Peptide]:
            length = end - start
            if length < _min or length > _max:
                return []
            pep_seq = _vseq[start:end]

            # PEFF mods + user variable mods as (0-based peptide position, tag) pairs
            applicable = _mods_in_span(_variant, start, end, all_mods, use_mod_names)
            if variable_mods:
                for aa, tags in variable_mods.items():
                    for i, res in enumerate(pep_seq):
                        if res == aa:
                            for tag in tags:
                                applicable.append((i, tag))
            if terminal_mods:
                for tm in terminal_mods:
                    if tm.position == "nterm":
                        is_term = (start == 0) if tm.protein_terminus else True
                        term_aa = pep_seq[0]
                        sentinel = _NTERM_POS
                    else:
                        is_term = (end == _seq_len) if tm.protein_terminus else True
                        term_aa = pep_seq[-1]
                        sentinel = _CTERM_POS
                    if not is_term:
                        continue
                    if tm.residue is not None and term_aa not in tm.residue:
                        continue
                    if tm.mod_type == "variable":
                        applicable.append((sentinel, tm.modification))

            # Drop any PEFF/variable mods at positions overridden by a fixed mod
            if fixed_mods:
                fixed_positions = {
                    i for i, res in enumerate(pep_seq) if res in fixed_mods
                }
                applicable = [
                    (pos, tag) for pos, tag in applicable if pos not in fixed_positions
                ]

            name = _variant_name(_variant, start, end)
            try:
                mod_variants = _yield_mod_variants(
                    pep_seq, applicable, max_ptm_per_peptide
                )
            except ValueError:
                return []  # unparseable sequence (e.g. contains '*' or spaces) — skip span
            span_results: list[Peptide] = []
            for ann in mod_variants:
                if fixed_mods:
                    for aa, mod_str in fixed_mods.items():
                        for i, res in enumerate(pep_seq):
                            if res == aa:
                                ann.append_internal_mod_at_index(i, mod_str)
                if terminal_mods:
                    for tm in terminal_mods:
                        if tm.mod_type != "fixed":
                            continue
                        if tm.position == "nterm":
                            is_term = (start == 0) if tm.protein_terminus else True
                            term_aa = pep_seq[0]
                        else:
                            is_term = (end == _seq_len) if tm.protein_terminus else True
                            term_aa = pep_seq[-1]
                        if not is_term:
                            continue
                        if tm.residue is not None and term_aa not in tm.residue:
                            continue
                        if tm.position == "nterm":
                            ann.append_nterm_mod(tm.modification)
                        else:
                            ann.append_cterm_mod(tm.modification)
                if annotate_variants and name is not None:
                    ann.peptide_name = name
                span_results.append(
                    Peptide(proforma=ann, missed_cleavages=mc, semi_enzymatic=is_semi)
                )
            return span_results

        # Fully enzymatic peptides
        for i in range(n_cuts - 1):
            for j in range(i + 1, min(i + 2 + missed_cleavages, n_cuts)):
                yield from _process_span(
                    cut_sites[i], cut_sites[j], mc=j - i - 1, is_semi=False
                )

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
                        yield from _process_span(enz_pos, end, mc=mc, is_semi=True)

                # Left-open: non-enzymatic N-term, enzymatic C-term
                for start in range(enz_pos - _min, -1, -1):
                    if _max and (enz_pos - start) > _max:
                        break
                    if start in cut_set:
                        continue
                    mc = enz_idx - bisect_left(cut_sites, start + 1)
                    if mc <= missed_cleavages:
                        yield from _process_span(start, enz_pos, mc=mc, is_semi=True)


def ann_to_map(ann: pt.ProFormaAnnotation) -> tuple[str, dict[int, str]]:

    mod_map: dict[int, str] = {}

    # Process C-terminal modifications
    if ann._cterm_mods is not None:
        for mod_name, count in ann._cterm_mods.items():
            if count != 1:
                raise ValueError("format does not support modification multipliers.")
            mod_map[-2] = mod_name

    # Process N-terminal modifications
    if ann._nterm_mods is not None:
        for mod_name, count in ann._nterm_mods.items():
            if count != 1:
                raise ValueError("format does not support modification multipliers.")
            mod_map[-1] = mod_name

    # Process internal modifications
    if ann._internal_mods is not None:
        for index, mods_dict in ann._internal_mods.items():
            if len(mods_dict) > 1:
                raise ValueError(
                    "format does not support multiple modifications at the same site."
                )
            for mod_name, count in mods_dict.items():
                if count != 1:
                    raise ValueError(
                        "format does not support modification multipliers."
                    )
                mod_map[index] = mod_name

    unmod_sequence = ann.stripped_sequence

    return unmod_sequence, dict(sorted(mod_map.items()))
