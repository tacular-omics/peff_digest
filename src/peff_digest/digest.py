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

import copy
import itertools
import logging
from bisect import bisect_left, bisect_right
from collections.abc import Generator
from dataclasses import dataclass

import pefftacular as pf
import peptacular as pt
import psimodpy
from psimodpy import AminoAcid, PsiModDatabase, TermSpec

from peff_digest.config import DigestConfig

logger = logging.getLogger(__name__)


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


def _variant_in_span(variant: _Variant, pep_start: int, pep_end: int) -> bool:
    """Return True if the variant's mutation site overlaps the peptide span [pep_start, pep_end)."""
    if variant.source is None:
        return False
    if isinstance(variant.source, pf.VariantSimple):
        mapped = variant.pos_map.get(int(variant.source.position))
        if mapped is None:
            return False
        return pep_start <= mapped < pep_end
    # VariantComplex: the new sequence occupies [start0, start0+len(new_sequence))
    v = variant.source
    start0 = int(v.start_pos) - 1
    new_len = len(v.new_sequence)
    # After complex variant, the inserted region sits at [start0, start0+new_len)
    # in the variant sequence. Check overlap with peptide span.
    var_end = start0 + new_len
    return pep_start < var_end and start0 < pep_end


def _mods_in_span(
    variant: _Variant,
    span_start: int,  # 0-based in variant.sequence
    span_end: int,  # 0-based exclusive
    mod_entries: list[pf.ModResPsi | pf.ModResUnimod],
    use_mod_names: bool = False,
    psi_db: PsiModDatabase | None = None,
    protein_len: int = 0,
) -> list[tuple[int, str]]:
    """Return (0-based peptide-local position, ProForma tag) for each PEFF mod in the span.

    When *psi_db* is provided each mod is looked up to:
    - Validate the origin residue against the (possibly variant-mutated) sequence.
      Mods whose origin does not match the current residue are silently dropped.
    - Promote terminus-specific mods (TermSpec.N_TERM / C_TERM) to the appropriate
      terminal sentinel position (_NTERM_POS / _CTERM_POS) when the mod sits at the
      protein terminus and the peptide spans that terminus.
    """
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
            if not (0 <= pep_local < (span_end - span_start)):
                continue
            tag = f"M:{mod.name}" if use_mod_names else mod.accession

            if psi_db is not None:
                entry = psi_db.get_by_id(mod.accession)
                if entry is not None:
                    # Validate: origin residue must match the (possibly mutated) residue.
                    origin = entry.origin
                    if isinstance(origin, AminoAcid) and origin != AminoAcid.ANY:
                        if variant.sequence[local0] != str(origin):
                            continue  # residue changed by variant — drop this mod

                    # Terminal promotion: convert positional mod to terminal sentinel.
                    # If a mod is terminal-specific but the conditions aren't met, drop it —
                    # it cannot be validly applied as an internal modification.
                    if entry.term_spec == TermSpec.N_TERM:
                        if local0 == 0 and span_start == 0:
                            result.append((_NTERM_POS, tag))
                        continue
                    if entry.term_spec == TermSpec.C_TERM:
                        if local0 == protein_len - 1 and span_end == protein_len:
                            result.append((_CTERM_POS, tag))
                        continue

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


def _count_mods(ann: pt.ProFormaAnnotation) -> int:
    """Return the total number of modifications on an annotation."""
    count = 0
    if ann._nterm_mods:
        count += sum(ann._nterm_mods.values())
    if ann._cterm_mods:
        count += sum(ann._cterm_mods.values())
    if ann._internal_mods:
        for mods_dict in ann._internal_mods.values():
            count += sum(mods_dict.values())
    return count


def _get_occupied_positions(ann: pt.ProFormaAnnotation) -> set[int]:
    """Return the set of positions (including terminal sentinels) already carrying a mod."""
    positions: set[int] = set()
    if ann._nterm_mods:
        positions.add(_NTERM_POS)
    if ann._cterm_mods:
        positions.add(_CTERM_POS)
    if ann._internal_mods:
        positions.update(ann._internal_mods.keys())
    return positions


def _yield_user_mod_variants(
    base_ann: pt.ProFormaAnnotation,
    user_applicable: list[tuple[int, str]],
    remaining: int,
) -> list[pt.ProFormaAnnotation]:
    """Return annotation variants produced by layering user mods on top of *base_ann*.

    Works like :func:`_yield_mod_variants` but starts from an existing annotation
    (which may already carry PEFF mods) rather than a bare sequence string.
    Uses :func:`copy.deepcopy` to avoid mutating the shared base.
    """
    variants: list[pt.ProFormaAnnotation] = [copy.deepcopy(base_ann)]
    if not user_applicable or remaining <= 0:
        return variants

    seq_len = len(base_ann.stripped_sequence)
    for n in range(1, min(remaining, len(user_applicable)) + 1):
        for combo in itertools.combinations(user_applicable, n):
            positions = [pos for pos, _ in combo]
            if len(positions) != len(set(positions)):
                continue  # multiple mods at the same site — skip
            ann = copy.deepcopy(base_ann)
            for pos, tag in combo:
                _apply_mod(ann, pos, tag, seq_len)
            variants.append(ann)
    return variants


# ---------------------------------------------------------------------------
# PSI-MOD → UniMod conversion
# ---------------------------------------------------------------------------


def _try_convert_psimod_to_unimod(
    ann: pt.ProFormaAnnotation,
    psi_db: PsiModDatabase,
    uni_db,
) -> pt.ProFormaAnnotation | None:
    """Replace MOD:/M: tags with UNIMOD:/U: equivalents.

    Returns None if any PSI-MOD mod has no UniMod xref — the caller should drop the peptide.
    Non-PSI-MOD tags are passed through unchanged.
    """
    sequence, mod_map = ann_to_map(ann)
    new_mod_map: dict[int, str] = {}
    for pos, tag in mod_map.items():
        if tag.startswith("MOD:"):
            psi_entry = psi_db.get_by_id(tag)
            if psi_entry is None or not psi_entry.xref_unimod:
                return None
            try:
                unimod_id = int(psi_entry.xref_unimod.replace("Unimod:", "").split("#")[0])
            except ValueError:
                return None
            if uni_db.get_by_id(unimod_id) is None:
                return None
            new_mod_map[pos] = f"UNIMOD:{unimod_id}"
        elif tag.startswith("M:"):
            name = tag[2:]
            psi_entry = psi_db.get_by_name(name)
            if psi_entry is None or not psi_entry.xref_unimod:
                return None
            uni_entry = uni_db.get_by_id(psi_entry.xref_unimod)
            if uni_entry is None:
                return None
            new_mod_map[pos] = f"U:{uni_entry.name}"
        else:
            new_mod_map[pos] = tag
    new_ann = pt.parse(sequence)
    seq_len = len(sequence)
    for pos, tag in new_mod_map.items():
        _apply_mod(new_ann, pos, tag, seq_len)
    return new_ann


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


@dataclass
class Peptide:
    proforma: pt.ProFormaAnnotation
    missed_cleavages: int
    semi_enzymatic: bool
    is_protein_nterm: bool = False
    is_protein_cterm: bool = False
    variant: pf.VariantSimple | pf.VariantComplex | None = None

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
    config: DigestConfig,
    psi_db: PsiModDatabase | None = None,
    uni_db=None,
) -> Generator[Peptide, None, None]:
    """
    Digest a PEFF SequenceEntry and return all peptide variants as Peptide objects.

    Each PEFF VariantSimple / VariantComplex is applied independently (not combined).
    PEFF PTMs (ModResPsi) are applied in combinations of up to
    ``config.max_ptm_per_peptide`` per peptide.  Pass 0 to skip PEFF PTMs entirely.

    When ``config.use_unimod_output`` is True, PSI-MOD tags are converted to UniMod;
    peptides with no UniMod xref are dropped.  Mass filters (``min_mass``,
    ``max_mass``, ``drop_invalid_mass``) are applied before yielding.
    """

    if psi_db is None:
        psi_db = psimodpy.load()
    if config.use_unimod_output and uni_db is None:
        import unimodpy as _unimodpy

        uni_db = _unimodpy.load()

    sequence = peff_entry.sequence
    _min = config.min_length if config.min_length is not None else 0
    _max = config.max_length if config.max_length is not None else len(sequence)

    restrict_after: set[str] = set(config.restrict_after)
    restrict_before: set[str] = set(config.restrict_before)
    cleave_on: set[str] = set(config.cleave_on)

    internal_mods = config.internal_mods or []
    terminal_mods = config.terminal_mods or []

    fixed_mods: dict[str, str] = {}
    variable_mods: dict[str, list[str]] = {}
    for m in internal_mods:
        for aa in m.residue:
            if m.mod_type == "fixed":
                fixed_mods[aa] = m.modification
            else:
                variable_mods.setdefault(aa, []).append(m.modification)

    # Collect all PEFF mod annotations (PSI-MOD)
    all_mods: list[pf.ModResPsi | pf.ModResUnimod] = []
    if config.use_psi_mods:
        all_mods.extend(peff_entry.mod_res_psi)

    # One variant per PEFF event (canonical + each simple/complex independently)
    variants: list[_Variant] = [_canonical(sequence)]
    if config.include_simple_variants:
        for v in peff_entry.variant_simple:
            variants.append(_apply_simple(sequence, v))
    if config.include_complex_variants:
        for v in peff_entry.variant_complex:
            variants.append(_apply_complex(sequence, v))

    max_ptm_per_peptide = config.max_ptm_per_peptide
    use_mod_names = config.use_mod_names
    semi_enzymatic = config.semi_enzymatic
    missed_cleavages = config.missed_cleavages

    seen: set[str] = set()

    for variant in variants:
        vseq = variant.sequence
        seq_len = len(vseq)

        cut_sites = get_cut_sites(vseq, cleave_on, restrict_after, restrict_before, config.cterminal)
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
            is_protein_nterm = start == 0
            is_protein_cterm = end == _seq_len

            # ----------------------------------------------------------------
            # Phase 1: PEFF mods only
            # When psi_db is provided, _mods_in_span validates each mod against
            # the (possibly variant-mutated) residue and promotes terminus-specific
            # mods to terminal sentinels.
            # ----------------------------------------------------------------
            peff_applicable = _mods_in_span(_variant, start, end, all_mods, use_mod_names, psi_db, _seq_len)

            # Compute fixed-mod positions; when override mode is on, drop PEFF mods there.
            fixed_positions: set[int] = set()
            if fixed_mods:
                fixed_positions = {i for i, res in enumerate(pep_seq) if res in fixed_mods}
                if config.fixed_mod_overrides_peff:
                    peff_applicable = [(pos, tag) for pos, tag in peff_applicable if pos not in fixed_positions]

            span_variant = _variant.source if _variant_in_span(_variant, start, end) else None
            try:
                peff_variants = _yield_mod_variants(pep_seq, peff_applicable, max_ptm_per_peptide)
            except ValueError:
                return []  # unparseable sequence (e.g. contains '*' or spaces) — skip span

            # ----------------------------------------------------------------
            # Build user variable mod candidates once for this span.
            # These are filtered per-peff-variant below to exclude occupied positions.
            # ----------------------------------------------------------------
            user_variable: list[tuple[int, str]] = []
            if variable_mods:
                for aa, tags in variable_mods.items():
                    for i, res in enumerate(pep_seq):
                        if res == aa:
                            for tag in tags:
                                user_variable.append((i, tag))
            if terminal_mods:
                for tm in terminal_mods:
                    if tm.mod_type != "variable":
                        continue
                    if tm.position == "nterm":
                        is_term = is_protein_nterm if tm.protein_terminus else True
                        term_aa = pep_seq[0]
                        sentinel = _NTERM_POS
                    else:
                        is_term = is_protein_cterm if tm.protein_terminus else True
                        term_aa = pep_seq[-1]
                        sentinel = _CTERM_POS
                    if not is_term:
                        continue
                    if tm.residue is not None and term_aa not in tm.residue:
                        continue
                    user_variable.append((sentinel, tm.modification))

            # ----------------------------------------------------------------
            # Phase 2: for each validated PEFF variant, apply user mods.
            # ----------------------------------------------------------------
            span_results: list[Peptide] = []
            for peff_ann in peff_variants:
                existing = _count_mods(peff_ann)
                remaining = max_ptm_per_peptide - existing
                occupied = _get_occupied_positions(peff_ann)

                # Exclude positions already carrying a PEFF mod or reserved by fixed mods.
                filtered_user = [
                    (pos, tag)
                    for pos, tag in user_variable
                    if pos not in occupied and (not fixed_mods or pos not in fixed_positions)
                ]

                user_variants = _yield_user_mod_variants(peff_ann, filtered_user, remaining)

                for ann in user_variants:
                    # Apply fixed internal mods (skip occupied positions when override is off).
                    if fixed_mods:
                        for aa, mod_str in fixed_mods.items():
                            for i, res in enumerate(pep_seq):
                                if res == aa:
                                    if not config.fixed_mod_overrides_peff and i in occupied:
                                        continue
                                    ann.append_internal_mod_at_index(i, mod_str)
                    # Apply fixed terminal mods.
                    if terminal_mods:
                        for tm in terminal_mods:
                            if tm.mod_type != "fixed":
                                continue
                            if tm.position == "nterm":
                                is_term = is_protein_nterm if tm.protein_terminus else True
                                term_aa = pep_seq[0]
                            else:
                                is_term = is_protein_cterm if tm.protein_terminus else True
                                term_aa = pep_seq[-1]
                            if not is_term:
                                continue
                            if tm.residue is not None and term_aa not in tm.residue:
                                continue
                            if tm.position == "nterm":
                                if not config.fixed_mod_overrides_peff and _NTERM_POS in occupied:
                                    continue
                                ann.append_nterm_mod(tm.modification)
                            else:
                                if not config.fixed_mod_overrides_peff and _CTERM_POS in occupied:
                                    continue
                                ann.append_cterm_mod(tm.modification)
                    span_results.append(
                        Peptide(
                            proforma=ann,
                            missed_cleavages=mc,
                            semi_enzymatic=is_semi,
                            is_protein_nterm=is_protein_nterm,
                            is_protein_cterm=is_protein_cterm,
                            variant=span_variant,
                        )
                    )
            return span_results

        def _yield_deduped(peptides: list[Peptide]) -> Generator[Peptide, None, None]:
            for peptide in peptides:
                key = str(peptide.proforma)
                if key not in seen:
                    seen.add(key)
                    yield peptide

        def _postprocess(peptides: list[Peptide]) -> Generator[Peptide, None, None]:
            for peptide in _yield_deduped(peptides):
                ann = peptide.proforma
                if config.use_unimod_output:
                    ann = _try_convert_psimod_to_unimod(ann, psi_db, uni_db)
                    if ann is None:
                        continue
                    peptide = Peptide(
                        proforma=ann,
                        missed_cleavages=peptide.missed_cleavages,
                        semi_enzymatic=peptide.semi_enzymatic,
                        is_protein_nterm=peptide.is_protein_nterm,
                        is_protein_cterm=peptide.is_protein_cterm,
                        variant=peptide.variant,
                    )
                if config.min_mass is not None or config.max_mass is not None or config.drop_invalid_mass:
                    try:
                        mass = None if "X" in ann.stripped_sequence else ann.mass()
                    except Exception:
                        mass = None
                    if mass is None:
                        logger.debug("Could not compute mass for peptide %s", str(ann))
                    if mass is None and config.drop_invalid_mass:
                        continue
                    if mass is not None and config.min_mass is not None and mass < config.min_mass:
                        continue
                    if mass is not None and config.max_mass is not None and mass > config.max_mass:
                        continue
                yield peptide

        # Fully enzymatic peptides
        for i in range(n_cuts - 1):
            for j in range(i + 1, min(i + 2 + missed_cleavages, n_cuts)):
                yield from _postprocess(_process_span(cut_sites[i], cut_sites[j], mc=j - i - 1, is_semi=False))

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
                        yield from _postprocess(_process_span(enz_pos, end, mc=mc, is_semi=True))

                # Left-open: non-enzymatic N-term, enzymatic C-term
                for start in range(enz_pos - _min, -1, -1):
                    if _max and (enz_pos - start) > _max:
                        break
                    if start in cut_set:
                        continue
                    mc = enz_idx - bisect_left(cut_sites, start + 1)
                    if mc <= missed_cleavages:
                        yield from _postprocess(_process_span(start, enz_pos, mc=mc, is_semi=True))


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
                raise ValueError("format does not support multiple modifications at the same site.")
            for mod_name, count in mods_dict.items():
                if count != 1:
                    raise ValueError("format does not support modification multipliers.")
                mod_map[index] = mod_name

    unmod_sequence = ann.stripped_sequence

    return unmod_sequence, dict(sorted(mod_map.items()))
