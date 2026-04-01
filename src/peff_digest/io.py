"""File I/O, database caching, and multiprocessing worker functions."""

from __future__ import annotations

import logging
from collections import Counter
from pathlib import Path

import pefftacular as pf
import peptacular as pt
from psimodpy import PsiModDatabase

from peff_digest.config import DigestConfig
from peff_digest.digest import ann_to_map, digest_peff_sequence

logger = logging.getLogger(__name__)

AA_RESIDUE_MASSES: dict[str, float] = {
    aa: m for aa in "ACDEFGHIKLMNPQRSTVWY" if (m := pt.AA_LOOKUP[aa].monoisotopic_mass) is not None
}


def _proforma_to_mass_array(ann: pt.ProFormaAnnotation, mod_mass_cache: dict[str, float]) -> list[float]:
    """Build a per-residue mass array from a ProFormaAnnotation."""
    masses = [AA_RESIDUE_MASSES[aa] for aa in ann.stripped_sequence]
    if ann._internal_mods is not None:
        for pos, mod_dict in ann._internal_mods.items():
            for mod, cnt in mod_dict.items():
                if mod not in mod_mass_cache:
                    mod_mass_cache[mod] = pt.ModificationTags.from_string(mod).get_mass(monoisotopic=True)
                masses[pos] += mod_mass_cache[mod] * cnt
    if ann._nterm_mods is not None:
        for mod, cnt in ann._nterm_mods.items():
            if mod not in mod_mass_cache:
                mod_mass_cache[mod] = pt.ModificationTags.from_string(mod).get_mass(monoisotopic=True)
            masses[0] += mod_mass_cache[mod] * cnt
    if ann._cterm_mods is not None:
        for mod, cnt in ann._cterm_mods.items():
            if mod not in mod_mass_cache:
                mod_mass_cache[mod] = pt.ModificationTags.from_string(mod).get_mass(monoisotopic=True)
            masses[-1] += mod_mass_cache[mod] * cnt
    return masses


# ---------------------------------------------------------------------------
# Lazy database loaders (cached per process for multiprocessing)
# ---------------------------------------------------------------------------

_PSI_DB: PsiModDatabase | None = None
_UNI_DB = None


def _get_psi_db() -> PsiModDatabase:
    """Load the bundled PSI-MOD database once per process and cache it."""
    global _PSI_DB
    if _PSI_DB is None:
        import psimodpy as _psimodpy

        _PSI_DB = _psimodpy.load()
    return _PSI_DB


def _get_uni_db():
    """Load the bundled UniMod database once per process and cache it."""
    global _UNI_DB
    if _UNI_DB is None:
        import unimodpy as _unimodpy

        _UNI_DB = _unimodpy.load()
    return _UNI_DB


# ---------------------------------------------------------------------------
# File reading (PEFF and FASTA)
# ---------------------------------------------------------------------------

_FASTA_EXTENSIONS = {".fasta", ".fa", ".faa", ".fas"}


def _iter_fasta(path: str):
    """Yield (header, sequence) pairs from a FASTA file."""
    header = None
    chunks: list[str] = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:]
                chunks = []
            elif header is not None:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def _fasta_to_entry(header: str, sequence: str) -> pf.SequenceEntry:
    """Convert a FASTA record to a minimal SequenceEntry."""
    first_token = header.split()[0]
    parts = first_token.split("|")
    if len(parts) >= 3:
        prefix, db_unique_id = parts[0], parts[1]
    else:
        prefix, db_unique_id = "", first_token
    return pf.SequenceEntry(prefix=prefix, db_unique_id=db_unique_id, sequence=sequence)


def read_sequences(path: str) -> tuple[list[pf.SequenceEntry], int]:
    """Read sequences from a PEFF or FASTA file. Returns (entries, n_malformed)."""
    sequences: list[pf.SequenceEntry] = []
    n_malformed = 0
    if Path(path).suffix.lower() in _FASTA_EXTENSIONS:
        for header, seq in _iter_fasta(path):
            try:
                sequences.append(_fasta_to_entry(header, seq))
            except Exception:
                n_malformed += 1
                logger.warning("Skipping malformed FASTA entry: %s", header.split()[0] if header else "<unknown>")
    else:
        reader = iter(pf.PeffReader(path))
        while True:
            try:
                sequences.append(next(reader))
            except StopIteration:
                break
            except Exception:
                n_malformed += 1
                logger.warning("Skipping malformed PEFF entry at position %d", len(sequences) + n_malformed)
    return sequences, n_malformed


# ---------------------------------------------------------------------------
# Multiprocessing worker functions
# ---------------------------------------------------------------------------

_Row = tuple[str, str, str, str | None, int, float | None, int, bool, list[float] | None, int]


def _dedup_rows(rows: list[tuple]) -> list[_Row]:
    """Deduplicate rows by peptidoform_id, keeping the best per unique form.

    Selection: lowest missed_cleavages wins; ties broken by preferring
    non-semi-enzymatic (False < True). First-appearance order is preserved.
    The final field n_peptidoforms records how many raw rows shared the key.
    """
    # Dedup on the full ProForma sequence string so positional mod variants are kept
    # as separate rows.  Only exact duplicates (same bare seq + same mod positions)
    # are collapsed, keeping the one with lowest missed_cleavages (non-semi preferred).
    seq_best: dict[str, tuple] = {}
    seq_order: list[str] = []
    for row in rows:
        seq = row[1]
        if seq not in seq_best:
            seq_best[seq] = row
            seq_order.append(seq)
        elif (row[6], row[7]) < (seq_best[seq][6], seq_best[seq][7]):
            seq_best[seq] = row

    # n_peptidoforms: count of distinct ProForma sequences per peptidoform_id group
    pform_count: dict[str, int] = {}
    for row in seq_best.values():
        pform_count[row[2]] = pform_count.get(row[2], 0) + 1

    return [(*seq_best[s], pform_count[seq_best[s][2]]) for s in seq_order]


def _format_protein_id(entry: pf.SequenceEntry) -> str:
    """Return the fullest available protein identifier (e.g. sp|Q9Y2X3|NOP5_HUMAN)."""
    parts = [p for p in (entry.prefix, entry.db_unique_id, entry.id) if p]
    return "|".join(parts) if len(parts) > 1 else entry.db_unique_id


def _format_variant(variant: pf.VariantSimple | pf.VariantComplex | None) -> str | None:
    """Format a variant object as a PEFF-notation string for output."""
    if variant is None:
        return None
    if isinstance(variant, pf.VariantSimple):
        return f"({variant.position}|{variant.new_amino_acid})"
    return f"({variant.start_pos}|{variant.end_pos}|{variant.new_sequence})"


def _digest_worker(
    sequence: pf.SequenceEntry,
    config: DigestConfig,
) -> list[_Row]:
    psi_db = _get_psi_db() if config.use_psi_mods else None
    uni_db = _get_uni_db() if config.use_unimod_output else None
    protein_id = _format_protein_id(sequence)
    _mod_mass_cache: dict[str, float] = {}
    rows = []
    for peptide in digest_peff_sequence(sequence, config, psi_db=psi_db, uni_db=uni_db):
        ann = peptide.proforma
        bare_seq, mod_map = ann_to_map(ann)
        mod_counter = Counter()
        for pos in mod_map:
            mod_counter[mod_map[pos]] += 1

        mod_str = "".join(f"[{k}]^{v}" if v != 1 else f"[{k}]" for k, v in sorted(mod_counter.items()))
        peptidoform_id = f"{mod_str}?{bare_seq}" if mod_str else bare_seq
        variant_str = _format_variant(peptide.variant)
        try:
            mass = None if "X" in ann.stripped_sequence else ann.mass()
        except Exception:
            mass = None
        if mass is None:
            logger.debug("Could not compute mass for peptide %s (protein %s)", str(ann), protein_id)
        try:
            mass_array: list[float] | None = _proforma_to_mass_array(ann, _mod_mass_cache)
        except Exception:
            mass_array = None
        rows.append(
            (
                protein_id,
                str(ann),
                peptidoform_id,
                variant_str,
                len(ann),
                mass,
                peptide.missed_cleavages,
                peptide.semi_enzymatic,
                mass_array,
            )
        )
    return _dedup_rows(rows)


def _digest_batch_worker(
    batch: list[pf.SequenceEntry],
    config: DigestConfig,
) -> list[tuple[str, str, str | None, int, float | None]]:
    rows = []
    for sequence in batch:
        rows.extend(_digest_worker(sequence, config))
    return rows
