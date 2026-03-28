"""File I/O, database caching, and multiprocessing worker functions."""

from __future__ import annotations

import logging
from pathlib import Path

import pefftacular as pf
import peptacular as pt
from psimodpy import PsiModDatabase

from peff_digest.config import DigestConfig
from peff_digest.digest import _apply_mod, ann_to_map, digest_peff_sequence

logger = logging.getLogger(__name__)

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
# PSI-MOD → UniMod conversion
# ---------------------------------------------------------------------------


def _try_convert_psimod_to_unimod(
    ann: pt.ProFormaAnnotation,
    psi_db: PsiModDatabase,
    uni_db,
) -> pt.ProFormaAnnotation | None:
    """Replace MOD:NNNNN tags with UNIMOD:N accessions.

    Returns None if any PSI-MOD mod has no UniMod xref — the caller should drop the peptide.
    Non-PSI-MOD tags (user-added mods) are passed through unchanged.
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
        else:
            new_mod_map[pos] = tag
    new_ann = pt.parse(sequence)
    seq_len = len(sequence)
    for pos, tag in new_mod_map.items():
        _apply_mod(new_ann, pos, tag, seq_len)
    return new_ann


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


def _digest_worker(
    sequence: pf.SequenceEntry,
    config: DigestConfig,
) -> list[tuple[str, str, str | None, int, float | None]]:
    psi_db = _get_psi_db() if config.use_psi_mods else None
    protein_id = sequence.db_unique_id
    peptides = digest_peff_sequence(sequence, config, psi_db=psi_db)
    rows = []
    for peptide in peptides:
        ann = peptide.proforma
        name = ann.peptide_name
        ann.peptide_name = None
        if config.use_unimod_output:
            ann = _try_convert_psimod_to_unimod(ann, _get_psi_db(), _get_uni_db())
            if ann is None:
                continue
        try:
            mass = ann.mass()
        except Exception:
            mass = None
        if mass is None:
            logger.debug("Could not compute mass for peptide %s (protein %s)", str(ann), protein_id)
        if mass is None and config.drop_invalid_mass:
            continue
        if mass is not None and config.min_mass is not None and mass < config.min_mass:
            continue
        if mass is not None and config.max_mass is not None and mass > config.max_mass:
            continue
        rows.append((protein_id, str(ann), name, len(ann), mass))
    return rows


def _digest_batch_worker(
    batch: list[pf.SequenceEntry],
    config: DigestConfig,
) -> list[tuple[str, str, str | None, int, float | None]]:
    rows = []
    for sequence in batch:
        rows.extend(_digest_worker(sequence, config))
    return rows
