from __future__ import annotations

from pathlib import Path

import pefftacular as pf
import psimodpy
import pytest

from peff_digest import DigestConfig, InternalMod

DATA_DIR = Path(__file__).parent / "data"
PEFF_FILES = sorted(DATA_DIR.glob("*.peff"))
FASTA_FILES = sorted(DATA_DIR.glob("*.fasta"))


@pytest.fixture(scope="session")
def psi_db() -> psimodpy.PsiModDatabase:
    return psimodpy.load()


def _make_entry(sequence: str, **kwargs) -> pf.SequenceEntry:
    return pf.SequenceEntry(prefix="sp", db_unique_id="TEST_ID", sequence=sequence, **kwargs)


def _cfg(**kwargs) -> DigestConfig:
    """Create a DigestConfig with test-friendly defaults."""
    defaults = dict(min_length=1, max_length=40, missed_cleavages=0, max_ptm_per_peptide=0)
    return DigestConfig(**{**defaults, **kwargs})


@pytest.fixture()
def tiny_fasta(tmp_path: Path) -> Path:
    fasta = tmp_path / "tiny.fasta"
    fasta.write_text(">sp|P12345|TEST\nACDEFGHIKLMNPQRSTVWYK\n")
    return fasta
