from __future__ import annotations

from pathlib import Path

import pefftacular as pf
import pytest
from conftest import FASTA_FILES, PEFF_FILES, _cfg

from peff_digest import digest_peff_sequence
from peff_digest.io import read_sequences


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
            _cfg(cleave_on="KR", missed_cleavages=1, min_length=4, max_length=50, max_ptm_per_peptide=2),
        )
        assert hasattr(result, "__iter__")
        n_digested += 1


@pytest.mark.parametrize("fasta_path", FASTA_FILES, ids=[f.name for f in FASTA_FILES])
def test_digest_fasta_file_no_crash(fasta_path: Path) -> None:
    """Every sequence in each example FASTA file should digest without error."""
    sequences, _ = read_sequences(str(fasta_path))
    assert len(sequences) > 0
    for entry in sequences:
        result = digest_peff_sequence(
            entry,
            _cfg(cleave_on="KR", missed_cleavages=1, min_length=4, max_length=50, max_ptm_per_peptide=2),
        )
        assert hasattr(result, "__iter__")
