from __future__ import annotations

import multiprocessing as mp
from functools import partial

import polars as pl

from peff_digest.cli import _digest_worker, read_sequences
from peff_digest.config import DigestConfig


def digest(config: DigestConfig) -> pl.DataFrame:
    """
    Run a full digest and return results as a Polars DataFrame.

    Columns: protein_id (str), sequence (str), variant (str | null),
             length (i64), mass (f64 | null).
    """
    sequences, _ = read_sequences(config.input_file)

    worker = partial(_digest_worker, config=config)
    with mp.Pool(config.workers) as pool:
        results = list(pool.imap(worker, sequences))

    protein_ids: list[str] = []
    seqs: list[str] = []
    variants: list[str | None] = []
    lengths: list[int] = []
    masses: list[float | None] = []

    for rows in results:
        for protein_id, seq, variant, length, mass in rows:
            protein_ids.append(protein_id)
            seqs.append(seq)
            variants.append(variant)
            lengths.append(length)
            masses.append(mass)

    return pl.DataFrame(
        {
            "protein_id": protein_ids,
            "sequence": seqs,
            "variant": variants,
            "length": lengths,
            "mass": masses,
        },
        schema={
            "protein_id": pl.String,
            "sequence": pl.String,
            "variant": pl.String,
            "length": pl.Int64,
            "mass": pl.Float64,
        },
    )
