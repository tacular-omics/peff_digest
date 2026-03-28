from __future__ import annotations

import csv
import itertools
import logging
import multiprocessing as mp
import tempfile
from functools import partial
from pathlib import Path

import polars as pl
from tqdm import tqdm

from peff_digest.cli import _digest_batch_worker, read_sequences
from peff_digest.config import DigestConfig

logger = logging.getLogger(__name__)


def digest(
    config: DigestConfig,
    show_progress: bool = False,
    write_file: bool = False,
) -> pl.DataFrame | None:
    """
    Run a full digest.

    If write_file is False (default), returns a Polars DataFrame with columns:
        protein_id (str), sequence (str), variant (str | null), length (i64), mass (f64 | null).
    If write_file is True, writes results to config.output_file as CSV and returns None.
    """
    logger.info("Reading sequences from %s", config.input_file)
    sequences, n_malformed = read_sequences(config.input_file)
    skipped = f", {n_malformed} malformed entries skipped" if n_malformed else ""
    logger.info("%d sequences loaded%s", len(sequences), skipped)

    batches = list(itertools.batched(sequences, config.batch_size))
    logger.info(
        "Digesting with %d worker(s), batch size %d (%d batches total)",
        config.workers,
        config.batch_size,
        len(batches),
    )
    worker = partial(_digest_batch_worker, config=config)
    with mp.Pool(config.workers) as pool:
        imap = pool.imap(worker, batches)
        results = list(tqdm(imap, total=len(batches), disable=not show_progress, desc="Digesting", unit="batch"))

    if write_file:
        output_path = Path(config.output_file)
        tmp_fd, tmp_path = tempfile.mkstemp(dir=output_path.parent, suffix=".tmp")
        try:
            with open(tmp_fd, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["protein_id", "sequence", "variant", "length", "mass"])
                for rows in results:
                    writer.writerows(rows)
            Path(tmp_path).replace(output_path)
        except Exception:
            Path(tmp_path).unlink(missing_ok=True)
            raise
        n_peptides = sum(len(r) for r in results)
        logger.info("%d peptides written to %s", n_peptides, config.output_file)
        return None

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

    logger.info("%d peptides returned", len(protein_ids))
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
