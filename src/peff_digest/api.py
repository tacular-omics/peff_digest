"""High-level digest API — single orchestration entry point."""

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

from peff_digest.config import DigestConfig
from peff_digest.io import _digest_batch_worker, read_sequences

logger = logging.getLogger(__name__)


def _write_csv(
    results: list[list[_Row]],
    output_file: str,
) -> int:
    """Write digest results to CSV atomically. Returns the number of peptides written."""
    output_path = Path(output_file)
    tmp_fd, tmp_path = tempfile.mkstemp(dir=output_path.parent, suffix=".tmp")
    try:
        with open(tmp_fd, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(
                [
                    "protein_id",
                    "sequence",
                    "peptidoform_id",
                    "variant",
                    "length",
                    "neutral_mass",
                    "missed_cleavages",
                    "semi_enzymatic",
                    "mass_array",
                    "n_peptidoforms",
                ]
            )
            for rows in results:
                for row in rows:
                    protein_id, seq, pform_id, variant, length, mass, mc, semi, mass_array, n_pf = row
                    mass_array_str = str(mass_array) if mass_array is not None else ""
                    writer.writerow([protein_id, seq, pform_id, variant, length, mass, mc, semi, mass_array_str, n_pf])
        Path(tmp_path).replace(output_path)
    except Exception:
        Path(tmp_path).unlink(missing_ok=True)
        raise
    return sum(len(r) for r in results)


_Row = tuple[str, str, str, str | None, int, float | None, int, bool, list[float] | None, int]


def _to_dataframe(results: list[list[_Row]]) -> pl.DataFrame:
    """Convert digest results to a Polars DataFrame."""
    protein_ids: list[str] = []
    seqs: list[str] = []
    peptidoform_ids: list[str] = []
    variants: list[str | None] = []
    lengths: list[int] = []
    masses: list[float | None] = []
    missed_cleavages: list[int] = []
    semi_enzymatic: list[bool] = []
    n_peptidoforms: list[int] = []

    mass_arrays: list[list[float] | None] = []

    for rows in results:
        for protein_id, seq, peptidoform_id, variant, length, mass, mc, semi, mass_array, n_pf in rows:
            protein_ids.append(protein_id)
            seqs.append(seq)
            peptidoform_ids.append(peptidoform_id)
            variants.append(variant)
            lengths.append(length)
            masses.append(mass)
            missed_cleavages.append(mc)
            semi_enzymatic.append(semi)
            mass_arrays.append(mass_array)
            n_peptidoforms.append(n_pf)

    return pl.DataFrame(
        {
            "protein_id": protein_ids,
            "sequence": seqs,
            "peptidoform_id": peptidoform_ids,
            "variant": variants,
            "length": lengths,
            "neutral_mass": masses,
            "missed_cleavages": missed_cleavages,
            "semi_enzymatic": semi_enzymatic,
            "mass_array": mass_arrays,
            "n_peptidoforms": n_peptidoforms,
        },
        schema={
            "protein_id": pl.String,
            "sequence": pl.String,
            "peptidoform_id": pl.String,
            "variant": pl.String,
            "length": pl.Int64,
            "neutral_mass": pl.Float64,
            "missed_cleavages": pl.Int32,
            "semi_enzymatic": pl.Boolean,
            "mass_array": pl.List(pl.Float64),
            "n_peptidoforms": pl.Int32,
        },
    )


def _consolidate_proteins(results: list[list[_Row]]) -> list[list[_Row]]:
    """Merge rows sharing a peptidoform_id across proteins.

    For each peptidoform_id the best row is selected (lowest missed_cleavages,
    then prefer non-semi-enzymatic). protein_id is replaced with all contributing
    protein IDs joined by ';', in order of first appearance. n_peptidoforms is
    summed across all contributing proteins.
    Returns a single-element list[list[row]] compatible with _write_csv / _to_dataframe.
    """
    # Consolidate on the full ProForma sequence so distinct positional variants
    # remain as separate rows; only the same exact sequence across proteins is merged.
    best: dict[str, _Row] = {}
    pids: dict[str, list[str]] = {}
    order: list[str] = []

    for batch_rows in results:
        for row in batch_rows:
            protein_id, seq = row[0], row[1]
            if seq not in best:
                best[seq] = row
                pids[seq] = [protein_id]
                order.append(seq)
            else:
                if protein_id not in pids[seq]:
                    pids[seq].append(protein_id)
                if (row[6], row[7]) < (best[seq][6], best[seq][7]):
                    best[seq] = row

    # Rebuild merged rows (without n_peptidoforms), then recompute the count
    merged = [(";".join(pids[s]), *best[s][1:9]) for s in order]
    pform_count: dict[str, int] = {}
    for row in merged:
        pform_count[row[2]] = pform_count.get(row[2], 0) + 1

    consolidated: list[_Row] = [(*row, pform_count[row[2]]) for row in merged]
    return [consolidated]


def digest(
    config: DigestConfig,
    show_progress: bool = False,
    write_file: bool = False,
) -> pl.DataFrame | None:
    """
    Run a full digest pipeline.

    If write_file is False (default), returns a Polars DataFrame with columns:
        protein_id, sequence, peptidoform_id, variant, length, neutral_mass, missed_cleavages, semi_enzymatic.
    If write_file is True, writes results to config.output_file as CSV and returns None.
    """
    logger.info("Reading sequences from %s", config.input_file)
    sequences, n_malformed = read_sequences(config.input_file)
    skipped = f", {n_malformed} malformed entries skipped" if n_malformed else ""
    logger.info("%d sequences loaded%s", len(sequences), skipped)

    batches = list(itertools.batched(sequences, config.batch_size))
    n_workers = config.workers or mp.cpu_count()
    logger.info(
        "Digesting with %d worker(s), batch size %d (%d batches total)",
        n_workers,
        config.batch_size,
        len(batches),
    )
    worker = partial(_digest_batch_worker, config=config)
    with mp.Pool(config.workers) as pool:
        imap = pool.imap(worker, batches)
        results = list(tqdm(imap, total=len(batches), disable=not show_progress, desc="Digesting", unit="batch"))

    if config.consolidate_proteins:
        results = _consolidate_proteins(results)
        logger.info("Consolidated to %d unique peptidoforms", len(results[0]))

    results = [
        sorted(
            (row for batch in results for row in batch),
            key=lambda r: (r[5] is None, r[5] or 0.0),
        )
    ]

    if write_file:
        n_peptides = _write_csv(results, config.output_file)
        logger.info("%d peptides written to %s", n_peptides, config.output_file)
        return None

    logger.info("%d peptides returned", sum(len(r) for r in results))
    return _to_dataframe(results)
