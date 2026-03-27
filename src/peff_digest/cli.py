from __future__ import annotations

import argparse
import csv
import itertools
import logging
import multiprocessing as mp
import sys
import tempfile
from collections.abc import Generator
from functools import partial
from pathlib import Path

import pefftacular as pf
import pydantic
from tqdm import tqdm

from peff_digest.config import DigestConfig
from peff_digest.digest import Peptide, digest_peff_sequence

logger = logging.getLogger(__name__)

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


def digest_sequences(
    sequences: list[pf.SequenceEntry],
    config: DigestConfig,
    show_progress: bool = False,
) -> Generator[Peptide, None, None]:
    """Digest all sequences in a PEFF/FASTA file, optionally showing a progress bar."""
    for seq in tqdm(sequences, disable=not show_progress, desc="Digesting", unit="protein"):
        yield from digest_peff_sequence(
            seq,
            cleave_on=config.cleave_on,
            missed_cleavages=config.missed_cleavages,
            semi_enzymatic=config.semi_enzymatic,
            max_ptm_per_peptide=config.max_ptm_per_peptide,
            min_length=config.min_length,
            max_length=config.max_length,
            restrict_after=config.restrict_after,
            restrict_before=config.restrict_before,
            cterminal=config.cterminal,
            internal_mods=config.internal_mods or None,
            terminal_mods=config.terminal_mods or None,
            annotate_variants=config.annotate_variants,
            use_mod_names=config.use_mod_names,
            use_psi_mods=config.use_psi_mods,
            include_simple_variants=config.include_simple_variants,
            include_complex_variants=config.include_complex_variants,
        )


def digest_sequence(
    sequence: pf.SequenceEntry,
    config: DigestConfig,
) -> Generator[Peptide, None, None]:
    """
    Digest a single PEFF sequence entry and return a set of Peptide objects.
    """
    return digest_peff_sequence(
        sequence,
        cleave_on=config.cleave_on,
        missed_cleavages=config.missed_cleavages,
        semi_enzymatic=config.semi_enzymatic,
        max_ptm_per_peptide=config.max_ptm_per_peptide,
        min_length=config.min_length,
        max_length=config.max_length,
        restrict_after=config.restrict_after,
        restrict_before=config.restrict_before,
        cterminal=config.cterminal,
        internal_mods=config.internal_mods or None,
        terminal_mods=config.terminal_mods or None,
        annotate_variants=config.annotate_variants,
        use_mod_names=config.use_mod_names,
        use_psi_mods=config.use_psi_mods,
        include_simple_variants=config.include_simple_variants,
        include_complex_variants=config.include_complex_variants,
    )


def _digest_worker(
    sequence: pf.SequenceEntry,
    config: DigestConfig,
) -> list[tuple[str, str, str | None, int, float | None]]:
    protein_id = sequence.db_unique_id
    peptides = digest_sequence(sequence, config)
    rows = []
    for peptide in peptides:
        ann = peptide.proforma
        name = ann.peptide_name
        ann.peptide_name = None
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



def main() -> None:
    mp.freeze_support()
    logging.basicConfig(level=logging.INFO, format="%(message)s", stream=sys.stderr)

    parser = argparse.ArgumentParser(
        prog="peff-digest",
        description="Digest a PEFF or FASTA file and write peptides to CSV.",
    )
    parser.add_argument("--config", metavar="FILE", help="JSON or TOML config file")
    parser.add_argument("--input-file", metavar="FILE", help="Input PEFF or FASTA file")
    parser.add_argument(
        "--output-file", metavar="FILE", default=argparse.SUPPRESS, help="Output CSV file (default: peptides.csv)"
    )
    parser.add_argument(
        "--cleave-on", metavar="AAS", default=argparse.SUPPRESS, help="Amino acids to cleave on (default: KR)"
    )
    parser.add_argument("--missed-cleavages", type=int, metavar="N", default=argparse.SUPPRESS)
    parser.add_argument("--semi-enzymatic", action="store_true", default=argparse.SUPPRESS)
    parser.add_argument("--max-ptm-per-peptide", type=int, metavar="N", default=argparse.SUPPRESS)
    parser.add_argument("--min-length", type=int, metavar="N", default=argparse.SUPPRESS)
    parser.add_argument("--max-length", type=int, metavar="N", default=argparse.SUPPRESS)
    parser.add_argument("--restrict-after", metavar="AAS", default=argparse.SUPPRESS)
    parser.add_argument("--restrict-before", metavar="AAS", default=argparse.SUPPRESS)
    parser.add_argument(
        "--workers", type=int, metavar="N", default=argparse.SUPPRESS, help="Worker processes (default: all CPUs)"
    )
    parser.add_argument(
        "--batch-size", type=int, metavar="N", default=argparse.SUPPRESS, help="Sequences per worker batch (default: 1)"
    )
    parser.add_argument("--min-mass", type=float, metavar="DA", default=argparse.SUPPRESS)
    parser.add_argument("--max-mass", type=float, metavar="DA", default=argparse.SUPPRESS)
    parser.add_argument(
        "--drop-invalid-mass",
        action="store_true",
        default=argparse.SUPPRESS,
        help="Exclude peptides whose mass cannot be computed",
    )
    parser.add_argument(
        "--no-annotate-variants",
        dest="annotate_variants",
        action="store_false",
        default=argparse.SUPPRESS,
        help="Do not set peptide_name on variant peptides",
    )
    parser.add_argument(
        "--no-psi-mods",
        dest="use_psi_mods",
        action="store_false",
        default=argparse.SUPPRESS,
        help="Do not apply PSI-MOD (ModResPsi) annotations from the PEFF file",
    )
    parser.add_argument(
        "--no-simple-variants",
        dest="include_simple_variants",
        action="store_false",
        default=argparse.SUPPRESS,
        help="Skip VariantSimple (single AA substitution) entries",
    )
    parser.add_argument(
        "--no-complex-variants",
        dest="include_complex_variants",
        action="store_false",
        default=argparse.SUPPRESS,
        help="Skip VariantComplex (multi-residue insertion/deletion/substitution) entries",
    )
    args = parser.parse_args()

    # Build config: file first, then CLI overrides
    cli_overrides = {
        k: v for k, v in vars(args).items() if k != "config" and k in DigestConfig.model_fields and v is not None
    }
    cli_overrides = {k.replace("-", "_"): v for k, v in cli_overrides.items()}

    if args.input_file:
        cli_overrides["input_file"] = args.input_file

    try:
        if args.config:
            config = DigestConfig.from_file(args.config, **cli_overrides)
            logger.info("Config loaded from %s", args.config)
        else:
            config = DigestConfig(**cli_overrides)
    except pydantic.ValidationError as exc:
        for error in exc.errors():
            loc = ".".join(str(x) for x in error["loc"])
            logger.error("  %s: %s", loc, error["msg"])
        sys.exit(1)

    logger.info("Reading sequences from %s", config.input_file)
    sequences, n_malformed = read_sequences(config.input_file)
    skipped = f", {n_malformed} malformed entries skipped" if n_malformed else ""
    logger.info("%d sequences loaded%s", len(sequences), skipped)

    batches = list(itertools.batched(sequences, config.batch_size))
    logger.info(
        "Digesting with %d worker(s), batch size %d (%d batches total)",
        config.workers, config.batch_size, len(batches),
    )
    worker = partial(_digest_batch_worker, config=config)
    with mp.Pool(config.workers) as pool:
        results = list(tqdm(pool.imap(worker, batches), total=len(batches), desc="Digesting", unit="batch"))

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


if __name__ == "__main__":
    main()
