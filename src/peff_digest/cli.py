from __future__ import annotations

import argparse
import csv
import json
import multiprocessing as mp
import sys
import tempfile
import tomllib
from collections.abc import Generator
from functools import partial
from pathlib import Path

import pefftacular as pf
import pydantic
from tqdm import tqdm

from peff_digest.config import DigestConfig
from peff_digest.digest import Peptide, digest_peff_sequence

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
    else:
        reader = iter(pf.PeffReader(path))
        while True:
            try:
                sequences.append(next(reader))
            except StopIteration:
                break
            except Exception:
                n_malformed += 1
    return sequences, n_malformed


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
        if mass is None and config.drop_invalid_mass:
            continue
        if mass is not None and config.min_mass is not None and mass < config.min_mass:
            continue
        if mass is not None and config.max_mass is not None and mass > config.max_mass:
            continue
        rows.append((protein_id, str(ann), name, len(ann), mass))
    return rows


def _load_config_file(path: str) -> dict:
    p = Path(path)
    if p.suffix == ".toml":
        with open(p, "rb") as f:
            return tomllib.load(f)
    with open(p) as f:
        return json.load(f)


def main() -> None:
    mp.freeze_support()

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
    args = parser.parse_args()

    # Build config: file first, then CLI overrides
    file_data: dict = {}
    if args.config:
        file_data = _load_config_file(args.config)

    cli_overrides = {
        k: v for k, v in vars(args).items() if k != "config" and k in DigestConfig.model_fields and v is not None
    }
    cli_overrides = {k.replace("-", "_"): v for k, v in cli_overrides.items()}

    if args.input_file:
        cli_overrides["input_file"] = args.input_file

    try:
        config = DigestConfig(**{**file_data, **cli_overrides})
    except pydantic.ValidationError as exc:
        for error in exc.errors():
            loc = ".".join(str(x) for x in error["loc"])
            print(f"  {loc}: {error['msg']}", file=sys.stderr)
        sys.exit(1)

    sequences, n_malformed = read_sequences(config.input_file)
    print(f"{len(sequences)} sequences loaded, {n_malformed} malformed entries skipped")

    worker = partial(_digest_worker, config=config)
    with mp.Pool(config.workers) as pool:
        results = list(tqdm(pool.imap(worker, sequences), total=len(sequences)))

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

    print(f"Written to {config.output_file}")


if __name__ == "__main__":
    main()
