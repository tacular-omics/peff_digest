"""Command-line interface for peff-digest."""

from __future__ import annotations

import argparse
import logging
import multiprocessing as mp
import sys

import pydantic

from peff_digest.api import digest
from peff_digest.config import DigestConfig


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
        "--no-psi-mods",
        dest="use_psi_mods",
        action="store_false",
        default=argparse.SUPPRESS,
        help="Do not apply PSI-MOD (ModResPsi) annotations from the PEFF file",
    )
    parser.add_argument(
        "--use-unimod-output",
        action="store_true",
        default=argparse.SUPPRESS,
        help="Convert PSI-MOD mod tags to UniMod accessions in output; peptides with unmapped mods are dropped",
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

    logger = logging.getLogger(__name__)
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

    digest(config, show_progress=True, write_file=True)


if __name__ == "__main__":
    main()
