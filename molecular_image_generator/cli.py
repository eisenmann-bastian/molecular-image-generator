"""Generate 2D PNG images from molecular structures or identifiers."""

import argparse
import csv
from pathlib import Path

from molecular_selector import MolecularSelector
from molecular_utilities.logging import configure_logging, get_logger

from .generator import ImageGenerator

NAMESPACE = "molecular_image_generator"
logger = get_logger("cli", namespace=NAMESPACE)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--sdf", type=Path, help="input SDF file (all records)")
    inputs.add_argument("--smiles", help="SMILES string")
    inputs.add_argument("--inchi", help="InChI string")
    inputs.add_argument("--name", help="molecule name resolved using PubChem")
    inputs.add_argument("--csv", type=Path, help="CSV file with molecular input columns")
    inputs.add_argument("--list-selectors", action="store_true", help="list named highlight selectors")
    parser.add_argument("-o", "--output-dir", type=Path, default=Path("images"))
    parser.add_argument("--width", type=int, default=600)
    parser.add_argument("--height", type=int, default=600)
    parser.add_argument("--delimiter", default=",", help="CSV delimiter (default: comma)")
    parser.add_argument(
        "--highlight", action="append", default=[], metavar="NAME",
        help="highlight a molecular-selector name; repeat for multiple selectors",
    )
    parser.add_argument(
        "--highlight-smarts", action="append", default=[], metavar="SMARTS",
        help="highlight a SMARTS pattern; repeat for multiple patterns",
    )
    parser.add_argument(
        "--highlight-common", action="store_true",
        help="highlight common functional groups and aromatic atoms and bonds",
    )
    args = parser.parse_args(argv)
    configure_logging(namespace=NAMESPACE)
    try:
        if args.list_selectors:
            if args.highlight or args.highlight_smarts or args.highlight_common:
                raise ValueError("--list-selectors cannot be combined with highlight options")
            for name, definition in MolecularSelector().registry.definitions.items():
                aliases = f" (aliases: {', '.join(definition.aliases)})" if definition.aliases else ""
                print(f"{name}{aliases}")
            return 0
        generator = ImageGenerator(
            width=args.width, height=args.height,
            highlight_names=args.highlight, highlight_smarts=args.highlight_smarts,
            highlight_common=args.highlight_common,
        )
        if args.csv is not None:
            paths, errors = generator.generate_csv(args.csv, args.output_dir, delimiter=args.delimiter)
        else:
            kind = next(key for key in ("sdf", "smiles", "inchi", "name")
                        if getattr(args, key) is not None)
            paths = generator.generate(kind, getattr(args, kind), args.output_dir)
            errors = []
        for path in paths:
            logger.info("Saved %s", path)
        for error in errors:
            logger.error("%s", error)
        return 1 if errors else 0
    except (ValueError, OSError, RuntimeError, csv.Error) as exc:
        logger.error("%s", exc)
        return 1
