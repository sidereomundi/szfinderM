"""Python translation of ``src/main_mknull.c``."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import numpy as np

from fits_utils import save_imagef
from simulation import DEFAULT_PARAMETER_FILE, load_map_geometry


def parse_args(argv: Iterable[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate an SZ-null map (all zeros) that matches the simulation geometry.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("sznull.fits"),
        help="Destination FITS filename for the null map",
    )
    parser.add_argument(
        "--param-file",
        type=Path,
        default=DEFAULT_PARAMETER_FILE,
        help="Parameter header defining NDIM and related geometry settings",
    )
    parser.add_argument("--ndim", type=int, help="Override NDIM from the parameter file")
    parser.add_argument(
        "--fieldsize",
        type=float,
        help="Override FIELDSIZE (degrees) from the parameter file",
    )
    parser.add_argument(
        "--extendeddim",
        type=int,
        help="Override EXTENDEDDIM from the parameter file",
    )
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    geometry = load_map_geometry(
        args.param_file,
        ndim=args.ndim,
        fieldsize=args.fieldsize,
        extendeddim=args.extendeddim,
    )

    data = np.zeros((geometry.ndim, geometry.ndim), dtype=np.float32)
    save_imagef(args.output, data)
    print(f"wrote SZ-null map to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
