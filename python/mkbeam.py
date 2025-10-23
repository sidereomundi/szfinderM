"""Python translation of ``src/main_mkbeam.c``.

The original C program builds a Gaussian beam map for the SZ finder pipeline,
writes it as a FITS image, and stores the forward Fourier transform in the
binary format produced by :func:`fftsave`.  This module mirrors that behaviour
using NumPy for the numerical work and ``astropy`` for FITS I/O.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Iterable

import numpy as np

from fft_io import save_fft
from fits_utils import save_imagef
from simulation import (
    DEFAULT_PARAMETER_FILE,
    MapGeometry,
    SimulationMap,
    initialise_map,
    load_map_geometry,
)

# Constant used by the original C implementation to convert FWHM to sigma.
FWHM_TO_SIGMA: float = 2.354


def build_gaussian_beam(
    fwhm_arcmin: float,
    normalization: float,
    *,
    geometry: MapGeometry | None = None,
    sim_map: SimulationMap | None = None,
) -> SimulationMap:
    """Populate ``sim_map`` with a rotationally symmetric Gaussian beam."""

    if fwhm_arcmin <= 0:
        raise ValueError("The FWHM must be positive.")
    if normalization <= 0:
        raise ValueError("The normalization must be positive.")

    if sim_map is None:
        if geometry is None:
            raise ValueError("A geometry must be supplied when no simulation map is provided.")
        sim_map = initialise_map(geometry)
    elif geometry is not None and geometry != sim_map.geometry:
        msg = "The provided geometry does not match the existing simulation map."
        raise ValueError(msg)

    fwhm_pixels = fwhm_arcmin / sim_map.pixsize
    sigma_pixels = fwhm_pixels / FWHM_TO_SIGMA
    sigma_squared = sigma_pixels * sigma_pixels
    if sigma_squared == 0:
        raise ValueError("The sigma derived from the provided FWHM is zero.")

    scale = normalization / (2.0 * math.pi * sigma_squared)

    x_indices = np.arange(sim_map.nsidepix, dtype=np.float64)
    x_indices = np.minimum(x_indices, sim_map.nsidepix - x_indices)
    y_indices = np.arange(sim_map.nsidepix, dtype=np.float64)
    y_indices = np.minimum(y_indices, sim_map.nsidepix - y_indices)
    radius_squared = x_indices[np.newaxis, :] ** 2 + y_indices[:, np.newaxis] ** 2
    gaussian = np.exp(-0.5 * radius_squared / sigma_squared)

    np.multiply(scale, gaussian, out=sim_map.image, casting="unsafe")
    return sim_map


def parse_args(argv: Iterable[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate a Gaussian beam map.")
    parser.add_argument("fwhm", type=float, help="Full width at half maximum in arcminutes")
    parser.add_argument("output", help="Base name for the generated files")
    parser.add_argument(
        "normalization",
        type=float,
        nargs="?",
        default=1.0,
        help="Optional normalization factor applied to the Gaussian profile",
    )
    parser.add_argument(
        "--param-file",
        type=Path,
        default=DEFAULT_PARAMETER_FILE,
        help="Path to a parameter file containing NDIM/FIELDSIZE/EXTENDEDDIM definitions",
    )
    parser.add_argument("--ndim", type=int, help="Override the NDIM value from the parameter file")
    parser.add_argument(
        "--fieldsize",
        type=float,
        help="Override the FIELDSIZE value (degrees) from the parameter file",
    )
    parser.add_argument(
        "--extendeddim",
        type=int,
        help="Override the EXTENDEDDIM value from the parameter file",
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
    sim_map = build_gaussian_beam(args.fwhm, args.normalization, geometry=geometry)

    fits_name = Path(f"{args.output}.fits")
    save_imagef(fits_name, sim_map.image)
    save_fft(sim_map.image, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
