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
from typing import Iterable, Sequence

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


def _coerce_to_list(value, *, name: str) -> list:
    """Return ``value`` as a list while preserving scalar inputs."""
    if isinstance(value, (str, bytes, Path)):
        return [value]
    if isinstance(value, Sequence):
        return list(value)
    return [value]


def generate_gaussian_beams(
    fwhm: float | Sequence[float],
    output: str | Path | Sequence[str | Path],
    normalization: float | Sequence[float] | None = None,
    *,
    geometry: MapGeometry | None = None,
    param_file: Path = DEFAULT_PARAMETER_FILE,
    ndim: int | None = None,
    fieldsize: float | None = None,
    extendeddim: int | None = None,
) -> list[SimulationMap]:
    """Create one or more Gaussian beams and write their outputs to disk."""
    if geometry is None:
        geometry = load_map_geometry(
            param_file,
            ndim=ndim,
            fieldsize=fieldsize,
            extendeddim=extendeddim,
        )

    fwhm_values = [float(value) for value in _coerce_to_list(fwhm, name="fwhm")]
    output_values = [Path(value) for value in _coerce_to_list(output, name="output")]
    if len(output_values) != len(fwhm_values):
        raise ValueError("The number of outputs must match the number of FWHM values.")

    if normalization is None:
        norm_values = [1.0] * len(fwhm_values)
    else:
        norm_entries = _coerce_to_list(normalization, name="normalization")
        if len(norm_entries) not in (1, len(fwhm_values)):
            raise ValueError(
                "The number of normalizations must be either one or match the number of FWHM values."
            )
        if len(norm_entries) == 1 and len(fwhm_values) > 1:
            norm_entries = norm_entries * len(fwhm_values)
        norm_values = [float(value) for value in norm_entries]

    maps: list[SimulationMap] = []
    for fwhm_value, output_name, norm_value in zip(fwhm_values, output_values, norm_values):
        sim_map = build_gaussian_beam(fwhm_value, norm_value, geometry=geometry)
        if output_name.suffix.lower() == ".fits":
            fits_name = output_name
            fft_base = output_name.with_suffix("")
        else:
            fits_name = output_name.with_suffix(".fits")
            fft_base = output_name
        save_imagef(fits_name, sim_map.image)
        save_fft(sim_map.image, str(fft_base))
        maps.append(sim_map)
    return maps


def mkbeam(
    fwhm: float | Sequence[float],
    output: str | Path | Sequence[str | Path],
    *,
    normalization: float | Sequence[float] | None = None,
    geometry: MapGeometry | None = None,
    param_file: Path = DEFAULT_PARAMETER_FILE,
    ndim: int | None = None,
    fieldsize: float | None = None,
    extendeddim: int | None = None,
) -> list[SimulationMap]:
    """Convenience wrapper matching the historical function-style interface."""
    return generate_gaussian_beams(
        fwhm,
        output,
        normalization,
        geometry=geometry,
        param_file=param_file,
        ndim=ndim,
        fieldsize=fieldsize,
        extendeddim=extendeddim,
    )


def parse_args(argv: Iterable[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate a Gaussian beam map.")
    parser.add_argument("fwhm", type=float, nargs="?", help="Full width at half maximum in arcminutes")
    parser.add_argument("output", nargs="?", help="Base name for the generated files")
    parser.add_argument(
        "normalization",
        type=float,
        nargs="?",
        default=1.0,
        help="Optional normalization factor applied to the Gaussian profile",
    )
    parser.add_argument(
        "--fwhm-list",
        dest="fwhm_list",
        type=float,
        nargs="+",
        help="Generate multiple beams by providing a list of FWHM values (arcminutes)",
    )
    parser.add_argument(
        "--output-names",
        dest="output_list",
        nargs="+",
        help="Output base names corresponding to the provided FWHM values",
    )
    parser.add_argument(
        "--normalizations",
        dest="normalization_list",
        type=float,
        nargs="+",
        help=(
            "Optional normalization factors for each FWHM value. "
            "Omit to use a default normalization of 1.0."
        ),
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

    if args.fwhm_list is not None:
        outputs = args.output_list
        if outputs is None:
            raise ValueError("--output-names must be supplied when using --fwhm-list")
        generate_gaussian_beams(
            args.fwhm_list,
            outputs,
            args.normalization_list,
            geometry=geometry,
        )
    else:
        if args.fwhm is None or args.output is None:
            raise ValueError("Both FWHM and output name must be provided for single-beam generation.")
        generate_gaussian_beams(
            args.fwhm,
            args.output,
            args.normalization,
            geometry=geometry,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
