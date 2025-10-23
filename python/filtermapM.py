"""Python translation of ``src/main_filtermapM.c``."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Iterable

import numpy as np

from fft_io import load_fft
from fits_utils import read_imagef, save_imagef
from simulation import (
    DEFAULT_PARAMETER_FILE,
    initialise_map,
    load_map_geometry,
    parse_parameter_file,
)


def parse_args(argv: Iterable[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Filter a pair of observed maps with the multi-frequency filters.",
    )
    parser.add_argument(
        "prefix",
        nargs="?",
        help=(
            "Prefix for the observed maps. The script expects '<prefix>150.fits' and "
            "'<prefix>95.fits' unless explicit filenames are provided."
        ),
    )
    parser.add_argument(
        "--map-150",
        dest="map_150",
        type=Path,
        help="FITS image containing the 150 GHz observation",
    )
    parser.add_argument(
        "--map-95",
        dest="map_95",
        type=Path,
        help="FITS image containing the 95 GHz observation",
    )
    parser.add_argument(
        "--param-file",
        type=Path,
        default=DEFAULT_PARAMETER_FILE,
        help="Parameter header containing the geometry and filter configuration",
    )
    parser.add_argument(
        "--ndim",
        type=int,
        help="Override the NDIM value from the parameter file",
    )
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
    parser.add_argument(
        "--filter-prefix",
        help="Base filename used for the filters (defaults to FILTERFILE)",
    )
    parser.add_argument(
        "--sigma-file",
        type=Path,
        help="File that stores the per-filter variance (defaults to FILTERMSIGFILE)",
    )
    parser.add_argument(
        "--num-filters",
        type=int,
        help="Number of filters to apply (defaults to NFILTER)",
    )
    parser.add_argument(
        "--output-template",
        default="filteredM{index:02d}.fits",
        help=(
            "Pattern used for the output FITS maps. Include '{index}' to insert the "
            "filter number."
        ),
    )
    return parser.parse_args(argv)


def _load_parameters(param_file: Path) -> dict[str, object]:
    values = parse_parameter_file(param_file)
    required = ["FILTERFILE", "FILTERMSIGFILE", "NFILTER"]
    missing = [name for name in required if name not in values]
    if missing:
        joined = ", ".join(missing)
        raise ValueError(f"Missing parameters in {param_file}: {joined}")
    return values


def _resolve_map_path(prefix: str | None, override: Path | None, suffix: str) -> Path:
    if override is not None:
        return override
    if prefix:
        return Path(f"{prefix}{suffix}")
    default = "szmock150.fits" if suffix == "150.fits" else "szmock95.fits"
    return Path(default)


def _load_sigma_values(path: Path, count: int) -> list[tuple[float, float]]:
    results: list[tuple[float, float]] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                theta = float(parts[0])
                sigma = float(parts[1])
            except ValueError:
                continue
            results.append((theta, sigma))
            if len(results) == count:
                break
    if len(results) < count:
        raise ValueError(
            f"Only {len(results)} entries were found in {path}, expected {count} filter variances.",
        )
    return results


def _fft_filter_path(prefix: str, band: str, index: int) -> Path:
    return Path(f"{prefix}{band}_{index}_fft.dat")


def _prepare_map(path: Path, expected_shape: tuple[int, int]) -> np.ndarray:
    image = read_imagef(path)
    if image.shape != expected_shape:
        raise ValueError(
            f"Image {path} has shape {image.shape}, expected {expected_shape} to match the simulation geometry.",
        )
    array = image.astype(np.float64, copy=False)
    mean = float(np.mean(array, dtype=np.float64))
    print(f"mean of map {path} is {mean:.6e}")
    return array - mean


def _normalise_map(filtered: np.ndarray) -> tuple[np.ndarray, float, float]:
    mean_value = float(np.mean(filtered, dtype=np.float64))
    square_mean = float(np.mean(np.square(filtered, dtype=np.float64), dtype=np.float64))
    variance = square_mean - mean_value * mean_value
    if variance <= 0.0:
        raise ValueError("Filtered map variance is non-positive; cannot normalise.")
    stddev = math.sqrt(variance)
    normalised = (filtered / stddev).astype(np.float32)
    return normalised, mean_value, stddev


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    parameters = _load_parameters(args.param_file)

    geometry = load_map_geometry(
        args.param_file,
        ndim=args.ndim,
        fieldsize=args.fieldsize,
        extendeddim=args.extendeddim,
    )
    sim_map = initialise_map(geometry)
    expected_shape = (sim_map.nsidepix, sim_map.nsidepix)

    map150_path = _resolve_map_path(args.prefix, args.map_150, "150.fits")
    map95_path = _resolve_map_path(args.prefix, args.map_95, "95.fits")

    print(" # preparing maps")
    image150 = _prepare_map(map150_path, expected_shape)
    image95 = _prepare_map(map95_path, expected_shape)

    print(" # computing FFTs of the observed maps")
    fft150 = np.fft.rfft2(image150)
    fft95 = np.fft.rfft2(image95)

    filter_prefix = str(args.filter_prefix or parameters["FILTERFILE"])
    sigma_file = Path(args.sigma_file or parameters["FILTERMSIGFILE"])
    num_filters = int(args.num_filters or parameters["NFILTER"])
    if num_filters <= 0:
        raise ValueError("The number of filters must be positive.")

    sigma_values = _load_sigma_values(sigma_file, num_filters)

    template = args.output_template
    if "{index" not in template:
        raise ValueError("The output template must contain '{index}' to insert the filter number.")

    print(f" # applying {num_filters} filters using prefix {filter_prefix}")
    for index in range(num_filters):
        _theta, sigma0 = sigma_values[index]
        print(f"# Filtering case {index}")

        filter150_path = _fft_filter_path(filter_prefix, "M150", index)
        filter95_path = _fft_filter_path(filter_prefix, "M90", index)

        filter150 = load_fft(filter150_path, sim_map.nsidepix)
        filter95 = load_fft(filter95_path, sim_map.nsidepix)

        filtered150 = np.fft.irfft2(fft150 * np.conjugate(filter150), s=expected_shape)
        filtered95 = np.fft.irfft2(fft95 * np.conjugate(filter95), s=expected_shape)

        combined = -(filtered150 + filtered95)

        normalised, mean_value, stddev = _normalise_map(combined)
        print(f"case {index}:: {mean_value:e} {stddev:e} {sigma0:e}")

        output_path = Path(template.format(index=index))
        save_imagef(output_path, normalised)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
