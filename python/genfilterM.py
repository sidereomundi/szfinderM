"""Python translation of ``src/main_genfilterM.c``."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Iterable

import numpy as np

from cmb import load_cmb_power_spectrum
from fft_io import save_fft
from fits_utils import read_imagef, save_imagef
from simulation import (
    DEFAULT_PARAMETER_FILE,
    initialise_map,
    load_map_geometry,
    parse_parameter_file,
)


def _resolve_beam_path(name: str | Path) -> Path:
    path = Path(name)
    if path.suffix.lower() == ".fits":
        return path
    return path.with_suffix(".fits")


def convert_noise(noise_k_arcmin: float, pixsize_arcmin: float) -> float:
    """Convert a white-noise level from K/arcmin to Fourier-space power density."""

    d_t_pixel = noise_k_arcmin * 1.0e-6 / pixsize_arcmin
    scale = pixsize_arcmin * math.pi / (180.0 * 60.0)
    return (d_t_pixel * scale) ** 2


def thetasize(index: int, theta_start: float, theta_step: float) -> float:
    """Mirror ``thetasize`` from ``sub_stdroutines.c`` (result in arcminutes)."""

    return theta_start + index * theta_step


def parse_args(argv: Iterable[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate the multi-frequency matched filters used by the pipeline.",
    )
    parser.add_argument(
        "--param-file",
        type=Path,
        default=DEFAULT_PARAMETER_FILE,
        help="Parameter header that defines the geometry, beams, and physics constants",
    )
    parser.add_argument(
        "--cl-file",
        type=Path,
        help="Override the CMB power spectrum (defaults to the clfile entry in the parameter file)",
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
    parser.add_argument(
        "--noise-150",
        dest="noise_150",
        type=float,
        help="Override the 150 GHz white-noise level (K/arcmin)",
    )
    parser.add_argument(
        "--noise-95",
        dest="noise_95",
        type=float,
        help="Override the 95 GHz white-noise level (K/arcmin)",
    )
    parser.add_argument(
        "--beam-150",
        dest="beam_150",
        type=Path,
        help="Override the 150 GHz beam FITS image",
    )
    parser.add_argument(
        "--beam-95",
        dest="beam_95",
        type=Path,
        help="Override the 95 GHz beam FITS image",
    )
    parser.add_argument(
        "--filter-prefix",
        help="Base filename used for the generated filters (defaults to FILTERFILE)",
    )
    parser.add_argument(
        "--sigma-output",
        type=Path,
        help="File that records the filter variance (defaults to FILTERMSIGFILE)",
    )
    parser.add_argument(
        "--num-filters",
        type=int,
        help="Number of filter scales to generate (defaults to NFILTER)",
    )
    return parser.parse_args(argv)


def _load_parameters(param_file: Path) -> dict[str, object]:
    values = parse_parameter_file(param_file)
    required = [
        "FREQUENCY150",
        "FREQUENCY90",
        "DT150",
        "DT90",
        "BEAM150FILE",
        "BEAM90FILE",
        "FILTERFILE",
        "FILTERMSIGFILE",
        "NFILTER",
        "THETA_START",
        "THETA_STEP",
        "BETA",
        "T_CMB",
        "KB",
        "H",
        "clfile",
    ]
    missing = [name for name in required if name not in values]
    if missing:
        joined = ", ".join(missing)
        raise ValueError(f"Missing parameters in {param_file}: {joined}")
    return values


def _load_beam(path: Path, nside: int) -> np.ndarray:
    beam = read_imagef(path)
    if beam.shape != (nside, nside):
        raise ValueError(
            f"Beam image {path} has shape {beam.shape}, expected {(nside, nside)} for the simulation geometry.",
        )
    return beam.astype(np.float64, copy=False)


def _prepare_cmb_power(values: dict[str, object], cl_override: Path | None) -> tuple[Path, object]:
    cl_path = cl_override or Path(str(values["clfile"]))
    print(f" # Reading CMB power spectrum from {cl_path}")
    spectrum = load_cmb_power_spectrum(cl_path)
    return cl_path, spectrum


def _compute_cmb_power_density(
    spectrum,
    radius_squared: np.ndarray,
    templ: float,
) -> np.ndarray:
    ell_values = np.sqrt(radius_squared, dtype=np.float64) / templ
    flat = ell_values.ravel()
    cmb_rms = spectrum.interpolate_array(flat)
    zero_mask = flat == 0.0
    if np.any(zero_mask):
        cmb_rms[zero_mask] = 0.0
    return np.square(cmb_rms, dtype=np.float64).reshape(radius_squared.shape)


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
    pixsize_arcmin = sim_map.pixsize
    nside = sim_map.nsidepix
    ntot = sim_map.npix

    noise150 = float(args.noise_150 if args.noise_150 is not None else parameters["DT150"])
    noise95 = float(args.noise_95 if args.noise_95 is not None else parameters["DT90"])
    frequency150 = float(parameters["FREQUENCY150"])
    frequency95 = float(parameters["FREQUENCY90"])
    t_cmb = float(parameters["T_CMB"])
    kb = float(parameters["KB"])
    planck_h = float(parameters["H"])
    beta = float(parameters["BETA"])

    beam150_path = _resolve_beam_path(args.beam_150 or parameters["BEAM150FILE"])
    beam95_path = _resolve_beam_path(args.beam_95 or parameters["BEAM90FILE"])
    filter_prefix = args.filter_prefix or str(parameters["FILTERFILE"])
    sigma_output = Path(args.sigma_output or parameters["FILTERMSIGFILE"])

    num_filters = int(args.num_filters or parameters["NFILTER"])
    if num_filters <= 0:
        raise ValueError("The number of filters must be positive.")
    theta_start = float(parameters["THETA_START"])
    theta_step = float(parameters["THETA_STEP"])

    print(" # Importing Beam profile")
    beam150 = _load_beam(beam150_path, nside)
    beam95 = _load_beam(beam95_path, nside)

    print(" # Generating filter profile")
    fbeam150 = np.fft.fft2(beam150)
    fbeam95 = np.fft.fft2(beam95)

    beam_power_150 = np.square(fbeam150.real) + np.square(fbeam150.imag)
    beam_power_95 = np.square(fbeam95.real) + np.square(fbeam95.imag)
    b1b2real = fbeam150.real * fbeam95.real + fbeam150.imag * fbeam95.imag
    b1b2imag = -fbeam150.real * fbeam95.imag + fbeam150.imag * fbeam95.real

    indices = np.arange(nside, dtype=np.float64)
    folded = np.minimum(indices, nside - indices)
    x = folded[np.newaxis, :]
    y = folded[:, np.newaxis]
    radius_squared = x**2 + y**2

    templ = nside * pixsize_arcmin * math.pi / (180.0 * 60.0) / (2.0 * math.pi)
    if templ == 0.0:
        raise ValueError("Derived Fourier-space scaling factor is zero; check the geometry settings.")

    _, spectrum = _prepare_cmb_power(parameters, args.cl_file)
    cmb_power = _compute_cmb_power_density(spectrum, radius_squared, templ)

    noise_power_150 = convert_noise(noise150, pixsize_arcmin)
    noise_power_95 = convert_noise(noise95, pixsize_arcmin)

    tempx = planck_h * frequency150 * 1.0e9 / kb / t_cmb
    tempexpx = math.exp(tempx)
    j150 = tempx * (tempexpx + 1.0) / (tempexpx - 1.0) - 4.0
    j150 = -j150

    tempx = planck_h * frequency95 * 1.0e9 / kb / t_cmb
    tempexpx = math.exp(tempx)
    j95 = tempx * (tempexpx + 1.0) / (tempexpx - 1.0) - 4.0
    j95 = -j95

    sigma_output.parent.mkdir(parents=True, exist_ok=True)
    with sigma_output.open("w", encoding="utf-8") as sig_file:
        for k in range(num_filters):
            theta_arcmin = thetasize(k, theta_start, theta_step)
            theta_pixels = theta_arcmin / pixsize_arcmin
            print(f"# Calculating {theta_arcmin}(arcmin)/{theta_pixels}(pixel) case")
            if theta_pixels <= 0:
                raise ValueError("Theta converted to pixel units must be positive.")

            theta_sq = theta_pixels * theta_pixels
            base = 1.0 + radius_squared / theta_sq
            tau_map = np.power(base, (1.0 - 3.0 * beta) / 2.0, dtype=np.float64)

            idata = np.fft.fft2(tau_map) / ntot
            idata_real = idata.real
            idata_imag = idata.imag

            cmb_j95 = cmb_power * j95
            cmb_j150 = cmb_power * j150
            p90 = (noise_power_95 + beam_power_95 * cmb_power) * j150
            p150 = (noise_power_150 + beam_power_150 * cmb_power) * j95
            ptot = (
                noise_power_150 * noise_power_95
                + beam_power_150 * cmb_power * noise_power_95
                + beam_power_95 * cmb_power * noise_power_150
            )
            if np.any(ptot == 0.0):
                raise ValueError("Encountered zero total power while building the filter profile.")
            sigma_factor = j150 * p90 + j95 * p150 - 2.0 * cmb_power * j150 * j95 * b1b2real

            f150_real = (
                idata_real * (p90 - b1b2real * cmb_j95)
                + b1b2imag * cmb_j95 * idata_imag
            ) / ptot
            f150_imag = (
                idata_imag * p90
                - cmb_j95 * (b1b2real * idata_imag + b1b2imag * idata_real)
            ) / ptot

            f95_real = (
                idata_real * (p150 - b1b2real * cmb_j150)
                - b1b2imag * cmb_j150 * idata_imag
            ) / ptot
            f95_imag = (
                idata_imag * p150
                - cmb_j150 * (b1b2real * idata_imag - b1b2imag * idata_real)
            ) / ptot

            sigma_terms = (idata_real + np.square(idata_imag)) ** 2 * sigma_factor / ptot
            sigma_value = float(np.sum(sigma_terms))
            if sigma_value <= 0.0:
                raise ValueError("Filter normalisation sigma became non-positive.")

            scale = sigma_value * ntot
            f150 = (f150_real + 1j * f150_imag) / scale
            f95 = (f95_real + 1j * f95_imag) / scale

            sigma_output_value = math.sqrt(1.0 / sigma_value)
            sig_file.write(f"{theta_pixels:f} {sigma_output_value:e}\n")

            image150 = (np.fft.ifft2(f150) * ntot).real.astype(np.float32)
            image95 = (np.fft.ifft2(f95) * ntot).real.astype(np.float32)

            output150 = Path(f"{filter_prefix}M150_{k}")
            output95 = Path(f"{filter_prefix}M90_{k}")

            save_imagef(output150.with_suffix(".fits"), image150)
            save_fft(image150, output150)

            save_imagef(output95.with_suffix(".fits"), image95)
            save_fft(image95, output95)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
