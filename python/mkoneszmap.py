"""Python translation of ``src/main_mkoneszmap.c``."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

from cmb import generate_cmb_map, load_cmb_power_spectrum
from fits_utils import read_imagef, save_imagef
from nr_random import NRRandom
from simulation import (
    DEFAULT_PARAMETER_FILE,
    initialise_map,
    load_map_geometry,
    map_sze_image,
    parse_parameter_file,
)
from sz_map_utils import apply_beam_and_noise, load_beam_fft, mkbolo, remove_mean, validate_sze_image


def _load_parameters(param_file: Path) -> dict[str, object]:
    values = parse_parameter_file(param_file)
    required = ["T_CMB", "KB", "H", "clfile"]
    missing = [name for name in required if name not in values]
    if missing:
        joined = ", ".join(missing)
        raise ValueError(f"Missing parameters in {param_file}: {joined}")
    return values


def parse_args(argv: Iterable[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate a single-frequency mock SZ map.")
    parser.add_argument("rseed", type=int, help="Random seed used for the CMB realisation")
    parser.add_argument("frequency", type=float, help="Band centre in GHz")
    parser.add_argument("white_noise", type=float, help="White noise level in K/arcmin")
    parser.add_argument("beam", help="Beam name (expects a '<beam>_fft.dat' companion file)")
    parser.add_argument("y_map", help="Input Compton-y map")
    parser.add_argument("output", help="Output FITS filename")
    parser.add_argument("noise_seed", type=int, help="Random seed dedicated to the instrumental noise")
    parser.add_argument(
        "--param-file",
        type=Path,
        default=DEFAULT_PARAMETER_FILE,
        help="Parameter header containing the geometry and physical constants",
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
        "--cmb-output",
        type=Path,
        default=Path("cmb.fits"),
        help="Filename used for the intermediate CMB map",
    )
    parser.add_argument(
        "--noise-output",
        type=Path,
        default=Path("noise.fits"),
        help="Filename used for the generated white noise map",
    )
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    parameters = _load_parameters(args.param_file)
    geometry = load_map_geometry(
        args.param_file,
        ndim=args.ndim,
        fieldsize=args.fieldsize,
        extendeddim=args.extendeddim,
    )

    y_image = validate_sze_image(read_imagef(args.y_map), geometry)
    sim_map = initialise_map(geometry)
    map_sze_image(sim_map, y_image)

    cl_path = args.cl_file or Path(str(parameters["clfile"]))
    spectrum = load_cmb_power_spectrum(cl_path)

    cmb_rng = NRRandom(args.rseed)
    noise_rng = NRRandom(args.noise_seed)

    print(f"# Read the SZE image {args.y_map}")
    print("# Input CMB")
    cmb_map = generate_cmb_map(sim_map.nsidepix, sim_map.pixsize, spectrum, cmb_rng, float(parameters["T_CMB"]))
    save_imagef(args.cmb_output, cmb_map)

    print(f"# Bolometry @ {args.frequency} GHz")
    noise = mkbolo(
        sim_map,
        args.frequency,
        cmb_map,
        args.white_noise,
        noise_rng,
        kb=float(parameters["KB"]),
        planck_h=float(parameters["H"]),
    )
    save_imagef(args.noise_output, noise)

    print("# Adding beam and white noise")
    beam_fft = load_beam_fft(args.beam, sim_map.nsidepix)
    apply_beam_and_noise(sim_map, noise, beam_fft)

    print("# Subtracting mean temperature")
    remove_mean(sim_map)

    print(f"# Output to {args.output}")
    save_imagef(args.output, sim_map.image)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
