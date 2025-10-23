"""Python translation of ``src/main_mkszmap.c``."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import numpy as np

from cmb import generate_cmb_map, load_cmb_power_spectrum
from fits_utils import read_imagef, save_imagef
from nr_random import NRRandom
from simulation import (
    DEFAULT_PARAMETER_FILE,
    MapGeometry,
    SimulationMap,
    initialise_map,
    load_map_geometry,
    map_sze_image,
    parse_parameter_file,
)
from sz_map_utils import (
    apply_beam_and_noise,
    load_beam_fft,
    mkbolo,
    remove_mean,
    validate_sze_image,
)


def _load_parameters(param_file: Path) -> dict[str, object]:
    values = parse_parameter_file(param_file)
    required = [
        "FREQUENCY150",
        "DT150",
        "BEAM150FILE",
        "FREQUENCY90",
        "DT90",
        "BEAM90FILE",
        "clfile",
        "T_CMB",
        "KB",
        "H",
    ]
    missing = [name for name in required if name not in values]
    if missing:
        joined = ", ".join(missing)
        raise ValueError(f"Missing parameters in {param_file}: {joined}")
    return values


def _prepare_simulation_map(geometry: MapGeometry, y_map: np.ndarray) -> SimulationMap:
    sim_map = initialise_map(geometry)
    map_sze_image(sim_map, y_map)
    return sim_map


def _default_cl_path(values: dict[str, object]) -> Path:
    return Path(str(values["clfile"]))


def parse_args(argv: Iterable[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate mock SZ maps at 150 and 95 GHz.")
    parser.add_argument("rseed", type=int, help="Random seed used for the CMB and noise generation")
    parser.add_argument(
        "y_map",
        nargs="?",
        default="sznull.fits",
        help="Input Compton-y map (defaults to sznull.fits)",
    )
    parser.add_argument(
        "--param-file",
        type=Path,
        default=DEFAULT_PARAMETER_FILE,
        help="Parameter header containing frequencies, noise levels, and geometry",
    )
    parser.add_argument(
        "--cl-file",
        type=Path,
        help="Override the CMB power spectrum (defaults to the clfile entry in the parameter file)",
    )
    parser.add_argument("--output-150", type=Path, default=Path("spt150.fits"))
    parser.add_argument("--output-95", type=Path, default=Path("spt95.fits"))
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
    parameters = _load_parameters(args.param_file)
    geometry = load_map_geometry(
        args.param_file,
        ndim=args.ndim,
        fieldsize=args.fieldsize,
        extendeddim=args.extendeddim,
    )

    y_image = validate_sze_image(read_imagef(args.y_map), geometry)
    sim_map = _prepare_simulation_map(geometry, y_image)

    cl_path = args.cl_file or _default_cl_path(parameters)
    spectrum = load_cmb_power_spectrum(cl_path)

    rng = NRRandom(args.rseed)
    print(f"# Read the SZE image {args.y_map}")
    cmb_map = generate_cmb_map(sim_map.nsidepix, sim_map.pixsize, spectrum, rng, float(parameters["T_CMB"]))

    print("# Input CMB")
    print("# Bolometry")
    noise150 = mkbolo(
        sim_map,
        float(parameters["FREQUENCY150"]),
        cmb_map,
        float(parameters["DT150"]),
        rng,
        kb=float(parameters["KB"]),
        planck_h=float(parameters["H"]),
    )

    beam150 = load_beam_fft(parameters["BEAM150FILE"], sim_map.nsidepix)
    print("# Adding beam and white noise")
    apply_beam_and_noise(sim_map, noise150, beam150)
    remove_mean(sim_map)
    print("# Output")
    save_imagef(args.output_150, sim_map.image)

    print("Addition for 90 GHz")
    map_sze_image(sim_map, y_image)
    noise95 = mkbolo(
        sim_map,
        float(parameters["FREQUENCY90"]),
        cmb_map,
        float(parameters["DT90"]),
        rng,
        kb=float(parameters["KB"]),
        planck_h=float(parameters["H"]),
    )

    beam95 = load_beam_fft(parameters["BEAM90FILE"], sim_map.nsidepix)
    apply_beam_and_noise(sim_map, noise95, beam95)
    remove_mean(sim_map)
    save_imagef(args.output_95, sim_map.image)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
