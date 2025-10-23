"""Python translation of ``src/main_szpeakfinder.cpp``."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List

import numpy as np

from fits_utils import read_imagef
from genfilterM import thetasize
from simulation import (
    DEFAULT_PARAMETER_FILE,
    MapGeometry,
    load_map_geometry,
    parse_parameter_file,
)


@dataclass(slots=True)
class Peak:
    """Container for detected SZ peaks."""

    xpix: int
    ypix: int
    sn: float
    coresize: float


def parse_args(argv: Iterable[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Identify SZ peaks within the filtered maps produced by the pipeline.",
    )
    parser.add_argument(
        "prefix",
        nargs="?",
        help=(
            "Prefix of the filtered FITS maps. Defaults to 'filtered' which matches the "
            "naming used by the C tools."
        ),
    )
    parser.add_argument(
        "identifier",
        nargs="?",
        help=(
            "Optional identifier appended to the output list name as 'szpeaks_<ID>.dat'. "
            "If omitted the SZPEAKLIST entry from the parameter file is used."
        ),
    )
    parser.add_argument(
        "--param-file",
        type=Path,
        default=DEFAULT_PARAMETER_FILE,
        help="Parameter header that defines the geometry, filter count, and thresholds",
    )
    parser.add_argument("--ndim", type=int, help="Override the NDIM value from the parameter file")
    parser.add_argument(
        "--fieldsize",
        type=float,
        help="Override the FIELDSIZE (degrees) value from the parameter file",
    )
    parser.add_argument(
        "--extendeddim",
        type=int,
        help="Override the EXTENDEDDIM value from the parameter file",
    )
    parser.add_argument(
        "--num-filters",
        type=int,
        help="Override the NFILTER value from the parameter file",
    )
    parser.add_argument(
        "--theta-start",
        type=float,
        help="Override the THETA_START (arcmin) value from the parameter file",
    )
    parser.add_argument(
        "--theta-step",
        type=float,
        help="Override the THETA_STEP (arcmin) value from the parameter file",
    )
    parser.add_argument(
        "--sigma-threshold",
        type=float,
        help="Override the SIGMA_THRESHOD value from the parameter file",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Write the peak catalogue to this explicit filename",
    )
    parser.add_argument(
        "--image-template",
        default="{prefix}{index:02d}.fits",
        help=(
            "Template used to locate the filtered FITS maps. The placeholders {prefix} and "
            "{index} are substituted for each filter case."
        ),
    )
    return parser.parse_args(argv)


def _load_parameters(param_file: Path) -> dict[str, object]:
    values = parse_parameter_file(param_file)
    required = [
        "NFILTER",
        "THETA_START",
        "THETA_STEP",
        "SIGMA_THRESHOD",
        "SZPEAKLIST",
    ]
    missing = [name for name in required if name not in values]
    if missing:
        joined = ", ".join(missing)
        raise ValueError(f"Missing parameters in {param_file}: {joined}")
    return values


def _resolve_output(args: argparse.Namespace, parameters: dict[str, object]) -> Path:
    if args.output is not None:
        return args.output
    if args.identifier:
        return Path(f"szpeaks_{args.identifier}.dat")
    return Path(str(parameters["SZPEAKLIST"]))


def _resolve_num_filters(args: argparse.Namespace, parameters: dict[str, object]) -> int:
    count = int(args.num_filters or parameters["NFILTER"])
    if count <= 0:
        raise ValueError("The number of filters must be positive.")
    return count


def _load_image(path: Path, expected_shape: tuple[int, int]) -> np.ndarray:
    image = read_imagef(path)
    if image.shape != expected_shape:
        raise ValueError(
            f"Image {path} has shape {image.shape}, expected {expected_shape} for the simulation geometry.",
        )
    return np.asarray(image, dtype=np.float32)


def _bounds(geometry: MapGeometry) -> tuple[int, int, int]:
    extended = geometry.extendeddim
    if extended == 0:
        lowbound = 3
        highbound = geometry.ndim - 3
    else:
        lowbound = extended
        highbound = geometry.ndim + extended
    nside = geometry.ndim + 2 * extended
    if highbound <= lowbound:
        raise ValueError(
            "Invalid geometry bounds: ensure NDIM is larger than the padding so peaks can be identified.",
        )
    return lowbound, highbound, nside


def _find_filter_peaks(
    image: np.ndarray,
    geometry: MapGeometry,
    index: int,
    sigma_threshold: float,
    theta_start: float,
    theta_step: float,
) -> List[Peak]:
    lowbound, highbound, _ = _bounds(geometry)
    peaks: List[Peak] = []
    coresize = thetasize(index, theta_start, theta_step)
    extended = geometry.extendeddim

    for y in range(lowbound, highbound):
        for x in range(lowbound, highbound):
            value = float(image[y, x])
            if value <= sigma_threshold:
                continue
            if (
                value > float(image[y, x + 1])
                and value > float(image[y, x - 1])
                and value > float(image[y + 1, x])
                and value > float(image[y - 1, x])
            ):
                peaks.append(Peak(x - extended, y - extended, value, coresize))
    return peaks


def _combine_peaks(peaks: List[Peak], geometry: MapGeometry) -> None:
    if not peaks:
        return
    pixelsize = geometry.fieldsize_deg / geometry.ndim * 60.0
    if pixelsize <= 0.0:
        raise ValueError("Pixel size must be positive to combine peaks.")

    i = 0
    while i < len(peaks):
        j = i + 1
        removed_current = False
        while j < len(peaks):
            thetalimit = max(peaks[i].coresize, peaks[j].coresize, 1.0)
            pixel_limit = thetalimit / pixelsize
            limit_squared = pixel_limit * pixel_limit
            dx = peaks[i].xpix - peaks[j].xpix
            dy = peaks[i].ypix - peaks[j].ypix
            if (dx * dx + dy * dy) < limit_squared:
                if peaks[i].sn > peaks[j].sn:
                    del peaks[j]
                    continue
                else:
                    del peaks[i]
                    removed_current = True
                    if i > 0:
                        i -= 1
                    break
            j += 1
        if removed_current:
            continue
        i += 1


def _write_output(peaks: Iterable[Peak], destination: Path) -> None:
    with destination.open("w", encoding="utf-8") as handle:
        for peak in peaks:
            handle.write(f"{peak.xpix} {peak.ypix} {peak.sn:.6f} {peak.coresize:.6f}\n")


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    parameters = _load_parameters(args.param_file)

    geometry = load_map_geometry(
        args.param_file,
        ndim=args.ndim,
        fieldsize=args.fieldsize,
        extendeddim=args.extendeddim,
    )
    nside = geometry.ndim + geometry.extendeddim * 2
    expected_shape = (nside, nside)

    num_filters = _resolve_num_filters(args, parameters)
    theta_start = float(args.theta_start or parameters["THETA_START"])
    theta_step = float(args.theta_step or parameters["THETA_STEP"])
    sigma_threshold = float(args.sigma_threshold or parameters["SIGMA_THRESHOD"])

    prefix = args.prefix or "filtered"
    template = args.image_template
    if "{index" not in template:
        raise ValueError("The image template must contain '{index}' to insert the filter number.")

    peaks: List[Peak] = []
    for index in range(num_filters):
        image_name = template.format(prefix=prefix, index=index)
        image_path = Path(image_name)
        print(f" # reading {image_path}")
        image = _load_image(image_path, expected_shape)
        peaks.extend(
            _find_filter_peaks(
                image,
                geometry,
                index,
                sigma_threshold,
                theta_start,
                theta_step,
            )
        )

    print(f"total peaks {len(peaks)}")
    _combine_peaks(peaks, geometry)
    print(f"total peaks after merge {len(peaks)}")

    output_path = _resolve_output(args, parameters)
    _write_output(peaks, output_path)
    print(f"peaks written to {output_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
