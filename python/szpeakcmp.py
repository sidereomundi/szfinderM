"""Python translation of ``src/main_szpeakcmp.cpp``."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List

from simulation import (
    DEFAULT_PARAMETER_FILE,
    MapGeometry,
    load_map_geometry,
    parse_parameter_file,
)


@dataclass(slots=True)
class Cluster:
    """Representation of a simulated cluster entry."""

    xpos: float
    ypos: float
    mass: float
    redshift: float
    identifier: int


@dataclass(slots=True)
class Peak:
    """Entry within the detected SZ peak catalogue."""

    xpix: int
    ypix: int
    sn: float
    coresize: float


def parse_args(argv: Iterable[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare detected SZ peaks with the underlying cluster catalogue.",
    )
    parser.add_argument(
        "identifier",
        nargs="?",
        help=(
            "Optional identifier appended to the default input/output filenames. If "
            "provided the script reads 'clusterlist_<ID>.dat' and 'szpeaks_<ID>.dat' "
            "and writes 'matchpeak_<ID>.dat'."
        ),
    )
    parser.add_argument(
        "--cluster-file",
        type=Path,
        help="Explicit path to the formatted cluster catalogue",
    )
    parser.add_argument(
        "--peaks-file",
        type=Path,
        help="Explicit path to the detected peaks file",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Destination filename for the match catalogue",
    )
    parser.add_argument(
        "--param-file",
        type=Path,
        default=DEFAULT_PARAMETER_FILE,
        help="Parameter header defining NDIM, FIELDSIZE, and related settings",
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


def _load_parameters(param_file: Path) -> dict[str, object]:
    values = parse_parameter_file(param_file)
    required = ["CLUSTERLIST", "SZPEAKLIST"]
    missing = [name for name in required if name not in values]
    if missing:
        joined = ", ".join(missing)
        raise ValueError(f"Missing parameters in {param_file}: {joined}")
    return values


def _resolve_cluster_path(args: argparse.Namespace, parameters: dict[str, object]) -> Path:
    if args.cluster_file is not None:
        return args.cluster_file
    if args.identifier:
        return Path(f"clusterlist_{args.identifier}.dat")
    return Path(str(parameters["CLUSTERLIST"]))


def _resolve_peaks_path(args: argparse.Namespace, parameters: dict[str, object]) -> Path:
    if args.peaks_file is not None:
        return args.peaks_file
    if args.identifier:
        return Path(f"szpeaks_{args.identifier}.dat")
    return Path(str(parameters["SZPEAKLIST"]))


def _resolve_output_path(args: argparse.Namespace, parameters: dict[str, object]) -> Path:
    if args.output is not None:
        return args.output
    if args.identifier:
        return Path(f"matchpeak_{args.identifier}.dat")
    return Path("matchpeak.dat")


def _load_clusters(path: Path, geometry: MapGeometry) -> List[Cluster]:
    clusters: List[Cluster] = []
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except FileNotFoundError as exc:
        msg = f"Unable to open cluster catalogue: {path}"
        raise FileNotFoundError(msg) from exc

    for raw in lines:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 5:
            continue
        identifier = int(parts[0])
        xpos = float(parts[1]) * geometry.ndim
        ypos = float(parts[2]) * geometry.ndim
        mass = float(parts[3])
        redshift = float(parts[4])
        clusters.append(Cluster(xpos, ypos, mass, redshift, identifier))
    return clusters


def _load_peaks(path: Path) -> List[Peak]:
    peaks: List[Peak] = []
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except FileNotFoundError as exc:
        msg = f"Unable to open peak catalogue: {path}"
        raise FileNotFoundError(msg) from exc

    for raw in lines:
        line = raw.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) < 4:
            continue
        xpix = int(parts[0])
        ypix = int(parts[1])
        sn = float(parts[2])
        coresize = float(parts[3])
        peaks.append(Peak(xpix, ypix, sn, coresize))
    return peaks


def _match_peaks(
    clusters: List[Cluster],
    peaks: List[Peak],
    geometry: MapGeometry,
) -> list[tuple[Peak, Cluster | None, int]]:
    pixelsize = geometry.fieldsize_deg / geometry.ndim * 60.0
    if pixelsize <= 0.0:
        raise ValueError("Pixel size must be positive to perform the matching.")

    matches: list[tuple[Peak, Cluster | None, int]] = []
    for peak in peaks:
        match_count = 0
        thetalimit = max(1.0, peak.coresize) / pixelsize
        limit_squared = thetalimit * thetalimit

        for cluster in clusters:
            dx = peak.xpix - cluster.xpos
            dy = peak.ypix - cluster.ypos
            if dx * dx + dy * dy <= limit_squared:
                match_count += 1
                matches.append((peak, cluster, match_count))

        if match_count == 0:
            matches.append((peak, None, 0))

    return matches


def _write_matches(path: Path, matches: Iterable[tuple[Peak, Cluster | None, int]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for peak, cluster, count in matches:
            if cluster is None:
                handle.write(
                    f"{peak.xpix} {peak.ypix} 0 0 {peak.sn:.6f} {peak.coresize:.6f} 0 0 0 0\n"
                )
            else:
                handle.write(
                    (
                        f"{peak.xpix} {peak.ypix} {cluster.xpos:.6f} {cluster.ypos:.6f} "
                        f"{peak.sn:.6f} {peak.coresize:.6f} {cluster.mass:.6g} "
                        f"{cluster.redshift:.6f} {cluster.identifier} {count}\n"
                    )
                )


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    parameters = _load_parameters(args.param_file)
    geometry = load_map_geometry(
        args.param_file,
        ndim=args.ndim,
        fieldsize=args.fieldsize,
        extendeddim=args.extendeddim,
    )

    cluster_path = _resolve_cluster_path(args, parameters)
    peaks_path = _resolve_peaks_path(args, parameters)
    output_path = _resolve_output_path(args, parameters)

    clusters = _load_clusters(cluster_path, geometry)
    peaks = _load_peaks(peaks_path)

    matches = _match_peaks(clusters, peaks, geometry)
    _write_matches(output_path, matches)
    print(f"wrote {len(matches)} entries to {output_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
