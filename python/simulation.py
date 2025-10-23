"""Core data structures and geometry helpers for the SZ finder pipeline."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
from typing import Any

import numpy as np

DEFAULT_PARAMETER_FILE = Path(__file__).resolve().parents[1] / "include" / "parameters.h"


@dataclass(frozen=True, slots=True)
class MapGeometry:
    """Geometry parameters that describe the simulation grid."""

    ndim: int
    fieldsize_deg: float
    extendeddim: int = 0


@dataclass(slots=True)
class SimulationMap:
    """Python representation of the ``struct simmap`` used in the C code."""

    geometry: MapGeometry
    pixsize: float
    nsidepix: int
    npix: int
    image: np.ndarray


def parse_parameter_file(path: Path) -> dict[str, Any]:
    """Parse a C-style parameter header and return the macro definitions."""

    try:
        content = path.read_text(encoding="utf-8")
    except FileNotFoundError as exc:  # pragma: no cover - depends on environment.
        msg = f"Unable to open parameter file: {path}"
        raise FileNotFoundError(msg) from exc

    pattern = re.compile(r"^#define\s+([A-Za-z0-9_]+)\s+([^/]+?)\s*$")
    values: dict[str, Any] = {}
    for raw_line in content.splitlines():
        line = raw_line.split("//", 1)[0].strip()
        if not line:
            continue
        match = pattern.match(line)
        if not match:
            continue
        name, raw_value = match.groups()
        value = raw_value.strip()
        if value.startswith('"') and value.endswith('"'):
            values[name] = value.strip('"')
            continue
        try:
            values[name] = int(value, 0)
            continue
        except ValueError:
            pass
        try:
            values[name] = float(value)
        except ValueError:
            continue
    return values


def _lookup_parameter(
    values: dict[str, Any], name: str, override: Any | None, *, default: Any | None = None
) -> Any:
    if override is not None:
        return override
    if name in values:
        return values[name]
    if default is not None:
        return default
    msg = f"Parameter '{name}' was not found and no override was supplied."
    raise ValueError(msg)


def load_map_geometry(
    param_file: Path,
    *,
    ndim: int | None = None,
    fieldsize: float | None = None,
    extendeddim: int | None = None,
) -> MapGeometry:
    """Load geometry settings from ``param_file`` with optional overrides."""

    try:
        file_values = parse_parameter_file(param_file)
    except FileNotFoundError:
        file_values = {}
        if ndim is None or fieldsize is None:
            raise

    resolved_ndim = int(_lookup_parameter(file_values, "NDIM", ndim))
    resolved_fieldsize = float(_lookup_parameter(file_values, "FIELDSIZE", fieldsize))
    resolved_extendeddim = int(_lookup_parameter(file_values, "EXTENDEDDIM", extendeddim, default=0))
    return MapGeometry(resolved_ndim, resolved_fieldsize, resolved_extendeddim)


def initialise_map(geometry: MapGeometry) -> SimulationMap:
    """Construct a :class:`SimulationMap` for the provided geometry."""

    pixsize = geometry.fieldsize_deg / geometry.ndim * 60.0  # Convert degrees to arcminutes.
    nsidepix = geometry.ndim + geometry.extendeddim * 2
    image = np.zeros((nsidepix, nsidepix), dtype=np.float32)
    return SimulationMap(
        geometry=geometry,
        pixsize=pixsize,
        nsidepix=nsidepix,
        npix=nsidepix**2,
        image=image,
    )


def map_sze_image(sim_map: SimulationMap, sze_image: np.ndarray) -> None:
    """Embed a Compton ``y`` map into the simulation map with optional padding."""

    geom = sim_map.geometry
    if sze_image.shape != (geom.ndim, geom.ndim):
        raise ValueError(
            f"Expected an input map of shape {(geom.ndim, geom.ndim)}, got {sze_image.shape}."
        )

    mapped = np.zeros_like(sim_map.image)
    extended = geom.extendeddim
    mapped[
        extended : extended + geom.ndim,
        extended : extended + geom.ndim,
    ] = sze_image.astype(np.float32, copy=False)
    sim_map.image[:, :] = mapped
