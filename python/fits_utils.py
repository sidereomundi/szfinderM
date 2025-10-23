"""Shared FITS helpers for the Python translations."""

from __future__ import annotations

from pathlib import Path

import numpy as np

try:  # ``astropy`` is the easiest way to emit FITS files from Python.
    from astropy.io import fits
except ImportError as exc:  # pragma: no cover - depends on the runtime env.
    raise SystemExit(
        "The Python translations require the `astropy` package to write FITS "
        "files. Install it via `pip install astropy`."
    ) from exc


def save_imagef(name: str | Path, data: np.ndarray) -> None:
    """Replicate ``saveimagef`` from ``sub_stdroutines.c`` in Python."""

    path = Path(name)
    if path.suffix.lower() != ".fits":
        msg = f"Only able to write FITS images, current name: {path.name}"
        raise ValueError(msg)

    array32 = np.asarray(data, dtype=np.float32)
    fits.PrimaryHDU(array32).writeto(path, overwrite=True)


def read_imagef(path: str | Path) -> np.ndarray:
    """Replicate ``readimage`` from ``sub_stdroutines.c`` in Python."""

    location = Path(path)
    if location.suffix.lower() != ".fits":
        msg = f"Only able to read FITS images, current name: {location.name}"
        raise ValueError(msg)

    data = fits.getdata(location)
    if data.ndim != 2:
        msg = f"Expected a 2D FITS image, found {data.ndim} dimensions"
        raise ValueError(msg)
    return np.asarray(data, dtype=np.float32)
