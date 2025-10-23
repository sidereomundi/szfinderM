"""Shared helpers for the SZ map generation scripts."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from fft_io import load_fft
from nr_random import NRRandom
from simulation import MapGeometry, SimulationMap


def validate_sze_image(image: np.ndarray, geometry: MapGeometry) -> np.ndarray:
    """Ensure the loaded SZE map matches the expected geometry."""

    array = np.asarray(image, dtype=np.float32)
    if array.shape != (geometry.ndim, geometry.ndim):
        msg = f"Expected an SZE map of shape {(geometry.ndim, geometry.ndim)}, found {array.shape}."
        raise ValueError(msg)
    return array


def resolve_beam_fft_path(name: str | Path) -> Path:
    """Mirror the ``sprintf("%s_fft.dat")`` pattern from the C binaries."""

    base = Path(name)
    if base.suffix.lower() == ".dat" and base.name.endswith("_fft.dat"):
        return base
    if base.suffix:
        return base
    return base.with_name(f"{base.name}_fft.dat")


def load_beam_fft(name: str | Path, nsidepix: int) -> np.ndarray:
    """Load the FFT beam associated with ``name``."""

    path = resolve_beam_fft_path(name)
    return load_fft(path, nsidepix)


def mkbolo(
    sim_map: SimulationMap,
    frequency_ghz: float,
    cmb_map: np.ndarray,
    white_noise: float,
    rng: NRRandom,
    *,
    kb: float,
    planck_h: float,
) -> np.ndarray:
    """Replicate ``mkbolo`` from ``sub_mkbolo.c`` using NumPy."""

    freq_hz = frequency_ghz * 1.0e9
    tsig_per_pixel = white_noise * 1.0e-6 / sim_map.pixsize
    hvk = planck_h * freq_hz / kb

    temperature = np.asarray(cmb_map, dtype=np.float64)
    x = hvk / temperature
    expx = np.exp(x)
    denom = np.where(np.abs(expx - 1.0) < 1.0e-12, 1.0e-12, expx - 1.0)
    ftemp = x * (expx + 1.0) / denom - 4.0

    sim_image = sim_map.image.astype(np.float64, copy=False)
    updated = (sim_image * ftemp + 1.0) * temperature
    sim_map.image[:, :] = updated.astype(np.float32)

    noise_values = rng.gasdev_array(sim_map.image.shape).astype(np.float64, copy=False)
    noise_values *= tsig_per_pixel
    return noise_values.astype(np.float32)


def apply_beam_and_noise(
    sim_map: SimulationMap, noise: np.ndarray, beam_fft: np.ndarray
) -> None:
    """Convolve the map with ``beam_fft`` and add ``noise``."""

    if beam_fft.shape != (sim_map.nsidepix, sim_map.nsidepix // 2 + 1):
        msg = (
            f"Beam FFT shape {beam_fft.shape} does not match the simulation map "
            f"dimensions ({sim_map.nsidepix}, {sim_map.nsidepix // 2 + 1})."
        )
        raise ValueError(msg)

    image_fft = np.fft.rfft2(sim_map.image.astype(np.float64, copy=False))
    dirty = np.fft.irfft2(image_fft * beam_fft)
    sim_map.image[:, :] = (dirty.real + noise).astype(np.float32)


def remove_mean(sim_map: SimulationMap) -> None:
    """Subtract the mean temperature from the simulated map."""

    mean_value = float(np.mean(sim_map.image, dtype=np.float64))
    sim_map.image[:, :] = (sim_map.image - mean_value).astype(np.float32)
