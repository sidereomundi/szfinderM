"""CMB power spectrum helpers and map synthesis."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from nr_random import NRRandom


@dataclass(slots=True)
class CMBPowerSpectrum:
    """Cubic-spline representation of the input CMB power spectrum."""

    ell: np.ndarray
    delta_t: np.ndarray
    second_derivatives: np.ndarray

    def interpolate(self, ell_value: float) -> float:
        """Evaluate the spline at the requested multipole."""

        if ell_value <= self.ell[0]:
            return float(self.delta_t[0])
        if ell_value >= self.ell[-1]:
            return float(self.delta_t[-1])

        klo = 0
        khi = self.ell.size - 1
        while khi - klo > 1:
            k = (khi + klo) // 2
            if self.ell[k] > ell_value:
                khi = k
            else:
                klo = k
        h = self.ell[khi] - self.ell[klo]
        if h == 0:
            raise ValueError("CMB power spectrum contains duplicate ell values.")
        a = (self.ell[khi] - ell_value) / h
        b = (ell_value - self.ell[klo]) / h
        term = (
            (a**3 - a) * self.second_derivatives[klo]
            + (b**3 - b) * self.second_derivatives[khi]
        ) * (h * h) / 6.0
        return float(a * self.delta_t[klo] + b * self.delta_t[khi] + term)

    def interpolate_array(self, ell_values: np.ndarray) -> np.ndarray:
        """Vectorised spline evaluation for an array of multipoles."""

        ell_array = np.asarray(ell_values, dtype=np.float64)
        result = np.empty_like(ell_array)

        low_mask = ell_array <= self.ell[0]
        high_mask = ell_array >= self.ell[-1]
        mid_mask = ~(low_mask | high_mask)

        result[low_mask] = self.delta_t[0]
        result[high_mask] = self.delta_t[-1]

        if np.any(mid_mask):
            mid_values = ell_array[mid_mask]
            indices = np.searchsorted(self.ell, mid_values, side="right")
            klo = indices - 1
            khi = indices
            h = self.ell[khi] - self.ell[klo]
            if np.any(h == 0):
                raise ValueError("CMB power spectrum contains duplicate ell values.")
            a = (self.ell[khi] - mid_values) / h
            b = (mid_values - self.ell[klo]) / h
            term = (
                (a**3 - a) * self.second_derivatives[klo]
                + (b**3 - b) * self.second_derivatives[khi]
            ) * (h * h) / 6.0
            result[mid_mask] = a * self.delta_t[klo] + b * self.delta_t[khi] + term

        return result


def _compute_spline_second_derivatives(
    ell: np.ndarray, values: np.ndarray, yp1: float, ypn: float
) -> np.ndarray:
    n = ell.size
    y2 = np.zeros(n, dtype=np.float64)
    u = np.zeros(n - 1, dtype=np.float64)

    if yp1 > 0.99e30:
        y2[0] = 0.0
    else:
        y2[0] = -0.5
        u[0] = (
            (3.0 / (ell[1] - ell[0]))
            * ((values[1] - values[0]) / (ell[1] - ell[0]) - yp1)
        )

    for i in range(1, n - 1):
        sig = (ell[i] - ell[i - 1]) / (ell[i + 1] - ell[i - 1])
        p = sig * y2[i - 1] + 2.0
        y2[i] = (sig - 1.0) / p
        u[i] = (
            (values[i + 1] - values[i]) / (ell[i + 1] - ell[i])
            - (values[i] - values[i - 1]) / (ell[i] - ell[i - 1])
        )
        u[i] = (6.0 * u[i] / (ell[i + 1] - ell[i - 1]) - sig * u[i - 1]) / p

    if ypn > 0.99e30:
        qn = 0.0
        un = 0.0
    else:
        qn = 0.5
        un = (
            (3.0 / (ell[-1] - ell[-2]))
            * (ypn - (values[-1] - values[-2]) / (ell[-1] - ell[-2]))
        )
    y2[-1] = (un - qn * u[-2]) / (qn * y2[-2] + 1.0)

    for k in range(n - 2, -1, -1):
        y2[k] = y2[k] * y2[k + 1] + u[k]
    return y2


def load_cmb_power_spectrum(path: Path) -> CMBPowerSpectrum:
    """Load the CAMB-format power spectrum used by the original pipeline."""

    data = np.loadtxt(path, dtype=np.float64)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError("CMB power spectrum must contain at least two columns.")

    ell = data[:, 0]
    # Convert the CAMB output to the RMS temperature per mode, mirroring ``clinput``.
    cl = data[:, 1] 
    cl = cl * (2.0 * math.pi / (ell * (ell + 1.0)))
    cl = np.sqrt(np.clip(cl * 1.0e-12, a_min=0.0, a_max=None))
    cl = cl*1.e6
    
    yp1 = (cl[1] - cl[0]) / (ell[1] - ell[0])
    ypn = (cl[-1] - cl[-2]) / (ell[-1] - ell[-2])
    second_derivatives = _compute_spline_second_derivatives(ell, cl, yp1, ypn)
    return CMBPowerSpectrum(ell=ell, delta_t=cl, second_derivatives=second_derivatives)


def generate_cmb_map(
    sim_map_nside: int,
    pixsize_arcmin: float,
    spectrum: CMBPowerSpectrum,
    rng: NRRandom,
    t_cmb: float,
) -> np.ndarray:
    """Synthesize a CMB temperature map using the provided power spectrum."""

    imsize = sim_map_nside
    pixsize_rad = pixsize_arcmin * math.pi / (180.0 * 60.0)
    total_scale = imsize * pixsize_rad
    if total_scale == 0.0:
        raise ValueError("Pixel size must be positive to generate a CMB map.")

    fourier = np.zeros((imsize, imsize), dtype=np.complex128)
    half = imsize // 2

    for u in range(half + 1):
        if u != 0 and u != half:
            v_values = range(-half + 1, half + 1)
        else:
            v_values = range(0, half + 1)
        for v in v_values:
            ell_value = math.hypot(u, v) / total_scale * 2.0 * math.pi 
            d_t = spectrum.interpolate(ell_value)
            scale = d_t / math.sqrt(2.0) / imsize / pixsize_rad
            amplitude_r = scale * rng.gasdev()
            amplitude_i = scale * rng.gasdev()

            x = u % imsize
            y = v % imsize
            value = amplitude_r + 1j * amplitude_i
            fourier[y, x] = value

            x2 = (-u) % imsize
            y2 = (-v) % imsize
            fourier[y2, x2] = amplitude_r - 1j * amplitude_i

    fourier[0, 0] = 0.0
    if imsize % 2 == 0:
        fourier[half, 0] = fourier[half, 0].real + 0.0j
        fourier[0, half] = fourier[0, half].real + 0.0j
        fourier[half, half] = fourier[half, half].real + 0.0j

    cmb = np.fft.ifft2(fourier).real + t_cmb
    return np.asarray(cmb, dtype=np.float32)
