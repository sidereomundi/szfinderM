"""Numerical Recipes compatible random deviates used by the C pipeline."""

from __future__ import annotations

import math
from typing import Iterable

import numpy as np


class NRRandom:
    """Reproduce ``ran1``/``gasdev`` from Numerical Recipes."""

    # Constants copied from Numerical Recipes ``ran1`` implementation.
    _IM1 = 2147483563
    _IM2 = 2147483399
    _AM = 1.0 / _IM1
    _IMM1 = _IM1 - 1
    _IA1 = 40014
    _IA2 = 40692
    _IQ1 = 53668
    _IQ2 = 52774
    _IR1 = 12211
    _IR2 = 3791
    _NTAB = 32
    _NDIV = 1 + _IMM1 // _NTAB
    _EPS = 1.2e-7
    _RNMX = 1.0 - _EPS

    def __init__(self, seed: int) -> None:
        if seed == 0:
            seed = 1
        self._idum = -abs(seed)
        self._idum2 = 123456789
        self._iy = 0
        self._iv = [0] * self._NTAB
        self._stored_gaussian: float | None = None

    def _ran1(self) -> float:
        if self._idum <= 0:
            if -self._idum < 1:
                self._idum = 1
            else:
                self._idum = -self._idum
            self._idum2 = 123456789
            for j in range(self._NTAB + 7, -1, -1):
                k = self._idum // self._IQ1
                self._idum = self._IA1 * (self._idum - k * self._IQ1) - k * self._IR1
                if self._idum < 0:
                    self._idum += self._IM1
                if j < self._NTAB:
                    self._iv[j] = self._idum
            self._iy = self._iv[0]
        k = self._idum // self._IQ1
        self._idum = self._IA1 * (self._idum - k * self._IQ1) - k * self._IR1
        if self._idum < 0:
            self._idum += self._IM1
        k = self._idum2 // self._IQ2
        self._idum2 = self._IA2 * (self._idum2 - k * self._IQ2) - k * self._IR2
        if self._idum2 < 0:
            self._idum2 += self._IM2
        j = self._iy // self._NDIV
        self._iy = self._iv[j] - self._idum2
        self._iv[j] = self._idum
        if self._iy < 1:
            self._iy += self._IMM1
        temp = self._AM * self._iy
        return temp if temp < self._RNMX else self._RNMX

    def gasdev(self) -> float:
        """Return a Gaussian deviate with unit variance."""

        if self._stored_gaussian is not None:
            value = self._stored_gaussian
            self._stored_gaussian = None
            return value

        while True:
            v1 = 2.0 * self._ran1() - 1.0
            v2 = 2.0 * self._ran1() - 1.0
            rsq = v1 * v1 + v2 * v2
            if rsq >= 1.0 or rsq == 0.0:
                continue
            fac = math.sqrt(-2.0 * math.log(rsq) / rsq)
            self._stored_gaussian = v1 * fac
            return v2 * fac

    def gasdev_array(self, shape: Iterable[int]) -> np.ndarray:
        """Return an array of Gaussian deviates with unit variance."""

        total = 1
        for size in shape:
            total *= int(size)
        values = np.empty(total, dtype=np.float64)
        for idx in range(total):
            values[idx] = self.gasdev()
        return values.reshape(tuple(int(size) for size in shape))
