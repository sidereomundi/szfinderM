"""Python translation of ``src/test.c``.

The original C program acted as a simple integration test-bed for various
utility routines in the SZ finder project.  Most of those tests were
commented-out and the only active code path generated a zero-valued image and
wrote it to a FITS file.  This module mirrors that behaviour in Python while
keeping the overall structure of the original program.
"""

from __future__ import annotations

from typing import Iterable

import numpy as np

from fits_utils import save_imagef

# ``NDIM`` is defined in ``include/parameters.h`` in the original code base.
# The value represents the number of pixels on one side of the generated map.
NDIM: int = 4800


def generate_test_image(ndim: int) -> np.ndarray:
    """Generate the zero-valued image used by the C test harness."""

    image = np.zeros((ndim, ndim), dtype=np.float32)

    # The C code contained a commented-out block that would draw a few bright
    # pixels in three consecutive rows.  Keeping the code here (disabled by
    # default) documents the original behaviour while making it easy to
    # re-enable if desired.
    #
    # for column in range(15, 30):
    #     image[1, column] = 1e-4
    #     image[2, column] = 1e-4
    #     image[3, column] = 1e-4

    return image


def main(_: Iterable[str] | None = None) -> int:
    """Entry point matching the signature of ``int main(void)`` in C."""

    image = generate_test_image(NDIM)
    save_imagef("sznull.fits", image)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
