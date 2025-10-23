"""High-level example for running the SZ finder Python pipeline.

This script demonstrates how the translated Python tools can be chained to
replicate the original C workflow.  It encodes the per-frequency inputs as
arrays so you can see, at a glance, which parameters are used at each stage of
pipeline execution.

By default the script **prints** the exact commands that would be executed.  Use
``--execute`` to run them instead.  Running the full pipeline requires that the
supporting data files from the original project (``include/parameters.h``,
``cl.dat``, cluster catalogue, etc.) are present in the working directory.

Example usage (from the repository root)::

    PYTHONPATH=python python python/examples/pipeline_example.py
    PYTHONPATH=python python python/examples/pipeline_example.py --execute --rseed 1234

The example assumes two observing bands (150 and 95 GHz) that match the default
``parameters.h`` file.  Adapt the ``BANDS`` configuration below for alternate
survey configurations.
"""

from __future__ import annotations

import argparse
import shlex
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, List, Sequence

REPO_ROOT = Path(__file__).resolve().parents[2]
PYTHON_DIR = REPO_ROOT / "python"
if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from genfilterM import main as genfilter_main  # type: ignore
from filtermapM import main as filtermap_main  # type: ignore
from mkbeam import main as mkbeam_main  # type: ignore
from mknull import main as mknull_main  # type: ignore
from mkoneszmap import main as mkoneszmap_main  # type: ignore
from mkszmap import main as mkszmap_main  # type: ignore
from simulation import DEFAULT_PARAMETER_FILE
from szpeakcmp import main as szpeakcmp_main  # type: ignore
from szpeakfinder import main as szpeakfinder_main  # type: ignore


@dataclass(frozen=True)
class BandConfig:
    """Inputs for a single observing band used throughout the pipeline."""

    label: str
    frequency_ghz: float
    fwhm_arcmin: float
    normalization: float
    white_noise_k_arcmin: float
    beam_output: str
    map_output: Path


BANDS: Sequence[BandConfig] = (
    BandConfig(
        label="150",
        frequency_ghz=154.477,
        fwhm_arcmin=1.1,
        normalization=1.0,
        white_noise_k_arcmin=18.0,
        beam_output="beam150T",
        map_output=Path("szmock150.fits"),
    ),
    BandConfig(
        label="95",
        frequency_ghz=97.6473,
        fwhm_arcmin=1.6,
        normalization=1.0,
        white_noise_k_arcmin=44.0,
        beam_output="beam95T",
        map_output=Path("szmock95.fits"),
    ),
)


def _format_command(program: str, args: Iterable[str]) -> str:
    parts = [program, *args]
    return " ".join(shlex.quote(part) for part in parts)


def _run_stage(
    description: str,
    func: Callable[[Iterable[str]], int],
    args: List[str],
    *,
    execute: bool,
) -> None:
    command = _format_command(description, args)
    print(f"\n# {description}")
    print(f"$ {command}")
    if execute:
        status = func(args)
        if status != 0:
            raise SystemExit(f"Stage '{description}' failed with status {status}")


def _mkbeam_args(param_file: Path) -> List[List[str]]:
    return [
        [
            f"{band.fwhm_arcmin}",
            band.beam_output,
            f"{band.normalization}",
            "--param-file",
            str(param_file),
        ]
        for band in BANDS
    ]


def _mkszmap_args(param_file: Path, rseed: int, y_map: Path, output_150: Path, output_95: Path) -> List[str]:
    return [
        str(rseed),
        str(y_map),
        "--param-file",
        str(param_file),
        "--output-150",
        str(output_150),
        "--output-95",
        str(output_95),
    ]


def _mkoneszmap_args(
    param_file: Path,
    rseed: int,
    noise_seed: int,
    band: BandConfig,
    y_map: Path,
) -> List[str]:
    return [
        str(rseed),
        f"{band.frequency_ghz}",
        f"{band.white_noise_k_arcmin}",
        band.beam_output,
        str(y_map),
        str(band.map_output),
        str(noise_seed),
        "--param-file",
        str(param_file),
    ]


def _filtermap_args(param_file: Path, map150: Path, map95: Path) -> List[str]:
    return [
        "--map-150",
        str(map150),
        "--map-95",
        str(map95),
        "--param-file",
        str(param_file),
    ]


def _szpeakfinder_args(param_file: Path, prefix: str) -> List[str]:
    return [prefix, "--param-file", str(param_file)]


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Demonstrate the SZ finder pipeline")
    parser.add_argument(
        "--param-file",
        type=Path,
        default=DEFAULT_PARAMETER_FILE,
        help="Parameter header that defines the survey geometry",
    )
    parser.add_argument(
        "--y-map",
        type=Path,
        default=Path("sznull.fits"),
        help="Input Compton-y map used by the SZ map generators",
    )
    parser.add_argument(
        "--rseed",
        type=int,
        default=12345,
        help="Random seed forwarded to the SZ map generators",
    )
    parser.add_argument(
        "--noise-seed",
        type=int,
        default=67890,
        help="Noise seed used when invoking mkoneszmap",
    )
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Run the commands instead of only printing them",
    )
    parser.add_argument(
        "--skip-single-band",
        action="store_true",
        help="Do not invoke mkoneszmap (only print the command)",
    )
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)

    print("Configured observing bands:")
    print("  frequencies (GHz):", [band.frequency_ghz for band in BANDS])
    print("  FWHM (arcmin):", [band.fwhm_arcmin for band in BANDS])
    print("  normalisations:", [band.normalization for band in BANDS])
    print("  white noise (K/arcmin):", [band.white_noise_k_arcmin for band in BANDS])
    print("  beam outputs:", [band.beam_output for band in BANDS])
    print("  map outputs:", [str(band.map_output) for band in BANDS])

    _run_stage(
        "mknull",
        mknull_main,
        ["--output", str(args.y_map), "--param-file", str(args.param_file)],
        execute=args.execute,
    )

    for beam_args in _mkbeam_args(args.param_file):
        band_label = beam_args[1]
        _run_stage(f"mkbeam ({band_label})", mkbeam_main, beam_args, execute=args.execute)

    mkszmap_args = _mkszmap_args(
        args.param_file,
        args.rseed,
        args.y_map,
        BANDS[0].map_output,
        BANDS[1].map_output,
    )
    _run_stage("mkszmap", mkszmap_main, mkszmap_args, execute=args.execute)

    if not args.skip_single_band:
        for offset, band in enumerate(BANDS):
            mkone_args = _mkoneszmap_args(
                args.param_file,
                args.rseed + offset,
                args.noise_seed + offset,
                band,
                args.y_map,
            )
            _run_stage(f"mkoneszmap ({band.label})", mkoneszmap_main, mkone_args, execute=args.execute)

    _run_stage("genfilterM", genfilter_main, ["--param-file", str(args.param_file)], execute=args.execute)

    filter_args = _filtermap_args(args.param_file, BANDS[0].map_output, BANDS[1].map_output)
    _run_stage("filtermapM", filtermap_main, filter_args, execute=args.execute)

    _run_stage("szpeakfinder", szpeakfinder_main, _szpeakfinder_args(args.param_file, "filteredM"), execute=args.execute)

    _run_stage("szpeakcmp", szpeakcmp_main, ["--param-file", str(args.param_file)], execute=args.execute)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
