#!/usr/bin/env python3
"""Export compact, portable SASE-input diagnostics from a large HDF5 file."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any

import h5py
import numpy as np


ATTRIBUTES = (
    "generator",
    "shot_noise_model",
    "represented_electrons",
    "requested_length_m",
    "effective_length_m",
    "resonant_wavelength_m",
    "electrons_per_wavelength_cell",
    "electrons_per_beamlet",
    "wavelength_cells",
    "beamlets_per_wavelength",
    "macroparticles_per_beamlet",
    "macroparticles",
    "random_seed",
    "shot_noise_harmonics",
    "mean_beta_z",
    "mean_current_A",
    "observed_global_bunching_abs",
    "expected_global_bunching_rms",
    "mean_electrons_per_cell_times_local_bunching_squared",
)


def serializable(value: Any) -> Any:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    array = np.asarray(value)
    if array.ndim == 0:
        return array.item()
    return array.tolist()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--input-sha256", required=True)
    parser.add_argument("--generator-sha256", required=True)
    args = parser.parse_args()

    with h5py.File(args.input, "r") as handle:
        group = handle["/particles"]
        complete = int(np.asarray(group.attrs.get("complete", 0)).reshape(-1)[0])
        if complete != 1:
            raise SystemExit(f"incomplete SASE input: {args.input}")
        missing = [name for name in ATTRIBUTES if name not in group.attrs]
        if missing:
            raise SystemExit("missing SASE attributes: " + ", ".join(missing))
        payload = {
            "schema": "aprl-sase-input-metadata-v1",
            "input_file_name": args.input.name,
            "input_sha256": args.input_sha256,
            "generator_sha256": args.generator_sha256,
            "attributes": {
                name: serializable(group.attrs[name]) for name in ATTRIBUTES
            },
        }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(f".{args.output.name}.tmp-{os.getpid()}")
    temporary.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(args.output)


if __name__ == "__main__":
    main()
