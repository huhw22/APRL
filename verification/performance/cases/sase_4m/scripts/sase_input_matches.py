#!/usr/bin/env python3
"""Return success only when a completed SASE input matches this campaign."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import h5py


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=Path)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--electrons", type=float, required=True)
    parser.add_argument("--length-m", type=float, required=True)
    parser.add_argument("--resonant-wavelength-m", type=float, required=True)
    parser.add_argument("--gamma", type=float, required=True)
    parser.add_argument("--sigma-gamma", type=float, required=True)
    parser.add_argument("--sigma-x-m", type=float, required=True)
    parser.add_argument("--sigma-y-m", type=float, required=True)
    parser.add_argument("--emittance-n-x-m", type=float, required=True)
    parser.add_argument("--emittance-n-y-m", type=float, required=True)
    parser.add_argument("--harmonics", type=int, required=True)
    parser.add_argument("--beamlets-per-wavelength", type=int, required=True)
    parser.add_argument("--boost-gamma", type=float, required=True)
    return parser.parse_args()


def close_float(actual: object, expected: float) -> bool:
    return math.isclose(float(actual), expected, rel_tol=2.0e-12, abs_tol=0.0)


def main() -> int:
    args = parse_args()
    expected_float = {
        "represented_electrons": args.electrons,
        "requested_length_m": args.length_m,
        "resonant_wavelength_m": args.resonant_wavelength_m,
        "gamma_mean": args.gamma,
        "sigma_gamma": args.sigma_gamma,
        "sigma_x_m": args.sigma_x_m,
        "sigma_y_m": args.sigma_y_m,
        "emittance_n_x_m": args.emittance_n_x_m,
        "emittance_n_y_m": args.emittance_n_y_m,
    }
    expected_integer = {
        "random_seed": args.seed,
        "shot_noise_harmonics": args.harmonics,
        "beamlets_per_wavelength": args.beamlets_per_wavelength,
    }
    mismatches: list[str] = []
    try:
        with h5py.File(args.path, "r") as handle:
            group = handle["/particles"]
            attributes = group.attrs
            if int(attributes.get("complete", 0)) != 1:
                mismatches.append("complete != 1")
            if int(attributes.get("format_version", -1)) != 2:
                mismatches.append("format_version != 2")
            if int(attributes.get("transverse_centering_version", -1)) != 1:
                mismatches.append("transverse_centering_version != 1")
            if ("transverse_centering_boost_gamma" not in attributes or
                    not close_float(attributes["transverse_centering_boost_gamma"],
                                    args.boost_gamma)):
                mismatches.append("transverse_centering_boost_gamma mismatch")
            if "transverse_centering_mean_beta_box" not in attributes:
                mismatches.append("missing transverse_centering_mean_beta_box")
            else:
                mean_beta = attributes["transverse_centering_mean_beta_box"]
                if math.hypot(float(mean_beta[0]), float(mean_beta[1])) > 1.0e-12:
                    mismatches.append("transverse mean beta is not centered")
            for name, expected in expected_float.items():
                if name not in attributes or not close_float(attributes[name], expected):
                    mismatches.append(
                        f"{name}={attributes.get(name, '<missing>')} expected {expected}"
                    )
            for name, expected in expected_integer.items():
                if name not in attributes or int(attributes[name]) != expected:
                    mismatches.append(
                        f"{name}={attributes.get(name, '<missing>')} expected {expected}"
                    )
    except (OSError, KeyError, IndexError, TypeError, ValueError) as error:
        mismatches.append(str(error))

    if mismatches:
        print(f"SASE input does not match campaign: {args.path}")
        for mismatch in mismatches:
            print("  - " + mismatch)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
