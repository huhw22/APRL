#!/usr/bin/env python3
"""Center a generated SASE snapshot without changing particle energies/z."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import h5py
import numpy as np


CENTERING_VERSION = 1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=Path)
    parser.add_argument("--boost-gamma", type=float, required=True)
    parser.add_argument("--tolerance", type=float, default=1.0e-13)
    return parser.parse_args()


def weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    return float(np.sum(weights * values, dtype=np.longdouble) /
                 np.sum(weights, dtype=np.longdouble))


def main() -> int:
    args = parse_args()
    if not args.boost_gamma >= 1.0 or not math.isfinite(args.boost_gamma):
        raise SystemExit("--boost-gamma must be finite and at least one")
    if not args.tolerance > 0.0 or not math.isfinite(args.tolerance):
        raise SystemExit("--tolerance must be finite and positive")

    with h5py.File(args.path, "r+") as handle:
        group = handle["/particles"]
        if int(group.attrs.get("complete", 0)) != 1:
            raise SystemExit("refusing to center an incomplete SASE input")
        dataset = group["records"]
        records = dataset[...]
        positions = records["position_m"]
        velocity = records["proper_velocity"]
        weights = np.asarray(records["macro_weight"], dtype=np.float64)
        if (not np.all(np.isfinite(weights))) or np.any(weights <= 0.0):
            raise SystemExit("SASE input contains invalid macro weights")

        original_gamma = np.sqrt(
            1.0 + np.sum(velocity * velocity, axis=1, dtype=np.float64)
        )
        transverse_position_removed = np.array(
            [weighted_mean(positions[:, axis], weights) for axis in (0, 1)]
        )
        positions[:, 0] -= transverse_position_removed[0]
        positions[:, 1] -= transverse_position_removed[1]

        gamma_b = args.boost_gamma
        beta_b = math.sqrt(max(0.0, 1.0 - 1.0 / (gamma_b * gamma_b)))
        proper_velocity_removed = np.zeros(2, dtype=np.float64)
        mean_beta_box = np.zeros(2, dtype=np.float64)
        for _ in range(12):
            longitudinal_squared = (
                original_gamma * original_gamma - 1.0
                - velocity[:, 0] * velocity[:, 0]
                - velocity[:, 1] * velocity[:, 1]
            )
            if np.any(longitudinal_squared <= 0.0):
                raise SystemExit(
                    "transverse centering would make a particle momentum unphysical"
                )
            velocity[:, 2] = np.sqrt(longitudinal_squared)
            gamma_box = gamma_b * (original_gamma - beta_b * velocity[:, 2])
            if np.any(gamma_box <= 0.0) or not np.all(np.isfinite(gamma_box)):
                raise SystemExit("invalid boosted gamma during transverse centering")
            inverse_gamma_box = 1.0 / gamma_box
            mean_inverse_gamma = weighted_mean(inverse_gamma_box, weights)
            for axis in (0, 1):
                mean_beta_box[axis] = weighted_mean(
                    velocity[:, axis] * inverse_gamma_box, weights
                )
            if float(np.hypot(*mean_beta_box)) <= args.tolerance:
                break
            correction = mean_beta_box / mean_inverse_gamma
            velocity[:, 0] -= correction[0]
            velocity[:, 1] -= correction[1]
            proper_velocity_removed += correction
        else:
            raise SystemExit("transverse centering did not converge")

        final_gamma = np.sqrt(
            1.0 + np.sum(velocity * velocity, axis=1, dtype=np.float64)
        )
        gamma_error = float(np.max(np.abs(final_gamma - original_gamma)))
        gamma_scale = max(1.0, float(np.max(original_gamma)))
        if gamma_error > 64.0 * np.finfo(np.float64).eps * gamma_scale:
            raise SystemExit("transverse centering did not preserve particle gamma")

        group.attrs.modify("complete", 0)
        handle.flush()
        dataset[...] = records
        group.attrs["transverse_centering_version"] = CENTERING_VERSION
        group.attrs["transverse_centering_boost_gamma"] = gamma_b
        group.attrs["transverse_centering_mean_beta_box"] = mean_beta_box
        group.attrs["transverse_position_mean_removed_m"] = (
            transverse_position_removed
        )
        group.attrs["transverse_proper_velocity_mean_removed"] = (
            proper_velocity_removed
        )
        group.attrs.modify("complete", 1)
        handle.flush()

    print(
        "Centered SASE input"
        f" mean_beta_box=({mean_beta_box[0]:.12g}, {mean_beta_box[1]:.12g})"
        f" removed_u=({proper_velocity_removed[0]:.12g},"
        f" {proper_velocity_removed[1]:.12g})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
