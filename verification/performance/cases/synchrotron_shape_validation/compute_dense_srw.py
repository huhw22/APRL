"""Compute a dense independent SRW axis spectrum for notch diagnosis."""
import csv
import json
from pathlib import Path
import sys

import numpy as np
from scipy.constants import c, e, hbar

from analyze_srw import srw_point_spectrum


def relative_l2(value, reference):
    return float(np.linalg.norm(value - reference) /
                 max(np.linalg.norm(reference), 1e-300))


def main(root, output):
    root, output = Path(root).resolve(), Path(output).resolve()
    if output.exists():
        raise RuntimeError("Dense SRW output already exists")
    output.mkdir(parents=True)
    campaign = json.loads((root / "campaign.json").read_text())
    common = campaign["common"]
    p = common["parameters"]
    energies = np.linspace(4., 12., 321)
    _, point = srw_point_spectrum(energies, p, 1.5e-3, 100000, 1.)
    _, check = srw_point_spectrum(energies, p, 1.5e-3, 20000, 1.)
    coarse = np.genfromtxt(root / "srw_analysis/gaussian_comparison.csv",
                           delimiter=",", names=True)
    correction = np.interp(
        energies, coarse["photon_eV"], coarse["LW_transverse_correction"])
    omega = energies * e / hbar
    form2 = np.exp(-(omega * common["sigma_z_m"] /
                     (p["beta"] * c))**2)
    gaussian = (point * common["represented_electrons"]**2 * form2 *
                correction)
    with (output / "dense_srw_axis.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["photon_eV", "SRW_point_per_electron_J_per_eV_m2",
                         "longitudinal_form_factor_squared",
                         "transverse_correction_interpolated",
                         "SRW_gaussian_J_per_eV_m2"])
        writer.writerows(zip(energies, point, form2, correction, gaussian))
    result = {
        "complete": True,
        "energies": len(energies),
        "energy_range_eV": [float(energies[0]), float(energies[-1])],
        "energy_step_eV": float(energies[1] - energies[0]),
        "trajectory_points": 100000,
        "check_trajectory_points": 20000,
        "point_spectrum_20000_vs_100000_relative_L2":
            relative_l2(check, point),
        "transverse_correction_source":
            str(root / "srw_analysis/gaussian_comparison.csv"),
        "normalization": "none",
    }
    (output / "VALIDATION.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise RuntimeError("usage: compute_dense_srw.py ROOT OUTPUT")
    main(*sys.argv[1:])
