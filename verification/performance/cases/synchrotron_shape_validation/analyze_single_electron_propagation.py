"""Compare the one-electron APRL spectrum with SRW at every screen."""
import csv
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analyze_single_electron_synchrotron import (PLANES, read_axis,
                                                  relative_l2, ripple)
from analyze_srw import srw_point_spectrum


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def main(root, srw_csv, output):
    root, srw_csv, output = map(lambda value: Path(value).resolve(),
                                (root, srw_csv, output))
    require(not output.exists(), "Output already exists")
    output.mkdir(parents=True)
    plan = json.loads((root / "bend" / "PLAN.json").read_text())
    parameters = dict(plan["parameters"])
    entrance = float(plan["magnet_entrance_z_m"])
    table = np.genfromtxt(srw_csv, delimiter=",", names=True)
    energies = np.asarray(table["photon_eV"])
    curves = {}
    rows = []
    for plane, z_lab in PLANES:
        spectrum, _, _ = read_axis(root, plane, energies)
        local = dict(parameters, screen=z_lab - entrance)
        require(local["screen"] > local["length"],
                f"Screen {plane} is inside the magnet")
        _, reference = srw_point_spectrum(energies, local, 1.5e-3,
                                          100000, 1.)
        selected_4_8 = (energies >= 4.) & (energies <= 8.)
        selected_8_12 = (energies >= 8.) & (energies <= 12.)
        ripple_rms, _ = ripple(spectrum, energies)
        row = {
            "plane": plane,
            "lab_z_m": z_lab,
            "srw_screen_from_magnet_entrance_m": local["screen"],
            "ripple_rms": ripple_rms,
        }
        for name, selected in (("4_8eV", selected_4_8),
                               ("8_12eV", selected_8_12)):
            row[name + "_relative_L2"] = relative_l2(
                spectrum, reference, selected)
            row[name + "_mean_ratio"] = float(np.mean(
                spectrum[selected] / reference[selected]))
            row[name + "_integral_ratio"] = float(
                np.trapezoid(spectrum[selected], energies[selected]) /
                np.trapezoid(reference[selected], energies[selected]))
        rows.append(row)
        curves[plane] = (spectrum, reference)

    names = list(rows[0])
    with (output / "propagation_metrics.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=names)
        writer.writeheader()
        writer.writerows(rows)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    for (plane, z_lab), axis in zip(PLANES, axes.flat):
        spectrum, reference = curves[plane]
        axis.plot(energies, reference * 1e16, "k-", lw=2,
                  label="SRW at matching distance")
        axis.plot(energies, spectrum * 1e16, color="#e45756", lw=1.6,
                  label="APRL bend - straight")
        axis.set_title(f"z={z_lab:g} m; ripple={100*next(row['ripple_rms'] for row in rows if row['plane'] == plane):.3g}%")
        axis.set_xlabel("photon energy [eV]")
        axis.set_ylabel(r"$10^{-16}$ J / (eV m$^2$)")
        axis.grid(alpha=.2)
        axis.legend(fontsize=8)
    fig.suptitle("One-electron propagation: absolute APRL versus distance-matched SRW")
    fig.savefig(output / "single_electron_propagation.png", dpi=200)
    fig.savefig(output / "single_electron_propagation.pdf")
    plt.close(fig)

    result = {
        "complete": True,
        "bunch_form_factor_applied": False,
        "angular_filter_applied": False,
        "fitted_normalization_applied": False,
        "rows": rows,
        "interpretation_scope": (
            "Every APRL detector is compared with a fresh native SRW "
            "calculation at the same distance from the magnet entrance."),
    }
    (output / "FINAL.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    require(len(sys.argv) == 4,
            "usage: analyze_single_electron_propagation.py ROOT SRW_CSV OUTPUT")
    main(*sys.argv[1:])
