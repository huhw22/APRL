"""Compare 8x8 and 12x12 transverse sampling in the notch campaign."""
import csv
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def load_spectra(path):
    table = np.genfromtxt(path, delimiter=",", names=True, dtype=None,
                          encoding="utf-8")
    return table


def per_plane(table):
    names = list(dict.fromkeys(np.atleast_1d(table["plane"]).tolist()))
    return {name: table[table["plane"] == name] for name in names}


def main(order8_analysis, order12_analysis, order12_root, output):
    order8_analysis, order12_analysis, order12_root, output = map(
        lambda value: Path(value).resolve(),
        (order8_analysis, order12_analysis, order12_root, output))
    require(not output.exists(), "Output already exists")
    require((order12_root / "PIC_COMPLETE.json").is_file(),
            "12x12 campaign is incomplete")
    output.mkdir(parents=True)
    campaign = json.loads((order12_root / "campaign.json").read_text())
    a = per_plane(load_spectra(order8_analysis / "propagation_spectra.csv"))
    b = per_plane(load_spectra(order12_analysis / "propagation_spectra.csv"))
    require(list(a) == list(b), "Plane lists differ")
    planes = campaign["planes"]
    names = [plane["name"] for plane in planes]
    require(names == list(a), "Unexpected plane order")
    rows = []
    for plane in planes:
        name, z = plane["name"], float(plane["z_m"])
        require(np.array_equal(a[name]["photon_eV"], b[name]["photon_eV"]),
                "Energy grids differ")
        energy = a[name]["photon_eV"]
        fixed = int(np.argmin(abs(energy - 9.825)))
        band = (energy >= 8.) & (energy <= 11.)
        ia = np.where(band)[0][np.argmin(
            a[name]["raw_axis_J_per_eV_m2"][band])]
        ib = np.where(band)[0][np.argmin(
            b[name]["raw_axis_J_per_eV_m2"][band])]
        rows.append({
            "plane": name, "z_m": z,
            "order8_near_nyquist_at_9p825":
                float(a[name]["near_nyquist_fraction"][fixed]),
            "order12_near_nyquist_at_9p825":
                float(b[name]["near_nyquist_fraction"][fixed]),
            "near_nyquist_order12_over_order8": float(
                b[name]["near_nyquist_fraction"][fixed] /
                a[name]["near_nyquist_fraction"][fixed]),
            "order8_deepest_notch_eV": float(energy[ia]),
            "order8_deepest_notch_over_curve_max": float(
                a[name]["raw_axis_J_per_eV_m2"][ia] /
                np.max(a[name]["raw_axis_J_per_eV_m2"])),
            "order12_deepest_notch_eV": float(energy[ib]),
            "order12_deepest_notch_over_curve_max": float(
                b[name]["raw_axis_J_per_eV_m2"][ib] /
                np.max(b[name]["raw_axis_J_per_eV_m2"])),
            "order8_low_angle_at_9p825_J_per_eV_m2": float(
                a[name]["low_angle_axis_J_per_eV_m2"][fixed]),
            "order12_low_angle_at_9p825_J_per_eV_m2": float(
                b[name]["low_angle_axis_J_per_eV_m2"][fixed]),
        })
    with (output / "transverse_sampling_summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    last = names[-1]
    energy = a[last]["photon_eV"]
    z = np.array([row["z_m"] for row in rows])
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    axes[0, 0].plot(energy, a[last]["raw_axis_J_per_eV_m2"] * 1e9,
                    color="#e45756", lw=2, label="8x8 transverse")
    axes[0, 0].plot(energy, b[last]["raw_axis_J_per_eV_m2"] * 1e9,
                    color="#4c78a8", lw=2, label="12x12 transverse")
    axes[0, 0].set(title="Final raw axis spectrum", ylabel=r"nJ / (eV m$^2$)")
    axes[0, 0].legend()
    axes[0, 1].plot(energy, a[last]["low_angle_axis_J_per_eV_m2"] * 1e9,
                    color="#e45756", lw=2, label="8x8 transverse")
    axes[0, 1].plot(energy, b[last]["low_angle_axis_J_per_eV_m2"] * 1e9,
                    color="#4c78a8", lw=2, label="12x12 transverse")
    axes[0, 1].set(title=r"Final low-angle field ($k_\perp/k<3/\gamma$)",
                   ylabel=r"nJ / (eV m$^2$)")
    axes[0, 1].legend()
    axes[1, 0].plot(z, [row["order8_near_nyquist_at_9p825"] for row in rows],
                    "o-", color="#e45756", label="8x8 transverse")
    axes[1, 0].plot(z, [row["order12_near_nyquist_at_9p825"] for row in rows],
                    "o-", color="#4c78a8", label="12x12 transverse")
    axes[1, 0].set(title="Propagated Ex near transverse Nyquist at 9.825 eV",
                   xlabel="detector z [m]", ylabel="power fraction")
    axes[1, 0].legend()
    alias = campaign["static_source_alias"]
    axes[1, 1].bar(["8x8", "12x12"],
                   [alias["order_8"], alias["order_12"]],
                   color=["#e45756", "#4c78a8"])
    axes[1, 1].set(title="Input transverse alias after tent response",
                   ylabel="near-Nyquist power fraction")
    for axis in axes.flat:
        axis.grid(alpha=.2)
    for axis in axes[0]:
        axis.set_xlabel("Photon energy [eV]")
        axis.set_xlim(4, 12)
    fig.suptitle("Does sparse transverse macro sampling create the synchrotron axis notch?")
    fig.savefig(output / "transverse_sampling_comparison.png", dpi=200)
    fig.savefig(output / "transverse_sampling_comparison.pdf")
    plt.close(fig)

    first, final = rows[0], rows[-1]
    source_reduction = alias["order_8"] / alias["order_12"]
    first_reduction = (first["order8_near_nyquist_at_9p825"] /
                       first["order12_near_nyquist_at_9p825"])
    notch_lift = (final["order12_deepest_notch_over_curve_max"] /
                  final["order8_deepest_notch_over_curve_max"])
    source_attributed = bool(source_reduction > 5 and first_reduction > 2 and
                             notch_lift > 2)
    result = {
        "complete": True,
        "status": "TRANSVERSE_SOURCE_SAMPLING_ATTRIBUTION",
        "rows": rows,
        "static_source_alias_reduction_order8_over_order12":
            float(source_reduction),
        "first_plane_near_nyquist_reduction_order8_over_order12":
            float(first_reduction),
        "final_deepest_notch_floor_lift_order12_over_order8":
            float(notch_lift),
        "source_sampling_attributed": source_attributed,
        "interpretation": (
            "Attribution requires the independently predicted source alias, the "
            "first-plane high-k field, and the final notch floor all to improve. "
            "The propagation solver remains unchanged; downstream phase evolution "
            "can still determine where a surviving high-k mode cancels the axis field."),
    }
    (output / "FINAL.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    require(len(sys.argv) == 5,
            "usage: compare_transverse_sampling.py ORDER8_ANALYSIS "
            "ORDER12_ANALYSIS ORDER12_ROOT OUTPUT")
    main(*sys.argv[1:])
