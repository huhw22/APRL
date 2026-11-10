"""Track the transverse high-k axis notch through four detector planes."""
import csv
import json
from pathlib import Path
import sys

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c, e, h, hbar, mu_0


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def poynting(values):
    return 2 * e / (h * mu_0) * np.real(
        values["Ex"] * np.conj(values["By"]) -
        values["Ey"] * np.conj(values["Bx"]))


def files_identical(left, right):
    with h5py.File(left, "r") as a, h5py.File(right, "r") as b:
        ga, gb = a["field_plane"], b["field_plane"]
        if not np.array_equal(ga["time_s"][:], gb["time_s"][:]):
            return False
        if ga["fields"].shape != gb["fields"].shape:
            return False
        for start in range(0, len(ga["fields"]), 8):
            stop = min(len(ga["fields"]), start + 8)
            if not np.array_equal(ga["fields"][start:stop],
                                  gb["fields"][start:stop]):
                return False
    return True


def closed_group(path):
    before = (path.stat().st_size, path.stat().st_mtime_ns)
    handle = h5py.File(path, "r")
    group = handle["field_plane"]
    require(int(group["complete"][()]) == 1, f"Incomplete field plane: {path}")
    require(int(group["committed_samples"][()]) == len(group["time_s"]),
            f"Uncommitted field tail: {path}")
    return handle, group, before


def close_unchanged(path, opened):
    handle, _, before = opened
    handle.close()
    require(before == (path.stat().st_size, path.stat().st_mtime_ns),
            f"Closed field changed while reading: {path}")


def plane_spectra(bend_path, straight_path, energies, gamma):
    paths = {"bend": bend_path, "straight": straight_path}
    opened = {kind: closed_group(path) for kind, path in paths.items()}
    groups = {kind: value[1] for kind, value in opened.items()}
    time = groups["bend"]["time_s"][:]
    require(np.array_equal(time, groups["straight"]["time_s"][:]),
            "Bend/straight time mismatch")
    dt = float((time[-1] - time[0]) / (len(time) - 1))
    phase = np.exp(1j * energies[:, None] * e / hbar *
                   (time - time[0])[None, :]) * dt
    ny, nx = int(groups["bend"].attrs["ny"]), int(groups["bend"].attrs["nx"])
    dx, dy = float(groups["bend"].attrs["dx_m"]), float(groups["bend"].attrs["dy_m"])
    components = {}
    for name, field, axis in [("Ex", "electric_V_per_m", 0),
                              ("Ey", "electric_V_per_m", 1),
                              ("Bx", "magnetic_T", 0),
                              ("By", "magnetic_T", 1)]:
        trace = (groups["bend"]["fields"][field][..., axis] -
                 groups["straight"]["fields"][field][..., axis])
        components[name] = (phase @ trace.reshape(len(time), -1)).reshape(
            len(energies), ny, nx)
    kx = 2 * np.pi * np.fft.fftfreq(nx, dx)
    ky = 2 * np.pi * np.fft.fftfreq(ny, dy)
    kx_grid, ky_grid = np.meshgrid(kx, ky)
    k_perp = np.sqrt(kx_grid**2 + ky_grid**2)
    wave_number = energies * e / (hbar * c)
    low_mask = k_perp[None, :, :] < (3 / gamma) * wave_number[:, None, None]
    nyquist_mask = ((abs(kx_grid) > .7 * np.pi / dx) |
                     (abs(ky_grid) > .7 * np.pi / dy))
    high_angle_mask = k_perp[None, :, :] > (5 / gamma) * wave_number[:, None, None]
    x = float(groups["bend"].attrs["x_first_m"]) + np.arange(nx) * dx
    y = float(groups["bend"].attrs["y_first_m"]) + np.arange(ny) * dy
    ix = np.sort(np.argsort(abs(x))[:2])
    iy = np.sort(np.argsort(abs(y))[:2])
    def axis(values):
        return {name: value[:, iy[:, None], ix].mean(axis=(1, 2))
                for name, value in values.items()}
    raw_axis = axis(components)
    low_components = {}
    ex_fourier = None
    for name, value in components.items():
        spatial = np.fft.fft2(value, axes=(1, 2))
        if name == "Ex":
            ex_fourier = spatial
        low_components[name] = np.fft.ifft2(
            spatial * low_mask, axes=(1, 2))
    low_axis = axis(low_components)
    ex_power = abs(ex_fourier)**2
    total = ex_power.sum(axis=(1, 2))
    result = {
        "raw": poynting(raw_axis),
        "low": poynting(low_axis),
        "nyquist_fraction": ex_power[:, nyquist_mask].sum(axis=1) / total,
        "high_angle_fraction":
            (ex_power * high_angle_mask).sum(axis=(1, 2)) / total,
        "low_angle_fraction":
            (ex_power * low_mask).sum(axis=(1, 2)) / total,
        "samples": len(time),
        "time_span_s": float(time[-1] - time[0]),
    }
    for kind, value in opened.items():
        close_unchanged(paths[kind], value)
    return result


def main(root, feedback_root, output):
    root, feedback_root, output = map(lambda value: Path(value).resolve(),
                                      (root, feedback_root, output))
    require(not output.exists(), "Output already exists")
    output.mkdir(parents=True)
    campaign = json.loads((root / "campaign.json").read_text())
    require((root / "PIC_COMPLETE.json").is_file(), "PIC campaign is incomplete")
    energies = np.linspace(4., 12., 321)
    gamma = 1174.171072271572
    results = {}
    for plane in campaign["planes"]:
        name = plane["name"]
        results[name] = plane_spectra(
            root / "bend/detectors" / (name + ".h5"),
            root / "straight/detectors" / (name + ".h5"), energies, gamma)

    fixed_energy = 9.825
    fixed = int(np.argmin(abs(energies - fixed_energy)))
    rows = []
    summary = []
    for plane in campaign["planes"]:
        name, z = plane["name"], float(plane["z_m"])
        value = results[name]
        search = (energies >= 8.) & (energies <= 11.)
        indices = np.where(search)[0]
        notch = int(indices[np.argmin(value["raw"][search])])
        summary.append({
            "name": name,
            "z_m": z,
            "notch_eV": float(energies[notch]),
            "notch_over_axis_max":
                float(value["raw"][notch] / np.max(value["raw"])),
            "raw_axis_at_9p825_J_per_eV_m2": float(value["raw"][fixed]),
            "low_angle_axis_at_9p825_J_per_eV_m2": float(value["low"][fixed]),
            "near_nyquist_fraction_at_9p825":
                float(value["nyquist_fraction"][fixed]),
            "outside_5_over_gamma_fraction_at_9p825":
                float(value["high_angle_fraction"][fixed]),
            "inside_3_over_gamma_fraction_at_9p825":
                float(value["low_angle_fraction"][fixed]),
            "samples": value["samples"],
            "time_span_s": value["time_span_s"],
        })
        for index, energy in enumerate(energies):
            rows.append([name, z, energy, value["raw"][index],
                         value["low"][index], value["nyquist_fraction"][index],
                         value["high_angle_fraction"][index],
                         value["low_angle_fraction"][index]])
    with (output / "propagation_spectra.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["plane", "z_m", "photon_eV",
                         "raw_axis_J_per_eV_m2", "low_angle_axis_J_per_eV_m2",
                         "near_nyquist_fraction", "outside_5_over_gamma_fraction",
                         "inside_3_over_gamma_fraction"])
        writer.writerows(rows)
    with (output / "propagation_summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)

    colors = ["#4c78a8", "#f58518", "#54a24b", "#b279a2"]
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    for plane, color in zip(campaign["planes"], colors):
        name, z = plane["name"], float(plane["z_m"])
        label = f"z={z:.3f} m"
        value = results[name]
        axes[0, 0].plot(energies, value["raw"] * 1e9, color=color, label=label)
        axes[0, 1].plot(energies, value["low"] * 1e9, color=color, label=label)
    axes[0, 0].set(title="Raw axis spectra", ylabel=r"nJ / (eV m$^2$)")
    axes[0, 1].set(title=r"Low-angle spectra ($k_\perp/k<3/\gamma$)",
                   ylabel=r"nJ / (eV m$^2$)")
    z = np.array([row["z_m"] for row in summary])
    axes[1, 0].plot(z, [row["notch_over_axis_max"] for row in summary],
                    "o-", color="#e45756")
    axes[1, 0].set(title="Deepest 8–11 eV axis notch",
                   xlabel="detector z [m]", ylabel="notch / axis-curve maximum")
    axes[1, 1].plot(z, [row["near_nyquist_fraction_at_9p825"]
                        for row in summary], "o-", color="#e45756",
                    label="near Nyquist")
    axes[1, 1].plot(z, [row["outside_5_over_gamma_fraction_at_9p825"]
                        for row in summary], "s-", color="#f58518",
                    label=r"outside $5/\gamma$")
    axes[1, 1].plot(z, [row["inside_3_over_gamma_fraction_at_9p825"]
                        for row in summary], "^-", color="#4c78a8",
                    label=r"inside $3/\gamma$")
    axes[1, 1].set(title="Transverse Ex angular content at 9.825 eV",
                   xlabel="detector z [m]", ylabel="power fraction")
    for axis in axes.flat:
        axis.grid(alpha=.2)
    for axis in axes[0]:
        axis.set_xlabel("Photon energy [eV]")
        axis.legend(fontsize=8)
    axes[1, 1].legend(fontsize=8)
    fig.suptitle("Where the synchrotron axis notch and transverse high-k mode appear")
    fig.savefig(output / "notch_propagation.png", dpi=200)
    fig.savefig(output / "notch_propagation.pdf")
    plt.close(fig)

    reproduced = all(files_identical(
        root / kind / "detectors/z03524.h5",
        feedback_root / "linear-notch/no-feedback" / kind /
        "detectors/screen.h5") for kind in ["bend", "straight"])
    first, last = summary[0], summary[-1]
    nyquist_growth = (last["near_nyquist_fraction_at_9p825"] /
                      max(first["near_nyquist_fraction_at_9p825"], 1e-300))
    if first["near_nyquist_fraction_at_9p825"] > .1:
        origin = "high-k mode already material at the first 0.060 m plane"
    else:
        origin = "high-k mode small at first plane and grows during propagation"
    final = {
        "complete": True,
        "status": "MULTIPLANE_TRANSVERSE_HIGH_K_ATTRIBUTION",
        "summary": summary,
        "original_screen_bitwise_reproduced": reproduced,
        "near_nyquist_fraction_growth_first_to_last": float(nyquist_growth),
        "origin_classification": origin,
        "interpretation_limit": (
            "A material first-plane high-k fraction locates generation upstream of that "
            "plane but does not by itself distinguish particle sampling from current "
            "deposition or the short 7.6 mm post-dipole propagation interval."),
    }
    (output / "FINAL.json").write_text(json.dumps(final, indent=2) + "\n")
    print(json.dumps(final, indent=2))


if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise RuntimeError("usage: analyze_notch_propagation.py ROOT FEEDBACK_ROOT OUTPUT")
    main(*sys.argv[1:])
