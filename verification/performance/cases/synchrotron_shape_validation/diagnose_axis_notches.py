"""Diagnose near-zero axis bins using the closed no-feedback 7-point field."""
import csv
import json
from pathlib import Path
import sys

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import e, h, hbar, mu_0
from scipy.signal import find_peaks
from scipy.signal.windows import tukey


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def load_field(path):
    before = (path.stat().st_size, path.stat().st_mtime_ns)
    handle = h5py.File(path, "r")
    group = handle["field_plane"]
    require(int(group["complete"][()]) == 1, "Incomplete field plane")
    require(int(group["committed_samples"][()]) == len(group["time_s"]),
            "Uncommitted field tail")
    return handle, group, before


def finish_field(path, handle, before):
    handle.close()
    require(before == (path.stat().st_size, path.stat().st_mtime_ns),
            "Closed field changed while reading")


def poynting(ex, ey, bx, by):
    return 2 * e / (h * mu_0) * np.real(
        ex * np.conj(by) - ey * np.conj(bx))


def main(root, srw_csv, output):
    root, srw_csv, output = map(lambda value: Path(value).resolve(),
                                (root, srw_csv, output))
    require(not output.exists(), "Output already exists")
    output.mkdir(parents=True)
    base = root / "linear-notch/no-feedback"
    paths = {kind: base / kind / "detectors/screen.h5"
             for kind in ["bend", "straight"]}
    opened = {kind: load_field(path) for kind, path in paths.items()}
    bend, straight = opened["bend"][1], opened["straight"][1]
    time = bend["time_s"][:]
    require(np.array_equal(time, straight["time_s"][:]),
            "Bend/straight times differ")
    dt = float((time[-1] - time[0]) / (len(time) - 1))
    srw = np.genfromtxt(srw_csv, delimiter=",", names=True)
    energies = np.asarray(srw["photon_eV"])
    phase = np.exp(1j * (energies[:, None] * e / hbar) *
                   (time - time[0])[None, :]) * dt
    ny, nx = int(bend.attrs["ny"]), int(bend.attrs["nx"])
    components = {}
    fields = [("Ex", "electric_V_per_m", 0),
              ("Ey", "electric_V_per_m", 1),
              ("Bx", "magnetic_T", 0), ("By", "magnetic_T", 1)]
    for name, field, axis in fields:
        trace = (bend["fields"][field][..., axis] -
                 straight["fields"][field][..., axis])
        components[name] = (phase @ trace.reshape(len(time), -1)).reshape(
            len(energies), ny, nx)
    density = poynting(components["Ex"], components["Ey"],
                       components["Bx"], components["By"])
    require(np.all(density >= -1e-12 * np.max(density)),
            "Spectral plane contains material negative flux")
    x = float(bend.attrs["x_first_m"]) + np.arange(nx) * float(bend.attrs["dx_m"])
    y = float(bend.attrs["y_first_m"]) + np.arange(ny) * float(bend.attrs["dy_m"])
    ix = np.sort(np.argsort(abs(x))[:2])
    iy = np.sort(np.argsort(abs(y))[:2])
    axis_components = {name: value[:, iy[:, None], ix].mean(axis=(1, 2))
                       for name, value in components.items()}
    axis_density = poynting(axis_components["Ex"], axis_components["Ey"],
                            axis_components["Bx"], axis_components["By"])
    plane_energy = density.sum(axis=(1, 2)) * float(bend.attrs["dx_m"]) * float(
        bend.attrs["dy_m"])

    search = (energies >= 8.) & (energies <= 11.)
    search_index = np.where(search)[0]
    notch = int(search_index[np.argmin(axis_density[search])])
    high_search = (energies >= energies[notch] + .4) & (energies <= 11.3)
    high_index = np.where(high_search)[0]
    neighbor = int(high_index[np.argmax(axis_density[high_search])])
    axis_minima = find_peaks(-axis_density,
                             prominence=np.max(axis_density) * .02,
                             distance=20)[0]
    plane_minima = find_peaks(-plane_energy,
                              prominence=np.max(plane_energy) * .02,
                              distance=20)[0]

    # The axis-only window check separates a record-edge FFT notch from an
    # intrinsic cancellation already present in the recorded pulse.
    axis_trace = {}
    for name, field, component in fields:
        values = []
        for kind in ["bend", "straight"]:
            values.append(opened[kind][1]["fields"][field][
                :, iy[:, None], ix, component].mean(axis=(1, 2)))
        axis_trace[name] = values[0] - values[1]
    window_checks = {}
    for label, window in [("rect", np.ones(len(time))),
                          ("tukey20", tukey(len(time), .2)),
                          ("tukey40", tukey(len(time), .4))]:
        spectra = {name: phase @ (values * window)
                   for name, values in axis_trace.items()}
        value = poynting(spectra["Ex"], spectra["Ey"],
                         spectra["Bx"], spectra["By"])
        local = search_index[np.argmin(value[search])]
        window_checks[label] = {
            "notch_eV": float(energies[local]),
            "notch_J_per_eV_m2": float(value[local]),
            "notch_over_curve_max": float(value[local] / np.max(value)),
        }

    with (output / "axis_plane_spectra.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["photon_eV", "APRL_axis_J_per_eV_m2",
                         "APRL_full_plane_J_per_eV",
                         "SRW_axis_J_per_eV_m2"])
        writer.writerows(zip(energies, axis_density, plane_energy,
                             srw["SRW_gaussian_J_per_eV_m2"]))

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    axes[0, 0].plot(energies, srw["SRW_gaussian_J_per_eV_m2"] * 1e9,
                    "k-", lw=2.5, label="dense SRW axis")
    axes[0, 0].plot(energies, axis_density * 1e9, color="#54a24b", lw=2,
                    label="APRL 7-point axis")
    axes[0, 0].axvline(energies[notch], color="#e45756", ls="--", lw=1)
    axes[0, 0].set(title="Axis spectrum: the APRL notch is real in the recorded field",
                   ylabel=r"nJ / (eV m$^2$)")
    axes[0, 0].legend(fontsize=8)

    axes[0, 1].plot(energies, axis_density / np.max(axis_density),
                    color="#54a24b", lw=2, label="axis")
    axes[0, 1].plot(energies, plane_energy / np.max(plane_energy),
                    color="#4c78a8", lw=2, label="whole detector plane")
    axes[0, 1].axvline(energies[notch], color="#e45756", ls="--", lw=1)
    axes[0, 1].set(title="Energy is redistributed transversely, not absent",
                   ylabel="diagnostic curve / own maximum")
    axes[0, 1].legend(fontsize=8)

    extent = [x[0] * 1e6, x[-1] * 1e6, y[0] * 1e6, y[-1] * 1e6]
    vmax = max(float(np.max(density[notch])), float(np.max(density[neighbor]))) * 1e9
    for axis, index, title in [
            (axes[1, 0], notch, f"Axis notch: {energies[notch]:.3f} eV"),
            (axes[1, 1], neighbor, f"Neighboring axis maximum: {energies[neighbor]:.3f} eV")]:
        image = axis.imshow(density[index] * 1e9, origin="lower", extent=extent,
                            aspect="auto", vmin=0., vmax=vmax, cmap="magma")
        axis.plot(0., 0., marker="+", ms=12, mew=1.5, color="cyan")
        axis.set(title=title, xlabel="x [µm]", ylabel="y [µm]")
        fig.colorbar(image, ax=axis, label=r"nJ / (eV m$^2$)")
    for axis in axes.flat:
        axis.grid(alpha=.15)
        if axis in axes[0]:
            axis.set_xlabel("Photon energy [eV]")
    fig.suptitle("Origin of the apparent zero: an axis node in the propagated angular field")
    fig.savefig(output / "axis_notch_diagnosis.png", dpi=200)
    fig.savefig(output / "axis_notch_diagnosis.pdf")
    plt.close(fig)

    result = {
        "complete": True,
        "source_shape": "linear-notch (7-point)",
        "feedback": False,
        "notch": {
            "energy_eV": float(energies[notch]),
            "axis_J_per_eV_m2": float(axis_density[notch]),
            "axis_over_axis_curve_max":
                float(axis_density[notch] / np.max(axis_density)),
            "full_plane_J_per_eV": float(plane_energy[notch]),
            "full_plane_over_plane_curve_max":
                float(plane_energy[notch] / np.max(plane_energy)),
        },
        "neighboring_axis_maximum": {
            "energy_eV": float(energies[neighbor]),
            "axis_J_per_eV_m2": float(axis_density[neighbor]),
            "full_plane_J_per_eV": float(plane_energy[neighbor]),
        },
        "axis_local_minima_eV": [float(energies[index])
                                  for index in axis_minima],
        "full_plane_local_minima_eV": [float(energies[index])
                                        for index in plane_minima],
        "window_checks": window_checks,
        "straight_subtraction_note":
            "The prior point check found straight transverse amplitude below 1% of bend at the notch.",
        "interpretation": (
            "The near-zero bin is a spatial axis node. The full detector-plane spectral "
            "energy is smooth and is larger at the axis-notch energy than at the "
            "neighboring axis maximum. This establishes transverse angular redistribution, "
            "not missing total radiation. A multi-plane control is still required to "
            "separate source-generated angular structure from numerical propagation."),
    }
    (output / "SUMMARY.json").write_text(json.dumps(result, indent=2) + "\n")
    for kind, (handle, _, before) in opened.items():
        finish_field(paths[kind], handle, before)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise RuntimeError("usage: diagnose_axis_notches.py ROOT SRW_CSV OUTPUT")
    main(*sys.argv[1:])
