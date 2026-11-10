"""Absolute one-electron APRL/SRW comparison without angular filtering."""
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


PLANES = [("z0060", .060), ("z0100", .100),
          ("z0200", .200), ("z03524", .3524)]


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def read_axis(root, plane, energies):
    values = {}
    time = None
    for kind in ("bend", "straight"):
        path = root / kind / "detectors" / f"{plane}.h5"
        before = (path.stat().st_size, path.stat().st_mtime_ns)
        with h5py.File(path, "r") as handle:
            group = handle["field_plane"]
            require(int(group["complete"][()]) == 1, f"Incomplete {path}")
            current_time = group["time_s"][:]
            if time is None:
                time = current_time
            else:
                require(np.array_equal(time, current_time), "Time mismatch")
            nx, ny = int(group.attrs["nx"]), int(group.attrs["ny"])
            dx, dy = float(group.attrs["dx_m"]), float(group.attrs["dy_m"])
            x = float(group.attrs["x_first_m"]) + np.arange(nx) * dx
            y = float(group.attrs["y_first_m"]) + np.arange(ny) * dy
            ix = np.sort(np.argsort(abs(x))[:2])
            iy = np.sort(np.argsort(abs(y))[:2])
            array = np.stack([group["fields"][:, int(j), int(i)]
                              for j in iy for i in ix], axis=1)
            values[kind] = {
                name: array[field][..., axis].mean(axis=1)
                for name, field, axis in
                [("Ex", "electric_V_per_m", 0),
                 ("Ey", "electric_V_per_m", 1),
                 ("Bx", "magnetic_T", 0),
                 ("By", "magnetic_T", 1)]
            }
        require(before == (path.stat().st_size, path.stat().st_mtime_ns),
                f"Input changed while reading {path}")
    difference = {name: values["bend"][name] - values["straight"][name]
                  for name in values["bend"]}
    dt = float((time[-1] - time[0]) / (len(time) - 1))
    phase = np.exp(1j * energies[:, None] * e / hbar *
                   (time - time[0])[None, :]) * dt
    field = {name: phase @ value for name, value in difference.items()}
    spectrum = 2 * e / (h * mu_0) * np.real(
        field["Ex"] * np.conj(field["By"]) -
        field["Ey"] * np.conj(field["Bx"]))
    amplitude = np.sqrt(difference["Ex"]**2 + difference["Ey"]**2 +
                        (c * difference["Bx"])**2 +
                        (c * difference["By"])**2)
    return spectrum, time, amplitude


def relative_l2(value, reference, selected):
    return float(np.linalg.norm(value[selected] - reference[selected]) /
                 np.linalg.norm(reference[selected]))


def ripple(values, energies):
    require(np.all(values > 0), "Nonpositive point-axis spectrum")
    trend = np.exp(np.polyval(np.polyfit(energies, np.log(values), 3), energies))
    ratio = values / trend
    return float(np.std(ratio - 1)), ratio


def main(root, srw_csv, output):
    root, srw_csv, output = map(lambda value: Path(value).resolve(),
                                (root, srw_csv, output))
    require(not output.exists(), "Output already exists")
    for kind in ("bend", "straight"):
        require((root / kind / "PIC_COMPLETE.json").is_file(),
                f"Missing completed {kind} case")
    output.mkdir(parents=True)
    srw = np.genfromtxt(srw_csv, delimiter=",", names=True)
    energies = np.asarray(srw["photon_eV"])
    reference = np.asarray(srw["SRW_point_per_electron_J_per_eV_m2"])
    curves = {}
    rows = []
    final_time = final_amplitude = None
    for plane, z in PLANES:
        spectrum, time, amplitude = read_axis(root, plane, energies)
        value, detrended = ripple(spectrum, energies)
        curves[plane] = (spectrum, detrended)
        rows.append({"plane": plane, "z_m": z, "detrended_rms": value})
        if plane == PLANES[-1][0]:
            final_time, final_amplitude = time, amplitude
    final = curves[PLANES[-1][0]][0]
    with (output / "axis_comparison.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["photon_eV", "APRL_single_electron_J_per_eV_m2",
                         "SRW_single_electron_J_per_eV_m2"])
        writer.writerows(zip(energies, final, reference))
    bands = {}
    for low, high in ((4., 8.), (8., 12.), (4., 12.)):
        selected = (energies >= low) & (energies <= high)
        bands[f"{low:g}_{high:g}eV"] = {
            "relative_L2": relative_l2(final, reference, selected),
            "mean_pointwise_ratio": float(np.mean(final[selected] /
                                                    reference[selected])),
            "integral_ratio": float(np.trapezoid(final[selected],
                                                   energies[selected]) /
                                    np.trapezoid(reference[selected],
                                                 energies[selected])),
        }
    span = float(final_time[-1] - final_time[0])
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    axes[0, 0].plot(energies, reference * 1e16, "k-", lw=2.2,
                    label="native single-electron SRW")
    axes[0, 0].plot(energies, final * 1e16, color="#e45756", lw=1.8,
                    label="APRL bend - straight")
    axes[0, 0].set(title="Absolute one-electron point-axis spectrum",
                   ylabel=r"$10^{-16}$ J / (eV m$^2$)")
    axes[0, 0].legend(fontsize=8)
    axes[0, 1].plot(energies, final / reference, color="#e45756")
    axes[0, 1].axhline(1, color="black", lw=1)
    axes[0, 1].set(title="Absolute ratio; no fitted normalization",
                   ylabel="APRL / SRW")
    axes[1, 0].plot([row["z_m"] for row in rows],
                    [100 * row["detrended_rms"] for row in rows], "o-")
    axes[1, 0].set(title="Point-axis ripple versus propagation",
                   xlabel="detector z [m]", ylabel="detrended RMS [%]")
    axes[1, 1].plot((final_time - final_time[0]) * 1e15,
                    final_amplitude / np.max(final_amplitude))
    axes[1, 1].set(title=f"Final trace; T={span*1e15:.3f} fs",
                   xlabel="time from first sample [fs]",
                   ylabel="field magnitude / maximum")
    axes[1, 1].set_yscale("symlog", linthresh=1e-3)
    for axis in axes.flat:
        axis.grid(alpha=.2)
    for axis in axes[0]:
        axis.set_xlabel("photon energy [eV]")
    fig.suptitle("One-electron synchrotron control; no bunch form factor or angular filter")
    fig.savefig(output / "single_electron_srw_comparison.png", dpi=200)
    fig.savefig(output / "single_electron_srw_comparison.pdf")
    plt.close(fig)
    result = {
        "complete": True,
        "particle_count": 1,
        "macro_weight": 1,
        "bunch_form_factor_applied": False,
        "angular_filter_applied": False,
        "fitted_normalization_applied": False,
        "bands": bands,
        "plane_ripple": rows,
        "record_span_fs": span * 1e15,
        "physical_accuracy_accepted": False,
        "acceptance_note": (
            "Acceptance requires inspecting the absolute APRL/SRW metrics and "
            "their stability under transverse-grid refinement."),
    }
    (output / "FINAL.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    require(len(sys.argv) == 4,
            "usage: analyze_single_electron_synchrotron.py ROOT SRW_CSV OUTPUT")
    main(*sys.argv[1:])
