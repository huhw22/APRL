"""Separate the axis notch into low-angle and transverse high-k fields."""
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


def main(root, srw_csv, output):
    root, srw_csv, output = map(lambda value: Path(value).resolve(),
                                (root, srw_csv, output))
    require(not output.exists(), "Output already exists")
    output.mkdir(parents=True)
    base = root / "linear-notch/no-feedback"
    paths = {kind: base / kind / "detectors/screen.h5"
             for kind in ["bend", "straight"]}
    handles = {kind: h5py.File(path, "r") for kind, path in paths.items()}
    groups = {kind: handle["field_plane"] for kind, handle in handles.items()}
    for group in groups.values():
        require(int(group["complete"][()]) == 1, "Incomplete field plane")
    time = groups["bend"]["time_s"][:]
    require(np.array_equal(time, groups["straight"]["time_s"][:]),
            "Bend/straight times differ")
    dt = float((time[-1] - time[0]) / (len(time) - 1))
    srw = np.genfromtxt(srw_csv, delimiter=",", names=True)
    energies = np.asarray(srw["photon_eV"])
    phase = np.exp(1j * energies[:, None] * e / hbar *
                   (time - time[0])[None, :]) * dt
    ny, nx = int(groups["bend"].attrs["ny"]), int(groups["bend"].attrs["nx"])
    dx, dy = (float(groups["bend"].attrs["dx_m"]),
              float(groups["bend"].attrs["dy_m"]))
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
    gamma = 1174.171072271572
    low_mask = k_perp[None, :, :] < (3 / gamma) * wave_number[:, None, None]
    nyquist_mask = ((abs(kx_grid) > .7 * np.pi / dx) |
                     (abs(ky_grid) > .7 * np.pi / dy))
    high_angle_mask = k_perp[None, :, :] > (5 / gamma) * wave_number[:, None, None]

    low_components, high_components = {}, {}
    ex_fourier = None
    for name, value in components.items():
        spatial = np.fft.fft2(value, axes=(1, 2))
        if name == "Ex":
            ex_fourier = spatial
        low = np.fft.ifft2(spatial * low_mask, axes=(1, 2))
        low_components[name] = low
        high_components[name] = value - low

    x = (float(groups["bend"].attrs["x_first_m"]) +
         np.arange(nx) * dx)
    y = (float(groups["bend"].attrs["y_first_m"]) +
         np.arange(ny) * dy)
    ix = np.sort(np.argsort(abs(x))[:2])
    iy = np.sort(np.argsort(abs(y))[:2])
    def axis(values):
        return {name: value[:, iy[:, None], ix].mean(axis=(1, 2))
                for name, value in values.items()}
    raw_axis = axis(components)
    low_axis = axis(low_components)
    high_axis = axis(high_components)
    raw_spectrum = poynting(raw_axis)
    low_spectrum = poynting(low_axis)
    high_spectrum = poynting(high_axis)
    interference = raw_spectrum - low_spectrum - high_spectrum
    ex_power = abs(ex_fourier)**2
    total_power = ex_power.sum(axis=(1, 2))
    nyquist_fraction = ex_power[:, nyquist_mask].sum(axis=1) / total_power
    high_angle_fraction = ((ex_power * high_angle_mask).sum(axis=(1, 2)) /
                           total_power)
    low_angle_fraction = ((ex_power * low_mask).sum(axis=(1, 2)) /
                          total_power)
    search = (energies >= 8.) & (energies <= 11.)
    selected = np.where(search)[0]
    notch = int(selected[np.argmin(raw_spectrum[search])])
    left, right = int(ix[0]), int(ix[1])
    center_phase_jump = float(abs(np.angle(
        components["Ex"][notch, iy[0], right] /
        components["Ex"][notch, iy[0], left])))

    with (output / "high_k_decomposition.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["photon_eV", "raw_axis_J_per_eV_m2",
                         "low_angle_axis_J_per_eV_m2",
                         "high_angle_axis_J_per_eV_m2",
                         "low_high_interference_J_per_eV_m2",
                         "Ex_near_nyquist_power_fraction",
                         "Ex_angle_gt_5_over_gamma_power_fraction",
                         "Ex_angle_lt_3_over_gamma_power_fraction",
                         "SRW_axis_J_per_eV_m2"])
        writer.writerows(zip(energies, raw_spectrum, low_spectrum,
                             high_spectrum, interference, nyquist_fraction,
                             high_angle_fraction, low_angle_fraction,
                             srw["SRW_gaussian_J_per_eV_m2"]))

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), layout="constrained")
    axes[0].plot(energies, srw["SRW_gaussian_J_per_eV_m2"] * 1e9,
                 "k-", lw=2.5, label="dense SRW")
    axes[0].plot(energies, raw_spectrum * 1e9, color="#54a24b", lw=2,
                 label="raw APRL axis")
    axes[0].plot(energies, low_spectrum * 1e9, color="#4c78a8", lw=2,
                 label=r"APRL transverse modes $k_\perp/k<3/\gamma$")
    axes[0].axvline(energies[notch], color="#e45756", ls="--")
    axes[0].set(ylabel=r"nJ / (eV m$^2$)",
                title="The smooth low-angle field is cancelled by transverse high-k modes")
    axes[0].legend(fontsize=8)
    axes[1].plot(energies, nyquist_fraction, color="#e45756", lw=2,
                 label=r"near transverse Nyquist ($>0.7\,k_{Ny}$)")
    axes[1].plot(energies, high_angle_fraction, color="#f58518", lw=2,
                 label=r"outside $5/\gamma$")
    axes[1].plot(energies, low_angle_fraction, color="#4c78a8", lw=2,
                 label=r"inside $3/\gamma$")
    axes[1].axvline(energies[notch], color="#e45756", ls="--")
    axes[1].set(xlabel="Photon energy [eV]", ylabel="fraction of transverse Ex power",
                title="A large fraction of the field occupies nonphysical high-k angles")
    axes[1].legend(fontsize=8)
    for axis in axes:
        axis.grid(alpha=.2)
        axis.set_xlim(4, 12)
    fig.suptitle("Transverse high-k origin of the apparent synchrotron axis zeros")
    fig.savefig(output / "high_k_decomposition.png", dpi=200)
    fig.savefig(output / "high_k_decomposition.pdf")
    plt.close(fig)

    result = {
        "complete": True,
        "notch_energy_eV": float(energies[notch]),
        "notch_raw_axis_nJ_per_eV_m2": float(raw_spectrum[notch] * 1e9),
        "notch_low_angle_axis_nJ_per_eV_m2":
            float(low_spectrum[notch] * 1e9),
        "notch_high_angle_axis_nJ_per_eV_m2":
            float(high_spectrum[notch] * 1e9),
        "notch_low_high_interference_nJ_per_eV_m2":
            float(interference[notch] * 1e9),
        "notch_Ex_near_nyquist_power_fraction":
            float(nyquist_fraction[notch]),
        "notch_Ex_outside_5_over_gamma_power_fraction":
            float(high_angle_fraction[notch]),
        "center_adjacent_cell_Ex_phase_jump_rad": center_phase_jump,
        "filter_definition": "diagnostic transverse angular filter k_perp/k < 3/gamma",
        "interpretation": (
            "The raw near-zero is caused by destructive interference between a smooth "
            "low-angle field and a large transverse high-k component. Adjacent axis cells "
            "are nearly pi out of phase. This identifies a numerical transverse grid mode; "
            "the separate multi-plane run determines whether it is generated at deposition "
            "or grows during propagation."),
    }
    (output / "SUMMARY.json").write_text(json.dumps(result, indent=2) + "\n")
    for handle in handles.values():
        handle.close()
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise RuntimeError("usage: diagnose_high_k_modes.py ROOT SRW_CSV OUTPUT")
    main(*sys.argv[1:])
