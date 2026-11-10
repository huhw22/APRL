"""Diagnose residual point-axis oscillations without angular filtering."""
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
from scipy.signal.windows import tukey


PLANES = [("z0060", 0.060), ("z0100", 0.100),
          ("z0200", 0.200), ("z03524", 0.3524)]
COMPONENTS = [("Ex", "electric_V_per_m", 0),
              ("Ey", "electric_V_per_m", 1),
              ("Bx", "magnetic_T", 0),
              ("By", "magnetic_T", 1)]


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def read_patch(root, plane, width=16):
    values = {}
    shared_time = None
    for kind in ("bend", "straight"):
        path = root / kind / "detectors" / f"{plane}.h5"
        before = (path.stat().st_size, path.stat().st_mtime_ns)
        with h5py.File(path, "r") as handle:
            group = handle["field_plane"]
            require(int(group["complete"][()]) == 1, f"Incomplete {path}")
            require(int(group["committed_samples"][()]) == len(group["time_s"]),
                    f"Uncommitted tail in {path}")
            time = group["time_s"][:]
            if shared_time is None:
                shared_time = time
            else:
                require(np.array_equal(shared_time, time),
                        f"Bend/straight time mismatch for {plane}")
            nx, ny = int(group.attrs["nx"]), int(group.attrs["ny"])
            dx, dy = float(group.attrs["dx_m"]), float(group.attrs["dy_m"])
            x = float(group.attrs["x_first_m"]) + np.arange(nx) * dx
            y = float(group.attrs["y_first_m"]) + np.arange(ny) * dy
            ix = np.sort(np.argsort(abs(x))[:width])
            iy = np.sort(np.argsort(abs(y))[:width])
            require(np.array_equal(ix, np.arange(ix[0], ix[-1] + 1)),
                    "Central x patch is not contiguous")
            require(np.array_equal(iy, np.arange(iy[0], iy[-1] + 1)),
                    "Central y patch is not contiguous")
            array = group["fields"][:, iy[0]:iy[-1] + 1,
                                    ix[0]:ix[-1] + 1]
            values[kind] = {
                name: array[field][..., axis]
                for name, field, axis in COMPONENTS
            }
        require(before == (path.stat().st_size, path.stat().st_mtime_ns),
                f"Closed input changed while reading: {path}")
    return shared_time, {
        name: values["bend"][name] - values["straight"][name]
        for name, _, _ in COMPONENTS
    }


def fourier(fields, time, energies, window=None):
    dt = float((time[-1] - time[0]) / (len(time) - 1))
    phase = np.exp(1j * energies[:, None] * e / hbar *
                   (time - time[0])[None, :]) * dt
    if window is None:
        window = np.ones(len(time))
    return {
        name: np.einsum("et,tij->eij", phase, value * window[:, None, None],
                        optimize=True)
        for name, value in fields.items()
    }


def poynting(values):
    return 2 * e / (h * mu_0) * np.real(
        values["Ex"] * np.conj(values["By"]) -
        values["Ey"] * np.conj(values["Bx"]))


def detrended(values, energies):
    require(np.all(values > 0), "Cannot log-detrend a non-positive spectrum")
    fit = np.exp(np.polyval(np.polyfit(energies, np.log(values), 3), energies))
    return values / fit


def ripple(values, energies):
    ratio = detrended(values, energies)
    extrema = np.where(np.diff(np.sign(np.diff(ratio))) != 0)[0]
    return {
        "detrended_rms": float(np.std(ratio - 1)),
        "detrended_min": float(np.min(ratio)),
        "detrended_max": float(np.max(ratio)),
        "interior_extrema": int(len(extrema)),
    }


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main(root, output):
    root, output = Path(root).resolve(), Path(output).resolve()
    require((root / "PIC_COMPLETE.json").is_file(), "PIC result is incomplete")
    require(not output.exists(), "Output already exists")
    output.mkdir(parents=True)
    axis = np.genfromtxt(root / "analysis/axis_comparison.csv",
                         delimiter=",", names=True)
    energies = np.asarray(axis["photon_eV"])
    plane_rows = []
    plane_curves = {}
    final_time = final_fields = final_fourier = None
    for plane, z in PLANES:
        time, fields = read_patch(root, plane)
        transformed = fourier(fields, time, energies)
        centre = {name: value[:, 7:9, 7:9].mean(axis=(1, 2))
                  for name, value in transformed.items()}
        spectrum = poynting(centre)
        row = {"plane": plane, "z_m": z, **ripple(spectrum, energies)}
        plane_rows.append(row)
        plane_curves[plane] = detrended(spectrum, energies)
        if plane == PLANES[-1][0]:
            final_time, final_fields, final_fourier = time, fields, transformed
    write_csv(output / "plane_ripple.csv", plane_rows)

    patch_rows = []
    patch_curves = {}
    for width in (2, 4, 6, 8, 12, 16):
        lower = (16 - width) // 2
        upper = lower + width
        patch = {name: value[:, lower:upper, lower:upper]
                 for name, value in final_fourier.items()}
        coherent = poynting({name: value.mean(axis=(1, 2))
                             for name, value in patch.items()})
        pixel_mean = poynting(patch).mean(axis=(1, 2))
        patch_rows.append({
            "width_cells": width,
            "width_um": width * 5.0,
            **{"coherent_" + key: value
               for key, value in ripple(coherent, energies).items()},
            **{"pixel_intensity_mean_" + key: value
               for key, value in ripple(pixel_mean, energies).items()},
        })
        if width in (2, 6):
            patch_curves[f"coherent_{width}"] = detrended(coherent, energies)
            patch_curves[f"pixel_{width}"] = detrended(pixel_mean, energies)
    write_csv(output / "patch_ripple.csv", patch_rows)

    audit = json.loads((root / "SOURCE_AUDIT.json").read_text())
    sample_ff = np.asarray(audit["sample_longitudinal_form_factor_power"])
    analytic_ff = np.asarray(audit["analytic_longitudinal_form_factor_power"])
    ff_ratio = sample_ff / analytic_ff
    ff_ripple = ripple(ff_ratio, energies)

    amplitude = np.sqrt(final_fields["Ex"][:, 7:9, 7:9].mean((1, 2))**2 +
                        final_fields["Ey"][:, 7:9, 7:9].mean((1, 2))**2 +
                        (c * final_fields["Bx"][:, 7:9, 7:9].mean((1, 2)))**2 +
                        (c * final_fields["By"][:, 7:9, 7:9].mean((1, 2)))**2)
    relative_time_fs = (final_time - final_time[0]) * 1e15
    span = float(final_time[-1] - final_time[0])

    rect_centre = {name: value[:, 7:9, 7:9].mean(axis=(1, 2))
                   for name, value in final_fourier.items()}
    rect_spectrum = poynting(rect_centre)
    taper_rows = []
    for name, alpha in (("rectangular", 0.0), ("tukey_0p10", .10),
                        ("tukey_0p25", .25)):
        window = np.ones(len(final_time)) if alpha == 0 else \
            tukey(len(final_time), alpha)
        transformed = fourier(final_fields, final_time, energies, window)
        centre = {key: value[:, 7:9, 7:9].mean(axis=(1, 2))
                  for key, value in transformed.items()}
        spectrum = poynting(centre)
        taper_rows.append({
            "window": name,
            **ripple(spectrum, energies),
            "relative_L2_to_rectangular": float(
                np.linalg.norm(spectrum - rect_spectrum) /
                np.linalg.norm(rect_spectrum)),
        })
    write_csv(output / "window_sensitivity.csv", taper_rows)

    selected = []
    central_ex = final_fourier["Ex"][:, 7:9, 7:9]
    for requested in (6.1, 8.55, 9.775, 11.05):
        index = int(np.argmin(abs(energies - requested)))
        values = central_ex[index].ravel()
        phases = np.angle(values / np.mean(values))
        pairwise = np.angle(values[:, None] / values[None, :])
        selected.append({
            "photon_eV": float(energies[index]),
            "maximum_absolute_pairwise_Ex_phase_rad":
                float(np.max(abs(pairwise))),
            "minimum_pixel_amplitude_over_mean_amplitude":
                float(np.min(abs(values)) / abs(np.mean(values))),
            "maximum_pixel_amplitude_over_mean_amplitude":
                float(np.max(abs(values)) / abs(np.mean(values))),
            "relative_phases_rad": [float(value) for value in phases],
        })
    (output / "central_pixel_phase.json").write_text(
        json.dumps(selected, indent=2) + "\n")

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.4),
                             layout="constrained")
    raw = np.asarray(axis["APRL_raw_axis_J_per_eV_m2"])
    srw = np.asarray(axis["SRW_continuum_J_per_eV_m2"])
    axes[0, 0].plot(energies, srw * 1e9, "k-", lw=2.2,
                    label="SRW coherent Gaussian")
    axes[0, 0].plot(energies, raw * 1e9, color="#e45756", lw=1.8,
                    label="APRL point-axis (central 2x2 fields)")
    axes[0, 0].set(title="Absolute point-axis comparison",
                   ylabel=r"nJ / (eV m$^2$)")
    axes[0, 0].legend(fontsize=8)

    axes[0, 1].plot(energies, patch_curves["coherent_2"],
                    color="#e45756", label="2x2 coherent field")
    axes[0, 1].plot(energies, patch_curves["pixel_2"],
                    color="#f58518", label="2x2 pixel intensity mean")
    axes[0, 1].plot(energies, patch_curves["pixel_6"],
                    color="#4c78a8", label="6x6 pixel intensity mean")
    axes[0, 1].plot(energies, detrended(ff_ratio, energies), "k--",
                    label="actual/analytic longitudinal FF")
    axes[0, 1].axhline(1, color="grey", lw=.8)
    axes[0, 1].set(title="Oscillation after removal of a cubic log trend",
                   ylabel="detrended ratio")
    axes[0, 1].legend(fontsize=8)

    axes[1, 0].plot([row["z_m"] for row in plane_rows],
                    [100 * row["detrended_rms"] for row in plane_rows],
                    "o-", color="#e45756")
    axes[1, 0].set(title="Point-axis ripple grows during propagation",
                   xlabel="detector z [m]", ylabel="detrended RMS [%]")

    axes[1, 1].plot(relative_time_fs, amplitude / np.max(amplitude),
                    color="#4c78a8")
    axes[1, 1].set(title=(f"Final point-axis time trace: T={span*1e15:.3f} fs, "
                         f"h/T={h/span/e:.3f} eV"),
                   xlabel="time from first sample [fs]",
                   ylabel="field-vector magnitude / maximum")
    axes[1, 1].set_yscale("symlog", linthresh=1e-3)
    for axis in axes.flat:
        axis.grid(alpha=.2)
    for axis in axes[0]:
        axis.set_xlabel("photon energy [eV]")
    fig.suptitle("Residual axis oscillation diagnosis; no angular filtering")
    fig.savefig(output / "axis_oscillation_diagnosis.png", dpi=200)
    fig.savefig(output / "axis_oscillation_diagnosis.pdf")
    plt.close(fig)

    result = {
        "complete": True,
        "angular_filter_used": False,
        "record_span_fs": span * 1e15,
        "rectangular_window_energy_scale_eV": h / span / e,
        "plane_ripple": plane_rows,
        "patch_ripple": patch_rows,
        "longitudinal_sample_FF_ripple": ff_ripple,
        "window_sensitivity": taper_rows,
        "central_pixel_phase": selected,
        "interpretation": (
            "The residual oscillation is mainly point-axis interference among "
            "transverse grid samples whose phases separate during propagation. "
            "Finite record length contributes, but mild endpoint tapers do not "
            "remove the ripple. The longitudinal finite-sample form factor is smooth."),
        "physical_accuracy_accepted": False,
        "next_discriminating_test": (
            "Repeat with transverse grid refinement and particle-count/seed scans, "
            "and compare a declared finite detector aperture with the same SRW aperture."),
    }
    (output / "FINAL.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    require(len(sys.argv) == 3,
            "usage: diagnose_independent_axis_oscillations.py ROOT OUTPUT")
    main(*sys.argv[1:])
