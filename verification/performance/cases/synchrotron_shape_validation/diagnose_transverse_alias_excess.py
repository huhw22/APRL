"""Attribute raw axis excess to modes outside the analytic Gaussian source band."""
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


SIGMA_XY = 20e-6
COMPONENTS = [("Ex", "electric_V_per_m", 0),
              ("Ey", "electric_V_per_m", 1),
              ("Bx", "magnetic_T", 0),
              ("By", "magnetic_T", 1)]


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def poynting(values):
    return 2 * e / (h * mu_0) * np.real(
        values["Ex"] * np.conj(values["By"]) -
        values["Ey"] * np.conj(values["Bx"]))


def axis_sample(values, iy, ix):
    return {name: value[:, iy[:, None], ix].mean(axis=(1, 2))
            for name, value in values.items()}


def ripple(values, energies):
    require(np.all(values > 0), "Nonpositive spectrum")
    trend = np.exp(np.polyval(np.polyfit(energies, np.log(values), 3),
                              energies))
    return float(np.std(values / trend - 1))


def band_metrics(energies, raw, supported, alias, cross, reference, low, high):
    selected = (energies >= low) & (energies <= high)
    ref_integral = np.trapz(reference[selected], energies[selected])
    def integral_ratio(values):
        return float(np.trapz(values[selected], energies[selected]) /
                     ref_integral)
    return {
        "raw_mean_pointwise_ratio": float(np.mean(
            raw[selected] / reference[selected])),
        "raw_integral_ratio": integral_ratio(raw),
        "source_band_integral_ratio": integral_ratio(supported),
        "outside_source_band_integral_ratio": integral_ratio(alias),
        "cross_integral_ratio": integral_ratio(cross),
    }


def main(root, srw_csv, output):
    root, srw_csv, output = map(lambda value: Path(value).resolve(),
                                (root, srw_csv, output))
    require((root / "PIC_COMPLETE.json").is_file(), "PIC result is incomplete")
    require(not output.exists(), "Output already exists")
    output.mkdir(parents=True)
    srw = np.genfromtxt(srw_csv, delimiter=",", names=True)
    energies = np.asarray(srw["photon_eV"])
    campaign = json.loads((root / "campaign.json").read_text())
    audit = json.loads((root / "SOURCE_AUDIT.json").read_text())
    represented = float(campaign["represented_electrons"])
    reference = (np.asarray(srw["SRW_point_per_electron_J_per_eV_m2"]) *
                 represented**2 *
                 np.asarray(srw["longitudinal_form_factor_squared"]) *
                 np.asarray(srw["transverse_correction_interpolated"]))
    paths = {kind: root / kind / "detectors/z03524.h5"
             for kind in ("bend", "straight")}
    handles = {kind: h5py.File(path, "r") for kind, path in paths.items()}
    groups = {kind: handle["field_plane"] for kind, handle in handles.items()}
    for group in groups.values():
        require(int(group["complete"][()]) == 1, "Incomplete field plane")
    time = groups["bend"]["time_s"][:]
    require(np.array_equal(time, groups["straight"]["time_s"][:]),
            "Bend/straight time mismatch")
    dt = float((time[-1] - time[0]) / (len(time) - 1))
    phase = np.exp(1j * energies[:, None] * e / hbar *
                   (time - time[0])[None, :]) * dt
    ny, nx = int(groups["bend"].attrs["ny"]), int(groups["bend"].attrs["nx"])
    dx, dy = float(groups["bend"].attrs["dx_m"]), float(groups["bend"].attrs["dy_m"])
    x = float(groups["bend"].attrs["x_first_m"]) + np.arange(nx) * dx
    y = float(groups["bend"].attrs["y_first_m"]) + np.arange(ny) * dy
    ix = np.sort(np.argsort(abs(x))[:2])
    iy = np.sort(np.argsort(abs(y))[:2])

    components = {}
    spatial = {}
    for name, field, component in COMPONENTS:
        trace = (groups["bend"]["fields"][field][..., component] -
                 groups["straight"]["fields"][field][..., component])
        value = (phase @ trace.reshape(len(time), -1)).reshape(
            len(energies), ny, nx)
        components[name] = value
        spatial[name] = np.fft.fft2(value, axes=(1, 2))

    kx = 2 * np.pi * np.fft.fftfreq(nx, dx)
    ky = 2 * np.pi * np.fft.fftfreq(ny, dy)
    kx_grid, ky_grid = np.meshgrid(kx, ky)
    k_sigma = np.sqrt(kx_grid**2 + ky_grid**2) * SIGMA_XY
    source_mask = k_sigma <= 3.
    supported = {name: np.fft.ifft2(value * source_mask,
                                    axes=(1, 2))
                 for name, value in spatial.items()}
    outside = {name: components[name] - supported[name]
               for name in components}
    raw_spectrum = poynting(axis_sample(components, iy, ix))
    supported_spectrum = poynting(axis_sample(supported, iy, ix))
    outside_spectrum = poynting(axis_sample(outside, iy, ix))
    cross_spectrum = raw_spectrum - supported_spectrum - outside_spectrum
    ex_power = abs(spatial["Ex"])**2
    outside_power_fraction = (ex_power[:, ~source_mask].sum(axis=1) /
                              ex_power.sum(axis=(1, 2)))

    with (output / "transverse_alias_decomposition.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["photon_eV", "raw_axis_J_per_eV_m2",
                         "source_band_axis_J_per_eV_m2",
                         "outside_source_band_axis_J_per_eV_m2",
                         "cross_axis_J_per_eV_m2",
                         "Ex_outside_source_band_power_fraction",
                         "SRW_continuum_J_per_eV_m2"])
        writer.writerows(zip(energies, raw_spectrum, supported_spectrum,
                             outside_spectrum, cross_spectrum,
                             outside_power_fraction, reference))

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), layout="constrained")
    axes[0].plot(energies, reference * 1e9, "k-", lw=2.2,
                 label="SRW Gaussian continuum")
    axes[0].plot(energies, raw_spectrum * 1e9, color="#e45756", lw=1.8,
                 label="raw APRL axis")
    axes[0].plot(energies, supported_spectrum * 1e9, color="#4c78a8", lw=1.8,
                 label=r"APRL source band $k_\perp\sigma_\perp\leq3$")
    axes[0].set(ylabel=r"nJ / (eV m$^2$)",
                title="Diagnostic decomposition; raw spectrum is not filtered")
    axes[0].legend(fontsize=8)
    axes[1].plot(energies, outside_power_fraction, color="#f58518", lw=2)
    axes[1].axhline(np.exp(-9), color="black", ls="--", lw=1,
                    label=r"continuous Gaussian power envelope $e^{-9}$")
    axes[1].set(xlabel="photon energy [eV]",
                ylabel="fraction of transverse Ex Fourier power",
                title=r"Power at $k_\perp\sigma_\perp>3$")
    axes[1].legend(fontsize=8)
    for axis in axes:
        axis.grid(alpha=.2)
        axis.set_xlim(4, 12)
    fig.suptitle("Transverse source-band attribution of the low-energy excess")
    fig.savefig(output / "transverse_alias_excess.png", dpi=200)
    fig.savefig(output / "transverse_alias_excess.pdf")
    plt.close(fig)

    result = {
        "complete": True,
        "raw_spectrum_modified": False,
        "diagnostic_source_band_definition": "k_perp * sigma_xy <= 3",
        "sigma_xy_m": SIGMA_XY,
        "represented_electrons": represented,
        "sample_longitudinal_form_factor_power_relative_L2_to_gaussian":
            audit["longitudinal_form_factor_power_relative_L2_to_gaussian_4_12eV"],
        "continuous_gaussian_amplitude_at_boundary": float(np.exp(-4.5)),
        "continuous_gaussian_power_at_boundary": float(np.exp(-9)),
        "bands": {
            "4_8eV": band_metrics(energies, raw_spectrum, supported_spectrum,
                                   outside_spectrum, cross_spectrum, reference,
                                   4., 8.),
            "8_12eV": band_metrics(energies, raw_spectrum, supported_spectrum,
                                    outside_spectrum, cross_spectrum, reference,
                                    8., 12.),
        },
        "ripple_rms": {
            "raw": ripple(raw_spectrum, energies),
            "source_band": ripple(supported_spectrum, energies),
        },
        "outside_source_band_Ex_power_fraction": {
            "minimum": float(np.min(outside_power_fraction)),
            "maximum": float(np.max(outside_power_fraction)),
            "mean_4_8eV": float(np.mean(outside_power_fraction[
                (energies >= 4.) & (energies <= 8.)])),
        },
        "interpretation_scope": (
            "This is a causal diagnostic decomposition based on the known analytic "
            "Gaussian transverse source bandwidth. It is not used as the reported "
            "physical spectrum or as an acceptance filter."),
    }
    (output / "FINAL.json").write_text(json.dumps(result, indent=2) + "\n")
    for handle in handles.values():
        handle.close()
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    require(len(sys.argv) == 4,
            "usage: diagnose_transverse_alias_excess.py ROOT SRW_CSV OUTPUT")
    main(*sys.argv[1:])
