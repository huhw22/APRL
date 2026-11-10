"""Compare the independent-electron synchrotron pair with dense SRW."""
import csv
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analyze_notch_propagation import plane_spectra, require


def relative_l2(value, reference, selected=None):
    if selected is not None:
        value, reference = value[selected], reference[selected]
    return float(np.linalg.norm(value - reference) /
                 max(np.linalg.norm(reference), 1e-300))


def main(root, srw_csv, output):
    root, srw_csv, output = map(lambda value: Path(value).resolve(),
                                (root, srw_csv, output))
    require((root / "PIC_COMPLETE.json").is_file(), "PIC campaign incomplete")
    require(not output.exists(), "Analysis output exists")
    output.mkdir(parents=True)
    campaign = json.loads((root / "campaign.json").read_text())
    audit = json.loads((root / "SOURCE_AUDIT.json").read_text())
    srw = np.genfromtxt(srw_csv, delimiter=",", names=True)
    energies = np.asarray(srw["photon_eV"])
    continuum = np.asarray(srw["SRW_gaussian_J_per_eV_m2"])
    sample_ff = np.asarray(audit["sample_longitudinal_form_factor_power"])
    analytic_ff = np.asarray(audit["analytic_longitudinal_form_factor_power"])
    require(len(energies) == len(sample_ff), "Source/SRW energy grids differ")
    sample_theory = continuum * sample_ff / analytic_ff
    values = {}
    for plane in campaign["planes"]:
        name = plane["name"]
        values[name] = plane_spectra(
            root / "bend/detectors" / f"{name}.h5",
            root / "straight/detectors" / f"{name}.h5",
            energies, 1174.171072271572)
    final_name = campaign["planes"][-1]["name"]
    final = values[final_name]
    band = (energies >= 8.0) & (energies <= 12.0)
    de = float(energies[1] - energies[0])
    rows = []
    for index, energy in enumerate(energies):
        rows.append({
            "photon_eV": energy,
            "APRL_raw_axis_J_per_eV_m2": final["raw"][index],
            "APRL_low_angle_axis_J_per_eV_m2": final["low"][index],
            "SRW_continuum_J_per_eV_m2": continuum[index],
            "SRW_sample_longitudinal_FF_J_per_eV_m2": sample_theory[index],
            "near_nyquist_fraction": final["nyquist_fraction"][index],
            "outside_5_over_gamma_fraction": final["high_angle_fraction"][index],
            "inside_3_over_gamma_fraction": final["low_angle_fraction"][index],
        })
    with (output / "axis_comparison.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fixed = int(np.argmin(abs(energies - 9.825)))
    propagation = []
    for plane in campaign["planes"]:
        value = values[plane["name"]]
        propagation.append({
            "plane": plane["name"], "z_m": plane["z_m"],
            "raw_axis_at_9p825_J_per_eV_m2": float(value["raw"][fixed]),
            "low_axis_at_9p825_J_per_eV_m2": float(value["low"][fixed]),
            "near_nyquist_at_9p825": float(value["nyquist_fraction"][fixed]),
            "outside_5gamma_at_9p825": float(value["high_angle_fraction"][fixed]),
            "inside_3gamma_at_9p825": float(value["low_angle_fraction"][fixed]),
        })
    with (output / "propagation_9p825.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(propagation[0]))
        writer.writeheader()
        writer.writerows(propagation)

    metrics = {
        "raw_relative_L2_to_SRW_continuum_4_12eV":
            relative_l2(final["raw"], continuum),
        "low_angle_relative_L2_to_SRW_continuum_4_12eV":
            relative_l2(final["low"], continuum),
        "raw_relative_L2_to_sample_FF_theory_4_12eV":
            relative_l2(final["raw"], sample_theory),
        "low_angle_relative_L2_to_sample_FF_theory_4_12eV":
            relative_l2(final["low"], sample_theory),
        "raw_relative_L2_to_sample_FF_theory_8_12eV":
            relative_l2(final["raw"], sample_theory, band),
        "low_angle_relative_L2_to_sample_FF_theory_8_12eV":
            relative_l2(final["low"], sample_theory, band),
        "raw_band_8_12_J_per_m2": float(np.sum(final["raw"][band]) * de),
        "low_angle_band_8_12_J_per_m2": float(np.sum(final["low"][band]) * de),
        "SRW_continuum_band_8_12_J_per_m2": float(np.sum(continuum[band]) * de),
        "sample_FF_theory_band_8_12_J_per_m2": float(np.sum(sample_theory[band]) * de),
        "maximum_near_nyquist_fraction_4_12eV":
            float(np.max(final["nyquist_fraction"])),
        "near_nyquist_fraction_at_9p825eV":
            float(final["nyquist_fraction"][fixed]),
    }

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.5), layout="constrained")
    axes[0, 0].plot(energies, continuum * 1e9, "k-", lw=2.4,
                    label="dense SRW continuum")
    axes[0, 0].plot(energies, sample_theory * 1e9, color="#b279a2", lw=1.8,
                    label="SRW × actual longitudinal sample FF")
    axes[0, 0].plot(energies, final["raw"] * 1e9, color="#e45756", lw=1.5,
                    label="APRL raw axis")
    axes[0, 0].plot(energies, final["low"] * 1e9, color="#4c78a8", lw=1.8,
                    label=r"APRL $k_\perp/k<3/\gamma$")
    axes[0, 0].set(title="Absolute axis synchrotron spectrum",
                   ylabel=r"nJ / (eV m$^2$)")
    axes[0, 0].legend(fontsize=7)

    axes[0, 1].plot(energies, final["raw"] / sample_theory,
                    color="#e45756", label="raw / sample theory")
    axes[0, 1].plot(energies, final["low"] / sample_theory,
                    color="#4c78a8", label="low-angle / sample theory")
    axes[0, 1].axhline(1, color="k", lw=1)
    axes[0, 1].set(title="Direct ratio; no fitted normalization", ylabel="ratio")
    axes[0, 1].legend(fontsize=8)

    axes[1, 0].semilogy(energies, final["nyquist_fraction"], color="#e45756",
                       label="near Nyquist")
    axes[1, 0].semilogy(energies, final["high_angle_fraction"], color="#f58518",
                       label=r"outside $5/\gamma$")
    axes[1, 0].semilogy(energies, final["low_angle_fraction"], color="#4c78a8",
                       label=r"inside $3/\gamma$")
    axes[1, 0].set(title="Final transverse Ex power", ylabel="power fraction")
    axes[1, 0].legend(fontsize=8)

    z = [row["z_m"] for row in propagation]
    axes[1, 1].plot(z, [row["near_nyquist_at_9p825"] for row in propagation],
                    "o-", label="near Nyquist")
    axes[1, 1].plot(z, [row["outside_5gamma_at_9p825"] for row in propagation],
                    "s-", label=r"outside $5/\gamma$")
    axes[1, 1].plot(z, [row["inside_3gamma_at_9p825"] for row in propagation],
                    "^-", label=r"inside $3/\gamma$")
    axes[1, 1].set(title="9.825 eV propagation", xlabel="detector z [m]",
                   ylabel="power fraction")
    axes[1, 1].legend(fontsize=8)
    for axis in axes.flat:
        axis.grid(alpha=.2)
    for axis in axes[0]:
        axis.set_xlabel("Photon energy [eV]")
    axes[1, 0].set_xlabel("Photon energy [eV]")
    fig.suptitle("Independent unit-electron synchrotron/SRW comparison")
    fig.savefig(output / "independent_electron_synchrotron.png", dpi=200)
    fig.savefig(output / "independent_electron_synchrotron.pdf")
    plt.close(fig)

    result = {
        "complete": True,
        "status": "INDEPENDENT_UNIT_ELECTRON_SYNCHROTRON_COMPARISON",
        "particle_count": campaign["particle_count"],
        "sampling": campaign["sampling"],
        "absolute_comparison_no_fitted_normalization": True,
        "metrics": metrics,
        "propagation_at_9p825eV": propagation,
        "source_longitudinal_form_factor_audit": {
            key: audit[key] for key in audit
            if key.startswith("longitudinal_form_factor_power_")
        },
        "physical_accuracy_accepted": False,
        "acceptance_note": (
            "One finite-N realization tests the mechanism and absolute scale but does "
            "not by itself establish seed convergence or a production macro-weight limit."),
    }
    (output / "FINAL.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    require(len(sys.argv) == 4,
            "usage: analyze_independent_electron_synchrotron.py ROOT SRW_CSV OUTPUT")
    main(*sys.argv[1:])
