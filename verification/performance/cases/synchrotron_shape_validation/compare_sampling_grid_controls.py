"""Compare closed raw-axis sampling and transverse-grid controls."""
import csv
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analyze_notch_propagation import plane_spectra, require


def relative_l2(value, reference, selected):
    return float(np.linalg.norm(value[selected] - reference[selected]) /
                 np.linalg.norm(reference[selected]))


def ripple(values, energies):
    require(np.all(values > 0), "Nonpositive raw-axis spectrum")
    trend = np.exp(np.polyval(np.polyfit(energies, np.log(values), 3),
                              energies))
    ratio = values / trend
    return float(np.std(ratio - 1)), ratio


def main(srw_csv, output, specifications):
    srw_csv, output = Path(srw_csv).resolve(), Path(output).resolve()
    require(not output.exists(), "Output already exists")
    output.mkdir(parents=True)
    srw = np.genfromtxt(srw_csv, delimiter=",", names=True)
    energies = np.asarray(srw["photon_eV"])
    curves = []
    rows = []
    for specification in specifications:
        require("=" in specification, "Each case must be LABEL=ROOT")
        label, root_text = specification.split("=", 1)
        root = Path(root_text).resolve()
        require((root / "PIC_COMPLETE.json").is_file(),
                f"Incomplete PIC campaign: {root}")
        campaign = json.loads((root / "campaign.json").read_text())
        audit = json.loads((root / "SOURCE_AUDIT.json").read_text())
        represented = float(campaign["represented_electrons"])
        analytic_ff = np.asarray(audit[
            "analytic_longitudinal_form_factor_power"])
        sample_ff = np.asarray(audit[
            "sample_longitudinal_form_factor_power"])
        continuum = (np.asarray(srw[
            "SRW_point_per_electron_J_per_eV_m2"]) * represented**2 *
            np.asarray(srw["longitudinal_form_factor_squared"]) *
            np.asarray(srw["transverse_correction_interpolated"]))
        sample_theory = continuum * sample_ff / analytic_ff
        plane = campaign["planes"][-1]["name"]
        values = plane_spectra(root / "bend/detectors" / f"{plane}.h5",
                               root / "straight/detectors" / f"{plane}.h5",
                               energies, 1174.171072271572)
        raw = values["raw"]
        ripple_rms, detrended = ripple(raw, energies)
        row = {"label": label, "root": str(root),
               "represented_electrons": represented,
               "ripple_rms_4_12eV": ripple_rms}
        for low, high in ((4., 8.), (8., 12.)):
            selected = (energies >= low) & (energies <= high)
            prefix = f"{low:g}_{high:g}eV"
            row[prefix + "_relative_L2_to_actual_sample_theory"] = \
                relative_l2(raw, sample_theory, selected)
            row[prefix + "_integral_ratio_to_actual_sample_theory"] = float(
                np.trapz(raw[selected], energies[selected]) /
                np.trapz(sample_theory[selected], energies[selected]))
            row[prefix + "_mean_pointwise_ratio_to_actual_sample_theory"] = \
                float(np.mean(raw[selected] / sample_theory[selected]))
        rows.append(row)
        curves.append({"label": label, "represented": represented,
                       "raw": raw, "continuum": continuum,
                       "sample_theory": sample_theory,
                       "detrended": detrended})

    with (output / "sampling_grid_metrics.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.5),
                             layout="constrained")
    theory_per_pair = curves[0]["continuum"] / curves[0]["represented"]**2
    axes[0, 0].plot(energies, theory_per_pair * 1e16, "k-", lw=2.3,
                    label=r"SRW Gaussian / $N^2$")
    colors = plt.cm.tab10(np.linspace(0, .8, len(curves)))
    for curve, color in zip(curves, colors):
        axes[0, 0].plot(energies,
                        curve["raw"] / curve["represented"]**2 * 1e16,
                        lw=1.5, color=color, label=curve["label"])
        axes[0, 1].plot(energies, curve["raw"] / curve["sample_theory"],
                        lw=1.5, color=color, label=curve["label"])
        axes[1, 0].plot(energies, curve["detrended"], lw=1.5,
                        color=color, label=curve["label"])
    axes[0, 0].set(title="Absolute coherent scale (explicit N² division)",
                   ylabel=r"$10^{-16}$ J / (eV m$^2$ electron$^2$)")
    axes[0, 1].axhline(1, color="black", lw=1)
    axes[0, 1].set(title="Raw APRL / theory using each actual sample FF",
                   ylabel="absolute ratio")
    axes[1, 0].axhline(1, color="black", lw=1)
    axes[1, 0].set(title="Raw-axis ripple after cubic log-trend removal",
                   ylabel="detrended ratio")
    indices = np.arange(len(rows))
    width = .36
    axes[1, 1].bar(indices - width / 2,
                   [row["4_8eV_integral_ratio_to_actual_sample_theory"]
                    for row in rows], width, label="4–8 eV")
    axes[1, 1].bar(indices + width / 2,
                   [row["8_12eV_integral_ratio_to_actual_sample_theory"]
                    for row in rows], width, label="8–12 eV")
    axes[1, 1].axhline(1, color="black", lw=1)
    axes[1, 1].set(title="Signed band integral ratio", ylabel="APRL / theory",
                   xticks=indices, xticklabels=[row["label"] for row in rows])
    axes[1, 1].tick_params(axis="x", labelrotation=20)
    axes[1, 1].legend(fontsize=8)
    for axis in axes.flat:
        axis.grid(alpha=.2)
    for axis in (axes[0, 0], axes[0, 1], axes[1, 0]):
        axis.set_xlabel("photon energy [eV]")
        axis.legend(fontsize=7)
    fig.suptitle("Raw-axis synchrotron sampling/grid controls; no angular filter or fitted normalization")
    fig.savefig(output / "sampling_grid_controls.png", dpi=200)
    fig.savefig(output / "sampling_grid_controls.pdf")
    plt.close(fig)

    result = {
        "complete": True,
        "raw_axis_spectrum": True,
        "angular_filter_applied": False,
        "aperture_average_applied": False,
        "fitted_normalization_applied": False,
        "N_squared_display_scaling_is_explicit": True,
        "physical_accuracy_accepted": False,
        "acceptance_note": (
            "Visual smoothness is insufficient: acceptance requires both seed "
            "stability and coupled particle-count/transverse-grid convergence."),
        "rows": rows,
    }
    (output / "FINAL.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    require(len(sys.argv) >= 5,
            "usage: compare_sampling_grid_controls.py SRW_CSV OUTPUT LABEL=ROOT ...")
    main(sys.argv[1], sys.argv[2], sys.argv[3:])
