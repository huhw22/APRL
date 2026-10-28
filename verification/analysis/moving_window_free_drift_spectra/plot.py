#!/usr/bin/env python3
"""Compare saved free-drift and optional undulator spectra from small NPZ/JSON.

The frozen shared_plot.py supplies the existing validation and metric formulas.
This wrapper changes case identity and presentation, never the spectral arrays.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CASES = ("old_free_drift", "new_free_drift")
EXPECTED = {"old_free_drift": (18938, .5), "new_free_drift": (18993, .5)}
LABELS = {"old_free_drift": "18938 | earlier fixed window",
          "new_free_drift": "18993 | moving window"}
COLORS = {"old_free_drift": "#657b9e", "new_free_drift": "#cc6727"}
WINDOWS = ("rect", "hann")


def shared_module(path):
    spec = importlib.util.spec_from_file_location("frozen_shared_plot", path)
    shared = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(shared)
    shared.EXPECTED = EXPECTED
    shared.LABELS = LABELS
    shared.COLORS = COLORS
    return shared


def curve_label(run, window):
    return LABELS[run["name"]] + (" | rect" if window == "rect" else " | Hann")


def curve_style(shared, run, window):
    style = shared.style(window)
    if window == "rect" and run["name"] == "old_free_drift":
        style["ls"] = "-."
    return style


def flux_figure(shared, plt, runs, metrics, output):
    three = len(runs) == 3
    four = len(runs) == 4
    global_legend = three or four
    fig, axes = plt.subplots(2, 2, figsize=(15, 10.8))
    fig.subplots_adjust(left=.075, right=.975, top=.80 if global_legend else .85, bottom=.205,
                        hspace=.43, wspace=.25)
    high = max(r["arrays"]["photon_energy_eV"][-1] for r in runs)
    panels = ((0., high, "(a) Entire sampled band"),
              (0., 12., "(b) Low-frequency field"),
              (50., 115., "(c) 50–115 eV band"),
              (88., 95., "(d) 91.5 eV neighbourhood"))
    for ax, (lo, hi, title) in zip(axes.flat, panels):
        minimum, maximum = np.inf, -np.inf
        for run in runs:
            energy = run["arrays"]["photon_energy_eV"]
            mask = (energy > 0) & (energy >= lo) & (energy <= hi)
            shared.require(mask.any(), "No native bins in plotted interval")
            for window in WINDOWS:
                values = run["arrays"][f"{window}_native_flux_J_per_eV"][mask] * 1e9
                minimum = min(minimum, float(values.min()))
                maximum = max(maximum, float(values.max()))
                ax.plot(energy[mask], values, color=COLORS[run["name"]],
                        label=curve_label(run, window), **curve_style(shared, run, window))
        ax.set(xlim=(lo, hi), title=title, xlabel=r"Photon-energy coordinate $h\nu$ [eV]",
               ylabel=r"Native signed $dW/d(h\nu)$ [nJ/eV]")
        if minimum > 0:
            ax.set_yscale("log")
            ax.margins(y=.07)
        else:
            threshold = max(max(abs(minimum), abs(maximum)) * 1e-7, 1e-20)
            ax.set_yscale("symlog", linthresh=threshold)
            lower, upper = min(0., minimum * 1.25), max(0., maximum * 1.25)
            ax.set_ylim((lower, upper) if lower < upper else (-threshold, threshold))
            ax.axhline(0, color="#8c98a5", lw=.55)
        if lo <= 91.5 <= hi:
            ax.axvline(91.5, color="#5a6572", lw=.75, ls=":")
        ax.grid(alpha=.18, lw=.6)
        if not global_legend:
            ax.legend(frameon=False, fontsize=8.1, loc="best")
    if global_legend:
        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(.52, .891),
                   ncol=4 if four else 3, frameon=False, fontsize=9)
    heading = ("Field energy spectra: old and new undulator / free drift" if four else
               "Field energy spectra: old undulator and two free-drift runs" if three else
               "Free-drift field energy spectra at 0.5 m")
    shared.title(fig, heading,
                 "Same 140 × 140 physical aperture (0.49 mm²) | complete records | no subtraction" if global_legend
                 else "18993 moving window vs 18938 earlier fixed window | 140 × 140 physical pixels (0.49 mm²) | no subtraction")
    dc = "DC-bin energy [nJ], rect / Hann: " + "  |  ".join(
        f"{r['summary']['simulation_job']}: "
        f"{metrics[r['name']]['windows']['rect']['dc_native_energy_J']*1e9:.4g} / "
        f"{metrics[r['name']]['windows']['hann']['dc_native_energy_J']*1e9:.4g}"
        for r in runs)
    shared.footer(fig, [shared.sampling_note(runs),
        "Each curve keeps its own complete record and native FFT bins; no time cropping, rebinning or amplitude matching.",
        dc + ". DC remains in totals; only the isolated DC point is omitted from curves.",
        "Log scale where every displayed value is positive; signed symlog otherwise. No sign or peak is removed.",
        "Native lab E/B mix stored Yee time layers: this is total-field flux, not simultaneous-field or radiation-only energy.",
        "Rect: solid (old free drift dash-dot). Thin dashed: normalized Hann. Dotted vertical line: 91.5 eV reference.",
        "The sampled upper limit is Nyquist; high-energy bins do not establish physical or numerical resolution."])
    return shared.save(plt, fig, output, "energy_spectra_comparison")


def field_figure(shared, plt, runs, output):
    three = len(runs) == 3
    four = len(runs) == 4
    global_legend = three or four
    fig, axes = plt.subplots(2, 2, figsize=(15, 10.5))
    fig.subplots_adjust(left=.075, right=.975, top=.80 if global_legend else .85, bottom=.19,
                        hspace=.43, wspace=.25)
    high = max(r["arrays"]["photon_energy_eV"][-1] for r in runs)
    for i, ax in enumerate(axes.flat):
        kind = ("E", "cB")[i % 2]
        lo, hi = (0., high) if i < 2 else (88., 95.)
        for run in runs:
            energy = run["arrays"]["photon_energy_eV"]
            for window in WINDOWS:
                values = run["arrays"][f"{window}_{kind}_psd_V2_m2_per_eV"].sum(axis=1)
                mask = (energy > 0) & (energy >= lo) & (energy <= hi) & (values > 0)
                ax.plot(energy[mask], values[mask], color=COLORS[run["name"]],
                        label=curve_label(run, window), **curve_style(shared, run, window))
        ax.set(xlim=(lo, hi), yscale="log", title=f"({'abcd'[i]}) {kind}: " +
               ("entire sampled band" if i < 2 else "91.5 eV neighbourhood"),
               xlabel=r"Photon-energy coordinate $h\nu$ [eV]",
               ylabel=r"Pixel-averaged squared-field density [(V/m)$^2$/eV]")
        ax.axvline(91.5, color="#5a6572", lw=.75, ls=":")
        ax.grid(alpha=.18, lw=.6)
        if not global_legend:
            ax.legend(frameon=False, fontsize=8.1, loc="best")
    if global_legend:
        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(.52, .891),
                   ncol=4 if four else 3, frameon=False, fontsize=9)
    heading = ("Electric and magnetic spectra: four simulation cases" if four else
               "Electric and magnetic spectra: three simulation cases" if three else
               "Free-drift electric and magnetic field spectra")
    subtitle = ("Old/new undulator at 4.5 m and old/new free drift at 0.5 m | same physical aperture" if four else
                "Old undulator at 4.5 m and two free-drift runs at 0.5 m | same physical aperture" if three else
                "18993 moving window vs 18938 earlier fixed window | sum of three components | same physical aperture")
    shared.title(fig, heading, subtitle)
    shared.footer(fig, [shared.sampling_note(runs),
        "E and cB densities are pixel- and record-time-averaged mean squares; they are not radiation energy spectra.",
        "Rect: solid (old free drift dash-dot). Thin dashed: Hann / sqrt(mean(window²)); window sensitivity is not an error bar.",
        "DC is retained in data and metrics. Only DC and exact zeros are omitted from these logarithmic curves.",
        "No mean removal, time interpolation, background subtraction, spectral notch or native time-layer correction.",
        "Nyquist bounds sampled bandwidth; it does not establish high-energy physical or numerical resolution."])
    return shared.save(plt, fig, output, "field_spectra_comparison")


def main():
    global CASES, EXPECTED, LABELS, COLORS
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", "--root", dest="input_dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--shared-plot", type=Path, default=HERE / "shared_plot.py")
    parser.add_argument("--include-old-undulator", action="store_true",
                        help="Add the saved 18932/4.5m spectrum under old_undulator/; preserve native records")
    parser.add_argument("--include-new-undulator", action="store_true",
                        help="Add 18992/4.5m under new_undulator/ and imply --include-old-undulator; four complete native records")
    args = parser.parse_args()
    args.include_old_undulator |= args.include_new_undulator
    if args.include_old_undulator:
        CASES = (*CASES, "old_undulator")
        EXPECTED = {**EXPECTED, "old_undulator": (18932, 4.5)}
        LABELS = {"old_free_drift": "18938 | old free drift | 0.5 m",
                  "new_free_drift": "18993 | new free drift | 0.5 m",
                  "old_undulator": "18932 | old undulator | 4.5 m"}
        COLORS = {**COLORS, "old_undulator": "#a13c72"}
    if args.include_new_undulator:
        CASES = (*CASES, "new_undulator")
        EXPECTED = {**EXPECTED, "new_undulator": (18992, 4.5)}
        LABELS = {**LABELS, "new_undulator": "18992 | new undulator | 4.5 m"}
        COLORS = {**COLORS, "new_undulator": "#173b8f"}
    root = args.input_dir.resolve()
    output = args.output_dir.resolve() if args.output_dir else root / (
        "figures_four_cases" if args.include_new_undulator else
        "figures_three_cases" if args.include_old_undulator else "figures")
    shared = shared_module(args.shared_plot.resolve())
    runs = [shared.load_run(root, name) for name in CASES]
    shared.require(all(np.isclose(runs[0]["summary"]["physical_area_m2"],
                             r["summary"]["physical_area_m2"], rtol=2e-13, atol=0) for r in runs),
                   "Different physical apertures")
    output.mkdir(parents=True, exist_ok=True)
    metrics = {r["name"]: shared.metrics_for(r) for r in runs}
    time_keys = ("samples", "sample_interval_s", "time_start_s", "time_end_s")
    time_metadata_equal = all(runs[0]["summary"][key] == runs[1]["summary"][key]
                              for key in time_keys)
    frequency_arrays_equal = np.array_equal(runs[0]["arrays"]["frequency_Hz"],
                                            runs[1]["arrays"]["frequency_Hz"])
    by_name = {r["name"]: r for r in runs}
    pair_time_comparisons = {}
    for pair in ("free_drift", "undulator"):
        old_name, new_name = "old_" + pair, "new_" + pair
        if old_name not in by_name or new_name not in by_name:
            continue
        old, new = by_name[old_name], by_name[new_name]
        pair_time_comparisons[pair] = {
            "old_case": old_name, "new_case": new_name,
            "old": {key: old["summary"][key] for key in time_keys},
            "new": {key: new["summary"][key] for key in time_keys},
            "old_record_period_Ndt_s": old["duration"],
            "new_record_period_Ndt_s": new["duration"],
            "old_frequency_bin_eV": old["de"], "new_frequency_bin_eV": new["de"],
            "full_record_time_metadata_equal": all(old["summary"][key] == new["summary"][key]
                                                   for key in time_keys),
            "frequency_arrays_equal": bool(np.array_equal(old["arrays"]["frequency_Hz"],
                                                            new["arrays"]["frequency_Hz"])),
            "scope": "Saved count/cadence/endpoints and native frequency arrays; no raw timestamp-vector check, cropping or rebinning by this plotter."}
    record = {"format": "aprl-moving-free-drift-spectral-comparison-v1", "runs": metrics,
        "scope": "Each complete recorded raw lab field, physical aperture only; native mixed time layers. Old free drift is a reference, not a subtracted background.",
        "comparison_time_scope": "Complete recorded series and native frequency bins, without cropping or interpolation. Time-metadata equality below checks saved summaries; raw timestamp identity belongs to the source audit.",
        "full_record_time_metadata_equal": time_metadata_equal,
        "frequency_arrays_equal": bool(frequency_arrays_equal),
        "time_equality_check_scope": "Top-level equality flags refer only to the old/new 0.5m free-drift pair; pair_time_comparisons reports each available old/new pair separately.",
        "pair_time_comparisons": pair_time_comparisons,
        "old_undulator_included": args.include_old_undulator,
        "new_undulator_included": args.include_new_undulator,
        "old_undulator_record_note": "18932 uses its own 150.795822fs record and native bins; free drift uses 135.179018fs. Same Hann definition, different sample count." if args.include_old_undulator else None,
        "integration_rule": "Sum density times native deltaE for lower <= bin centre < upper; no partial-bin interpolation.",
        "power_definition": "Band energy divided by that run's own Ndt; record-average, not instantaneous peak or repetition-rate average.",
        "hann_definition": "Symmetric Hann / sqrt(mean(w^2)); separate finite-record weighting, not missing-tail recovery or error bars.",
        "input_hdf5_opened_by_plotter": False, "script_sha256": shared.sha(__file__),
        "shared_plot_sha256": shared.sha(args.shared_plot), "figures": []}
    with (output / "BANDS.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("case", "job", "plane_m", "window", "band", "bin_count", "first_bin_eV", "last_bin_eV",
                         "native_energy_J", "record_average_power_W", "record_period_s"))
        for name, item in metrics.items():
            for window in WINDOWS:
                for band, value in item["windows"][window]["bands"].items():
                    writer.writerow((name, item["simulation_job"], item["plane_z_m"], window, band,
                        value["bins"], value["first_bin_center_eV"], value["last_bin_center_eV"],
                        value["native_energy_J"], value["record_average_power_W"], item["fft_period_fs"] * 1e-15))
    plt = shared.plotting()
    figures = flux_figure(shared, plt, runs, metrics, output) + field_figure(shared, plt, runs, output)
    record["figures"] = [{"path": str(path), "sha256": shared.sha(path)} for path in figures]
    (output / "METRICS.json").write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"figures": [str(path) for path in figures], "metrics": str(output / "METRICS.json")}))


if __name__ == "__main__":
    main()
