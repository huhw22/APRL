#!/usr/bin/env python3
"""Plot the small spectral reductions of 18932 and 18938; never open raw H5.

INPUT/{undulator,free_drift}/{spectra.npz,summary.json} are produced by reduce.py.
The stored lab E/B retain their native, mixed Yee time layers. Spectral flux is
signed, includes the total field, and is not labelled as radiation-only energy.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

H_EV_S = 4.135667696e-15
COLORS = {"undulator": "#2563a6", "free_drift": "#cc6727"}
LABELS = {"undulator": "Undulator | 4.5 m | 18932",
          "free_drift": "Free drift | 0.5 m | 18938"}
EXPECTED = {"undulator": (18932, 4.5), "free_drift": (18938, 0.5)}
TEXT, MUTED = "#203040", "#536373"
WINDOWS = ("rect", "hann")
BANDS = (("all", 0., None), ("below_1_eV", 0., 1.),
         ("90_to_93_eV", 90., 93.), ("91p3_to_91p7_eV", 91.3, 91.7))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def ratio(value, denominator):
    return float(value / denominator) if denominator != 0 else None


def load_run(root, name):
    directory = root / name
    sp, ap = directory / "summary.json", directory / "spectra.npz"
    summary = json.loads(sp.read_text())
    with np.load(ap, allow_pickle=False) as archive:
        arrays = {key: archive[key] for key in archive.files}
    n = int(summary["samples"])
    require(n >= 3, f"{name}: too few samples")
    dt = float(summary["sample_interval_s"])
    require(dt > 0 and np.isfinite(dt), f"{name}: invalid sample interval")
    duration = n * dt
    de = H_EV_S / duration
    energy, frequency = arrays["photon_energy_eV"], arrays["frequency_Hz"]
    require(energy.shape == frequency.shape == (n // 2 + 1,), f"{name}: FFT extent")
    require(np.isfinite(energy).all() and np.isfinite(frequency).all(), f"{name}: nonfinite axis")
    require(energy[0] == frequency[0] == 0 and np.all(np.diff(energy) > 0),
            f"{name}: non-increasing axis")
    require(np.allclose(energy, H_EV_S * frequency, rtol=2e-13, atol=0),
            f"{name}: energy/frequency conversion")
    require(np.allclose(np.diff(energy), de, rtol=3e-10, atol=0), f"{name}: nonuniform bins")
    require(np.isclose(de, summary["frequency_bin_eV"], rtol=2e-13), f"{name}: bin width")
    require(np.isclose(duration * 1e15, summary["fft_period_fs"], rtol=2e-13),
            f"{name}: N*dt period")
    require(summary["maximum_timestamp_deviation_samples"] <= 1e-4, f"{name}: irregular samples")
    require(summary["physical_shape_yx"] == [140, 140], f"{name}: expected physical aperture")
    attrs = summary["source_attributes"]
    job, plane = EXPECTED[name]
    require(summary["simulation_job"] == job, f"{name}: wrong simulation job")
    require(np.isclose(attrs["plane_z_m"], plane, rtol=0, atol=1e-14), f"{name}: wrong plane")
    require(attrs["frame"] == "laboratory" and attrs["external_background_subtracted"] == "none",
            f"{name}: expected raw lab fields")
    require(summary["input_transaction"]["complete"] == 1 and
            summary["input_transaction"]["committed_samples"] == n and
            summary["input_transaction"]["field_shape"] == [n, 160, 160] and
            summary["input_transaction"]["time_shape"] == [n], f"{name}: transaction mismatch")
    require(summary["input_unchanged_during_read"], f"{name}: source changed during scan")
    expected_sampling = "global Yee trilinear stencil; adjacent MPI half-cell planes included"
    require(summary["input_transaction"]["spatial_sampling"] == expected_sampling,
            f"{name}: wrong detector sampling implementation")
    require(np.isclose(summary["physical_area_m2"], 140 * 140 * attrs["dx_m"] * attrs["dy_m"],
                       rtol=2e-13, atol=0), f"{name}: aperture area")
    checks = {}
    for window in WINDOWS:
        review = summary["windows"][window]
        integrals = []
        for kind in ("E", "cB"):
            values = arrays[f"{window}_{kind}_psd_V2_m2_per_eV"]
            require(values.shape == (energy.size, 3) and np.isfinite(values).all() and
                    np.all(values >= 0), f"{name}/{window}/{kind}: invalid density")
            integrals.extend((values.sum(axis=0) * de).tolist())
        time_ms = np.asarray(review["mean_square_components_E_then_cB"])
        rel = np.abs(np.asarray(integrals) - time_ms) / np.maximum(np.abs(time_ms), 1e-300)
        require(np.all(rel < 2e-10), f"{name}/{window}: mean-square Parseval")
        flux = arrays[f"{window}_native_flux_J_per_eV"]
        require(flux.shape == energy.shape and np.isfinite(flux).all(),
                f"{name}/{window}: invalid signed flux")
        total = float(flux.sum() * de)
        scale = max(float(np.abs(flux).sum() * de), abs(review["native_flux_time_sum_J"]), 1e-300)
        err = abs(total - review["native_flux_time_sum_J"]) / scale
        require(err < 2e-10, f"{name}/{window}: native flux Parseval")
        dc = float(arrays[f"{window}_E_psd_V2_m2_per_eV"][0].sum() * de)
        require(np.isclose(dc, review["dc_E_mean_square_V2_m2"], rtol=2e-12, atol=0),
                f"{name}/{window}: DC bin")
        checks[window] = {"mean_square_relative_errors": rel.tolist(),
                          "native_flux_relative_error": err}
    return dict(name=name, summary=summary, arrays=arrays, duration=duration, de=de,
                source_hashes={"summary.json": sha(sp), "spectra.npz": sha(ap)},
                parseval_checks=checks)


def metrics_for(run):
    a, s = run["arrays"], run["summary"]
    energy, de, duration = a["photon_energy_eV"], run["de"], run["duration"]
    result = {"simulation_job": s["simulation_job"], "source_file": s["source_file"],
              "plane_z_m": s["source_attributes"]["plane_z_m"], "samples": s["samples"],
              "recorded_span_fs": s["recorded_span_fs"], "fft_period_fs": s["fft_period_fs"],
              "frequency_bin_eV": de, "nyquist_eV": s["theoretical_nyquist_eV"],
              "physical_area_m2": s["physical_area_m2"],
              "batch_exit_code": s["batch_exit_code"], "input_hashes": run["source_hashes"],
              "parseval_checks": run["parseval_checks"], "windows": {}}
    for window in WINDOWS:
        flux = a[f"{window}_native_flux_J_per_eV"]
        e = a[f"{window}_E_psd_V2_m2_per_eV"].sum(axis=1)
        b = a[f"{window}_cB_psd_V2_m2_per_eV"].sum(axis=1)
        total_flux, total_e, total_b = (float(v.sum() * de) for v in (flux, e, b))
        bands = {}
        for label, lo, hi in BANDS:
            mask = energy >= lo
            if hi is not None:
                mask &= energy < hi
            indices = np.flatnonzero(mask)
            require(indices.size > 0, f"{run['name']}: empty band {label}")
            integral = float(flux[mask].sum() * de)
            bands[label] = {"requested_lower_eV": lo, "requested_upper_eV_exclusive": hi,
                "bins": int(indices.size), "first_bin_center_eV": float(energy[indices[0]]),
                "last_bin_center_eV": float(energy[indices[-1]]), "native_energy_J": integral,
                "record_average_power_W": integral / duration,
                "signed_native_energy_fraction": ratio(integral, total_flux),
                "E_mean_square_fraction": ratio(float(e[mask].sum() * de), total_e),
                "cB_mean_square_fraction": ratio(float(b[mask].sum() * de), total_b)}
        peaks = {}
        for label, lo, hi in (("target_50_115_eV", 50., 115.), ("zoom_88_95_eV", 88., 95.)):
            indices = np.flatnonzero((energy >= lo) & (energy <= hi))
            for kind, values in (("native_flux", flux), ("E_density", e)):
                i = int(indices[np.argmax(values[indices])])
                peaks[f"{label}_{kind}"] = {"energy_eV": float(energy[i]), "bin": i,
                    "value": float(values[i]),
                    "unit": "J/eV" if kind == "native_flux" else "V^2/m^2/eV",
                    "meaning": "Maximum stored ordinate in stated interval; not an identified radiation line"}
        result["windows"][window] = {"bands": bands, "peaks": peaks,
            "dc_native_energy_J": float(flux[0] * de),
            "dc_native_energy_fraction": ratio(float(flux[0] * de), total_flux),
            "dc_E_mean_square_fraction": ratio(float(e[0] * de), total_e),
            "negative_native_flux_bins": int(np.count_nonzero(flux < 0)),
            "negative_native_energy_J": float(flux[flux < 0].sum() * de),
            "native_total_energy_J": total_flux, "record_average_total_power_W": total_flux / duration}
    result["hann_vs_rect"] = {label: {
        "native_energy_ratio": ratio(result["windows"]["hann"]["bands"][label]["native_energy_J"],
                                      result["windows"]["rect"]["bands"][label]["native_energy_J"]),
        "meaning": "Window-weighting sensitivity, not correction or uncertainty interval"}
        for label, _, _ in BANDS}
    return result


def plotting():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
        "axes.titlesize": 11.5, "axes.labelsize": 10, "text.color": TEXT,
        "axes.labelcolor": TEXT, "xtick.color": MUTED, "ytick.color": MUTED,
        "axes.edgecolor": "#aeb9c4", "axes.spines.top": False, "axes.spines.right": False,
        "figure.facecolor": "white", "savefig.facecolor": "white", "pdf.fonttype": 42,
        "path.simplify": False})
    return plt


def title(fig, main, sub):
    fig.text(.07, .97, main, va="top", fontsize=18, weight="bold")
    fig.text(.07, .928, sub, va="top", fontsize=10.2, color=MUTED)


def footer(fig, lines):
    fig.text(.07, .025, "\n".join(lines), va="bottom", color=MUTED,
             fontsize=8.7, linespacing=1.42)


def sampling_note(runs):
    return "  |  ".join(f"{r['summary']['simulation_job']}: T = {r['duration']*1e15:.3f} fs, "
                        f"ΔE = {r['de']:.5f} eV" for r in runs)


def style(window):
    return dict(ls="-" if window == "rect" else "--", lw=1.05 if window == "rect" else .8,
                alpha=1.0 if window == "rect" else .85)


def label(run, window):
    return LABELS[run["name"]] + (" | rect" if window == "rect" else " | Hann")


def save(plt, fig, output, basename):
    files = []
    for ext in ("png", "pdf"):
        p = output / f"{basename}.{ext}"
        fig.savefig(p, dpi=190, bbox_inches="tight", pad_inches=.12)
        files.append(p)
    plt.close(fig)
    return files


def flux_figure(plt, runs, metrics, output):
    fig, axes = plt.subplots(2, 2, figsize=(15, 10.8))
    fig.subplots_adjust(left=.075, right=.975, top=.85, bottom=.205, hspace=.43, wspace=.25)
    high = max(r["arrays"]["photon_energy_eV"][-1] for r in runs)
    ranges = ((0., high, "(a) Entire sampled band"), (0., 12., "(b) Low-frequency field"),
              (50., 115., "(c) Target band"), (88., 95., "(d) 91.5 eV neighbourhood"))
    for ax, (lo, hi, panel) in zip(axes.flat, ranges):
        minimum, maximum = np.inf, -np.inf
        for run in runs:
            energy = run["arrays"]["photon_energy_eV"]
            mask = (energy > 0) & (energy >= lo) & (energy <= hi)
            for window in WINDOWS:
                value = run["arrays"][f"{window}_native_flux_J_per_eV"][mask] * 1e9
                minimum = min(minimum, float(np.min(value)))
                maximum = max(maximum, float(np.max(value)))
                ax.plot(energy[mask], value, color=COLORS[run["name"]],
                        label=label(run, window), **style(window))
        ax.set(xlim=(lo, hi), title=panel, xlabel=r"Photon-energy coordinate $h\nu$ [eV]",
               ylabel=r"Native signed $dW/d(h\nu)$ [nJ/eV]")
        if minimum > 0:
            # Do not allocate an empty negative region for positive data.
            # This decision inspects every displayed run/window without
            # deleting, clipping or taking absolute values of any bin.
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
        ax.legend(frameon=False, fontsize=8.1, loc="best")
    title(fig, "Stored-field energy spectra: undulator vs free drift",
          "Same 140 × 140 physical aperture (0.49 mm²) | per-pixel FFT, then area sum | no field or spectral subtraction")
    dc = "DC-bin energy [nJ], rect / Hann: " + "  |  ".join(
        f"{EXPECTED[r['name']][0]}: {metrics[r['name']]['windows']['rect']['dc_native_energy_J']*1e9:.4g} / "
        f"{metrics[r['name']]['windows']['hann']['dc_native_energy_J']*1e9:.4g}" for r in runs)
    all_positive = all(np.all(r["arrays"][f"{w}_native_flux_J_per_eV"] > 0)
                       for r in runs for w in WINDOWS)
    scale_note = ("All stored flux bins are positive in both runs/windows; logarithmic axes."
                  if all_positive else
                  "Log axes where all displayed bins are positive; signed symlog otherwise. No sign is discarded.")
    footer(fig, [sampling_note(runs) + "; each curve retains its own native FFT grid.",
        dc + ". DC remains in all totals; only its isolated bin is omitted from the curves.",
        scale_note + " Solid: rectangle; dashed: normalized Hann.",
        "Native lab E/B mix stored Yee time layers: this is a total-field flux estimate, not simultaneous-field or radiation-only energy.",
        "Dotted 91.5 eV line is a comparison reference. Window agreement is not proof of radiation purity; no mean, trend or peak was removed."])
    return save(plt, fig, output, "energy_spectra_comparison")


def field_figure(plt, runs, output):
    fig, axes = plt.subplots(2, 2, figsize=(15, 10.5))
    fig.subplots_adjust(left=.075, right=.975, top=.85, bottom=.19, hspace=.43, wspace=.25)
    high = max(r["arrays"]["photon_energy_eV"][-1] for r in runs)
    for i, ax in enumerate(axes.flat):
        kind = ("E", "cB")[i % 2]
        lo, hi = (0., high) if i < 2 else (88., 95.)
        for run in runs:
            energy = run["arrays"]["photon_energy_eV"]
            mask = (energy > 0) & (energy >= lo) & (energy <= hi)
            for window in WINDOWS:
                value = run["arrays"][f"{window}_{kind}_psd_V2_m2_per_eV"].sum(axis=1)[mask] / 1e12
                ax.plot(energy[mask], np.ma.masked_less_equal(value, 0),
                        color=COLORS[run["name"]], label=label(run, window), **style(window))
        ax.set(xlim=(lo, hi), yscale="log", title=f"({'abcd'[i]}) Total {kind} mean-square spectrum",
               xlabel=r"Photon-energy coordinate $h\nu$ [eV]",
               ylabel=r"Field PSD [$(\mathrm{MV/m})^2$/eV]")
        ax.axvline(91.5, color="#5a6572", lw=.75, ls=":")
        ax.grid(alpha=.18, lw=.6)
        ax.legend(frameon=False, fontsize=8.1, loc="best")
    title(fig, "Electric and magnetic field spectra",
          "Each pixel FFT before squared-amplitude averaging | E and cB in the same field units | full physical aperture")
    footer(fig, [sampling_note(runs) + "; normalization uses each record's own T = NΔt.",
        "These densities integrate to space/time mean-square field, not pulse energy. E and cB spectra alone do not separate bound fields.",
        "Solid: raw rectangular window; thin dashed: energy-normalized Hann. DC is retained in the data and reported in METRICS.json.",
        "Only exact zeros/DC are omitted from logarithmic curves. Stored native time layers remain; no background or spectral peak subtraction."])
    return save(plt, fig, output, "field_spectra_comparison")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-dir", type=Path, required=True)
    p.add_argument("--output-dir", type=Path)
    args = p.parse_args()
    root = args.input_dir.resolve()
    output = args.output_dir.resolve() if args.output_dir else root / "figures"
    runs = [load_run(root, n) for n in EXPECTED]
    require(np.isclose(runs[0]["summary"]["physical_area_m2"], runs[1]["summary"]["physical_area_m2"],
                       rtol=2e-13, atol=0), "Different physical apertures")
    output.mkdir(parents=True, exist_ok=True)
    metrics = {r["name"]: metrics_for(r) for r in runs}
    record = {"format": "aprl-undulator-free-drift-spectral-comparison-v1", "runs": metrics,
        "scope": "Full recorded raw lab fields, physical aperture only, no field/spectrum subtraction; native mixed time layers.",
        "integration_rule": "Sum bin-centre-selected density times native ΔE; finite bands use lower <= centre < upper. No interpolation or edge-fraction model.",
        "power_definition": "Band energy / each full record's own T=N*dt; not instantaneous peak power or repetition-rate average.",
        "hann_definition": "Symmetric Hann / sqrt(mean(w^2)); separately weighted record, not missing-tail recovery or error bar.",
        "input_hdf5_opened_by_plotter": False, "script_sha256": sha(__file__),
        "font": "DejaVu Sans; no CJK font available on this host", "figures": []}
    with (output / "BANDS.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(("case", "job", "plane_m", "window", "band", "bin_count", "first_bin_eV", "last_bin_eV",
                         "native_energy_J", "record_average_power_W", "record_period_s"))
        for name, item in metrics.items():
            for window in WINDOWS:
                for band, v in item["windows"][window]["bands"].items():
                    writer.writerow((name, item["simulation_job"], item["plane_z_m"], window, band,
                        v["bins"], v["first_bin_center_eV"], v["last_bin_center_eV"], v["native_energy_J"],
                        v["record_average_power_W"], item["fft_period_fs"] * 1e-15))
    plt = plotting()
    paths = flux_figure(plt, runs, metrics, output) + field_figure(plt, runs, output)
    record["figures"] = [{"path": str(p), "sha256": sha(p)} for p in paths]
    (output / "METRICS.json").write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"figures": [str(p) for p in paths], "metrics": str(output / "METRICS.json")}, indent=2))


if __name__ == "__main__":
    main()
