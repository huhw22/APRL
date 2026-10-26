#!/usr/bin/env python3
"""Plot existing whole-aperture E/cB spectral reductions, without reading HDF5.

INPUT/{retirement,original}/{spectra.npz,summary.json} must be produced by
analyze_spectra.py. Raw rectangular-window spectra are primary; normalized
Hann spectra expose window sensitivity. No mean, spike, or band is removed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


H_EV_S = 4.135667696e-15
COLORS = {"retirement": "#2563a6", "original": "#cc6727"}
LABELS = {"retirement": "18893", "original": "Original"}
COMPONENT_COLORS = ("#2563a6", "#d07822", "#38815e")
TEXT, MUTED = "#203040", "#536373"
PSD_UNIT = r"Field PSD [$(\mathrm{MV/m})^2$/eV]"
NOMINAL_EV = 91.481643


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load_run(root, name):
    directory = root / name
    summary_path, array_path = directory / "summary.json", directory / "spectra.npz"
    summary = json.loads(summary_path.read_text())
    with np.load(array_path, allow_pickle=False) as archive:
        arrays = {key: archive[key] for key in archive.files}
    n = int(summary["samples"])
    nf = n // 2 + 1
    energy, frequency = (np.asarray(arrays[k], dtype=float)
                         for k in ("photon_energy_eV", "frequency_Hz"))
    require(energy.shape == frequency.shape == (nf,), f"{name}: inconsistent FFT extent")
    require(np.isfinite(energy).all() and np.isfinite(frequency).all(), f"{name}: nonfinite axis")
    require(energy[0] == frequency[0] == 0 and np.all(np.diff(energy) > 0),
            f"{name}: expected increasing nonnegative FFT axis")
    require(np.allclose(energy, H_EV_S * frequency, rtol=2e-13, atol=0),
            f"{name}: photon-energy/frequency mismatch")
    de = float(summary["frequency_bin_eV"])
    require(np.allclose(np.diff(energy), de, rtol=2e-10, atol=0), f"{name}: irregular frequency bins")
    require(summary["physical_shape_yx"] == [140, 140], f"{name}: expected physical 140x140 aperture")
    require(summary["source_attributes"]["frame"] == "laboratory", f"{name}: unexpected frame")
    require(summary["source_attributes"]["external_background_subtracted"] == "none",
            f"{name}: expected raw, unsubtracted fields")
    require(np.isclose(energy[-1], summary["last_frequency_eV"], rtol=2e-13, atol=0),
            f"{name}: last frequency mismatch")
    for window in ("rect", "hann"):
        integrals = []
        for kind in ("E", "cB"):
            key = f"{window}_{kind}_psd_V2_m2_per_eV"
            value = np.asarray(arrays[key], dtype=float)
            require(value.shape == (nf, 3), f"{name}/{key}: shape mismatch")
            require(np.isfinite(value).all() and np.all(value >= 0), f"{name}/{key}: invalid PSD")
            integrals.extend((value.sum(axis=0) * de).tolist())
        review = summary["windows"][window]
        require(np.allclose(integrals, review["mean_square_components_E_then_cB"],
                            rtol=2e-10, atol=0), f"{name}/{window}: PSD integral disagrees with time mean square")
        flux = arrays[f"{window}_native_flux_J_per_eV"]
        require(flux.shape == (nf,) and np.isfinite(flux).all(), f"{name}/{window}: invalid signed flux")
        scale = max(float(np.abs(flux).sum() * de), abs(review["native_flux_spectral_sum_J"]), 1e-300)
        require(abs(float(flux.sum() * de) - review["native_flux_spectral_sum_J"]) <= 2e-10 * scale,
                f"{name}/{window}: signed flux integral mismatch")
        total = arrays[f"{window}_E_psd_V2_m2_per_eV"].sum(axis=1)
        require(np.isclose(total[0] * de, review["dc_E_mean_square_V2_m2"], rtol=2e-12, atol=0),
                f"{name}/{window}: DC integral mismatch")
    return dict(name=name, summary=summary, arrays=arrays, directory=directory,
                input_hashes={path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                              for path in (summary_path, array_path)})


def plotting():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.titlesize": 12, "axes.labelsize": 10,
                         "text.color": TEXT, "axes.labelcolor": TEXT,
                         "xtick.color": MUTED, "ytick.color": MUTED,
                         "axes.edgecolor": "#afbac5", "axes.spines.top": False,
                         "axes.spines.right": False, "figure.facecolor": "white",
                         "savefig.facecolor": "white", "pdf.fonttype": 42,
                         "path.simplify": False})
    return plt


def page_title(fig, title, subtitle):
    fig.text(.065, .973, title, fontsize=18, weight="bold", va="top")
    fig.text(.065, .928, subtitle, fontsize=10.3, color=MUTED, va="top")


def footnote(fig, lines, y=.025):
    fig.text(.065, y, "\n".join(lines), fontsize=8.8, color=MUTED,
             va="bottom", linespacing=1.45)


def save(plt, fig, output, name):
    paths = []
    for suffix in ("png", "pdf"):
        path = output / f"{name}.{suffix}"
        fig.savefig(path, dpi=190, bbox_inches="tight", pad_inches=.12)
        paths.append(str(path.resolve()))
    plt.close(fig)
    return paths


def curve_label(run, window):
    return LABELS[run["name"]] + (" | rectangular" if window == "rect" else " | Hann / sqrt(U)")


def positive(values):
    # Log axes cannot display exact zeros; do not invent a positive noise floor.
    return np.ma.masked_less_equal(values, 0.)


def psd_axes(ax, low, high, title):
    ax.set(title=title, xlabel=r"Photon-energy coordinate $h\nu$ [eV]", ylabel=PSD_UNIT,
           xlim=(low, high), yscale="log")
    ax.grid(alpha=.18, linewidth=.6)


def total_curves(ax, runs, low, high, kind="E", windows=("rect", "hann"), legend=True):
    for run in runs:
        energy = run["arrays"]["photon_energy_eV"]
        selected = (energy > 0) & (energy >= low) & (energy <= high)
        for window in windows:
            values = run["arrays"][f"{window}_{kind}_psd_V2_m2_per_eV"].sum(axis=1) / 1e12
            ax.plot(energy[selected], positive(values[selected]),
                    color=COLORS[run["name"]], ls="-" if window == "rect" else "--",
                    lw=.8 if window == "rect" else 1.05, label=curve_label(run, window))
    if legend:
        ax.legend(frameon=False, fontsize=8.5, loc="upper right")


def sampling_lines(runs):
    text = []
    for run in runs:
        s = run["summary"]
        text.append(f"{LABELS[run['name']]}: {s['recorded_span_fs']:.2f} fs, "
                    f"ΔE = {s['frequency_bin_eV']:.5f} eV, "
                    f"Nyquist = {s['theoretical_nyquist_eV']:.2f} eV")
    return "  |  ".join(text)


def config_note(runs):
    z = [r["summary"]["source_attributes"]["plane_z_m"] for r in runs]
    return (f"18893 plane z = {z[0]:.6f} m; original plane z = {z[1]:.6f} m. "
            "Startup, boundaries and other settings also differ: this is not an isolated retirement comparison.")


def overview(plt, runs, output):
    fig, axes = plt.subplots(2, 2, figsize=(14.5, 10.5))
    fig.subplots_adjust(left=.075, right=.97, bottom=.22, top=.84, hspace=.47, wspace=.25)
    high = max(r["arrays"]["photon_energy_eV"][-1] for r in runs)
    for ax, kind, title in ((axes[0, 0], "E", r"(a) Total electric spectrum: $D_{E_x}+D_{E_y}+D_{E_z}$"),
                            (axes[0, 1], "cB", r"(b) Total magnetic spectrum in electric units: $D_{cB}$")):
        psd_axes(ax, 0., high, title)
        ax.axvspan(50., 114.3520538, color="#e2e8ef", alpha=.45, zorder=-5)
        total_curves(ax, runs, 0., high, kind)
    psd_axes(axes[1, 0], 0., 12., "(c) Low-frequency electric field | DC separated")
    total_curves(axes[1, 0], runs, 0., 12.)
    axes[1, 0].axvline(1., color=MUTED, ls=":", lw=.7)
    dc = axes[1, 1]
    heights, labels, colors, fractions = [], [], [], []
    for run in runs:
        for window in ("rect", "hann"):
            review = run["summary"]["windows"][window]
            heights.append(review["dc_E_mean_square_V2_m2"] / 1e12)
            labels.append(LABELS[run["name"]] + "\n" + ("Rectangular" if window == "rect" else "Hann / sqrt(U)"))
            colors.append(COLORS[run["name"]])
            fractions.append(review["dc_E_fraction"])
    bars = dc.bar(np.arange(4), heights, color=colors, width=.61, alpha=.85)
    for i, bar in enumerate(bars):
        if i % 2:
            bar.set_hatch("///")
        dc.annotate(f"{heights[i]:.4g}\n{fractions[i]:.2%} of mean square",
                    (bar.get_x() + bar.get_width() / 2, bar.get_height()),
                    xytext=(0, 5), textcoords="offset points", ha="center", va="bottom", fontsize=8)
    dc.set(title="(d) Electric DC bin | integrated bin, not density",
           ylabel=r"DC mean square [$(\mathrm{MV/m})^2$]", xticks=np.arange(4), xticklabels=labels)
    dc.tick_params(axis="x", labelsize=8.5)
    dc.set_ylim(0, max(heights) * 1.33 if max(heights) > 0 else 1.)
    dc.grid(axis="y", alpha=.18)
    page_title(fig, "Whole-aperture spectra of the stored raw fields",
               "140 × 140 physical cells per plane  |  each pixel FFT first, then mean-square averaging  |  no background subtraction")
    low_frequency_note = "Raw E mean-square fraction below 1 eV, including DC: " + "  |  ".join(
        f"{LABELS[r['name']]} {r['summary']['windows']['rect']['bands']['below_1_eV']['E_fraction']:.4%}" for r in runs)
    footnote(fig, [sampling_lines(runs),
                   "Spectral density integrates to spatial/time mean-square field, not separated radiation energy. E and cB share the displayed units.",
                   "DC is retained and shown as its bin integral in (d); zero bins alone cannot appear on logarithmic axes. No mean, trend or spike was removed.",
                   low_frequency_note,
                   "Hann uses identical windows on all components, divided by sqrt(U), U = mean(window²). It tests finite-record sensitivity; it does not restore a missing tail.",
                   "Shading: 50–114.352 eV. Slab-transfer artifacts remain; nonzero endpoints and Nyquist do not certify physical fidelity.",
                   config_note(runs)])
    return save(plt, fig, output, "raw_field_spectra_overview")


def comb_marks(ax, runs, low, high):
    for row, run in enumerate(runs):
        spacing = float(run["summary"]["sampling_comb_spacing_eV"])
        first, last = max(1, int(np.ceil(low / spacing))), int(np.floor(high / spacing))
        for harmonic in range(first, last + 1):
            at = harmonic * spacing
            ax.axvline(at, color=COLORS[run["name"]], lw=.5, alpha=.18,
                       ls=":" if row == 0 else "--", zorder=-3)
            ax.plot([at, at], [.015 + .045 * row, .045 + .045 * row],
                    transform=ax.get_xaxis_transform(), color=COLORS[run["name"]], lw=1.2)


def band_view(plt, runs, output, nominal, zoom):
    fig, axes = plt.subplots(2, 2, figsize=(14.5, 10.2))
    fig.subplots_adjust(left=.075, right=.97, bottom=.22, top=.84, hspace=.43, wspace=.26)
    low, high = 50., 114.3520538
    zoom_low, zoom_high = nominal - zoom, nominal + zoom
    for row, window in enumerate(("rect", "hann")):
        window_title = "Raw rectangular window" if window == "rect" else "Energy-normalized Hann window"
        for column, (left, right) in enumerate(((low, high), (zoom_low, zoom_high))):
            ax = axes[row, column]
            psd_axes(ax, left, right, f"({'abcd'[row*2+column]}) {window_title}" + (" | resonance region" if column else ""))
            comb_marks(ax, runs, left, right)
            total_curves(ax, runs, left, right, windows=(window,), legend=False)
            ax.axvline(nominal, color="#52616b", lw=1, ls="-.", zorder=1)
            ax.legend(frameon=False, fontsize=8.5, loc="upper right")
            if column:
                for i, run in enumerate(runs):
                    energy = run["arrays"]["photon_energy_eV"]
                    values = run["arrays"][f"{window}_E_psd_V2_m2_per_eV"].sum(axis=1) / 1e12
                    indices = np.flatnonzero((energy >= left) & (energy <= right))
                    require(len(indices) > 0, "Resonance zoom does not contain an FFT sample")
                    index = int(indices[np.argmax(values[indices])])
                    ax.plot(energy[index], values[index], marker="o", color=COLORS[run["name"]], ms=3)
                    ax.text(.03, .94 - .07*i, f"{LABELS[run['name']]} local max: {energy[index]:.5f} eV",
                            transform=ax.transAxes, color=COLORS[run["name"]], fontsize=8.5,
                            va="top", bbox=dict(facecolor="white", edgecolor="none", alpha=.8, pad=1))
    comb = "  |  ".join(f"{LABELS[r['name']]} comb spacing: {r['summary']['sampling_comb_spacing_eV']:.5f} eV "
                         f"({r['summary']['slab_cells']} steps)" for r in runs)
    page_title(fig, "Raw electric-field spectrum in the target band",
               "Total E mean-square spectral density  |  separate native FFT grids, no peak deletion or interpolation")
    footnote(fig, [comb,
                   "Colored vertical marks show harmonics of the slab-transfer cadence; they flag possible contamination, not identified radiation lines.",
                   f"Gray dash-dot: nominal undulator resonance {nominal:.6f} eV. Markers identify each stored spectrum's local maximum in the zoom, not its physical origin.",
                   "Rectangular and normalized Hann windows preserve distinct finite-record weightings. Their agreement alone does not establish radiation purity.",
                   "The known detector interface issue can contaminate this band. Native field values, including all spikes and DC, remain unchanged.",
                   config_note(runs)])
    return save(plt, fig, output, "raw_field_spectra_target_band")


def components(plt, runs, output):
    fig, axes = plt.subplots(2, 2, figsize=(14.5, 10.2))
    fig.subplots_adjust(left=.075, right=.97, bottom=.20, top=.84, hspace=.43, wspace=.25)
    for row, run in enumerate(runs):
        energy = run["arrays"]["photon_energy_eV"]
        high = float(energy[-1])
        for column, kind in enumerate(("E", "cB")):
            ax = axes[row, column]
            psd_axes(ax, 0., high, f"({'abcd'[row*2+column]}) {LABELS[run['name']]} | {kind} components")
            ax.axvspan(50., 114.3520538, color="#e2e8ef", alpha=.45, zorder=-5)
            for window in ("rect", "hann"):
                for k, letter in enumerate("xyz"):
                    values = run["arrays"][f"{window}_{kind}_psd_V2_m2_per_eV"][1:, k] / 1e12
                    ax.plot(energy[1:], positive(values), color=COMPONENT_COLORS[k],
                            ls="-" if window == "rect" else "--", lw=.75,
                            label=f"{kind}{letter} | {'rect' if window == 'rect' else 'Hann'}")
            ax.legend(frameon=False, fontsize=8.2, ncol=2, loc="upper right")
    page_title(fig, "Supplement: all six field-component spectra",
               "E [V/m] and cB [V/m] treated identically  |  solid: raw rectangle; dashed: energy-normalized Hann")
    footnote(fig, ["Each component is Fourier transformed per pixel before its squared amplitude is averaged over the physical aperture.",
                   "DC is included in the numerical products and shown separately in the overview. No mean subtraction, detrending or spike filtering.",
                   "MPI-interface sampling artifacts and native Yee timing remain. Similar E/cB spectra do not prove that bound fields have been removed.",
                   config_note(runs)])
    return save(plt, fig, output, "raw_field_spectra_components")


def native_flux(plt, runs, output):
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 6.3))
    fig.subplots_adjust(left=.075, right=.97, bottom=.27, top=.80, wspace=.27)
    high = max(r["arrays"]["photon_energy_eV"][-1] for r in runs)
    for ax, left, right, title in ((axes[0], 0., high, "(a) Entire sampled band | DC bin separated"),
                                   (axes[1], 50., 114.3520538, "(b) Target band")):
        magnitude = 0.
        minimum = 0.
        for run in runs:
            energy = run["arrays"]["photon_energy_eV"]
            mask = (energy > 0) & (energy >= left) & (energy <= right)
            for window in ("rect", "hann"):
                values = run["arrays"][f"{window}_native_flux_J_per_eV"][mask] * 1e6
                magnitude = max(magnitude, float(np.max(np.abs(values))))
                minimum = min(minimum, float(np.min(values)))
                ax.plot(energy[mask], values, color=COLORS[run["name"]],
                        ls="-" if window == "rect" else "--", lw=.85, label=curve_label(run, window))
        ax.set(title=title, xlim=(left, right), xlabel=r"Photon-energy coordinate $h\nu$ [eV]",
               ylabel=r"Native signed spectral flux [$\mu$J/eV]")
        ax.set_yscale("symlog", linthresh=max(magnitude * 1e-6, 1e-30))
        if minimum == 0.:
            ax.set_ylim(bottom=0.)
        ax.axhline(0., color=MUTED, lw=.6)
        ax.grid(alpha=.18, linewidth=.6)
        ax.legend(frameon=False, fontsize=8.5, loc="upper right")
    page_title(fig, "Supplement: native spectral Poynting-flux estimate",
               r"Signed $\mathrm{Re}(\widetilde E_x\widetilde B_y^*-\widetilde E_y\widetilde B_x^*)/\mu_0$  |  finite aperture and recorded window only")
    dc_notes = []
    for run in runs:
        de = run["summary"]["frequency_bin_eV"]
        dc_rect = run["arrays"]["rect_native_flux_J_per_eV"][0] * de * 1e6
        dc_hann = run["arrays"]["hann_native_flux_J_per_eV"][0] * de * 1e6
        dc_notes.append(f"{LABELS[run['name']]} DC-bin flux: rectangular {dc_rect:.6g} µJ; Hann {dc_hann:.6g} µJ")
    footnote(fig, ["  |  ".join(dc_notes),
                   "Signed values are retained; the vertical scale is logarithmic away from zero and linear very close to zero.",
                   "B retains a Yee half-step stagger; the Lorentz transform mixes time levels. This is not exact simultaneous-field flux or separated radiation energy.",
                   "Sampling defects and nonzero record endpoints remain. Hann is a window-sensitivity check, not a correction.",
                   config_note(runs)])
    return save(plt, fig, output, "raw_field_spectra_native_flux")


def peak_record(run, window, low, high):
    energy = run["arrays"]["photon_energy_eV"]
    values = run["arrays"][f"{window}_E_psd_V2_m2_per_eV"].sum(axis=1)
    indices = np.flatnonzero((energy >= low) & (energy <= high))
    require(indices.size > 0, "Peak report interval has no FFT sample")
    i = int(indices[np.argmax(values[indices])])
    return dict(energy_eV=float(energy[i]), PSD_V2_m2_per_eV=float(values[i]), bin_index=i)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-dir", type=Path, help="Default: INPUT/figures")
    parser.add_argument("--nominal-resonance-eV", type=float, default=NOMINAL_EV)
    parser.add_argument("--zoom-half-width-eV", type=float, default=1.5)
    parser.add_argument("--no-native-flux", action="store_true")
    args = parser.parse_args()
    require(args.zoom_half_width_eV > 0 and args.nominal_resonance_eV > args.zoom_half_width_eV,
            "Invalid nominal resonance/zoom interval")
    root = args.input_dir.resolve()
    output = args.output_dir.resolve() if args.output_dir else root / "figures"
    runs = [load_run(root, name) for name in ("retirement", "original")]
    require(np.isclose(runs[0]["summary"]["physical_area_m2"],
                       runs[1]["summary"]["physical_area_m2"], rtol=2e-13, atol=0),
            "The two physical aperture areas differ")
    output.mkdir(parents=True, exist_ok=True)
    plt = plotting()
    paths = overview(plt, runs, output)
    paths += band_view(plt, runs, output, args.nominal_resonance_eV, args.zoom_half_width_eV)
    paths += components(plt, runs, output)
    if not args.no_native_flux:
        paths += native_flux(plt, runs, output)
    audit = dict(source_hdf5_opened=False, values_filtered=False, zero_bins_on_log_axes="masked for display only",
                 nominal_resonance_eV=args.nominal_resonance_eV,
                 nominal_line_meaning="Card-based undulator resonance reference, not radiation identification",
                 script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), figures=paths, runs={})
    for run in runs:
        s = run["summary"]
        record = dict(input_directory=str(run["directory"]), input_sha256=run["input_hashes"],
                      source_file=s["source_file"], samples=s["samples"], physical_shape_yx=s["physical_shape_yx"],
                      frequency_bin_eV=s["frequency_bin_eV"], theoretical_nyquist_eV=s["theoretical_nyquist_eV"],
                      last_frequency_eV=s["last_frequency_eV"], sampling_comb_spacing_eV=s["sampling_comb_spacing_eV"], windows={})
        for window in ("rect", "hann"):
            record["windows"][window] = dict(dc_E_fraction=s["windows"][window]["dc_E_fraction"],
                peak_target_band=peak_record(run, window, 50., 114.3520538),
                peak_resonance_zoom=peak_record(run, window, args.nominal_resonance_eV-args.zoom_half_width_eV,
                                                args.nominal_resonance_eV+args.zoom_half_width_eV))
        audit["runs"][run["name"]] = record
    audit_path = output / "SPECTRA_PLOT_AUDIT.json"
    audit_path.write_text(json.dumps(audit, indent=2, allow_nan=False) + "\n")
    print(json.dumps(dict(figures=paths, audit=str(audit_path)), indent=2))


if __name__ == "__main__":
    main()
