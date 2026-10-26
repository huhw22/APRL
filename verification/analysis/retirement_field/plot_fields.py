#!/usr/bin/env python3
"""Plot small APRL field reductions; never reopen the detector HDF5 file.

Inputs are lab_field_summary.json, lab_field_timeseries.npz and
lab_field_maps.npz from analyze_fields.py. Each figure is saved as PNG and PDF.
The physical aperture excludes transverse CPML, and all time integrals cover
only the recorded interval. No near-field subtraction is implied.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


COLORS = ("#2563a6", "#d97825", "#347f58")
TEXT = "#202c38"
MUTED = "#526272"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def checked_array(source, key, shape):
    array = np.asarray(source[key], dtype=float)
    require(array.shape == shape, f"{key}: expected {shape}, got {array.shape}")
    require(np.isfinite(array).all(), f"{key}: non-finite values")
    return array


def load_reductions(directory):
    summary = json.loads((directory / "lab_field_summary.json").read_text())
    with np.load(directory / "lab_field_timeseries.npz", allow_pickle=False) as archive:
        series = {key: archive[key] for key in archive.files}
    with np.load(directory / "lab_field_maps.npz", allow_pickle=False) as archive:
        maps = {key: archive[key] for key in archive.files}
    n, nx, ny = (int(summary[key]) for key in ("committed_samples", "nx", "ny"))
    require(n >= 2, "At least two detector samples are required")
    require(summary["frame"] == "laboratory", "Expected laboratory-frame fields")
    require(summary["field_is_raw_total"] is True, "Expected raw total-field reductions")
    require(summary["background_subtraction"] == "none", "Unexpected background subtraction")
    x = checked_array(maps, "x_m", (nx,))
    y = checked_array(maps, "y_m", (ny,))
    dx, dy = (float(summary[key]) for key in ("dx_m", "dy_m"))
    require(dx > 0 and dy > 0, "Pixel spacing must be positive")
    require(np.allclose(np.diff(x), dx, rtol=1e-10, atol=dx * 1e-12), "Invalid x spacing")
    require(np.allclose(np.diff(y), dy, rtol=1e-10, atol=dy * 1e-12), "Invalid y spacing")
    cx, cy = (int(v) for v in summary["transverse_cpml_cells"])
    require(0 <= 2 * cx < nx and 0 <= 2 * cy < ny, "Invalid CPML exclusion")
    view = (slice(cy, ny - cy), slice(cx, nx - cx))
    xp, yp = x[cx:nx-cx], y[cy:ny-cy]
    extent = [xp[0] - dx / 2, xp[-1] + dx / 2,
              yp[0] - dy / 2, yp[-1] + dy / 2]
    extent = (np.asarray(extent) * 1e6).tolist()
    times = checked_array(series, "time_lab_s", (n,))
    require(np.all(np.diff(times) > 0), "Detector times must increase strictly")
    relative = checked_array(series, "time_from_first_fs", (n,))
    require(np.allclose(relative, (times - times[0]) * 1e15, rtol=1e-12, atol=1e-10),
            "Relative times disagree with source timestamps")
    for key in ("physical_peak_E_V_per_m", "physical_rms_E_V_per_m",
                "physical_net_power_W", "physical_positive_power_W", "physical_negative_power_W"):
        checked_array(series, key, (n,))
    for letter in "xyz":
        checked_array(series, f"physical_rms_E{letter}_V_per_m", (n,))
        checked_array(series, f"axis_E{letter}_V_per_m", (n,))
    peak_map = checked_array(maps, "peak_E_magnitude_V_per_m", (ny, nx))
    time_rms = checked_array(maps, "time_rms_E_components_V_per_m", (ny, nx, 3))
    frame = checked_array(maps, "peak_frame_E_V_per_m", (ny, nx, 3))
    checked_array(maps, "native_signed_fluence_J_per_m2", (ny, nx))
    require(np.all(peak_map >= 0) and np.all(time_rms >= 0), "Negative magnitude/RMS")
    peak = float(peak_map[view].max())
    require(np.isclose(peak, series["physical_peak_E_V_per_m"].max(), rtol=2e-12, atol=0),
            "Physical peak map and timeseries disagree")
    require(np.isclose(peak, summary["physical_aperture"]["peak_E_magnitude_V_per_m"],
                       rtol=2e-12, atol=0), "Physical peak and summary disagree")
    component_squared = sum(series[f"physical_rms_E{letter}_V_per_m"] ** 2 for letter in "xyz")
    require(np.allclose(component_squared, series["physical_rms_E_V_per_m"] ** 2,
                        rtol=2e-12, atol=0), "Vector/component RMS disagreement")
    require(np.isclose(np.sqrt(np.sum(time_rms[view] ** 2, axis=2).mean()),
                       summary["physical_aperture"]["space_time_rms_E_V_per_m"],
                       rtol=2e-12, atol=0), "Time-RMS map and summary disagree")
    peak_index = int(summary["peak_E_sample_index"])
    require(0 <= peak_index < n, "Peak sample index outside recorded interval")
    peak_y, peak_x = np.unravel_index(np.argmax(np.sum(frame * frame, axis=2)), (ny, nx))
    require(np.allclose([x[peak_x], y[peak_y]], summary["peak_E_xy_m"],
                        rtol=1e-12, atol=min(dx, dy) * 1e-9), "Peak frame location mismatch")
    in_physical = cx <= peak_x < nx - cx and cy <= peak_y < ny - cy
    physical_frame_peak = float(np.sqrt(np.sum(frame[view] ** 2, axis=2).max()))
    frame_is_physical_peak = bool(in_physical and np.isclose(physical_frame_peak, peak, rtol=2e-12, atol=0))
    dominant = int(np.argmax(np.sum(frame[view] ** 2, axis=(0, 1))))
    central = np.stack([series[f"axis_E{letter}_V_per_m"] for letter in "xyz"], axis=1)
    central_peak_index = int(np.argmax(np.sum(central ** 2, axis=1)))
    central_cells = len(summary["central_sample_x_m"]) * len(summary["central_sample_y_m"])
    require(central_cells > 0, "Missing central-cell coordinates")
    geometry = dict(view=view, extent=extent, nx=len(xp), ny=len(yp),
                    peak_frame_is_physical_peak=frame_is_physical_peak,
                    full_peak_in_physical=in_physical, dominant_component=dominant,
                    central_peak_index=central_peak_index, central_cells=central_cells)
    return summary, series, maps, geometry


def setup_plotting():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 12,
        "axes.labelsize": 10, "axes.edgecolor": "#acb7c2", "axes.labelcolor": TEXT,
        "text.color": TEXT, "xtick.color": MUTED, "ytick.color": MUTED,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": False, "figure.facecolor": "white", "axes.facecolor": "white",
        "savefig.facecolor": "white", "pdf.fonttype": 42,
        "path.simplify": False, "image.interpolation": "nearest",
    })
    return plt


def page_title(fig, title, subtitle):
    fig.text(.065, .973, title, fontsize=18, weight="bold", va="top")
    fig.text(.065, .928, subtitle, fontsize=10.5, color=MUTED, va="top")


def caption(fig, text, y=.025, sampling_note=""):
    if sampling_note:
        text += "\n" + sampling_note
    fig.text(.065, y, text, fontsize=9, color=MUTED, va="bottom", linespacing=1.5)


def heatmap(fig, ax, values, extent, title, unit, signed=False):
    limit = max(float(np.max(np.abs(values))), np.finfo(float).tiny)
    options = dict(vmin=-limit if signed else 0, vmax=limit,
                   cmap="RdBu_r" if signed else "magma")
    im = ax.imshow(values, origin="lower", extent=extent, aspect="equal",
                   interpolation="nearest", rasterized=True, **options)
    ax.set(title=title, xlabel=r"$x$ [$\mu$m]", ylabel=r"$y$ [$\mu$m]")
    bar = fig.colorbar(im, ax=ax, fraction=.046, pad=.04)
    bar.set_label(unit)
    return im


def time_axis(ax, start, end):
    ax.set_xlim(start, end)
    ax.set_xlabel("Laboratory time from first recorded sample [fs]")
    ax.grid(axis="both", alpha=.19, linewidth=.6)


def save_figure(plt, fig, out, name):
    paths = []
    for suffix in ("png", "pdf"):
        path = out / f"{name}.{suffix}"
        fig.savefig(path, dpi=190, bbox_inches="tight", pad_inches=.12)
        paths.append(str(path.resolve()))
    plt.close(fig)
    return paths


def plot_maps(plt, summary, series, maps, geometry, output, sampling_note=""):
    fig, axes = plt.subplots(1, 3, figsize=(15, 6.1))
    fig.subplots_adjust(left=.065, right=.965, bottom=.23, top=.79, wspace=.44)
    view, extent = geometry["view"], geometry["extent"]
    peak = maps["peak_E_magnitude_V_per_m"][view] / 1e6
    rms = np.sqrt(np.sum(maps["time_rms_E_components_V_per_m"][view] ** 2, axis=2)) / 1e6
    component = geometry["dominant_component"]
    signed = maps["peak_frame_E_V_per_m"][view][:, :, component] / 1e6
    peak_time = float(series["time_from_first_fs"][summary["peak_E_sample_index"]])
    heatmap(fig, axes[0], peak, extent, r"(a) Peak over time: $\max_t |\mathbf{E}|$", "MV/m")
    heatmap(fig, axes[1], rms, extent,
            r"(b) Time RMS: $\sqrt{T^{-1}\int |\mathbf{E}|^2\,dt}$", "MV/m")
    heatmap(fig, axes[2], signed, extent,
            f"(c) Signed $E_{'xyz'[component]}$ at {peak_time:.3f} fs", "MV/m", signed=True)
    page_title(fig, "Electric field across the physical aperture",
               f"Laboratory plane z = {summary['attributes']['plane_z_m']:.6f} m  |  "
               f"{geometry['nx']} × {geometry['ny']} cells, CPML excluded  |  "
               f"recorded interval: {summary['window_duration_fs']:.2f} fs")
    selected = ("The full-aperture peak lies inside this aperture; panel (c) is also its physical-aperture peak frame."
                if geometry["peak_frame_is_physical_peak"] else
                "Panel (c) uses the full recorded aperture's peak frame; it is not labelled as the physical-aperture peak.")
    caption(fig,
            "(a) Each pixel may peak at a different time.  (b) RMS includes the entire recorded interval, including quiet time.\n"
            + selected + "\n"
            "Raw Maxwell fields; complete near-field removal is unproven. Native Lorentz transforms mix Yee-staggered E/B.\n"
            "Pixel cells are shown without interpolation. The recorded interval does not establish complete pulse capture.",
            sampling_note=sampling_note)
    return save_figure(plt, fig, output, "lab_electric_field_physical_maps")


def plot_pulse(plt, summary, series, maps, geometry, output, sampling_note=""):
    fig = plt.figure(figsize=(13.8, 10.2))
    grid = fig.add_gridspec(2, 2, left=.075, right=.965, bottom=.20, top=.85,
                           hspace=.46, wspace=.28, height_ratios=[1, 1])
    top = fig.add_subplot(grid[0, :])
    components = fig.add_subplot(grid[1, 0])
    waveform = fig.add_subplot(grid[1, 1])
    t = series["time_from_first_fs"]
    emax = series["physical_peak_E_V_per_m"] / 1e6
    erms = series["physical_rms_E_V_per_m"] / 1e6
    top.plot(t, emax, color=COLORS[0], lw=.85, label=r"Spatial maximum $\max_A |\mathbf{E}|$")
    top.plot(t, erms, color=COLORS[1], lw=1, label=r"Spatial RMS $\sqrt{\langle|\mathbf{E}|^2\rangle_A}$")
    top.set(title="(a) Field strength over the full recorded interval", ylabel="Electric field [MV/m]")
    time_axis(top, t[0], t[-1])
    top.set_ylim(bottom=0)
    top.axvline(t[-1], color=MUTED, lw=1, ls="--", clip_on=False)
    top.plot(t[-1], emax[-1], "o", color=COLORS[0], ms=4, clip_on=False)
    fraction = float(emax[-1] / emax.max()) if emax.max() else 0.
    top.annotate(f"Recorded end: {t[-1]:.2f} fs\nLast max |E| = {fraction:.1%} of peak",
                 xy=(t[-1], emax[-1]), xytext=(-10, 26), textcoords="offset points",
                 ha="right", fontsize=9, color=MUTED,
                 bbox=dict(facecolor="white", edgecolor="none", alpha=.9, pad=2))
    top.legend(loc="upper left", frameon=False, fontsize=9)
    for letter, color in zip("xyz", COLORS):
        components.plot(t, series[f"physical_rms_E{letter}_V_per_m"] / 1e6,
                        label=f"$E_{letter}$ RMS", color=color, lw=.8)
    components.set(title="(b) Spatial RMS of each component", ylabel="Component RMS [MV/m]")
    time_axis(components, t[0], t[-1])
    components.set_ylim(bottom=0)
    components.axvline(t[-1], color=MUTED, lw=1, ls="--", clip_on=False)
    components.legend(frameon=False, fontsize=9)
    centre = float(t[geometry["central_peak_index"]])
    width = min(1., float(t[-1] - t[0]))
    left = max(float(t[0]), min(centre - width / 2, float(t[-1]) - width))
    right = left + width
    selection = (t >= left) & (t <= right)
    require(selection.sum() >= 3, "Fewer than three samples in central waveform zoom")
    for letter, color in zip("xyz", COLORS):
        waveform.plot(t[selection], series[f"axis_E{letter}_V_per_m"][selection] / 1e6,
                      label=f"$E_{letter}$", color=color, lw=1)
    waveform.axhline(0, lw=.55, color=MUTED)
    waveform.set(title=f"(c) Central {geometry['central_cells']}-cell signed waveform | {width:g} fs",
                 ylabel="Mean central electric field [MV/m]")
    time_axis(waveform, left, right)
    waveform.ticklabel_format(axis="x", useOffset=False)
    waveform.legend(frameon=False, fontsize=9)
    page_title(fig, "Electric-field evolution at the detector",
               f"Physical aperture: {geometry['nx']} × {geometry['ny']} cells  |  "
               f"z = {summary['attributes']['plane_z_m']:.6f} m  |  "
               f"{summary['committed_samples']:,} native samples, no time averaging")
    caption(fig,
            f"Time zero is the first detector sample (absolute lab time {summary['time_start_s']:.16e} s).\n"
            "(c) Averages the nearest cell-centre vectors, then zooms around that averaged vector's largest magnitude; it is not an axis sample.\n"
            "Raw Maxwell fields; complete near-field removal is unproven. Native Lorentz transforms mix Yee-staggered E/B.\n"
            "The signal remains at the recorded end. These traces cover the saved interval, not a demonstrated complete pulse.",
            y=.035, sampling_note=sampling_note)
    return save_figure(plt, fig, output, "lab_electric_field_physical_traces")


def plot_flux(plt, summary, series, maps, geometry, output, sampling_note=""):
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 6.5), gridspec_kw={"width_ratios": [1.35, 1]})
    fig.subplots_adjust(left=.075, right=.965, bottom=.25, top=.80, wspace=.35)
    t = series["time_from_first_fs"]
    power = series["physical_net_power_W"] / 1e6
    axes[0].plot(t, power, color=COLORS[0], lw=.8, label="Native signed net power")
    axes[0].set(title=r"(a) $P_z(t)=\int_A (E_x B_y-E_y B_x)/\mu_0\,dA$", ylabel="Native net power [MW]")
    time_axis(axes[0], t[0], t[-1])
    axes[0].axhline(0, color=MUTED, lw=.6)
    axes[0].axvline(t[-1], color=MUTED, lw=1, ls="--", clip_on=False)
    axes[0].plot(t[-1], power[-1], "o", color=COLORS[0], ms=4, clip_on=False)
    axes[0].annotate(f"Recorded end\n{power[-1]:.3f} MW", xy=(t[-1], power[-1]),
                     xytext=(-8, 24), textcoords="offset points", ha="right", fontsize=9, color=MUTED,
                     bbox=dict(facecolor="white", edgecolor="none", alpha=.9, pad=2))
    fluence = maps["native_signed_fluence_J_per_m2"][geometry["view"]]
    heatmap(fig, axes[1], fluence, geometry["extent"],
            r"(b) Recorded-window fluence: $\int S_z\,dt$", r"J/m$^2$", signed=True)
    energy = float(summary["physical_aperture"]["net_energy_J"]) * 1e6
    page_title(fig, "Supplement: native Poynting-flux estimate",
               f"Physical aperture only  |  signed, no negative-value clipping  |  "
               f"recorded-window integral: {energy:.6f} µJ")
    caption(fig,
            "Native means E and B exactly as stored. B has Yee half-step staggering; the laboratory transform mixes those time levels.\n"
            "This is not an exact simultaneous-field flux or a measurement of separated radiation. Near-field removal is unproven.\n"
            "The integral covers only this finite aperture and the recorded time window. Power is nonzero at the recorded end.\n"
            "No interpolation or temporal smoothing is applied; a positive signed flux alone does not establish radiation purity.",
            sampling_note=sampling_note)
    return save_figure(plt, fig, output, "lab_native_flux_supplement")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True, help="Directory containing the three small reduction files")
    parser.add_argument("--output-dir", type=Path, help="PNG/PDF destination; default: INPUT/figures")
    parser.add_argument("--no-flux", action="store_true", help="Omit the separately labelled native-flux supplement")
    parser.add_argument("--sampling-note", default="", help="Optional sampling caveat shown on every figure and saved in PLOT_AUDIT.json")
    args = parser.parse_args()
    input_dir = args.input_dir.resolve()
    output = args.output_dir.resolve() if args.output_dir else input_dir / "figures"
    summary, series, maps, geometry = load_reductions(input_dir)
    output.mkdir(parents=True, exist_ok=True)
    plt = setup_plotting()
    paths = plot_maps(plt, summary, series, maps, geometry, output, args.sampling_note)
    paths += plot_pulse(plt, summary, series, maps, geometry, output, args.sampling_note)
    if not args.no_flux:
        paths += plot_flux(plt, summary, series, maps, geometry, output, args.sampling_note)
    inputs = {}
    for name in ("lab_field_summary.json", "lab_field_timeseries.npz", "lab_field_maps.npz"):
        path = input_dir / name
        inputs[name] = dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    audit = dict(inputs=inputs, source_hdf5_opened=False,
                 sampling_note=args.sampling_note,
                 plotting_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                 physical_aperture_shape_yx=[geometry["ny"], geometry["nx"]],
                 physical_aperture_edges_um=geometry["extent"],
                 full_aperture_peak_inside_physical=bool(geometry["full_peak_in_physical"]),
                 selected_frame_is_physical_peak=geometry["peak_frame_is_physical_peak"],
                 selected_signed_component="xyz"[geometry["dominant_component"]],
                 central_cells=geometry["central_cells"],
                 central_zoom_peak_sample=geometry["central_peak_index"],
                 interpolation="nearest cells, native time samples, no temporal averaging",
                 time_rms_definition="sqrt(trapezoidal integral of |E|^2 over recorded time / recorded duration)",
                 peak_map_definition="max over recorded time of sqrt(Ex^2+Ey^2+Ez^2) at each cell",
                 spatial_rms_definition="sqrt(mean over physical cells of Ex^2+Ey^2+Ez^2)",
                 notes=["Raw total laboratory Maxwell field, no background subtraction.",
                        "Native laboratory E/B inherit Yee time staggering through the Lorentz transform.",
                        "Recorded interval and finite aperture do not establish complete pulse capture."],
                 figures=paths)
    (output / "PLOT_AUDIT.json").write_text(json.dumps(audit, indent=2, allow_nan=False) + "\n")
    print(json.dumps(dict(figures=paths, audit=str(output / "PLOT_AUDIT.json")), indent=2))


if __name__ == "__main__":
    main()
