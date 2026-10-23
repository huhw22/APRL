#!/usr/bin/env python3
"""Stream a completed APRL background-subtracted field and plot its Poynting flux.

Reads only the reconstructed field file, never the original large field cube.
The raw and cleaned full-aperture power diagnostics are independently checked.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import tempfile

import h5py
import numpy as np

MU0 = 1.25663706212e-6
C = 299792458.0
POWER_NAMES = ("raw_signed_power_W", "cleaned_signed_power_W",
               "raw_forward_power_W", "cleaned_forward_power_W")
NOTES = [
    "Cleaned fields are the saved raw fields minus the prescribed uniform-velocity particle-background model, not an exact replay of the PIC self-field history.",
    "All fluxes use S=(E_clean cross B_clean)/mu0 at stored sample times and cell centres; the original Yee B time stagger is preserved by reconstruction.",
    "Positive and negative Sz are split pointwise before area integration; this is not a forward/backward modal decomposition.",
    "Raw-minus-cleaned power contains field cross terms and is not the independent energy of the subtracted particle near field.",
    "Integrals cover only the saved time interval and aperture. Padding/smoothing convergence and complete-pulse coverage are separate questions.",
    "The full recorded aperture includes the configured CPML edge cells. Physical-aperture diagnostics exclude --cpml-cells on every transverse edge.",
    "The displayed peak snapshot is selected at the maximum full-aperture cleaned signed power.",
]


def plain(value):
    if isinstance(value, bytes):
        return value.decode(errors="replace")
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value


def trapezoid(y, t):
    return np.sum((np.asarray(y)[1:] + np.asarray(y)[:-1]) * np.diff(t) / 2)


def save_csv(path, columns):
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(columns)
        writer.writerows(zip(*columns.values()))


def assert_close(actual, expected, name, absolute_scale):
    error = float(np.max(np.abs(np.asarray(actual) - np.asarray(expected))))
    tolerance = 2e-10 * max(float(absolute_scale), 1e-300)
    if error > tolerance:
        raise ValueError(f"{name} disagrees: max error {error:g}, tolerance {tolerance:g}")
    return {"max_absolute_error": error, "absolute_tolerance": tolerance}


def analyze(path, output, cpml=10, block_frames=8, plots=True):
    output.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "r") as f:
        g = f["reconstructed_field"]
        if int(g["complete"][()]) != 1:
            raise ValueError("Reconstructed field is incomplete")
        attrs = {key: plain(value) for key, value in g.attrs.items()}
        ed, bd = g["electric_V_per_m"], g["magnetic_T"]
        if ed.shape != bd.shape or len(ed.shape) != 4 or ed.shape[-1] != 3:
            raise ValueError("Expected matching [time,y,x,3] electric and magnetic fields")
        nt, ny, nx, _ = ed.shape
        if nt < 2 or min(nx, ny) <= 2 * cpml:
            raise ValueError("Insufficient samples or physical aperture after CPML exclusion")
        t = g["time_s"][:]
        if t.shape != (nt,) or not np.all(np.isfinite(t)) or not np.all(np.diff(t) > 0):
            raise ValueError("Time axis must match fields and be finite and strictly increasing")
        for key, value in (("samples", nt), ("nx", nx), ("ny", ny)):
            if key in attrs and int(attrs[key]) != value:
                raise ValueError(f"Geometry attribute {key} disagrees with dataset shape")
        dx, dy = float(attrs["dx_m"]), float(attrs["dy_m"])
        if min(dx, dy) <= 0:
            raise ValueError("Cell sizes must be positive")
        x = float(attrs["x_first_m"]) + np.arange(nx) * dx
        y = float(attrs["y_first_m"]) + np.arange(ny) * dy
        area = dx * dy
        dt = np.diff(t)
        weights = np.empty(nt)
        weights[0], weights[-1] = dt[0]/2, dt[-1]/2
        weights[1:-1] = (dt[:-1] + dt[1:])/2
        duration = t[-1] - t[0]
        physical = (slice(cpml, ny-cpml), slice(cpml, nx-cpml))
        traces = {"time_lab_s": t, "time_from_first_fs": (t-t[0])*1e15}
        for name in POWER_NAMES:
            value = g[name][:]
            if value.shape != (nt,) or not np.isfinite(value).all():
                raise ValueError(f"Invalid stored power diagnostic: {name}")
            traces[name] = value
        for aperture in ("full", "physical"):
            for sign in ("signed", "positive", "negative"):
                traces[f"{aperture}_cleaned_{sign}_power_W"] = np.empty(nt)
        traces["cleaned_peak_E_V_per_m"] = np.empty(nt)
        traces["cleaned_peak_B_T"] = np.empty(nt)
        flux_integral = np.zeros((ny, nx, 3))
        positive_fluence = np.zeros((ny, nx))
        negative_fluence = np.zeros((ny, nx))
        peak_abs_flux = np.zeros((ny, nx, 3))
        snapshot_index = -1
        snapshot_power = -np.inf
        snapshot_e = snapshot_b = snapshot_s = None
        for start in range(0, nt, block_frames):
            stop = min(start + block_frames, nt)
            electric, magnetic = ed[start:stop], bd[start:stop]
            if not np.isfinite(electric).all() or not np.isfinite(magnetic).all():
                raise ValueError(f"Nonfinite reconstructed fields in frames {start}:{stop}")
            flux = np.cross(electric, magnetic) / MU0
            sz = flux[..., 2]
            for aperture, view in (("full", (slice(None), slice(None))), ("physical", physical)):
                selected = sz[:, view[0], view[1]]
                traces[f"{aperture}_cleaned_signed_power_W"][start:stop] = selected.sum(axis=(1, 2)) * area
                traces[f"{aperture}_cleaned_positive_power_W"][start:stop] = np.maximum(selected, 0).sum(axis=(1, 2)) * area
                traces[f"{aperture}_cleaned_negative_power_W"][start:stop] = np.minimum(selected, 0).sum(axis=(1, 2)) * area
            traces["cleaned_peak_E_V_per_m"][start:stop] = np.sqrt(np.sum(electric*electric, axis=3).max(axis=(1, 2)))
            traces["cleaned_peak_B_T"][start:stop] = np.sqrt(np.sum(magnetic*magnetic, axis=3).max(axis=(1, 2)))
            flux_integral += np.einsum("t,tyxc->yxc", weights[start:stop], flux)
            positive_fluence += np.einsum("t,tyx->yx", weights[start:stop], np.maximum(sz, 0))
            negative_fluence += np.einsum("t,tyx->yx", weights[start:stop], np.minimum(sz, 0))
            peak_abs_flux = np.maximum(peak_abs_flux, np.abs(flux).max(axis=0))
            local = int(np.argmax(traces["full_cleaned_signed_power_W"][start:stop]))
            this_power = float(traces["full_cleaned_signed_power_W"][start+local])
            if this_power > snapshot_power:
                snapshot_index, snapshot_power = start+local, this_power
                snapshot_e, snapshot_b, snapshot_s = electric[local].copy(), magnetic[local].copy(), flux[local].copy()
            if start % (block_frames*128) == 0:
                print(f"cleaned flux: {stop}/{nt} frames", flush=True)
        power_scale = np.max(np.abs(traces["full_cleaned_positive_power_W"]) + np.abs(traces["full_cleaned_negative_power_W"]))
        checks = {
            "cleaned_signed_power": assert_close(traces["full_cleaned_signed_power_W"], traces["cleaned_signed_power_W"], "cleaned signed power", power_scale),
            "cleaned_forward_power": assert_close(traces["full_cleaned_positive_power_W"], traces["cleaned_forward_power_W"], "cleaned forward power", power_scale),
        }
        energy = {}
        for aperture, view in (("full", (slice(None), slice(None))), ("physical", physical)):
            values = {sign+"_energy_J": float(trapezoid(traces[f"{aperture}_cleaned_{sign}_power_W"], t))
                      for sign in ("signed", "positive", "negative")}
            absolute_energy = abs(values["positive_energy_J"]) + abs(values["negative_energy_J"])
            checks[f"{aperture}_space_time_integral"] = assert_close(values["signed_energy_J"], flux_integral[view][..., 2].sum()*area, "space/time integration", absolute_energy)
            checks[f"{aperture}_signed_partition"] = assert_close(values["signed_energy_J"], values["positive_energy_J"]+values["negative_energy_J"], "positive/negative partition", absolute_energy)
            power = traces[f"{aperture}_cleaned_signed_power_W"]
            values.update(peak_signed_power_W=float(power.max()), minimum_signed_power_W=float(power.min()),
                          initial_signed_power_W=float(power[0]), final_signed_power_W=float(power[-1]),
                          aperture_area_m2=float(flux_integral[view].shape[0]*flux_integral[view].shape[1]*area))
            energy[aperture] = values
        for prefix in ("raw", "cleaned"):
            for kind in ("signed", "forward"):
                key = f"{prefix}_{kind}_energy_J"
                measured = float(trapezoid(traces[f"{prefix}_{kind}_power_W"], t))
                scale = float(trapezoid(np.abs(traces[f"{prefix}_{kind}_power_W"]), t))
                if key not in attrs:
                    raise ValueError(f"Missing reconstruction energy attribute: {key}")
                checks[key] = assert_close(measured, attrs[key], key, scale)
        raw_energy = float(attrs["raw_signed_energy_J"])
        clean_energy = energy["full"]["signed_energy_J"]
        traces["raw_minus_cleaned_signed_power_W"] = traces["raw_signed_power_W"] - traces["cleaned_signed_power_W"]
        summary = dict(source_file=str(path.resolve()), frame="laboratory", input_complete=True,
                       samples=nt, ny=ny, nx=nx, dx_m=dx, dy_m=dy, cpml_cells=cpml,
                       window_duration_fs=float(duration*1e15), attributes=attrs,
                       peak_snapshot_sample=snapshot_index,
                       peak_snapshot_time_lab_s=float(t[snapshot_index]),
                       peak_snapshot_time_from_first_fs=float((t[snapshot_index]-t[0])*1e15),
                       raw_full_signed_energy_J=raw_energy,
                       cleaned_full_signed_energy_J=clean_energy,
                       raw_minus_cleaned_full_signed_energy_J=raw_energy-clean_energy,
                       relative_cleaned_minus_raw_signed_energy=(clean_energy-raw_energy)/raw_energy if raw_energy else None,
                       global_cleaned_peak_E_V_per_m=float(traces["cleaned_peak_E_V_per_m"].max()),
                       global_cleaned_peak_B_T=float(traces["cleaned_peak_B_T"].max()),
                       last_peak_E_over_global_peak=float(traces["cleaned_peak_E_V_per_m"][-1]/traces["cleaned_peak_E_V_per_m"].max()) if traces["cleaned_peak_E_V_per_m"].max() else None,
                       apertures=energy, checks=checks, notes=NOTES)
    save_csv(output/"reconstructed_flux_timeseries.csv", traces)
    np.savez_compressed(output/"reconstructed_flux_timeseries.npz", **traces)
    maps = dict(x_m=x, y_m=y, cleaned_signed_fluence_J_per_m2=flux_integral[..., 2],
                cleaned_vector_time_integral_J_per_m2=flux_integral,
                cleaned_positive_fluence_J_per_m2=positive_fluence,
                cleaned_negative_fluence_J_per_m2=negative_fluence,
                cleaned_time_mean_Sz_W_per_m2=flux_integral[..., 2]/duration,
                cleaned_peak_abs_S_components_W_per_m2=peak_abs_flux,
                peak_power_snapshot_S_W_per_m2=snapshot_s,
                peak_power_snapshot_E_V_per_m=snapshot_e,
                peak_power_snapshot_B_T=snapshot_b)
    np.savez_compressed(output/"reconstructed_flux_maps.npz", **maps)
    (output/"reconstructed_flux_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False, default=plain, allow_nan=False)+"\n")
    if plots:
        make_plots(output, traces, maps, summary)
    print(f"Cleaned-flux analysis complete: {output}", flush=True)
    return summary


def make_plots(output, traces, maps, summary):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": .22})
    x, y = maps["x_m"]*1e6, maps["y_m"]*1e6
    dx, dy = summary["dx_m"]*1e6, summary["dy_m"]*1e6
    extent = [x[0]-dx/2, x[-1]+dx/2, y[0]-dy/2, y[-1]+dy/2]
    cpml = summary["cpml_cells"]
    peak_time = summary["peak_snapshot_time_from_first_fs"]
    def save(fig, stem):
        for ext in ("png", "pdf"):
            fig.savefig(output/f"{stem}.{ext}", dpi=170, bbox_inches="tight")
        plt.close(fig)
    def heat(ax, data, title, unit):
        negative = bool(np.min(data) < 0)
        lim = max(float(np.abs(data).max()), 1e-300)
        options = {"cmap": "RdBu_r", "vmin": -lim, "vmax": lim} if negative else {"cmap": "magma", "vmin": 0}
        im = ax.imshow(data, origin="lower", extent=extent, aspect="equal", **options)
        if cpml:
            ax.add_patch(Rectangle((extent[0]+cpml*dx, extent[2]+cpml*dy),
                                  (summary["nx"]-2*cpml)*dx, (summary["ny"]-2*cpml)*dy,
                                  fill=False, edgecolor="cyan", lw=.8))
        ax.set(xlabel="x [um]", ylabel="y [um]", title=title); ax.grid(False)
        plt.colorbar(im, ax=ax, label=unit, shrink=.82)
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5), constrained_layout=True)
    heat(axes[0], maps["cleaned_signed_fluence_J_per_m2"], "Time-integrated cleaned Sz", "J/m²")
    heat(axes[1], maps["peak_power_snapshot_S_W_per_m2"][..., 2]/1e15,
         f"Cleaned Sz at peak net power (t={peak_time:.3f} fs)", "10¹⁵ W/m²")
    fig.suptitle("Particle-background-subtracted laboratory energy flux | cyan: physical aperture\nFinite recorded time window; subtraction uses the specified uniform-velocity model", fontsize=12)
    save(fig, "nearfield_subtracted_energy_flux")
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True, constrained_layout=True)
    t = traces["time_from_first_fs"]
    axes[0].plot(t, traces["raw_signed_power_W"]/1e6, lw=.8, label="Raw, full aperture")
    axes[0].plot(t, traces["cleaned_signed_power_W"]/1e6, lw=.8, label="Background-subtracted, full aperture")
    axes[0].plot(t, traces["physical_cleaned_signed_power_W"]/1e6, lw=.7, ls="--", label="Background-subtracted, physical aperture")
    axes[0].set_ylabel("Signed net power [MW]"); axes[0].legend(fontsize=9)
    axes[1].plot(t, traces["raw_minus_cleaned_signed_power_W"]/1e6, color="#b2182b", lw=.8)
    axes[1].axhline(0, color="gray", lw=.6)
    axes[1].set(xlabel="Time from first recorded sample [fs]", ylabel="Raw minus cleaned power [MW]")
    fig.suptitle("Raw and particle-background-subtracted Poynting power\nTheir difference includes cross terms; it is not the separate near-field energy", fontsize=12)
    save(fig, "nearfield_subtracted_power")
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)
    for i, ax in enumerate(axes):
        heat(ax, maps["peak_power_snapshot_S_W_per_m2"][..., i]/1e15,
             f"Cleaned S{'xyz'[i]} at peak net power", "10¹⁵ W/m²")
    fig.suptitle(f"Poynting-vector components at t={peak_time:.3f} fs", fontsize=12)
    save(fig, "nearfield_subtracted_flux_components")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), constrained_layout=True)
    heat(axes[0], np.linalg.norm(maps["peak_power_snapshot_E_V_per_m"], axis=2)/1e9,
         "Cleaned |E| at peak net power", "GV/m")
    heat(axes[1], np.linalg.norm(maps["peak_power_snapshot_B_T"], axis=2),
         "Cleaned |B| at peak net power", "T")
    fig.suptitle("Particle-background-subtracted field amplitudes at the selected snapshot", fontsize=12)
    save(fig, "nearfield_subtracted_peak_fields")


def self_test(output=None):
    with tempfile.TemporaryDirectory(prefix="aprl-clean-flux-test-") as temp:
        base = Path(temp)
        path, out = base/"synthetic.h5", base/"products" if output is None else output
        nt, ny, nx = 5, 6, 8
        t = np.array([0., 1., 3., 4., 7.])*1e-15
        a = np.array([1., 2., -3., 4., 1.])*1e7
        e = np.zeros((nt, ny, nx, 3)); b = np.zeros_like(e)
        e[..., 0] = a[:, None, None]
        b[..., 1] = (np.abs(a)/C)[:, None, None]
        b[..., 2] = 0.2*b[..., 1]
        sz = (e[..., 0]*b[..., 1])/MU0
        area = 2e-6*3e-6
        signed = sz.sum(axis=(1, 2))*area
        forward = np.maximum(sz, 0).sum(axis=(1, 2))*area
        with h5py.File(path, "w") as f:
            g = f.create_group("reconstructed_field")
            g["complete"] = 1; g["time_s"] = t
            g["electric_V_per_m"] = e; g["magnetic_T"] = b
            g.attrs.update(dx_m=2e-6, dy_m=3e-6, x_first_m=-7e-6, y_first_m=-7.5e-6,
                           nx=nx, ny=ny, samples=nt, plane_z_m=1.)
            for prefix, factor in (("raw", 1.2), ("cleaned", 1.)):
                for kind, values in (("signed", signed), ("forward", forward)):
                    g[f"{prefix}_{kind}_power_W"] = factor*values
                    g.attrs[f"{prefix}_{kind}_energy_J"] = factor*trapezoid(values, t)
        result = analyze(path, out, cpml=1, block_frames=2, plots=output is not None)
        assert np.isclose(result["cleaned_full_signed_energy_J"], trapezoid(signed, t), rtol=1e-13)
        assert np.isclose(result["apertures"]["physical"]["signed_energy_J"], trapezoid(signed, t)*24/48, rtol=1e-13)
        assert result["apertures"]["full"]["negative_energy_J"] < 0
        assert result["peak_snapshot_sample"] == 3
        with h5py.File(path, "r+") as f:
            f["reconstructed_field/complete"][()] = 0
        try:
            analyze(path, base/"invalid", cpml=1, plots=False)
        except ValueError as exc:
            assert "incomplete" in str(exc)
        else:
            raise AssertionError("Incomplete input was accepted")
    print("SELF-TEST PASS: signed/reversed flux, nonuniform time, CPML aperture, stored power/energy checks, incomplete-file rejection")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reconstructed-file", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--cpml-cells", type=int, default=10)
    parser.add_argument("--block-frames", type=int, default=8)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test(args.output_dir)
        return
    if args.reconstructed_file is None or args.output_dir is None:
        parser.error("--reconstructed-file and --output-dir are required")
    if args.cpml_cells < 0 or args.block_frames < 1:
        parser.error("--cpml-cells must be nonnegative and --block-frames must be positive")
    analyze(args.reconstructed_file, args.output_dir, args.cpml_cells, args.block_frames)


if __name__ == "__main__":
    main()
